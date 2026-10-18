use super::*;
use gammalooprs::{
    graph::edge::PossibleParticle,
    processes::{CutId, ProcessCollection},
};
use std::{collections::BTreeMap, f64::consts::PI};

#[test]
fn thermal_scalar_bubble_matches_literature_in_all_cut_domains() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let state_path =
        get_tests_workspace_path().join(format!("thermal_scalar_bubble_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/scalar_bubble/run.toml",
        &["generate", "rest_decay"],
        Some(state_path),
        None,
        true,
    )?;
    cli.state =
        gammaloop_api::state::State::load(cli.cli_settings.state.folder.clone(), None, None)?;
    cli.state.activate_loaded_integrand_backends(false)?;

    // Independent integrated two-body phase space. See Ho--Scherrer,
    // arXiv:1503.03534, (A17--A18), and Nishikawa--Morimatsu--Hidaka,
    // hep-ph/0302098, (III.35--42). All entries are Gamma_rho=-Im Pi_R/E
    // for -g X B1 B2 with distinct daughters and g=1, not forward loss.
    for (block, [energy, momentum, mass1, mass2, temperature, target]) in [
        (
            "rest_decay",
            [2.0_f64, 0.0, 0.4, 0.7, 0.6, 0.012104501639179106],
        ),
        (
            "boosted_decay",
            [2.5, 1.5, 0.4, 0.7, 0.6, 0.009094016657615578],
        ),
        ("closed", [1.0, 0.0, 0.6, 0.8, 0.6, 0.0]),
        ("absorption", [1.0, 0.0, 0.3, 1.8, 0.6, 0.00688050878313337]),
        ("spacelike", [1.0, 2.0, 0.5, 0.5, 1.0, 0.010364773879410214]),
    ] {
        cli.run_command(&format!("run {block}"))?;
        let points = arr2(&[[0.6, 0.0, 0.8], [-0.6, 0.0, -0.8]]);
        for process in ["equilibrium", "thermal_excess", "vacuum"] {
            let (process_id, integrand_name) = cli.state.find_integrand_ref(
                Some(&ProcessRef::Unqualified(process.to_owned())),
                Some(&"LO".to_owned()),
            )?;
            let ProcessCollection::CrossSections(cross_sections) =
                &cli.state.process_list.processes[process_id].collection
            else {
                unreachable!()
            };
            let graph = &cross_sections["LO"].supergraphs[0];
            // Read only the graph's physical momentum routing. Roots, radial
            // derivatives, residues and thermal weights below are independent.
            let routes = graph
                .graph
                .iter_edges_of(&graph.cuts[CutId(0)].cut)
                .map(|(_, edge, data)| {
                    let signature = &graph.graph.loop_momentum_basis.edge_signatures[edge];
                    let pdg = match &data.data.particle {
                        PossibleParticle::Particle(particle) => particle.pdg_code,
                        _ => panic!("expected a model bath particle"),
                    };
                    let mass = match pdg {
                        9000002 => mass1,
                        9000003 => mass2,
                        _ => panic!("unexpected bath field"),
                    };
                    (
                        *signature.internal.iter().next().unwrap() as i8 as f64,
                        signature
                            .external
                            .iter()
                            .map(|sign| *sign as i8 as f64)
                            .sum::<f64>()
                            * momentum,
                        mass,
                        pdg,
                    )
                })
                .collect_vec();
            assert_eq!(routes.len(), 2);
            let discrete = ndarray::Array2::<usize>::zeros((points.nrows(), 0));
            let batch = EvaluateSamples {
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
            for (point, sample) in points.rows().into_iter().zip(batch.samples) {
                let energies = |t: f64| {
                    routes
                        .iter()
                        .map(|(sign, shift, mass, _)| {
                            ((t * sign * point[0]).powi(2)
                                + (t * sign * point[1]).powi(2)
                                + (t * sign * point[2] + shift).powi(2)
                                + mass.powi(2))
                            .sqrt()
                        })
                        .collect_vec()
                };
                let mut expected_channels = BTreeMap::new();
                for signs in [[1.0_f64, 1.0], [1.0, -1.0], [-1.0, 1.0]] {
                    let equation = |t| {
                        energies(t)
                            .iter()
                            .zip(signs)
                            .map(|(e, s)| e * s)
                            .sum::<f64>()
                            - energy
                    };
                    let mut low = 0.0;
                    let mut high = 1.0;
                    // These specified rays have at most one root for each
                    // sign assignment. A finite sign search plus bisection
                    // supplies the oracle without the production root solver.
                    for _ in 0..32 {
                        if equation(low) * equation(high) < 0.0 {
                            break;
                        }
                        low = high;
                        high *= 2.0;
                    }
                    if equation(low) * equation(high) >= 0.0 {
                        continue;
                    }
                    for _ in 0..80 {
                        let mid = (low + high) / 2.0;
                        if equation(low) * equation(mid) <= 0.0 {
                            high = mid;
                        } else {
                            low = mid;
                        }
                    }
                    let t = (low + high) / 2.0;
                    let es = energies(t);
                    let derivative = routes
                        .iter()
                        .zip(&es)
                        .zip(signs)
                        .map(|(((sign, shift, _, _), e), sigma)| {
                            sigma
                                * (t * point.iter().map(|x| x * x).sum::<f64>()
                                    + sign * shift * point[2])
                                / e
                        })
                        .sum::<f64>();
                    let n = es
                        .iter()
                        .map(|e| 1.0 / (e / temperature).exp_m1())
                        .collect_vec();
                    let (vacuum, thermal) = if signs == [1.0, 1.0] {
                        (1.0, n[0] + n[1])
                    } else if signs[0] > 0.0 {
                        (0.0, n[1] - n[0])
                    } else {
                        (0.0, n[0] - n[1])
                    };
                    let weight = match process {
                        "vacuum" => vacuum,
                        "thermal_excess" => thermal,
                        _ => vacuum + thermal,
                    };
                    if weight == 0.0 {
                        continue;
                    }
                    let h = 2.0 / PI.sqrt() * (2.0 - t * t - 1.0 / (t * t)).exp();
                    let density = 2.0 * PI * t.powi(3) * h * weight
                        / (derivative.abs()
                            * 4.0
                            * es[0]
                            * es[1]
                            * (2.0 * PI).powi(3)
                            * 2.0
                            * energy);
                    let incoming = signs
                        .iter()
                        .position(|sign| *sign < 0.0)
                        .map_or(0, |i| routes[i].3);
                    *expected_channels.entry(incoming).or_insert(0.0) += density;
                }
                let evaluation = sample.evaluation;
                let expected = expected_channels.values().sum::<f64>();
                let value = evaluation.integrand_result;
                assert!(value.re.0.is_finite() && value.im.0.is_finite());
                assert!(
                    (value.re.0 - expected).abs() < 1e-12 + 1e-8 * expected.abs(),
                    "{block}/{process}/{point:?}: {value:?}, phase space={expected:e}"
                );
                assert!(value.im.0.abs() < 1e-12);
                let mut actual_channels = BTreeMap::new();
                for event in evaluation
                    .event_groups
                    .iter()
                    .flat_map(|group| group.iter())
                {
                    let incoming = event
                        .cut_info
                        .particle_pdgs
                        .0
                        .iter()
                        .copied()
                        .find(|pdg| *pdg != 9000001)
                        .unwrap_or(0);
                    *actual_channels.entry(incoming).or_insert(0.0) += event.weight.re.0;
                    let mut residual = [0.0_f64; 4];
                    for (sign, momenta, pdgs) in [
                        (
                            1.0,
                            event.kinematic_configuration.0.as_slice(),
                            &event.cut_info.particle_pdgs.0[..],
                        ),
                        (
                            -1.0,
                            event.kinematic_configuration.1.as_slice(),
                            &event.cut_info.particle_pdgs.1[..],
                        ),
                    ] {
                        for (p, pdg) in momenta.iter().zip(pdgs) {
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
                            if *pdg != 9000001 {
                                let mass = if *pdg == 9000002 { mass1 } else { mass2 };
                                assert!(
                                    (components[0].powi(2)
                                        - components[1..].iter().map(|x| x * x).sum::<f64>()
                                        - mass * mass)
                                        .abs()
                                        < 1e-10
                                );
                            }
                        }
                    }
                    assert!(residual.iter().all(|x| x.abs() < 1e-10));
                }
                assert_eq!(
                    actual_channels.len(),
                    expected_channels.len(),
                    "{block}/{process}: channel completeness"
                );
                for (channel, expected) in expected_channels {
                    let actual = actual_channels[&channel];
                    assert!(
                        (actual - expected).abs() < 1e-12 + 1e-8 * expected.abs(),
                        "{block}/{process}/{channel}: event phase space"
                    );
                }
                assert!(
                    (actual_channels.values().sum::<f64>() - value.re.0).abs()
                        < 1e-12 + 1e-10 * value.re.0.abs()
                );
            }
        }
        // Deterministic bath-frame integral over d^3k. Rotations about P
        // remove the azimuth; r=x/(1-x) covers the full half-line. The two
        // midpoint grids and Richardson extrapolation are independent of
        // the native parameterization, h normalization and thermal residues.
        let (process_id, integrand_name) = cli.state.find_integrand_ref(
            Some(&ProcessRef::Unqualified("equilibrium".into())),
            Some(&"LO".to_owned()),
        )?;
        let mut integrals = Vec::new();
        for count in [32, 64] {
            let angular_count = if momentum == 0.0 { 1 } else { count };
            let mut coordinates = Vec::new();
            let mut weights = Vec::new();
            for i in 0..count {
                let x = (i as f64 + 0.5) / count as f64;
                let radius = x / (1.0 - x);
                for j in 0..angular_count {
                    let cosine = 2.0 * (j as f64 + 0.5) / angular_count as f64 - 1.0;
                    coordinates.extend([
                        radius * (1.0 - cosine * cosine).sqrt(),
                        0.0,
                        radius * cosine,
                    ]);
                    weights.push(
                        4.0 * PI * radius.powi(2)
                            / ((1.0 - x).powi(2) * (count * angular_count) as f64),
                    );
                }
            }
            let points = ndarray::Array2::from_shape_vec((weights.len(), 3), coordinates)?;
            let discrete = ndarray::Array2::<usize>::zeros((points.nrows(), 0));
            let batch = EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(integrand_name.clone()),
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
            let integral = batch
                .samples
                .iter()
                .zip(weights)
                .map(|(sample, weight)| {
                    assert!(sample.evaluation.integrand_result.re.0.is_finite());
                    sample.evaluation.integrand_result.re.0 * weight
                })
                .sum::<f64>();
            integrals.push(integral);
        }
        let extrapolated = (4.0 * integrals[1] - integrals[0]) / 3.0;
        let relative_error = (target != 0.0).then(|| (extrapolated - target).abs() / target.abs());
        gammalooprs::debug_tags!(#integration, #thermal, #benchmark, #summary;
            stage = "thermal_rate_reference_comparison",
            fixture = "scalar_bubble", case = block, observable = "spectral_rate",
            energy, momentum, mass1, mass2, temperature,
            coarse = integrals[0], fine = integrals[1], extrapolated, target,
            absolute_error = (extrapolated - target).abs(),
            relative_error = ?relative_error,
            "Native scalar rate compared with independent two-body phase space"
        );
        assert!(
            (extrapolated - target).abs() < 1e-12 + 1e-4 * target.abs(),
            "{block}: native quadrature={integrals:?}, extrapolated={extrapolated:e}, literature={target:e}"
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn thermal_raised_scalar_bubble_matches_mass_derivative_on_both_root_orientations() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let state_path = get_tests_workspace_path().join(format!(
        "thermal_raised_scalar_bubble_{}",
        std::process::id()
    ));
    let mut cli = get_example_cli(
        "thermal_rates/scalar_bubble/run.toml",
        &["generate", "rest_decay"],
        Some(state_path),
        None,
        true,
    )?;
    // The serial insertion gives (i/D)(-i)(i/D)=i/D^2, exactly
    // d(i/D)/d(m_1^2). Model vertices retain their -ig factors.
    const DOT: &str = r#"digraph raised_thermal_bubble {
        a[int_id="V_XBB"]; b[int_id="V_XBB"]; insertion[num="-1𝑖"];
        ext0[style=invis,is_cut=0];
        ext0->a[id=0,particle=X]; b->ext0[id=4,particle=X];
        a->insertion[id=1,particle=B1,lmb_id=0];
        insertion->b[id=2,particle=B1]; a->b[id=3,particle=B2];
    }"#;
    cli.run_command("set global kv global.generation.medium.mode=thermodynamic_equilibrium global.generation.medium.vacuum_subtraction=false")?;
    cli.run_command(&format!(
        "import graphs --inline-dot \"\"\"{DOT}\"\"\" --process-spec 'X > B1 B2' -p raised -i LO -o"
    ))?;
    cli.run_command("generate existing -p raised -i LO")?;
    let (raised_id, raised_integrand) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified("raised".into())),
        Some(&"LO".to_owned()),
    )?;
    for (block, energy, mass1, mass2, target, derivative_target) in [
        (
            "rest_decay",
            2.0_f64,
            0.4_f64,
            0.7_f64,
            0.012104501639179106_f64,
            -0.005206172268741891_f64,
        ),
        (
            "absorption",
            1.0,
            0.3,
            1.8,
            0.00688050878313337,
            0.00036810356346982456,
        ),
    ] {
        cli.run_command(&format!("run {block}"))?;
        cli.run_command(&format!("set process -p raised -i LO kv kinematics.e_cm={energy} kinematics.externals.data.momenta=[[{energy},0.0,0.0,0.0]] general.inverse_temperature=1.6666666666666667"))?;
        let discriminant = ((energy.powi(2) - (mass1 + mass2).powi(2))
            * (energy.powi(2) - (mass1 - mass2).powi(2)))
        .sqrt();
        let physical_radius = discriminant / (2.0 * energy);
        let radius_derivative =
            (mass1 * mass1 - energy * energy - mass2 * mass2) / (2.0 * energy * discriminant);
        for radius in [0.4_f64, 1.1, 2.0] {
            let t = physical_radius / radius;
            let h = 2.0 / PI.sqrt() * (2.0 - t * t - 1.0 / (t * t)).exp();
            let dh = h * (-2.0 * t + 2.0 / t.powi(3));
            let expected = ((derivative_target * physical_radius + target * radius_derivative) * h
                + target * physical_radius * dh * radius_derivative / radius)
                / (4.0 * PI * radius.powi(4));
            let (_, actual) = Inspect {
                process: Some(ProcessRef::Unqualified("raised".into())),
                integrand_name: Some("LO".into()),
                point: vec![radius, 0.0, 0.0],
                momentum_space: true,
                ..Default::default()
            }
            .run(&mut cli.state)?;
            assert!(
                actual.im.abs() < 1e-12
                    && (actual.re - expected).abs() < 1e-12 + 1e-8 * expected.abs(),
                "{block}/r={radius}: raised={actual:?}, analytic derivative={expected:e}"
            );
            let mut shifted = Vec::new();
            let step = 1e-5;
            for sign in [-1.0, 1.0] {
                let mass = (mass1 * mass1 + sign * step).sqrt();
                cli.run_command(&format!("set model M1={mass}"))?;
                let (_, value) = Inspect {
                    process: Some(ProcessRef::Unqualified("equilibrium".into())),
                    integrand_name: Some("LO".into()),
                    point: vec![radius, 0.0, 0.0],
                    momentum_space: true,
                    ..Default::default()
                }
                .run(&mut cli.state)?;
                shifted.push(value.re);
            }
            cli.run_command(&format!("set model M1={mass1}"))?;
            let finite_difference = (shifted[1] - shifted[0]) / (2.0 * step);
            assert!(
                (actual.re - finite_difference).abs() < 1e-10 + 2e-6 * expected.abs(),
                "{block}/r={radius}: raised={actual:?}, native mass difference={finite_difference:e}"
            );
        }
        let mut integrals = Vec::new();
        for count in [64, 128] {
            let mut coordinates = Vec::new();
            let mut weights = Vec::new();
            for i in 0..count {
                let x = (i as f64 + 0.5) / count as f64;
                let radius = x / (1.0 - x);
                coordinates.extend([radius, 0.0, 0.0]);
                weights.push(4.0 * PI * radius.powi(2) / (count as f64 * (1.0 - x).powi(2)));
            }
            let points = ndarray::Array2::from_shape_vec((count, 3), coordinates)?;
            let discrete = ndarray::Array2::<usize>::zeros((count, 0));
            let batch = EvaluateSamples {
                process_id: Some(raised_id),
                integrand_name: Some(raised_integrand.clone()),
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
            integrals.push(
                batch
                    .samples
                    .iter()
                    .zip(weights)
                    .map(|(sample, weight)| sample.evaluation.integrand_result.re.0 * weight)
                    .sum::<f64>(),
            );
        }
        let extrapolated = (4.0 * integrals[1] - integrals[0]) / 3.0;
        gammalooprs::debug_tags!(#integration, #thermal, #benchmark, #summary;
            stage = "thermal_rate_reference_comparison",
            fixture = "scalar_bubble", case = block, observable = "mass1_squared_derivative",
            energy, mass1, mass2,
            coarse = integrals[0], fine = integrals[1], extrapolated,
            target = derivative_target,
            absolute_error = (extrapolated - derivative_target).abs(),
            relative_error = (extrapolated - derivative_target).abs() / derivative_target.abs(),
            "Native raised scalar rate compared with analytic mass-squared derivative"
        );
        assert!(
            (extrapolated - derivative_target).abs() < 1e-8 * derivative_target.abs(),
            "{block}: raised quadrature={integrals:?}, analytic mass derivative={derivative_target:e}"
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
