use super::*;
use std::{
    collections::BTreeMap,
    f64::consts::{FRAC_1_SQRT_2, PI},
    fs,
};

#[test]
fn cold_yukawa_rate_and_mass_insertions_match_pauli_blocked_phase_space() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let root = get_tests_workspace_path().join(format!("cold_yukawa_rate_{}", std::process::id()));
    clean_test(&root);
    let mut cli = get_test_cli(None, root.join("state"), None, true)?;

    // No bundled model has a two-point Dirac mass vertex. Reuse Hbb's
    // native spin and color tensors, removing only its scalar slot. A
    // scalar DOT num=-i would replace the vertex and lose both identities.
    let mut model =
        serde_json::to_value(gammalooprs::utils::load_generic_model("sm").to_serializable())?;
    let mut lorentz = model["lorentz_structures"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["name"] == "FFS4")
        .unwrap()
        .clone();
    lorentz["name"] = serde_json::json!("FF_mass_b");
    lorentz["spins"] = serde_json::json!([2, 2]);
    model["lorentz_structures"]
        .as_array_mut()
        .unwrap()
        .push(lorentz);
    let mut vertex = model["vertex_rules"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["name"] == "V_78")
        .unwrap()
        .clone();
    vertex["name"] = serde_json::json!("V_mass_b");
    vertex["particles"] = serde_json::json!(["b~", "b"]);
    vertex["lorentz_structures"] = serde_json::json!(["FF_mass_b"]);
    vertex["couplings"] = serde_json::json!([["GC_mass_b"]]);
    model["vertex_rules"].as_array_mut().unwrap().push(vertex);
    model["couplings"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({
            "name": "GC_mass_b", "expression": "-1𝑖", "orders": [], "value": [0.0, -1.0]
        }));
    for parameter in model["parameters"].as_array_mut().unwrap() {
        let value = match parameter["name"].as_str().unwrap() {
            "MB" => Some(1.5),
            "MH" => Some(5.0),
            "muB" => Some(7.5),
            "muQ" => Some(0.0),
            "Gf" | "ymb" => Some(FRAC_1_SQRT_2),
            _ => None,
        };
        if let Some(value) = value {
            parameter["value"] = serde_json::json!([value, 0.0]);
        }
    }
    let model_path = root.join("sm_with_mass_vertex.json");
    fs::write(&model_path, serde_json::to_vec_pretty(&model)?)?;
    cli.run_command(&format!("import model {}", model_path.display()))?;
    // Nonzero muB at import prevents simplification from freezing it to zero.
    // MB and ymb are independent: varying the physical mass keeps yb fixed.
    let yb = cli.state.model.get_parameter("yb").value.unwrap().re.0;
    // Resolving the electroweak parameters uses a difference of nearly
    // equal masses here; normalize to the resolved coupling below.
    assert!(yb.is_finite() && (yb - 1.0).abs() < 1e-8);
    run_commands(
        &mut cli,
        &[
            "set global kv global.n_cores.feyngen=1 global.n_cores.generate=1 global.n_cores.compile=1 global.n_cores.integrate=1",
            "set global kv global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.uv.softct=false global.generation.threshold_subtraction.enable_thresholds=false global.generation.tropical_subgraph_table.disable_tropical_generation=true",
            "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true global.generation.evaluator.iterative_orientation_optimization=false",
            r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
integral_unit = "none"
generate_events = true
[kinematics]
e_cm = 6.25
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[6.25, 0.0, 0.0, 3.75]]
helicities = [0]
[subtraction]
disable_threshold_subtraction = true
[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
[h_function]
function = "poly_exponential"
sigma = 1.0
'"#,
        ],
    )?;

    // Every numerator comes from the model. The minus sign is the closed
    // fermion-loop statistics factor, present also in the generated Hbb
    // forward graph. Serial insertion retains S(-i 1)S = dS/dm, including
    // the derivative of the Dirac numerator, with no extra Yukawa factor.
    for (process, mode, subtract_vacuum, fermion_edges) in [
        (
            "cold",
            "zero_temperature_equilibrium",
            false,
            "a->b[id=1,particle=b,lmb_id=0]; b->a[id=2,particle=b]; b->ext0[id=3,particle=H];",
        ),
        (
            "excess",
            "zero_temperature_equilibrium",
            true,
            "a->b[id=1,particle=b,lmb_id=0]; b->a[id=2,particle=b]; b->ext0[id=3,particle=H];",
        ),
        (
            "vacuum",
            "vacuum",
            false,
            "a->b[id=1,particle=b,lmb_id=0]; b->a[id=2,particle=b]; b->ext0[id=3,particle=H];",
        ),
        (
            "mass_b",
            "zero_temperature_equilibrium",
            false,
            "insertion[int_id=V_mass_b]; a->insertion[id=1,particle=b,lmb_id=0]; insertion->b[id=2,particle=b]; b->a[id=3,particle=b]; b->ext0[id=4,particle=H];",
        ),
        (
            "mass_antib",
            "zero_temperature_equilibrium",
            false,
            "insertion[int_id=V_mass_b]; a->b[id=1,particle=b,lmb_id=0]; b->insertion[id=2,particle=b]; insertion->a[id=3,particle=b]; b->ext0[id=4,particle=H];",
        ),
    ] {
        cli.run_command(&format!("set global kv global.generation.medium.mode={mode} global.generation.medium.vacuum_subtraction={subtract_vacuum}"))?;
        let dot = format!(
            r#"digraph {process} {{
            overall_factor=-1;
            a[int_id=V_78]; b[int_id=V_78]; ext0[style=invis,is_cut=0];
            ext0->a[id=0,particle=H]; {fermion_edges}
        }}"#
        );
        cli.run_command(&format!("import graphs --inline-dot \"\"\"{dot}\"\"\" --process-spec 'H > b b~' -p {process} -i LO -o"))?;
        cli.run_command(&format!("generate existing -p {process} -i LO"))?;
    }
    cli.run_command("save state -o")?;
    cli.state =
        gammaloop_api::state::State::load(cli.cli_settings.state.folder.clone(), None, None)?;
    cli.state.activate_loaded_integrand_backends(false)?;

    // Ho--Scherrer arXiv:1503.03534, (3.15--16), and the finite-mu
    // fermion bubble of Lundberg--Pasechnik arXiv:2007.01224, (VIII.3,9).
    // At Nc=3 and Hbb vertex -iyb/sqrt(2), Gamma_vac=3 yb^2 M^2 v^3/(16pi E).
    // E=25/4,p=15/4,M=5,m=3/2 gives v=4/5 and e in [13/8,37/8].
    // For |mu|=5/2 the unblocked fraction is 17/24. Differentiating the
    // full propagator on each line gives unequal partial derivatives;
    // their sum is dGamma/dm at fixed yb,P,mu. Reversing mu exchanges them.
    let processes = ["cold", "excess", "vacuum", "mass_b", "mass_antib"];
    let points = arr2(&[[1.0, 0.0, 0.0], [0.4, 0.0, 1.2], [1.0, 0.0, -0.6]]);
    for mu_sign in [1.0_f64, -1.0] {
        cli.run_command(&format!("set model muB={}", mu_sign * 7.5))?;
        assert!((cli.state.model.get_parameter("yb").value.unwrap().re.0 - yb).abs() < 1e-14);
        let mut point_values = BTreeMap::new();
        let mut fermi_event_count = 0;
        for process in processes {
            let (process_id, integrand_name) = cli.state.find_integrand_ref(
                Some(&ProcessRef::Unqualified(process.into())),
                Some(&"LO".into()),
            )?;
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
            let mut values = Vec::new();
            for sample in batch.samples {
                let evaluation = sample.evaluation;
                let value = evaluation.integrand_result;
                assert!(value.re.0.is_finite() && value.im.0.abs() < 1e-10);
                values.push(value.re.0);
                let mut event_weight = 0.0;
                for event in evaluation
                    .event_groups
                    .iter()
                    .flat_map(|group| group.iter())
                {
                    event_weight += event.weight.re.0;
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
                            for (total, component) in residual.iter_mut().zip(components) {
                                *total += sign * component;
                            }
                            if pdg.abs() == 5 {
                                assert!(
                                    (components[0].powi(2)
                                        - components[1..].iter().map(|x| x * x).sum::<f64>()
                                        - 2.25)
                                        .abs()
                                        < 1e-9
                                );
                                if process.starts_with("mass_")
                                    && *pdg == (mu_sign * 5.0) as isize
                                    && (components[0] - 2.5).abs() < 1e-9
                                {
                                    // The occupied daughter lies simultaneously on the
                                    // Fermi sphere k_F=2 and the physical two-body cut.
                                    assert!((components[3] / 2.0 - 5.0 / 12.0).abs() < 1e-9);
                                    fermi_event_count += 1;
                                }
                            }
                        }
                    }
                    assert!(residual.iter().all(|x| x.abs() < 1e-9));
                }
                assert!(
                    (event_weight - value.re.0).abs() < 1e-11 + 1e-9 * value.re.0.abs(),
                    "mu_sign={mu_sign}/{process}: event closure {event_weight:e} != {value:?}"
                );
            }
            point_values.insert(process, values);
        }
        assert!(
            fermi_event_count > 0,
            "mu_sign={mu_sign}: missing joint Fermi/cut events"
        );
        for ((cold, vacuum), excess) in point_values["cold"]
            .iter()
            .zip(&point_values["vacuum"])
            .zip(&point_values["excess"])
        {
            let expected = cold - vacuum;
            assert!((excess - expected).abs() < 1e-11 + 1e-9 * expected.abs());
        }

        // Integrate in the bath frame, splitting at both exact Pauli-step
        // directions and their reflected routes. This avoids estimating a
        // distributional mass derivative by finite differences of a grid
        // that cannot resolve the moving Fermi boundary.
        let boundary = 35.0 / (9.0 * 21.0_f64.sqrt());
        let angular_bounds = [-1.0, -boundary, -5.0 / 12.0, 5.0 / 12.0, boundary, 1.0];
        // Gauss--Legendre integration on each smooth angular interval cuts
        // the native evaluation budget by four. Orders 8 and 16 reproduce
        // the independent translated Fermi-profile integral within 9e-6
        // after extrapolation, below the unchanged physical-rate tolerance.
        let angular_rules: [&[(f64, f64)]; 2] = [
            &[
                (-0.9602898564975362, 0.10122853629037706),
                (-0.7966664774136267, 0.22238103445337443),
                (-0.525532409916329, 0.3137066458778869),
                (-0.18343464249564978, 0.36268378337836166),
                (0.18343464249564978, 0.36268378337836166),
                (0.525532409916329, 0.3137066458778869),
                (0.7966664774136267, 0.22238103445337443),
                (0.9602898564975362, 0.10122853629037706),
            ],
            &[
                (-0.9894009349916499, 0.027152459411754176),
                (-0.9445750230732326, 0.062253523938647456),
                (-0.8656312023878318, 0.0951585116824926),
                (-0.755404408355003, 0.12462897125553407),
                (-0.6178762444026438, 0.1495959888165767),
                (-0.45801677765722737, 0.16915651939500265),
                (-0.2816035507792589, 0.18260341504492364),
                (-0.09501250983763744, 0.18945061045506864),
                (0.09501250983763744, 0.18945061045506864),
                (0.2816035507792589, 0.18260341504492364),
                (0.45801677765722737, 0.16915651939500265),
                (0.6178762444026438, 0.1495959888165767),
                (0.755404408355003, 0.12462897125553407),
                (0.8656312023878318, 0.0951585116824926),
                (0.9445750230732326, 0.062253523938647456),
                (0.9894009349916499, 0.027152459411754176),
            ],
        ];
        let mut integrals = BTreeMap::new();
        for process in processes {
            let (process_id, integrand_name) = cli.state.find_integrand_ref(
                Some(&ProcessRef::Unqualified(process.into())),
                Some(&"LO".into()),
            )?;
            let mut estimates = Vec::new();
            for (count, angular_rule) in [32, 64].into_iter().zip(angular_rules) {
                // At negative mu the occupied antiparticle has momentum
                // P-k, so its Fermi profile is translated in this input LMB.
                // Its independent normalized-kernel integral needs twice the
                // radial resolution; both axes still refine by two below.
                let radial_count = if mu_sign < 0.0 && process.starts_with("mass_") {
                    2 * count
                } else {
                    count
                };
                let mut coordinates = Vec::new();
                let mut weights = Vec::new();
                for i in 0..radial_count {
                    let x = (i as f64 + 0.5) / radial_count as f64;
                    let radius = 2.0 * x / (1.0 - x);
                    for interval in angular_bounds.windows(2) {
                        let width = interval[1] - interval[0];
                        for &(node, angular_weight) in angular_rule {
                            let cosine = interval[0] + width * (node + 1.0) / 2.0;
                            coordinates.extend([
                                radius * (1.0 - cosine * cosine).sqrt(),
                                0.0,
                                radius * cosine,
                            ]);
                            weights.push(
                                2.0 * PI * width * angular_weight * radius.powi(2)
                                    / (1.0 - x).powi(2)
                                    / radial_count as f64,
                            );
                        }
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
                        let value = sample.evaluation.integrand_result;
                        assert!(value.re.0.is_finite() && value.im.0.abs() < 1e-9);
                        value.re.0 * weight
                    })
                    .sum::<f64>();
                estimates.push(integral);
            }
            let extrapolated = (4.0 * estimates[1] - estimates[0]) / 3.0;
            let target = yb * yb / PI
                * match process {
                    "cold" => 34.0 / 125.0,
                    "excess" => -14.0 / 125.0,
                    "vacuum" => 48.0 / 125.0,
                    "mass_b" if mu_sign > 0.0 => -9.0 / 100.0,
                    "mass_antib" if mu_sign < 0.0 => -9.0 / 100.0,
                    _ => -93.0 / 500.0,
                };
            gammalooprs::debug_tags!(#integration, #thermal, #benchmark, #summary;
                stage = "thermal_rate_reference_comparison",
                fixture = "cold_yukawa", process, chemical_potential = mu_sign * 2.5, yb,
                coarse = estimates[0], fine = estimates[1], extrapolated, target,
                absolute_error = (extrapolated - target).abs(),
                relative_error = (extrapolated - target).abs() / target.abs(),
                "Native cold Yukawa contribution compared with independent Pauli-blocked phase space"
            );
            assert!(
                (extrapolated - target).abs() < 1e-4 * target.abs(),
                "mu_sign={mu_sign}/{process}: native quadrature={estimates:?}, extrapolated={extrapolated:e}, literature={target:e}"
            );
            integrals.insert(process, extrapolated);
        }
        let derivative = integrals["mass_b"] + integrals["mass_antib"];
        let target = -69.0 * yb * yb / (250.0 * PI);
        gammalooprs::debug_tags!(#integration, #thermal, #benchmark, #summary;
            stage = "thermal_rate_reference_comparison",
            fixture = "cold_yukawa", observable = "common_mass_derivative",
            chemical_potential = mu_sign * 2.5, yb, derivative, target,
            absolute_error = (derivative - target).abs(),
            relative_error = (derivative - target).abs() / target.abs(),
            "Sum of native fermion mass insertions compared with common-mass derivative"
        );
        assert!((derivative - target).abs() < 1e-4 * target.abs());
        assert!((integrals["cold"] - integrals["vacuum"] - integrals["excess"]).abs() < 1e-10);
    }
    clean_test(&root);
    Ok(())
}
