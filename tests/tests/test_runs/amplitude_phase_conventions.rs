use super::utils::*;
use super::*;
use gammaloop_api::commands::integrate::RendererOption;
use gammalooprs::{
    feyngen::diagram_generator::evaluate_overall_factor, graph::FeynmanGraph,
    processes::ProcessCollection, settings::global::OrientationPattern,
};
use std::f64::consts::PI;
use symbolica::atom::Atom;

#[test]
#[serial]
fn scalar_tree_amplitudes_match_textbook_ufo_rules() -> Result<()> {
    // With lambda=2, the contact graph is -2i. Distinct external species
    // and the allowed cubic vertices select only s-channel scalar_2 exchange:
    // (-i lambda)^2 i/(s-m_2^2) = +4i/3 at s=1 and m_2=2.
    for (name, final_states, vertices, expected_imaginary) in [
        ("contact", "scalar_0 scalar_0", "V_4_SCALAR_0000", -2.0),
        (
            "exchange",
            "scalar_1 scalar_1",
            "V_3_SCALAR_002 V_3_SCALAR_112",
            4.0 / 3.0,
        ),
    ] {
        let test_name = format!("scalar_tree_amplitude_{name}");
        let root = get_tests_workspace_path().join(&test_name);
        clean_test(&root);
        let mut cli = get_test_cli(None, root, Some(test_name), true)?;
        run_commands(
            &mut cli,
            &[
                "import model scalars-default.json",
                "set model mass_scalar_1=0 mass_scalar_2=2 lam=2",
                "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=false",
                r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
integral_unit = "none"
disable_flux_factor = true
[subtraction]
disable_threshold_subtraction = true
[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
[kinematics]
e_cm = 1.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[0.5,0.0,0.0,0.5],[0.5,0.0,0.0,-0.5],[0.5,0.5,0.0,0.0],"dependent"]
helicities = [0,0,0,0]
'"#,
                &format!(
                    "generate amp scalar_0 scalar_0 > {final_states} [{{0}}] --allowed-vertex-interactions {vertices} -p tree -i scalar",
                ),
            ],
        )?;
        let ProcessCollection::Amplitudes(amplitudes) =
            &cli.state.process_list.processes[0].collection
        else {
            panic!("expected a scalar amplitude")
        };
        let [amplitude_graph] = amplitudes["scalar"].graphs.as_slice() else {
            panic!("the tree oracle must select exactly one graph")
        };
        assert_eq!(amplitude_graph.graph.get_loop_number(), 0);
        let (_, actual) = Inspect {
            process: Some(ProcessRef::Unqualified("tree".to_string())),
            integrand_name: Some("scalar".to_string()),
            point: vec![],
            momentum_space: true,
            graph_id: Some(0),
            ..Default::default()
        }
        .run(&mut cli.state)?;
        assert!(
            actual.re.abs() <= 1.0e-12 && (actual.im - expected_imaginary).abs() <= 1.0e-12,
            "{name}: the full UFO tree amplitude is {actual:e}, expected {expected_imaginary:e}i",
        );
        clean_test(&cli.cli_settings.state.folder);
    }
    Ok(())
}

#[test]
#[serial]
fn scalar_unit_amplitudes_match_convergent_mc_integrals_at_one_through_four_loops() -> Result<()> {
    // Each triangle carries one independent momentum and three identical massive
    // scalar denominators. The petals share a vertex, so the graph is connected,
    // but the full integral factorizes. Every energy and UV integral converges:
    // integral d^4k/(2π)^4 (k²-m²+i0)^(-3) = -i/(32π²m²).
    // At zero external momentum, each of the six acyclic triangle orientations
    // contributes equally. Selecting (+,+,-) gives -i/(192π²m²) per petal.
    // This is a literal scalar topology with unit numerator, not a UFO amplitude.
    for (index, phase) in [(0.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 0.0)]
        .into_iter()
        .enumerate()
    {
        let loops = index + 1;
        let name = format!("scalar_unit_amplitude_phase_{loops}");
        let root = get_tests_workspace_path().join(&name);
        clean_test(&root);
        let mut cli = get_test_cli(None, root, Some(name), true)?;
        let orientation = (0..loops)
            .flat_map(|_| ["+", "+", "-"])
            .chain(["0", "0"])
            .join(",");
        cli.cli_settings.global.generation.orientation_pattern =
            OrientationPattern::from_user_pattern(&format!("({orientation})"))?;
        let mut dot = format!(
            "digraph vacuum_triangles {{ edge [mass=1 num=1]; node [num=1]; \
             incoming [style=invis]; outgoing [style=invis]; \
             incoming -> a [id={} particle=scalar_0 mass=0]; \
             a -> outgoing [id={} particle=scalar_0 mass=0];",
            3 * loops,
            3 * loops + 1,
        );
        for petal in 0..loops {
            write!(
                dot,
                "a -> b{petal} [id={} lmb_id={petal}]; \
                 b{petal} -> c{petal} [id={}]; c{petal} -> a [id={}];",
                3 * petal,
                3 * petal + 1,
                3 * petal + 2,
            )?;
        }
        dot.push('}');
        // This analytically known axis only trains the sampling grid; the
        // integrator accumulates and the assertions below check both components.
        let training_phase = if loops.is_multiple_of(2) {
            "real"
        } else {
            "imag"
        };
        run_commands(
            &mut cli,
            &[
                "import model scalars-default.json",
                "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=false",
                &format!(
                    r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
integral_unit = "none"
disable_flux_factor = true
[subtraction]
disable_threshold_subtraction = true
[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0
[kinematics]
e_cm = 1.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[0.0,0.0,0.0,0.0],"dependent"]
helicities = [0,0]
[integrator]
integrated_phase = "{training_phase}"
n_start = 20000
n_increase = 10000
n_max = 200000
min_samples_for_update = 1000
seed = 90210
target_relative_accuracy = 0.01
'"#,
                ),
                &format!("import graphs --inline-dot \"\"\"{dot}\"\"\" -p vacuum -i scalar -o"),
                "generate existing -p vacuum -i scalar",
            ],
        )?;
        let ProcessCollection::Amplitudes(amplitudes) =
            &cli.state.process_list.processes[0].collection
        else {
            panic!("expected a scalar amplitude")
        };
        let [amplitude_graph] = amplitudes["scalar"].graphs.as_slice() else {
            panic!("the scalar integral oracle must select exactly one graph")
        };
        assert_eq!(amplitude_graph.graph.get_loop_number(), loops);
        assert_eq!(
            cli.state.process_list.processes[0]
                .get_integrand("scalar")?
                .require_generated()?
                .graph_orientation_count(0),
            Some(1),
        );

        let point: Vec<f64> = (0..loops)
            .flat_map(|petal| [0.2 + 0.1 * petal as f64, -0.3, 0.1])
            .collect();
        // The full energy residue is -3i/(16E^5), so the selected orientation
        // contributes -i/(32E^5). Its remaining spatial normalization is
        // (2π)^(-3), with no momentum-space Jacobian here.
        let density = point
            .chunks_exact(3)
            .map(|momentum| {
                let e_squared = 1.0 + momentum.iter().map(|k| k * k).sum::<f64>();
                1.0 / (256.0 * PI.powi(3) * e_squared.powf(2.5))
            })
            .product::<f64>();
        let (_, point_value) = Inspect {
            process: Some(ProcessRef::Unqualified("vacuum".to_string())),
            integrand_name: Some("scalar".to_string()),
            point,
            momentum_space: true,
            graph_id: Some(0),
            ..Default::default()
        }
        .run(&mut cli.state)?;
        for (actual, expected) in [
            (point_value.re, phase.0 * density),
            (point_value.im, phase.1 * density),
        ] {
            assert!(
                (actual - expected).abs() <= 1.0e-10 * density,
                "L={loops}: scalar energy residue {point_value:e} differs from {phase:?} * {density:e}",
            );
        }

        // The pointwise check certifies the selected orientation. MC integrates
        // that same contribution, whose independent target includes 1/6 per
        // petal relative to the full integral.
        let output = Integrate {
            process: vec![ProcessRef::Unqualified("vacuum".to_string())],
            integrand_name: vec!["scalar".to_string()],
            workspace_path: Some(cli.cli_settings.state.folder.join("integration_workspace")),
            n_cores: Some(1),
            restart: true,
            renderer: RendererOption::Tabled,
            show_max_weight_info: false,
            no_stream_iterations: true,
            no_stream_updates: true,
            ..Default::default()
        }
        .run(&mut cli.state, &cli.cli_settings)?;
        let integral = single_slot_integral(&output);
        let magnitude = (1.0 / (192.0 * PI * PI)).powi(loops as i32);
        for (value, error, sign) in [
            (integral.result.re.0, integral.error.re.0, phase.0),
            (integral.result.im.0, integral.error.im.0, phase.1),
        ] {
            assert!(value.is_finite() && error.is_finite());
            if sign == 0.0 {
                assert!(value.abs() <= 1.0e-12 * magnitude && error.abs() <= 1.0e-12 * magnitude);
            } else {
                assert!(
                    sign * value > 0.0,
                    "L={loops}: incorrect scalar phase: {integral:?}"
                );
                assert!(
                    error / magnitude < 0.025,
                    "L={loops}: insufficient MC precision: {integral:?}"
                );
                assert!(
                    (value - sign * magnitude).abs() <= (5.0 * error).max(1.0e-10 * magnitude),
                    "L={loops}: scalar MC value {value:e} ± {error:e}, oracle {}",
                    sign * magnitude,
                );
            }
        }
        clean_test(&cli.cli_settings.state.folder);
    }
    Ok(())
}

#[test]
#[serial]
fn scalar_bubble_absorptive_part_matches_the_causal_discontinuity() -> Result<()> {
    // For equal internal masses m=1 and p²=9, the Feynman-parameter logarithm
    // has imaginary part +π sqrt(1-4/9). The physical loop factor i/(16π²)
    // makes the corresponding part of A=-iM negative real. Two -i vertices,
    // two propagator i factors and symmetry 1/2 give Re A=-sqrt(5/9)/(32π).
    // Real local and integrated UV counterterms only affect Im A here.
    let expected = -(5.0_f64 / 9.0).sqrt() / (32.0 * PI);
    let name = "scalar_bubble_causal_discontinuity";
    let root = get_tests_workspace_path().join(name);
    clean_test(&root);
    let mut cli = get_test_cli(None, root, Some(name.to_string()), true)?;
    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "set model mass_scalar_1=1 mass_scalar_2=3 lam=1",
            "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.uv.subtract_uv=true global.generation.uv.generate_integrated=true global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.disable_integrated_ct=false",
            r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
integral_unit = "none"
disable_flux_factor = true
m_uv = 5.0
mu_r = 3.0
[subtraction]
disable_threshold_subtraction = false
[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0
[kinematics]
e_cm = 3.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[3.0,0.0,0.0,0.0],"dependent"]
helicities = [0,0]
[integrator]
integrated_phase = "real"
n_start = 10000
n_increase = 10000
n_max = 50000
min_samples_for_update = 1000
seed = 7331
target_relative_accuracy = 0.01
'"#,
            "generate amp scalar_2 > scalar_2 [{1}] --allowed-vertex-interactions V_3_SCALAR_112 -p bubble -i causal",
        ],
    )?;
    let ProcessCollection::Amplitudes(amplitudes) = &cli.state.process_list.processes[0].collection
    else {
        panic!("expected a scalar amplitude")
    };
    let [bubble] = amplitudes["causal"].graphs.as_slice() else {
        panic!("the discontinuity oracle requires exactly one bubble")
    };
    assert_eq!(bubble.graph.get_loop_number(), 1);
    assert_eq!(
        evaluate_overall_factor(bubble.graph.overall_factor.as_view()),
        Atom::num((1, 2))
    );
    let output = Integrate {
        process: vec![ProcessRef::Unqualified("bubble".to_string())],
        integrand_name: vec!["causal".to_string()],
        workspace_path: Some(cli.cli_settings.state.folder.join("integration_workspace")),
        n_cores: Some(1),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let integral = single_slot_integral(&output);
    let (value, error) = (integral.result.re.0, integral.error.re.0);
    assert!(value < 0.0 && value.is_finite() && error.is_finite());
    assert!(error < 0.025 * expected.abs(), "{integral:?}");
    assert!(
        (value - expected).abs() <= (5.0 * error).max(1.0e-10 * expected.abs()),
        "causal bubble: {value:e} ± {error:e}, discontinuity oracle {expected:e}",
    );
    assert!(integral.result.im.0.is_finite() && integral.error.im.0.is_finite());
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
