use super::utils::*;
use super::*;
use gammaloop_api::commands::integrate::RendererOption;
use gammalooprs::{
    graph::FeynmanGraph,
    processes::{CutId, ProcessCollection},
};
use std::f64::consts::PI;

fn scalar_contact_cli(
    multiplicity: usize,
    incoming: usize,
    final_mass: f64,
) -> Result<gammaloop_integration_tests::CLIState> {
    let name = format!("scalar_born_phase_n{multiplicity}_{incoming}_{final_mass}");
    let root = get_tests_workspace_path().join(&name);
    clean_test(&root);
    let mut cli = get_test_cli(None, root, Some(name), true)?;
    let (initial, vertex_suffix, e_cm, external_momenta, helicities) = if incoming == 1 {
        (
            "scalar_2",
            format!("{}2", "1".repeat(multiplicity)),
            1.0,
            "[[1.0,0.0,0.0,0.0]]",
            "[0]",
        )
    } else {
        assert_eq!(incoming, 2);
        (
            "scalar_0 scalar_0",
            format!("00{}", "1".repeat(multiplicity)),
            2.0,
            "[[1.0,0.0,0.0,1.0],[1.0,0.0,0.0,-1.0]]",
            "[0,0]",
        )
    };
    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            &format!("set model mass_scalar_1={final_mass} mass_scalar_2=1.0 lam=1.0"),
            "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=false",
            &format!(
                r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
integral_unit = "none"
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
[h_function]
function = "poly_exponential"
sigma = 1.0
[kinematics]
e_cm = {e_cm}
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = {external_momenta}
helicities = {helicities}
[integrator]
integrated_phase = "real"
n_start = 10000
n_increase = 10000
n_max = 250000
min_samples_for_update = 1000
seed = 7331
target_relative_accuracy = 0.025
'"#
            ),
            &format!(
                "generate xs {initial} > {} [{{{{{}}}}}] --allowed-vertex-interactions V_{}_SCALAR_{vertex_suffix} -p born -i LO",
                vec!["scalar_1"; multiplicity].join(" "),
                multiplicity - 1,
                multiplicity + incoming,
            ),
        ],
    )?;
    let ProcessCollection::CrossSections(cross_sections) =
        &cli.state.process_list.processes[0].collection
    else {
        panic!("expected a scalar cross section")
    };
    let [supergraph] = cross_sections["LO"].supergraphs.as_slice() else {
        panic!("contact acceptance must select exactly one graph")
    };
    assert_eq!(supergraph.cuts.len(), 1, "expected one physical Born cut");
    assert_eq!(supergraph.graph.get_loop_number(), multiplicity - 1);
    assert_eq!(
        supergraph
            .graph
            .iter_edges_of(&supergraph.cuts[CutId(0)].cut)
            .count(),
        multiplicity,
    );
    Ok(cli)
}

fn check_scalar_contact_point(
    cli: &mut gammaloop_integration_tests::CLIState,
    multiplicity: usize,
    incoming: usize,
    final_mass: f64,
    seed: usize,
) -> Result<()> {
    let point: Vec<f64> = (0..multiplicity - 1)
        .flat_map(|i| {
            [
                0.11 + 0.013 * (i + seed) as f64,
                -0.07 + 0.021 * i as f64,
                0.09 - 0.017 * (i + seed) as f64,
            ]
        })
        .collect();
    let ProcessCollection::CrossSections(cross_sections) =
        &cli.state.process_list.processes[0].collection
    else {
        unreachable!()
    };
    let supergraph = &cross_sections["LO"].supergraphs[0];
    let squared_momenta: Vec<f64> = supergraph
        .graph
        .iter_edges_of(&supergraph.cuts[CutId(0)].cut)
        .map(|(_, edge, _)| {
            let signature = &supergraph.graph.loop_momentum_basis.edge_signatures[edge];
            // In the contact COM frame, external shifts on a final line cancel.
            let shift = signature
                .external
                .iter()
                .zip(if incoming == 1 {
                    vec![0.0]
                } else {
                    vec![1.0, -1.0]
                })
                .map(|(sign, pz)| *sign as i8 as f64 * pz)
                .sum::<f64>();
            assert_eq!(shift, 0.0);
            (0..3)
                .map(|axis| {
                    signature
                        .internal
                        .iter()
                        .zip(point.chunks_exact(3))
                        .map(|(sign, momentum)| *sign as i8 as f64 * momentum[axis])
                        .sum::<f64>()
                        .powi(2)
                })
                .sum()
        })
        .collect();
    let e_cm: f64 = if incoming == 1 { 1.0 } else { 2.0 };
    let energy_sum = |t: f64| {
        squared_momenta
            .iter()
            .map(|q2| (t * t * q2 + final_mass * final_mass).sqrt())
            .sum::<f64>()
    };
    let (mut low, mut high) = (0.0, 1.0);
    while energy_sum(high) < e_cm {
        high *= 2.0;
    }
    for _ in 0..80 {
        let mid = (low + high) / 2.0;
        if energy_sum(mid) < e_cm {
            low = mid;
        } else {
            high = mid;
        }
    }
    let t = (low + high) / 2.0;
    let energies: Vec<f64> = squared_momenta
        .iter()
        .map(|q2| (t * t * q2 + final_mass * final_mass).sqrt())
        .collect();
    let eta_derivative: f64 = squared_momenta
        .iter()
        .zip(&energies)
        .map(|(q2, energy)| t * q2 / energy)
        .sum();
    // This is the phase-space measure itself, without CFF residues, a generated
    // numerator, or a fitted phase. The normalized h integrates to one on t>0.
    let h = 2.0 / PI.sqrt() * (2.0 - t * t - 1.0 / (t * t)).exp();
    let symmetry = (1..=multiplicity).product::<usize>() as f64;
    let flux = if incoming == 1 {
        2.0 * e_cm
    } else {
        2.0 * e_cm * e_cm
    };
    let expected = (2.0 * PI) * t.powi(3 * (multiplicity as i32 - 1)) * h
        / (eta_derivative
            * energies.iter().map(|energy| 2.0 * energy).product::<f64>()
            * (2.0 * PI).powi(3 * (multiplicity as i32 - 1))
            * symmetry
            * flux);
    assert!(expected.is_finite() && expected > 0.0);
    let (_, actual) = Inspect {
        process: Some(ProcessRef::Unqualified("born".to_string())),
        integrand_name: Some("LO".to_string()),
        point,
        momentum_space: true,
        graph_id: Some(0),
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(actual.re > 0.0, "Born cut is not positive: {actual:e}");
    assert!(
        actual.im.abs() <= 1.0e-12 * expected,
        "Born cut is not real: {actual:e}, expected {expected:e}",
    );
    assert!(
        (actual.re - expected).abs() <= 1.0e-10 * expected,
        "Born phase-space density mismatch: {actual:e}, expected {expected:e}",
    );
    Ok(())
}

#[test]
#[serial]
fn scalar_born_phases_match_phase_space_at_one_through_five_loops() -> Result<()> {
    for multiplicity in 2..=6 {
        let mut cli = scalar_contact_cli(multiplicity, 1, 0.0)?;
        for seed in 0..3 {
            check_scalar_contact_point(&mut cli, multiplicity, 1, 0.0, seed)?;
        }
        clean_test(&cli.cli_settings.state.folder);
    }
    for (incoming, final_mass) in [(1, 0.15), (2, 0.0)] {
        let mut cli = scalar_contact_cli(2, incoming, final_mass)?;
        check_scalar_contact_point(&mut cli, 2, incoming, final_mass, 1)?;
        clean_test(&cli.cli_settings.state.folder);
    }
    Ok(())
}

fn integrate_scalar_contact(multiplicity: usize, incoming: usize, final_mass: f64) -> Result<()> {
    let e_cm: f64 = if incoming == 1 { 1.0 } else { 2.0 };
    let factorial = |n: usize| (1..=n).product::<usize>() as f64;
    // Lorentz-invariant n-body phase space, including identical final scalars
    // and physical decay/scattering flux. The coupling is fixed to lambda=1.
    let phase_space = if final_mass == 0.0 {
        e_cm.powi(2 * multiplicity as i32 - 4)
            / (2.0
                * (4.0 * PI).powi(2 * multiplicity as i32 - 3)
                * factorial(multiplicity - 1)
                * factorial(multiplicity - 2))
    } else {
        assert_eq!(multiplicity, 2);
        (1.0 - 4.0 * final_mass * final_mass / (e_cm * e_cm)).sqrt() / (8.0 * PI)
    };
    let flux = if incoming == 1 {
        2.0 * e_cm
    } else {
        2.0 * e_cm * e_cm
    };
    let expected = phase_space / (factorial(multiplicity) * flux);
    let mut cli = scalar_contact_cli(multiplicity, incoming, final_mass)?;
    let output = Integrate {
        process: vec![ProcessRef::Unqualified("born".to_string())],
        integrand_name: vec!["LO".to_string()],
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
    assert!(value > 0.0 && error.is_finite() && error / expected < 0.05);
    assert!(
        (value - expected).abs() <= (5.0 * error).max(1.0e-10 * expected),
        "n={multiplicity}: {value:e} ± {error:e}, phase-space oracle={expected:e}",
    );
    assert!(
        integral.result.im.0.abs() <= 1.0e-12 * expected
            && integral.error.im.0.abs() <= 1.0e-12 * expected,
        "Born MC integral contains an imaginary component: {integral:?}",
    );
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

macro_rules! scalar_born_mc_test {
    ($name:ident, $multiplicity:expr) => {
        #[test]
        #[serial]
        fn $name() -> Result<()> {
            integrate_scalar_contact($multiplicity, 1, 0.0)
        }
    };
}

scalar_born_mc_test!(scalar_born_two_body_mc_matches_phase_space, 2);
scalar_born_mc_test!(scalar_born_three_body_mc_matches_phase_space, 3);
scalar_born_mc_test!(scalar_born_four_body_mc_matches_phase_space, 4);
scalar_born_mc_test!(scalar_born_five_body_mc_matches_phase_space, 5);
scalar_born_mc_test!(scalar_born_six_body_mc_matches_phase_space, 6);

#[test]
#[serial]
fn scalar_born_massive_decay_and_scattering_mc_match_phase_space() -> Result<()> {
    integrate_scalar_contact(2, 1, 0.15)?;
    integrate_scalar_contact(2, 2, 0.0)
}

#[test]
#[serial]
fn scalar_born_exchange_tree_matches_squared_amplitude() -> Result<()> {
    // This labelled Born topology has one internal scalar on each side. The
    // explicit identical-state factor is 1/3!, with no sum over the other
    // possible spectator choices or their interference graphs. All vertex and
    // propagator numerators are supplied by the model, including their i's.
    const DOT: &str = r#"digraph born_tree {
        overall_factor="1/6";
        a[int_id="V_3_SCALAR_012"]; b[int_id="V_3_SCALAR_001"];
        c[int_id="V_3_SCALAR_012"]; d[int_id="V_3_SCALAR_001"];
        ext0[style=invis,is_cut=0];
        ext0->a[particle=scalar_2]; c->ext0[particle=scalar_2];
        a->c[id=1,particle=scalar_0,lmb_id=0];
        b->d[id=2,particle=scalar_0,lmb_id=1];
        b->d[id=3,particle=scalar_0];
        a->b[id=4,particle=scalar_1]; d->c[id=5,particle=scalar_1];
    }"#;
    let mut cli = scalar_contact_cli(3, 1, 0.0)?;
    run_commands(
        &mut cli,
        &[
            "remove processes",
            "set model mass_scalar_1=2.0",
            &format!(
                "import graphs --inline-dot \"\"\"{DOT}\"\"\" --process-spec 'scalar_2 > scalar_0 scalar_0 scalar_0' -p born -i LO -o"
            ),
            "generate existing -p born -i LO",
        ],
    )?;
    let ProcessCollection::CrossSections(cross_sections) =
        &cli.state.process_list.processes[0].collection
    else {
        unreachable!()
    };
    assert_eq!(cross_sections["LO"].supergraphs.len(), 1);
    assert_eq!(cross_sections["LO"].supergraphs[0].cuts.len(), 1);
    for point in [
        vec![0.11, 0.07, 0.09, 0.123, 0.049, 0.073],
        vec![0.15, -0.04, 0.08, -0.08, 0.16, 0.12],
    ] {
        let spectator_radius = point[..3].iter().map(|x| x * x).sum::<f64>().sqrt();
        let second_radius = point[3..].iter().map(|x| x * x).sum::<f64>().sqrt();
        let third_radius = (0..3)
            .map(|axis| (point[axis] + point[3 + axis]).powi(2))
            .sum::<f64>()
            .sqrt();
        let radius_sum = spectator_radius + second_radius + third_radius;
        let t = 1.0 / radius_sum;
        let energies = [spectator_radius * t, second_radius * t, third_radius * t];
        let pair_mass_squared = 1.0 - 2.0 * energies[0];
        // A=(-i lambda)^2 i/(s_pair-m_med^2), with lambda=1 and
        // m_med=2 above the whole physical decay phase space.
        let amplitude_squared = 1.0 / (pair_mass_squared - 4.0).powi(2);
        let h = 2.0 / PI.sqrt() * (2.0 - t * t - 1.0 / (t * t)).exp();
        let expected = amplitude_squared * 2.0 * PI * t.powi(6) * h
            / (radius_sum
                * energies
                    .into_iter()
                    .map(|energy| 2.0 * energy)
                    .product::<f64>()
                * (2.0 * PI).powi(6)
                * 6.0
                * 2.0);
        let (_, actual) = Inspect {
            process: Some(ProcessRef::Unqualified("born".to_string())),
            integrand_name: Some("LO".to_string()),
            point,
            momentum_space: true,
            graph_id: Some(0),
            ..Default::default()
        }
        .run(&mut cli.state)?;
        assert!(actual.re > 0.0);
        assert!(actual.im.abs() <= 1.0e-12 * expected);
        assert!(
            (actual.re - expected).abs() <= 1.0e-10 * expected,
            "exchange-tree cut={actual:e}, |A|^2 dPhi oracle={expected:e}",
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
