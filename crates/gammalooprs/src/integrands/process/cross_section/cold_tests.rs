use std::io::Cursor;

use symbolica::state::State;

use super::*;
use crate::{
    feyngen::GenerationType,
    initialisation::test_initialise,
    integrands::process::{MomentumSpaceEvaluationInput, ProcessIntegrand},
    model::InputParamCard,
    processes::{CrossSection, ProcessDefinition},
    utils::load_generic_model,
};

#[test]
fn cold_rate_localizes_virtual_fermi_boundary_with_shared_lu_events() -> Result<()> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    let mut card = InputParamCard::from_file(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../assets/models/json/sm/restrict_thermal.json"),
    )?;
    card.insert("MB".into(), Complex::new_re(F(3.0)));
    card.insert("MH".into(), Complex::new_re(F(1.0)));
    card.insert("muB".into(), Complex::new_re(F(15.0)));
    model.simplify(&mut card)?;
    let mut global: GlobalSettings = toml::from_str(
        r#"
[generation.medium]
mode = "zero_temperature_equilibrium"
[generation.uv]
subtract_uv = false
generate_integrated = false
softct = false
[generation.threshold_subtraction]
enable_thresholds = false
[generation.tropical_subgraph_table]
disable_tropical_generation = true
[generation.evaluator]
compile = false
summed = false
summed_function_map = false
iterative_orientation_optimization = false
"#,
    )?;
    let settings: RuntimeSettings = toml::from_str(
        r#"
[general]
evaluator_method = "SingleParametric"
integral_unit = "none"
generate_events = true
[kinematics]
e_cm = 5.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0]]
helicities = [0]
[subtraction]
disable_threshold_subtraction = true
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = false
[h_function]
function = "exponential"
sigma = 1.0
"#,
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    let mut runtimes = Vec::new();
    for spectator in [false, true] {
        let fermion_cycle = if spectator {
            "A -> C [id=3 lmb_id=1 particle=\"b\"]; C -> A [id=4 particle=\"b\"];"
        } else {
            ""
        };
        let final_edge_id = if spectator { 5 } else { 3 };
        let graphs = Graph::from_string(
            format!(
                r#"digraph cold_spectator_decay {{
                    num=1; node [num=1]; edge [num=1 particle="H"];
                    ext0 [style=invis is_cut=0];
                    ext0 -> A [id=0];
                    A -> B [id=1 lmb_id=0]; A -> B [id=2];
                    {fermion_cycle}
                    B -> ext0 [id={final_edge_id}];
                }}"#
            ),
            &model,
        )?;
        global.generation.medium.vacuum_subtraction = spectator;
        let mut definition =
            ProcessDefinition::from_graph_list(&graphs, GenerationType::CrossSection, &model)?;
        // Imported topology also permits cutting the spectator cycle. Select
        // the scalar channel whose virtual finite-density correction we test.
        definition.final_pdgs_lists = vec![vec![25, 25]];
        let mut process = CrossSection::from_graph_list("cold_rate".into(), graphs, &model)?;
        process.preprocess(
            &model,
            &definition,
            &global.generation,
            (&settings).into(),
            &pool,
        )?;
        process.build_integrand(&model, &definition, &global, (&settings).into(), &pool)?;
        let mut runtime = process.integrand.take().unwrap();
        runtime.warm_up(&model)?;
        let ProcessIntegrand::CrossSection(generated) = &runtime else {
            unreachable!()
        };
        assert_eq!(
            generated.data.graph_terms[0]
                .fermi_surfaces
                .iter()
                .any(|sectors| !sectors.is_empty()),
            spectator,
        );
        runtimes.push(runtime);
    }
    let mut input = MomentumSpaceEvaluationInput {
        loop_momenta: vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0))],
        integrator_weight: F(1.0),
        graph_id: Some(0),
        group_id: None,
        orientation: None,
        channel_id: None,
    };
    let baseline = runtimes[0]
        .evaluate_momentum_configuration(&model, &input, false)?
        .integrand_result;
    assert!(baseline.re.0.abs() > 0.0 && baseline.im.0.abs() < 1e-14);
    let root = (2.5_f64.powi(2) - 1.0).sqrt();
    let mut last = None;
    for radius in [1.0_f64, 3.0] {
        input.loop_momenta.truncate(1);
        input
            .loop_momenta
            .push(ThreeMomentum::new(F(0.0), F(radius), F(0.0)));
        let result = runtimes[1].evaluate_momentum_configuration(&model, &input, false)?;
        // Gorda--Österman--Säppi, Phys.Rev.D106 105026, Eqs.(12--14):
        // the doubled fermion contour contains both the occupied-ball term
        // -Theta(mu-E)/(4E^3) and the boundary -delta(E-mu)/(4E^2).
        // The latter acts on the LU-rescaled spectator momentum; t^3 is
        // the extra loop's original LU spatial measure.
        let r = root * radius;
        let energy = (r * r + 9.0).sqrt();
        let shell_scale = 4.0 / r;
        let profile = 2.0 / std::f64::consts::PI.sqrt() * (-shell_scale.powi(2)).exp();
        let bulk = if energy < 5.0 {
            -1.0 / (4.0 * energy.powi(3))
        } else {
            0.0
        };
        let boundary = -profile * 16.0 / (20.0 * r.powi(4));
        let expected =
            baseline.re.0 * root.powi(3) * (bulk + boundary) / std::f64::consts::TAU.powi(3);
        let value = result.integrand_result;
        assert!(
            value.re.0.abs() < 1e-11 * expected.abs()
                && (value.im.0 / expected - 1.0).abs() < 1e-10,
            "radius={radius}: {value:?} != i*{expected}"
        );
        let events = result
            .event_groups
            .iter()
            .flat_map(|group| group.iter())
            .collect_vec();
        assert!(!events.is_empty());
        let sum = events
            .iter()
            .fold(Complex::new_re(F(0.0)), |sum, event| sum + event.weight);
        assert!((sum.im.0 / expected - 1.0).abs() < 1e-10);
        for event in events {
            let residual = event.kinematic_configuration.0.iter().fold(
                FourMomentum::from_args(F(0.0), F(0.0), F(0.0), F(0.0)),
                |sum, p| sum + *p,
            ) - event.kinematic_configuration.1.iter().fold(
                FourMomentum::from_args(F(0.0), F(0.0), F(0.0), F(0.0)),
                |sum, p| sum + *p,
            );
            assert!(residual.temporal.value.0.abs() < 1e-12);
            assert!(residual.spatial.norm_squared().0 < 1e-24);
        }
        last = Some((input.clone(), value));
    }
    // Integrating the spectator at fixed scalar-cut coordinates gives the
    // mass derivative of the occupied-ball tadpole,
    // -log((mu+p_F)/m)/(8*pi^2) = -log(3)/(8*pi^2).
    // The compactification puts the step at x=1/2 on both midpoint grids;
    // Richardson extrapolation removes their leading quadrature error.
    let mut integrals = Vec::new();
    for count in [128, 256] {
        let mut integral = 0.0;
        for index in 0..count {
            let x = (index as f64 + 0.5) / count as f64;
            let radius = 4.0 / root * x / (1.0 - x);
            input.loop_momenta[1] = ThreeMomentum::new(F(0.0), F(radius), F(0.0));
            let value = runtimes[1]
                .evaluate_momentum_configuration(&model, &input, false)?
                .integrand_result;
            integral += value.im.0 / baseline.re.0
                * 4.0
                * std::f64::consts::PI
                * radius.powi(2)
                * (4.0 / root / (1.0 - x).powi(2))
                / count as f64;
        }
        integrals.push(integral);
    }
    let integral = (4.0 * integrals[1] - integrals[0]) / 3.0;
    let expected = -3.0_f64.ln() / (8.0 * std::f64::consts::PI.powi(2));
    assert!((integral / expected - 1.0).abs() < 1e-8);

    let bytes = bincode::encode_to_vec(&runtimes[1], bincode::config::standard())?;
    let mut symbols = Vec::new();
    State::export(&mut symbols)?;
    let state_map = State::import(&mut Cursor::new(symbols), None)?;
    let context = GammaLoopContextContainer {
        model: &model,
        state_map: &state_map,
    };
    let (mut restored, consumed): (ProcessIntegrand, usize) =
        bincode::decode_from_slice_with_context(&bytes, bincode::config::standard(), context)?;
    assert_eq!(consumed, bytes.len());
    restored.activate_runtime_backends_after_load(false)?;
    restored.warm_up(&model)?;
    let (input, expected) = last.unwrap();
    assert_eq!(
        restored
            .evaluate_momentum_configuration(&model, &input, false)?
            .integrand_result,
        expected,
    );
    Ok(())
}
