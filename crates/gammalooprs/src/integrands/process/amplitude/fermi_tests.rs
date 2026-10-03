use std::io::Cursor;

use symbolica::state::State;

use super::*;
use crate::{
    initialisation::test_initialise,
    integrands::process::{MomentumSpaceEvaluationInput, ProcessIntegrand},
    processes::Amplitude,
    utils::load_generic_model,
};

#[test]
fn production_fermi_amplitude_preserves_orientation_sum_and_saved_state() -> Result<()> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    model.get_parameter_mut("mud")?.value = Some(Complex::new_re(F(5.0)));
    let graphs = Graph::from_string(
        r#"digraph fermi_vacuum_cycle {
            node [num=1]; edge [num=1 particle="d"];
            A -> B [id=0 lmb_id=0]; B -> A [id=1];
        }"#,
        &model,
    )?;
    let mut global: GlobalSettings = toml::from_str(
        r#"
[generation.medium]
mode = "zero_temperature_equilibrium"
vacuum_subtraction = true
[generation.uv]
subtract_uv = false
generate_integrated = false
[generation.threshold_subtraction]
enable_thresholds = false
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
[kinematics]
e_cm = 10.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = []
helicities = []
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = false
"#,
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    let mut reference = None;
    for explicit in [false, true] {
        global.generation.explicit_orientation_sum_only = explicit;
        let mut amplitude = Amplitude::from_graph_list("fermi_vacuum_cycle", graphs.clone())?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(
            &model,
            "fermi_vacuum_cycle",
            &global,
            (&settings).into(),
            &pool,
        )?;
        let runtime = amplitude.integrand.as_mut().unwrap();
        runtime.warm_up(&model)?;
        let ProcessIntegrand::Amplitude(generated) = runtime else {
            unreachable!()
        };
        assert!(!generated.data.graph_terms[0].fermi_surfaces.is_empty());
        assert_eq!(
            generated.data.graph_terms[0]
                .threshold_counterterm
                .generic_evaluator_count(),
            0
        );
        let input = MomentumSpaceEvaluationInput {
            loop_momenta: vec![ThreeMomentum::new(F(2.0), F(-1.0), F(0.5))],
            integrator_weight: F(1.0),
            graph_id: Some(0),
            group_id: None,
            orientation: None,
            channel_id: None,
        };
        let result = runtime.evaluate_momentum_configuration(&model, &input, false)?;
        let value = result.integrand_result;
        assert!(value.re.0.is_finite() && value.im.0.is_finite());
        assert!(value.re.0 != 0.0 || value.im.0 != 0.0);
        if let Some(expected) = reference {
            let expected: Complex<F<f64>> = expected;
            assert!((value.re.0 - expected.re.0).abs() < 1e-12);
            assert!((value.im.0 - expected.im.0).abs() < 1e-12);
        } else {
            reference = Some(value);
        }

        // Round-trip the actual runtime object, including the new shell maps
        // and compiled coefficient definitions, without a filesystem fixture.
        let bytes = bincode::encode_to_vec(&*runtime, bincode::config::standard())?;
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
        restored.warm_up(&model)?;
        assert_eq!(
            restored
                .evaluate_momentum_configuration(&model, &input, false)?
                .integrand_result,
            value,
        );

        let mut shifted_model = model.clone();
        shifted_model.get_parameter_mut("mud")?.value = Some(Complex::new_re(F(7.0)));
        restored.warm_up(&shifted_model)?;
        let shifted = restored
            .evaluate_momentum_configuration(&shifted_model, &input, false)?
            .integrand_result;
        assert_ne!(
            shifted, value,
            "warm-up must refresh the Fermi shell chemical potential"
        );
    }
    Ok(())
}
