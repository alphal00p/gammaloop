use std::io::Cursor;

use symbolica::state::State;

use super::*;
use crate::{
    initialisation::test_initialise,
    integrands::process::{MomentumSpaceEvaluationInput, ProcessIntegrand},
    processes::Amplitude,
    settings::runtime::{HFunction, HFunctionSettings},
    utils::load_generic_model,
};

/// The simplified `sm-thermal` model with vacuum-subtracted zero-temperature
/// amplitude settings; `muB` is its only nonzero chemical potential.
fn cold_dense_setup() -> Result<(
    Model,
    crate::model::InputParamCard<F<f64>>,
    GlobalSettings,
    RuntimeSettings,
)> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    let mut card = crate::model::InputParamCard::from_file(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../assets/models/json/sm/restrict_thermal.json"),
    )?;
    model.simplify(&mut card)?;
    let global: GlobalSettings = toml::from_str(
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
    Ok((model, card, global, settings))
}

#[test]
fn production_fermi_amplitude_preserves_orientation_sum_and_saved_state() -> Result<()> {
    let (mut model, mut card, mut global, settings) = cold_dense_setup()?;
    let mass = 3.0_f64;
    let chemical_potential = 5.0_f64;
    card.insert("MB".into(), Complex::new_re(F(mass)));
    card.insert("muB".into(), Complex::new_re(F(3.0 * chemical_potential)));
    model.apply_param_card(&card)?;
    let graphs = Graph::from_string(
        r#"digraph fermi_vacuum_cycle {
            node [num=1]; edge [num=1 particle="b"];
            A -> B [id=0 lmb_id=0]; B -> A [id=1];
        }"#,
        &model,
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
        // The clockwise energy contour of D^-1, D=q0²-E²+i0, is
        // -Theta(E-mu)/(2E). Differentiating in m² gives D^-2 and
        // vacuum subtraction leaves
        // -Theta(mu-E)/(4E³) - delta(E-mu)/(4E²).
        // For pF=sqrt(mu²-m²), inserting the normalized radial profile
        // localizes the second term to -h(pF/r) pF²/(4 mu r⁴).
        // The unit-numerator graph contributes exactly i/(2pi)³ from
        // its single loop; no symmetry factor or spatial Jacobian is
        // added by fixed-momentum evaluation.
        let fermi_radius_squared = chemical_potential.powi(2) - mass.powi(2);
        let mut last_point = None;
        for (profile, sigma) in [
            (HFunction::Exponential, 1.0_f64),
            (HFunction::PolyExponential, 1.3_f64),
        ] {
            let ProcessIntegrand::Amplitude(generated) = &mut *runtime else {
                unreachable!()
            };
            generated.settings.lu_h_function.function = profile.clone();
            generated.settings.lu_h_function.sigma = sigma;
            generated.settings.lu_h_function.power = Some(0);
            for coordinates in [[2.0_f64, -1.0, 0.5], [3.0, 1.0, 1.0], [5.0, -2.0, 1.0]] {
                let input = MomentumSpaceEvaluationInput {
                    loop_momenta: vec![ThreeMomentum::new(
                        F(coordinates[0]),
                        F(coordinates[1]),
                        F(coordinates[2]),
                    )],
                    integrator_weight: F(1.0),
                    graph_id: Some(0),
                    group_id: None,
                    orientation: None,
                    channel_id: None,
                };
                let radius_squared = coordinates.into_iter().map(|x| x * x).sum::<f64>();
                let energy = (radius_squared + mass.powi(2)).sqrt();
                let scale = (fermi_radius_squared / radius_squared).sqrt();
                let scaled = scale / sigma;
                // Evaluate the two normalized profiles explicitly, independently
                // of h(), h_dual(), the sector extractor and the localizer.
                let exponent = match profile {
                    HFunction::Exponential => -scaled.powi(2),
                    HFunction::PolyExponential => 2.0 - scaled.powi(2) - scaled.powi(-2),
                    _ => unreachable!(),
                };
                let profile_value = 2.0 * exponent.exp() / (std::f64::consts::PI.sqrt() * sigma);
                let bulk = if energy < chemical_potential {
                    -1.0 / (4.0 * energy.powi(3))
                } else {
                    0.0
                };
                let surface = -profile_value * fermi_radius_squared
                    / (4.0 * chemical_potential * radius_squared.powi(2));
                let expected_imaginary = (bulk + surface) / std::f64::consts::TAU.powi(3);
                let value = runtime
                    .evaluate_momentum_configuration(&model, &input, false)?
                    .integrand_result;
                assert!(value.re.0.is_finite() && value.im.0.is_finite());
                assert!(
                    value.re.0.abs() < 1e-12 * expected_imaginary.abs()
                        && (value.im.0 / expected_imaginary - 1.0).abs() < 1e-11,
                    "explicit={explicit}, profile={profile:?}, sigma={sigma}, momentum={coordinates:?}: {value:?} != i*{expected_imaginary}"
                );
                last_point = Some((input, value));
            }
        }
        let (input, value) = last_point.unwrap();
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
        let mut shifted_card = card.clone();
        shifted_card.insert("muB".into(), Complex::new_re(F(21.0)));
        shifted_model.apply_param_card(&shifted_card)?;
        restored.warm_up(&shifted_model)?;
        let shifted = restored
            .evaluate_momentum_configuration(&shifted_model, &input, false)?
            .integrand_result;
        assert_ne!(
            shifted, value,
            "warm-up must refresh the Fermi shell chemical potential"
        );

        // Localization inserts the profile as a unit scale integral.
        for h_function in [
            HFunctionSettings {
                function: HFunction::ExponentialCT,
                ..Default::default()
            },
            HFunctionSettings {
                power: Some(2),
                ..Default::default()
            },
            HFunctionSettings {
                sigma: 0.0,
                ..Default::default()
            },
        ] {
            restored.get_mut_settings().lu_h_function = h_function;
            let error = restored.warm_up(&shifted_model).unwrap_err();
            assert!(error.to_string().contains("h_function"), "{error}");
        }
    }
    Ok(())
}

#[test]
fn zero_chemical_potentials_have_no_fermi_surface() -> Result<()> {
    // In sm-thermal, muQ = muLe = 0 make the chemical potentials of G+ and e-
    // vanish. Their raised lines previously failed generation (a boson) and
    // the first sample (a degenerate massless onset).
    let (model, _, mut global, settings) = cold_dense_setup()?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    let input = MomentumSpaceEvaluationInput {
        loop_momenta: vec![ThreeMomentum::new(F(2.0), F(-1.0), F(0.5))],
        integrator_weight: F(1.0),
        graph_id: Some(0),
        group_id: None,
        orientation: None,
        channel_id: None,
    };
    for (particle, fermi_sectors) in [("G+", false), ("e-", true)] {
        let graphs = Graph::from_string(
            format!(
                r#"digraph zero_mu_cycle {{
                    node [num=1]; edge [num=1 particle="{particle}"];
                    A -> B [id=0 lmb_id=0]; B -> A [id=1];
                }}"#
            ),
            &model,
        )?;
        let mut values = Vec::new();
        for vacuum_subtraction in [false, true] {
            global.generation.medium.vacuum_subtraction = vacuum_subtraction;
            let mut amplitude = Amplitude::from_graph_list("zero_mu_cycle", graphs.clone())?;
            amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
            amplitude.build_integrand(
                &model,
                "zero_mu_cycle",
                &global,
                (&settings).into(),
                &pool,
            )?;
            let runtime = amplitude.integrand.as_mut().unwrap();
            runtime.warm_up(&model)?;
            let ProcessIntegrand::Amplitude(generated) = &*runtime else {
                unreachable!()
            };
            assert_eq!(
                !generated.data.graph_terms[0].fermi_surfaces.is_empty(),
                fermi_sectors,
                "{particle}"
            );
            values.push(
                runtime
                    .evaluate_momentum_configuration(&model, &input, false)?
                    .integrand_result,
            );
        }
        // Without a Fermi sea, the medium contribution is only binary64
        // rounding of the full (vacuum) value.
        let [full, medium] = [&values[0], &values[1]].map(|value| value.re.0.hypot(value.im.0));
        assert!(full > 0.0, "{particle}: {values:?}");
        assert!(
            medium <= 8.0 * f64::EPSILON * full,
            "{particle}: {values:?}"
        );
    }
    Ok(())
}

#[test]
fn thermal_amplitude_roundtrip_revalidates_runtime_and_model() -> Result<()> {
    use three_dimensional_reps::MediumMode;

    test_initialise()?;
    let mut model = load_generic_model("sm");
    model.get_parameter_mut("MW")?.value = Some(Complex::new_re(F(2.0)));
    model.get_parameter_mut("muG")?.value = Some(Complex::new_re(F(1.0)));
    let graphs = Graph::from_string(
        r#"digraph thermal_boson {
            node [num=1]; edge [num=1 particle="G+"];
            A -> A [id=0];
        }"#,
        &model,
    )?;
    let global: GlobalSettings = toml::from_str(
        r#"
[generation.medium]
mode = "thermodynamic_equilibrium"
[generation.uv]
subtract_uv = false
generate_integrated = false
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
inverse_temperature = 1.0
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
    let mut amplitude = Amplitude::from_graph_list("thermal_boson", graphs)?;
    amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
    amplitude.build_integrand(&model, "thermal_boson", &global, (&settings).into(), &pool)?;
    let runtime = amplitude.integrand.as_mut().unwrap();
    runtime.warm_up(&model)?;
    let input = MomentumSpaceEvaluationInput {
        loop_momenta: vec![ThreeMomentum::new(F(1.0), F(0.25), F(0.5))],
        integrator_weight: F(1.0),
        graph_id: Some(0),
        group_id: None,
        orientation: None,
        channel_id: None,
    };
    let expected = runtime
        .evaluate_momentum_configuration(&model, &input, false)?
        .integrand_result;
    assert!(expected.re.0.is_finite() && expected.im.0.is_finite());
    assert_ne!(expected, Complex::new_re(F(0.0)));

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
    let ProcessIntegrand::Amplitude(generated) = &restored else {
        unreachable!()
    };
    let term = &generated.data.graph_terms[0];
    assert_eq!(
        term.graph.param_builder.medium_mode,
        MediumMode::ThermodynamicEquilibrium
    );
    assert_eq!(
        term.param_builder.medium_mode,
        MediumMode::ThermodynamicEquilibrium
    );
    restored.warm_up(&model)?;
    assert_eq!(
        restored
            .evaluate_momentum_configuration(&model, &input, false)?
            .integrand_result,
        expected,
    );

    restored.get_mut_settings().general.inverse_temperature = 0.0;
    let error = restored.warm_up(&model).unwrap_err();
    assert!(error.to_string().contains("inverse_temperature"), "{error}");
    restored.get_mut_settings().general.inverse_temperature = 1.0;
    // Both signs of bosonic saturation and a changed mass must be revalidated
    // from refreshed model slots after loading the compiled integrand.
    for (mass, mu) in [(2.0, 2.0), (2.0, -2.0), (2.0, 3.0), (0.5, 1.0)] {
        model.get_parameter_mut("MW")?.value = Some(Complex::new_re(F(mass)));
        model.get_parameter_mut("muG")?.value = Some(Complex::new_re(F(mu)));
        let error = restored.warm_up(&model).unwrap_err();
        assert!(error.to_string().contains("|mu| < m"), "{error}");
    }
    model.get_parameter_mut("MW")?.value = Some(Complex::new_re(F(2.0)));
    model.get_parameter_mut("muG")?.value = Some(Complex::new_re(F(1.0)));
    restored.warm_up(&model)?;
    assert_eq!(
        restored
            .evaluate_momentum_configuration(&model, &input, false)?
            .integrand_result,
        expected,
    );
    Ok(())
}
