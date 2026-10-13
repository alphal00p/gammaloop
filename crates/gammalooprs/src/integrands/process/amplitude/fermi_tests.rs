use std::io::Cursor;

use symbolica::state::State;

use super::*;
use crate::{
    initialisation::test_initialise,
    integrands::process::{
        MomentumSpaceEvaluationInput, ProcessIntegrand,
        evaluators::{ActiveF64Backend, EvaluatorMethod},
    },
    processes::Amplitude,
    settings::{
        global::OrientationPattern,
        runtime::{DiscreteGraphSamplingSettings, HFunction, HFunctionSettings},
    },
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

/// Vacuum and vacuum-subtracted zero-temperature energy contours of a doubled
/// fermion propagator at spatial momentum squared `radius_squared`.
///
/// The clockwise energy contour of D^-1, D=q0²-E²+i0, is -Theta(E-mu)/(2E).
/// Differentiating in m² gives D^-2 and vacuum subtraction leaves
/// -Theta(mu-E)/(4E³) - delta(E-mu)/(4E²), while the vacuum contour is 1/(4E³).
/// For pF=sqrt(mu²-m²), inserting the normalized radial profile localizes the
/// delta to -h(pF/r) pF²/(4 mu r⁴).
fn double_pole_contours(
    mass: f64,
    chemical_potential: f64,
    profile: &HFunction,
    sigma: f64,
    radius_squared: f64,
) -> (f64, f64) {
    let fermi_radius_squared = chemical_potential.powi(2) - mass.powi(2);
    let energy = (radius_squared + mass.powi(2)).sqrt();
    let scaled = (fermi_radius_squared / radius_squared).sqrt() / sigma;
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
    let surface =
        -profile_value * fermi_radius_squared / (4.0 * chemical_potential * radius_squared.powi(2));
    (1.0 / (4.0 * energy.powi(3)), bulk + surface)
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
    for (explicit, compile) in itertools::iproduct!([false, true], [false, true]) {
        global.generation.explicit_orientation_sum_only = explicit;
        // Summed and iterative evaluators require the orientation catalogue.
        let evaluator = &mut global.generation.evaluator;
        evaluator.compile = compile;
        evaluator.summed = !explicit;
        evaluator.summed_function_map = !explicit;
        evaluator.iterative_orientation_optimization = !explicit;
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
        let orientation_count = generated.data.graph_terms[0].orientations.len();
        assert_eq!(
            generated.data.graph_terms[0]
                .threshold_counterterm
                .generic_evaluator_count(),
            0
        );
        assert_eq!(
            generated.active_f64_backend(),
            if compile {
                ActiveF64Backend::Symjit
            } else {
                ActiveF64Backend::Eager
            }
        );
        let methods = if explicit {
            &[EvaluatorMethod::SingleParametric][..]
        } else {
            &[
                EvaluatorMethod::Iterative,
                EvaluatorMethod::Summed,
                EvaluatorMethod::SummedFunctionMap,
                EvaluatorMethod::SingleParametric,
            ][..]
        };
        let mut last_point = None;
        // Reversing mu moves the Fermi sea to the opposite orientation; the
        // orientation sum is unchanged. Finish at mu > 0 for the round trip.
        for (method, mu_sign) in itertools::iproduct!(methods, [-1.0, 1.0]) {
            let mut signed_model = model.clone();
            let mut signed_card = card.clone();
            signed_card.insert(
                "muB".into(),
                Complex::new_re(F(mu_sign * 3.0 * chemical_potential)),
            );
            signed_model.apply_param_card(&signed_card)?;
            runtime.warm_up(&signed_model)?;
            for (profile, sigma) in [
                (HFunction::Exponential, 1.0_f64),
                (HFunction::PolyExponential, 1.3_f64),
            ] {
                let ProcessIntegrand::Amplitude(generated) = &mut *runtime else {
                    unreachable!()
                };
                generated.settings.general.evaluator_method = method.clone();
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
                    let (_, medium) = double_pole_contours(
                        mass,
                        chemical_potential,
                        &profile,
                        sigma,
                        coordinates.into_iter().map(|x| x * x).sum(),
                    );
                    // The unit-numerator graph contributes exactly i/(2pi)³
                    // from its single loop; no symmetry factor or spatial
                    // Jacobian is added by fixed-momentum evaluation.
                    let expected_imaginary = medium / std::f64::consts::TAU.powi(3);
                    let value = runtime
                        .evaluate_momentum_configuration(&signed_model, &input, false)?
                        .integrand_result;
                    assert!(value.re.0.is_finite() && value.im.0.is_finite());
                    assert!(
                        value.re.0.abs() < 1e-12 * expected_imaginary.abs()
                            && (value.im.0 / expected_imaginary - 1.0).abs() < 1e-11,
                        "explicit={explicit}, compile={compile}, method={method:?}, mu_sign={mu_sign}, profile={profile:?}, sigma={sigma}, momentum={coordinates:?}: {value:?} != i*{expected_imaginary}"
                    );
                    if !explicit && *method == EvaluatorMethod::SingleParametric {
                        let mut selected_sum = Complex::new(F(0.0), F(0.0));
                        for orientation in 0..orientation_count {
                            let selected = MomentumSpaceEvaluationInput {
                                orientation: Some(orientation),
                                ..input.clone()
                            };
                            selected_sum += runtime
                                .evaluate_momentum_configuration(&signed_model, &selected, false)?
                                .integrand_result;
                        }
                        assert!(
                            selected_sum.re.0.abs() < 1e-12 * expected_imaginary.abs()
                                && (selected_sum.im.0 / expected_imaginary - 1.0).abs() < 1e-11,
                            "selected orientations must sum to the analytic contour: {selected_sum:?} != i*{expected_imaginary}"
                        );
                    }
                    last_point = Some((input, value));
                }
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
        // Reactivate the frozen backend exactly as a loaded process does.
        restored.activate_runtime_backends_after_load(false)?;
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
fn production_fermi_orientations_preserve_filters_and_mc_weights() -> Result<()> {
    let (mut model, mut card, mut global, settings) = cold_dense_setup()?;
    card.insert("MB".into(), Complex::new_re(F(3.0)));
    card.insert("muB".into(), Complex::new_re(F(15.0)));
    model.apply_param_card(&card)?;
    let graphs = Graph::from_string(
        r#"digraph filtered_fermi_cycle {
            node [num=1]; edge [num=1 particle="b"];
            A -> B [id=0 lmb_id=0]; B -> A [id=1];
        }"#,
        &model,
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    // Outside the Fermi sea the regular medium term vanishes: a nonzero
    // result here certifies that the selected Fermi surface contributes.
    let input = MomentumSpaceEvaluationInput {
        loop_momenta: vec![ThreeMomentum::new(F(5.0), F(-2.0), F(1.0))],
        integrator_weight: F(1.0),
        graph_id: Some(0),
        group_id: None,
        orientation: None,
        channel_id: None,
    };
    let continuous = Sample::Continuous(F(1.0), vec![F(0.31), F(0.27), F(0.61)]);
    for compile in [false, true] {
        global.generation.evaluator.compile = compile;
        global.generation.orientation_pattern = Default::default();
        let mut full = Amplitude::from_graph_list("filtered_fermi_cycle", graphs.clone())?;
        full.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        full.build_integrand(
            &model,
            "filtered_fermi_cycle",
            &global,
            (&settings).into(),
            &pool,
        )?;
        let runtime = full.integrand.as_mut().unwrap();
        runtime.warm_up(&model)?;
        let ProcessIntegrand::Amplitude(generated) = &*runtime else {
            unreachable!()
        };
        let term = &generated.data.graph_terms[0];
        assert!(!term.fermi_surfaces.is_empty());
        let count = term.orientations.len();
        assert_eq!(
            term.original_integrand
                .production_orientation_ids()
                .iter()
                .map(|id| id.0)
                .collect_vec(),
            (0..count).collect_vec()
        );
        let pattern = OrientationPattern::from_orientation(term.orientations.last().unwrap());
        let selected_ids = term
            .orientations
            .iter_enumerated()
            .filter_map(|(id, orientation)| pattern.filter(orientation).then_some(id.0))
            .collect_vec();
        assert!(!selected_ids.is_empty() && selected_ids.len() < count && selected_ids[0] > 0);
        global.generation.orientation_pattern = pattern.clone();
        let mut subset = Amplitude::from_graph_list("filtered_fermi_cycle", graphs.clone())?;
        subset.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        subset.build_integrand(
            &model,
            "filtered_fermi_cycle",
            &global,
            (&settings).into(),
            &pool,
        )?;
        let filtered = subset.integrand.as_mut().unwrap();
        filtered.warm_up(&model)?;
        let ProcessIntegrand::Amplitude(generated) = &*filtered else {
            unreachable!()
        };
        let term = &generated.data.graph_terms[0];
        assert!(!term.fermi_surfaces.is_empty());
        assert_eq!(
            term.original_integrand
                .production_orientation_ids()
                .iter()
                .map(|id| id.0)
                .collect_vec(),
            selected_ids
        );
        let mut nonzero_filtered_surface = false;
        for mu_sign in [-1.0, 1.0] {
            card.insert("muB".into(), Complex::new_re(F(mu_sign * 15.0)));
            model.apply_param_card(&card)?;
            *runtime.get_mut_settings() = settings.clone();
            runtime.warm_up(&model)?;
            filtered.warm_up(&model)?;
            let mut selected_values = Vec::new();
            for orientation in 0..count {
                selected_values.push(
                    runtime
                        .evaluate_momentum_configuration(
                            &model,
                            &MomentumSpaceEvaluationInput {
                                orientation: Some(orientation),
                                ..input.clone()
                            },
                            false,
                        )?
                        .integrand_result,
                );
            }
            let expected = selected_ids
                .iter()
                .fold(Complex::<F<f64>>::new_zero(), |sum, id| {
                    sum + selected_values[*id]
                });
            nonzero_filtered_surface |= expected.im.0.abs() > 1e-12;
            runtime.get_mut_settings().general.orientation_pat = pattern.clone();
            runtime.warm_up(&model)?;
            let ProcessIntegrand::Amplitude(generated) = &*runtime else {
                unreachable!()
            };
            assert_eq!(
                generated.data.graph_terms[0]
                    .orientation_filter
                    .included_iter()
                    .map(|id| id.0)
                    .collect_vec(),
                selected_ids
            );
            // Both filters renumber the visible catalogue. Its index zero
            // must still select the original nonzero production ID.
            for (orientation, expected) in [
                (None, expected),
                (Some(0), selected_values[selected_ids[0]]),
            ] {
                let selected = MomentumSpaceEvaluationInput {
                    orientation,
                    ..input.clone()
                };
                for value in [
                    runtime
                        .evaluate_momentum_configuration(&model, &selected, false)?
                        .integrand_result,
                    filtered
                        .evaluate_momentum_configuration(&model, &selected, false)?
                        .integrand_result,
                ] {
                    assert!(
                        (value.re.0 - expected.re.0).abs() < 1e-13
                            && (value.im.0 - expected.im.0).abs() < 1e-13,
                        "compile={compile}, mu_sign={mu_sign}, orientation={orientation:?}: filtered {value:?} != {expected:?}"
                    );
                }
            }

            *runtime.get_mut_settings() = settings.clone();
            runtime.warm_up(&model)?;
            let summed = runtime
                .evaluate_sample_precise(
                    &continuous,
                    &model,
                    continuous.get_weight(),
                    false,
                    Complex::new_zero(),
                )?
                .try_into_f64()?;
            let expected = summed.integrand_result
                * summed.parameterization_jacobian.unwrap()
                * summed.integrator_weight;
            assert!(expected.im.0.is_finite() && expected.im.0.abs() > 1e-12);
            runtime.get_mut_settings().sampling =
                SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                    sample_orientations: true,
                    sampling_type: DiscreteGraphSamplingType::Default(
                        settings.sampling.get_parameterization_settings().unwrap(),
                    ),
                    ..Default::default()
                });
            runtime.warm_up(&model)?;
            let mut average = Complex::<F<f64>>::new_zero();
            for orientation in 0..count {
                // Deliberately unequal probabilities expose missing or doubled
                // inverse-probability weights in the actual MC sample path.
                let probability = 2.0 * (orientation + 1) as f64 / (count * (count + 1)) as f64;
                let sample = Sample::Discrete(
                    F(1.0 / probability),
                    0,
                    Some(Box::new(Sample::Discrete(
                        F(1.0),
                        orientation,
                        Some(Box::new(continuous.clone())),
                    ))),
                );
                let selected = runtime
                    .evaluate_sample_precise(
                        &sample,
                        &model,
                        sample.get_weight(),
                        false,
                        Complex::new_zero(),
                    )?
                    .try_into_f64()?;
                assert_eq!(selected.integrator_weight, F(1.0 / probability));
                average += selected.integrand_result
                    * selected.parameterization_jacobian.unwrap()
                    * selected.integrator_weight
                    * F(probability);
            }
            assert!(
                (average.re.0 - expected.re.0).abs() < 1e-11 * expected.im.0.abs()
                    && (average.im.0 / expected.im.0 - 1.0).abs() < 1e-11,
                "compile={compile}, mu_sign={mu_sign}: MC expectation {average:?} != {expected:?}"
            );
        }
        assert!(nonzero_filtered_surface);
    }
    Ok(())
}

#[test]
fn production_fermi_products_localize_independent_surfaces() -> Result<()> {
    let (mut model, mut card, mut global, settings) = cold_dense_setup()?;
    let mass = 3.0_f64;
    let chemical_potential = 5.0_f64;
    card.insert("MB".into(), Complex::new_re(F(mass)));
    card.insert("muB".into(), Complex::new_re(F(3.0 * chemical_potential)));
    model.apply_param_card(&card)?;
    // Two doubled fermion cycles meet at one vertex, so each orientation of
    // the product carries both Fermi deltas.
    let graphs = Graph::from_string(
        r#"digraph fermi_figure_eight {
            node [num=1]; edge [num=1 particle="b"];
            A -> B [id=0 lmb_id=0]; B -> A [id=1];
            A -> C [id=2 lmb_id=1]; C -> A [id=3];
        }"#,
        &model,
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()?;
    for explicit in [false, true] {
        global.generation.explicit_orientation_sum_only = explicit;
        let mut amplitude = Amplitude::from_graph_list("fermi_figure_eight", graphs.clone())?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(
            &model,
            "fermi_figure_eight",
            &global,
            (&settings).into(),
            &pool,
        )?;
        let runtime = amplitude.integrand.as_mut().unwrap();
        runtime.warm_up(&model)?;
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
            // Inside, inside; inside, outside; and outside, inside the Fermi sea.
            for momenta in [
                [[2.0_f64, -1.0, 0.5], [3.0, 1.0, 1.0]],
                [[3.0, 1.0, 1.0], [5.0, -2.0, 1.0]],
                [[5.0, -2.0, 1.0], [2.0, -1.0, 0.5]],
            ] {
                let input = MomentumSpaceEvaluationInput {
                    loop_momenta: momenta
                        .iter()
                        .map(|k| ThreeMomentum::new(F(k[0]), F(k[1]), F(k[2])))
                        .collect(),
                    integrator_weight: F(1.0),
                    graph_id: Some(0),
                    group_id: None,
                    orientation: None,
                    channel_id: None,
                };
                let [(vacuum_0, medium_0), (vacuum_1, medium_1)] = momenta.map(|k| {
                    double_pole_contours(
                        mass,
                        chemical_potential,
                        &profile,
                        sigma,
                        k.into_iter().map(|x| x * x).sum(),
                    )
                });
                // The loops factorize, each contributing i/(2pi)³. Vacuum
                // subtraction removes only the product of the vacuum contours,
                // retaining single-surface terms on either vacuum loop and the
                // double-delta intersection of the two independent surfaces.
                let expected_real = -((vacuum_0 + medium_0) * (vacuum_1 + medium_1)
                    - vacuum_0 * vacuum_1)
                    / std::f64::consts::TAU.powi(6);
                let value = runtime
                    .evaluate_momentum_configuration(&model, &input, false)?
                    .integrand_result;
                assert!(
                    value.re.0.is_finite()
                        && value.im.0.abs() < 1e-12 * expected_real.abs()
                        && (value.re.0 / expected_real - 1.0).abs() < 1e-11,
                    "explicit={explicit}, profile={profile:?}, sigma={sigma}, momenta={momenta:?}: {value:?} != {expected_real}"
                );
            }
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

/// The finite-temperature amplitude of `dot`, with interpreted evaluators and
/// neither UV nor threshold subtraction, built for `inverse_temperature = 1`.
fn thermal_amplitude(model: &Model, name: &str, dot: &str) -> Result<Amplitude> {
    let graphs = Graph::from_string(dot, model)?;
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
    let mut amplitude = Amplitude::from_graph_list(name, graphs)?;
    amplitude.preprocess(model, &global.generation, &(&settings).into(), &pool)?;
    amplitude.build_integrand(model, name, &global, (&settings).into(), &pool)?;
    Ok(amplitude)
}

#[test]
fn thermal_amplitude_roundtrip_revalidates_runtime_and_model() -> Result<()> {
    use three_dimensional_reps::MediumMode;

    test_initialise()?;
    let mut model = load_generic_model("sm");
    model.get_parameter_mut("MW")?.value = Some(Complex::new_re(F(2.0)));
    model.get_parameter_mut("muG")?.value = Some(Complex::new_re(F(1.0)));
    let mut amplitude = thermal_amplitude(
        &model,
        "thermal_boson",
        r#"digraph thermal_boson {
            node [num=1]; edge [num=1 particle="G+"];
            A -> A [id=0];
        }"#,
    )?;
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

#[test]
fn thermal_amplitude_warm_up_rejects_cycle_poles_with_equal_masses() -> Result<()> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    let [top, w] =
        ["MT", "MW"].map(|name| Atom::var(model.get_parameter(name).name.0).to_canonical_string());
    model.get_parameter_mut("MT")?.value = Some(Complex::new_re(F(2.0)));
    model.get_parameter_mut("MW")?.value = Some(Complex::new_re(F(3.0)));
    // Distinct mass parameters generate the cycle as a divided difference.
    let mut amplitude = thermal_amplitude(
        &model,
        "thermal_cycle",
        &format!(
            r#"digraph thermal_cycle {{
                node [num=1]; edge [num=1];
                A -> B [particle="G+" mass="{top}"];
                B -> A [particle="G+" mass="{w}"];
            }}"#
        ),
    )?;
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
    let result = runtime
        .evaluate_momentum_configuration(&model, &input, false)?
        .integrand_result;
    assert!(result.re.0.is_finite() && result.im.0.is_finite());

    model.get_parameter_mut("MW")?.value = Some(Complex::new_re(F(2.0)));
    let error = runtime.warm_up(&model).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("vanishes for every loop momentum"),
        "{error}"
    );
    Ok(())
}
