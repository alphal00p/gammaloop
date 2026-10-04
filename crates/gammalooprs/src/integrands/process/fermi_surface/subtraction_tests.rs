use color_eyre::Result;
use spenso::algebra::complex::Complex;
use symbolica::atom::AtomCore;

use crate::{
    dot,
    graph::{Graph, parse::IntoGraph},
    initialisation::test_initialise,
    integrands::process::{
        MomentumSpaceEvaluationInput, ProcessIntegrand,
        param_builder::{ParamBuilderGraph, ThermalDistributionReplacement},
    },
    momentum::ThreeMomentum,
    numerator::symbolica_ext::NumeratorAtomExt,
    processes::{Amplitude, AmplitudeGraph},
    settings::{
        GlobalSettings, RuntimeSettings,
        global::{GenerationSettings, MediumMode},
        runtime::HFunction,
    },
    utils::{F, GS, load_generic_model},
};

use super::FermiSurfaceSector;

#[test]
fn fermi_sectors_retain_vacuum_subtraction_and_both_uv_contributions() -> Result<()> {
    test_initialise()?;
    // The repeated fermion pole supplies N'(E_d). A proper scalar tadpole
    // counterterm leaves that cycle in its thermal cograph, so both local and
    // integrated UV contributions must reach the same Fermi-sector boundary.
    // Every numerator is the constant on its own edge or vertex.
    let original: AmplitudeGraph = dot!(digraph fermi_cycle_with_uv_tadpole {
        node [num=1]
        edge [num=1]
        A -> B [id=0 particle="d" lmb_id=0]
        B -> A [id=1 particle="d"]
        A -> A [id=2 particle="H" lmb_id=1]
    })?;
    let mut generated = Vec::new();
    for (vacuum_subtraction, subtract_uv, generate_integrated) in [
        (false, false, false),
        (true, false, false),
        (true, true, false),
        (true, true, true),
    ] {
        let mut settings = GenerationSettings::default();
        settings.medium.mode = MediumMode::ZeroTemperatureEquilibrium;
        settings.medium.vacuum_subtraction = vacuum_subtraction;
        settings.threshold_subtraction.enable_thresholds = false;
        settings.explicit_orientation_sum_only = true;
        settings.uv.subtract_uv = subtract_uv;
        settings.uv.generate_integrated = generate_integrated;
        settings.uv.softct = false;
        let mut amplitude = original.clone();
        amplitude.generate_cff(&settings)?;
        amplitude.build_integrands(&settings, crate::utils::vakint()?)?;

        // Exercise the actual evaluator input before resolving any retained
        // numerator definitions. Sector collection keeps those bodies opaque.
        let expression = &amplitude.derived_data.all_mighty_integrand;
        let (bulk, sectors) = FermiSurfaceSector::extract(&amplitude.graph, expression)?;
        assert!(
            !sectors.is_empty(),
            "the generated cycle must have Fermi support"
        );
        assert!(!FermiSurfaceSector::is_present(&bulk)?);
        let reconstructed = sectors.iter().fold(bulk, |sum, sector| {
            sum + sector
                .product
                .factors()
                .iter()
                .zip(&sector.orientations)
                .fold(
                    sector.coefficient.clone(),
                    |coefficient, (factor, orientation)| {
                        coefficient
                            * GS.thermal_distribution(
                                factor.edge_id.0 as i64,
                                factor.derivative_order as i64,
                                0,
                                factor.sign,
                                *orientation,
                            )
                    },
                )
        });
        let difference = reconstructed - expression.unwrap_function(GS.thermal_weight_wrapper);
        assert!(
            difference.expand().is_zero(),
            "sector extraction changed UV={subtract_uv}, integrated={generate_integrated}: {difference}"
        );

        // These are scalar diagnostic copies. Only now resolve the definitions
        // to compare separately generated UV settings with distinct scope tags.
        generated.push(amplitude.derived_data.resolved_integrand()?);
    }

    let raw = &generated[0];
    let vacuum = original.graph.make_thermal_distributions_explicit(
        raw,
        MediumMode::Vacuum,
        original.graph.iter_edge_ids(),
        ThermalDistributionReplacement::All,
    )?;
    let difference = (&generated[1] - (raw - vacuum)).unwrap_function(GS.thermal_weight_wrapper);
    assert!(
        difference.expand().is_zero(),
        "vacuum subtraction must act on the complete distribution weight: {difference}"
    );

    for (name, contribution) in [
        ("local UV", &generated[2] - &generated[1]),
        ("integrated UV", &generated[3] - &generated[2]),
    ] {
        let (_, sectors) = FermiSurfaceSector::extract(&original.graph, &contribution)?;
        assert!(
            sectors
                .iter()
                .any(|sector| !sector.coefficient.expand().is_zero()),
            "the {name} contribution must retain the untouched fermion cycle's Fermi delta"
        );
    }
    Ok(())
}

#[test]
fn production_fermi_uv_contributions_match_factorized_pointwise_oracle() -> Result<()> {
    test_initialise()?;
    let mut model = load_generic_model("sm");
    let mut card = crate::model::InputParamCard::from_file(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../assets/models/json/sm/restrict_thermal.json"),
    )?;
    model.simplify(&mut card)?;
    let chemical_potential = 5.0;
    let scalar_mass = 3.0;
    card.insert("muB".into(), Complex::new_re(F(3.0 * chemical_potential)));
    card.insert("MH".into(), Complex::new_re(F(scalar_mass)));
    model.apply_param_card(&card)?;
    let mut global: GlobalSettings = toml::from_str(
        r#"
[generation]
explicit_orientation_sum_only = true
[generation.uv]
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
    let mut values = Vec::new();
    for (subtract_uv, generate_integrated) in [(false, false), (true, false), (true, true)] {
        global.generation.uv.subtract_uv = subtract_uv;
        global.generation.uv.generate_integrated = generate_integrated;
        let mut runtimes = Vec::new();
        for (name, source, medium) in [
            (
                "fermi_cycle_with_uv_tadpole",
                r#"digraph fermi_cycle_with_uv_tadpole {
                    node [num=1]; edge [num=1];
                    A -> B [id=0 particle="d" lmb_id=0];
                    B -> A [id=1 particle="d"];
                    A -> A [id=2 particle="H" lmb_id=1];
                }"#,
                MediumMode::ZeroTemperatureEquilibrium,
            ),
            (
                "independent_uv_tadpole",
                r#"digraph independent_uv_tadpole {
                    node [num=1]; edge [num=1];
                    A -> A [id=0 particle="H"];
                }"#,
                MediumMode::Vacuum,
            ),
        ] {
            global.generation.medium.mode = medium;
            global.generation.medium.vacuum_subtraction = medium != MediumMode::Vacuum;
            let graphs = Graph::from_string(source, &model)?;
            let mut amplitude = Amplitude::from_graph_list(name, graphs)?;
            amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
            amplitude.build_integrand(&model, name, &global, (&settings).into(), &pool)?;
            runtimes.push(amplitude.integrand.take().unwrap());
        }
        let [mut full, mut scalar]: [ProcessIntegrand; 2] = runtimes.try_into().ok().unwrap();
        let mut mode_values = Vec::new();
        for (m_uv, mu_r, localization_scale, sigma) in [(1.5, 2.0, 1.2, 0.8), (2.5, 3.0, 2.3, 1.4)]
        {
            for runtime in [&mut full, &mut scalar] {
                let runtime_settings = runtime.get_mut_settings();
                runtime_settings.general.m_uv = m_uv;
                runtime_settings.general.mu_r = mu_r;
                runtime_settings.general.renormalization_localization_scale = localization_scale;
                runtime_settings.lu_h_function.function = HFunction::Exponential;
                runtime_settings.lu_h_function.sigma = sigma;
                runtime.warm_up(&model)?;
            }
            // One point inside the Fermi sea and one outside: the latter has
            // only the localized delta contribution, so a correct bulk term
            // cannot mask a dropped UV contribution on the shell.
            for radius in [2.0_f64, 7.0] {
                let scalar_momentum = ThreeMomentum::new(F(0.5), F(-1.0), F(2.0));
                let full_input = MomentumSpaceEvaluationInput {
                    loop_momenta: vec![
                        ThreeMomentum::new(F(0.6 * radius), F(0.8 * radius), F(0.0)),
                        scalar_momentum,
                    ],
                    integrator_weight: F(1.0),
                    graph_id: Some(0),
                    group_id: None,
                    orientation: None,
                    channel_id: None,
                };
                let scalar_input = MomentumSpaceEvaluationInput {
                    loop_momenta: vec![scalar_momentum],
                    ..full_input.clone()
                };
                let scalar_value = scalar
                    .evaluate_momentum_configuration(&model, &scalar_input, false)?
                    .integrand_result;
                let normalization = (2.0 * std::f64::consts::PI).powi(3);
                if !subtract_uv {
                    let scalar_energy = (5.25 + scalar_mass * scalar_mass).sqrt();
                    let expected = -1.0 / (2.0 * scalar_energy * normalization);
                    assert!(scalar_value.re.0.abs() < 1.0e-15);
                    assert!((scalar_value.im.0 - expected).abs() < 1.0e-12 * expected.abs());
                }
                // Independently differentiate the one-propagator contour with
                // respect to m^2, then remove its vacuum part. At m=0 this is
                // -i/(2pi)^3 [theta(nu-r)/(4r^3) + delta(r-nu)/(4r^2)].
                // The normalized scale insertion maps the delta term to
                // h(nu/r) nu/(4r^4). No production CFF coefficient or Fermi
                // evaluator is used to construct this factor.
                let bulk = if radius < chemical_potential {
                    1.0 / (4.0 * radius.powi(3))
                } else {
                    0.0
                };
                let scale = chemical_potential / radius;
                let profile =
                    2.0 * (-(scale / sigma).powi(2)).exp() / (std::f64::consts::PI.sqrt() * sigma);
                let shell = profile * chemical_potential / (4.0 * radius.powi(4));
                let cycle = Complex::new(F(0.0), F(-(bulk + shell) / normalization));
                let expected = cycle * scalar_value;
                let actual = full
                    .evaluate_momentum_configuration(&model, &full_input, false)?
                    .integrand_result;
                assert_ne!(expected.re, F(0.0));
                let tolerance = 2.0e-11 * expected.re.0.abs();
                assert!(
                    (actual.re.0 - expected.re.0).abs() < tolerance
                        && (actual.im.0 - expected.im.0).abs() < tolerance,
                    "UV={subtract_uv}, integrated={generate_integrated}, m_uv={m_uv}, r={radius}: {actual:?} != {expected:?}"
                );
                mode_values.push(actual);
            }
        }
        values.push(mode_values);
    }
    for (name, before, after) in [("local", 0, 1), ("integrated", 1, 2)] {
        assert!(
            values[before]
                .iter()
                .zip(&values[after])
                .all(|(before, after)| (after.re.0 - before.re.0).abs() > 1.0e-14),
            "the {name} UV contribution must be resolved at every selected point"
        );
    }
    Ok(())
}
