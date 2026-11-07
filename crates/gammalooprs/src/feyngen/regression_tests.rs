//! Regressions retained when the generator moved to FeynKit.
use crate::integrands::evaluation::StabilityStatus;
use crate::integrands::process::MomentumSpaceEvaluationInput;
use crate::model::{InputParamCard, Model};
use crate::momentum::ThreeMomentum;
use crate::processes::Process;
use crate::settings::RuntimeSettings;
use crate::settings::runtime::{RotationSetting, StabilityLevelSetting};
use crate::utils::F;
use crate::{
    initialisation::test_initialise,
    model::ModelGammaLoopExt,
    processes::ProcessDefinition,
    settings::{GlobalSettings, global::Parallelisation},
    utils::load_generic_model,
};
use color_eyre::Result;
use feynkit_generator::{
    GenerationFilter, GenerationOptions, GenerationType, ParticleSelector,
    Process as GenerationProcess, SewnFilterOptions, SnailFilterOptions, TadpoleFilterOptions,
};
use std::collections::BTreeMap;
use symbolica::atom::{Atom, AtomCore};

#[test]
fn complex_ckm_generation_preserves_named_and_inline_couplings() -> color_eyre::Result<()> {
    use crate::{
        model::{InputParamCard, ParameterType},
        utils::F,
    };

    for inline in [false, true] {
        let mut model = load_generic_model("sm");
        // UFO-real inputs remain real. Their derived CKM value can nevertheless
        // be complex; ordinary Feynman i factors are not intrinsic parameters.
        // The JSON defaults to diagonal CKM, so set every independent input.
        model.apply_param_card(&InputParamCard::<F<f64>>::from_str(
            "lamWS = [0.2253, 0.0]\nAWS = [0.808, 0.0]\nrhoWS = [0.132, 0.0]\netaWS = [0.341, 0.0]"
                .into(),
            "toml",
        )?)?;
        let ckm = model.get_parameter("CKM1x3");
        assert_eq!(ckm.parameter_type, ParameterType::Complex);
        assert!(ckm.value.unwrap().re > 0.0);
        assert!(ckm.value.unwrap().im < 0.0);
        if inline {
            let symbol = Atom::from(crate::model::UFOSymbol::from(ckm.name.as_str()));
            let expression = ckm.expression.as_ref().unwrap().clone();
            let mut definition = serde_json::to_value(&model)?;
            for name in ["GC_43", "GC_102"] {
                let expression = model
                    .get_coupling(name)
                    .expression
                    .replace(symbol.to_pattern())
                    .with(expression.to_pattern());
                let coupling = definition["couplings"]
                    .as_array_mut()
                    .unwrap()
                    .iter_mut()
                    .find(|coupling| coupling["name"] == name)
                    .unwrap();
                coupling["expression"] = serde_json::to_value(expression.to_canonical_string())?;
            }
            model = serde_json::from_value(definition)?;
            // Coupling aliases and inline expressions retain their full physical
            // values. No CP decision is inferred from a scalar coefficient.
            model.recompute_dependents()?;
        }
        // The unrelated leptonic coupling retains its ordinary imaginary factor.
        assert_ne!(model.get_coupling("GC_40").value.unwrap().im, 0.0);
        let mut process = ProcessDefinition {
            generation_type: feynkit_generator::GenerationType::CrossSection,
            process: GenerationProcess::new([24_i64], [2_i64, -5]),
            generation_options: GenerationOptions::default()
                .with_loop_count(1, 1)?
                .with_graph_filter(GenerationFilter::VertexAllow(vec![
                    "V_95".into(),
                    "V_125".into(),
                ])),
            ..Default::default()
        };
        let settings = GlobalSettings {
            n_cores: Parallelisation {
                feyngen: 1,
                ..Default::default()
            },
            ..Default::default()
        };
        for symmetrize in [false, true] {
            // Acceptance certifies the opt-in contract, not CP validity of
            // this complex coupling point or a physical optimized rate.
            process.generation_options = process
                .generation_options
                .clone()
                .symmetrize_left_right(symmetrize);
            assert_eq!(process.generate(&model, &settings)?.len(), 1);
        }
        let left = model.get_coupling("GC_43").value.unwrap();
        let right = model.get_coupling("GC_102").value.unwrap();
        assert_ne!(left.re, 0.0);
        assert_eq!(left.re, -right.re);
        assert_eq!(left.im, right.im);
    }
    Ok(())
}

#[test]
fn complex_ckm_updates_preserve_direct_integrand_warm_up() -> color_eyre::Result<()> {
    use crate::{
        integrands::process::ProcessIntegrand, model::InputParamCard, processes::Process,
        settings::RuntimeSettings, utils::F,
    };

    test_initialise()?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
    for (inline, symmetrize) in [(false, false), (false, true), (true, false), (true, true)] {
        let mut model = load_generic_model("sm");
        model.apply_param_card(&InputParamCard::<F<f64>>::from_str(
            "lamWS = [0.2253, 0.0]\nAWS = [0.808, 0.0]\nrhoWS = [0.132, 0.0]\netaWS = [0.0, 0.0]"
                .into(),
            "toml",
        )?)?;
        assert!(model.get_parameter("CKM1x3").value.unwrap().re > 0.0);
        assert_eq!(model.get_parameter("CKM1x3").value.unwrap().im, 0.0);
        if inline {
            let ckm = model.get_parameter("CKM1x3");
            let symbol = Atom::from(crate::model::UFOSymbol::from(ckm.name.as_str()));
            let expression = ckm.expression.as_ref().unwrap().clone();
            let mut definition = serde_json::to_value(&model)?;
            for name in ["GC_43", "GC_102"] {
                let expression = model
                    .get_coupling(name)
                    .expression
                    .replace(symbol.to_pattern())
                    .with(expression.to_pattern());
                let coupling = definition["couplings"]
                    .as_array_mut()
                    .unwrap()
                    .iter_mut()
                    .find(|coupling| coupling["name"] == name)
                    .unwrap();
                coupling["expression"] = serde_json::to_value(expression.to_canonical_string())?;
            }
            model = serde_json::from_value(definition)?;
            model.recompute_dependents()?;
        }

        let definition = ProcessDefinition {
            generation_type: feynkit_generator::GenerationType::CrossSection,
            process: GenerationProcess::new([24_i64], [2_i64, -5]),
            generation_options: GenerationOptions::default()
                .with_loop_count(1, 1)?
                .symmetrize_left_right(symmetrize)
                .with_graph_filter(GenerationFilter::VertexAllow(vec![
                    "V_95".into(),
                    "V_125".into(),
                ])),
            ..Default::default()
        };
        let mut settings = GlobalSettings::default();
        settings.n_cores.feyngen = 1;
        settings.generation.evaluator.compile = false;
        settings.generation.uv.subtract_uv = false;
        settings.generation.uv.generate_integrated = false;
        settings.generation.threshold_subtraction.enable_thresholds = false;
        let mass = model.get_parameter("MW").value.unwrap().re;
        // The rest frame is the natural decay frame; its longitudinal state is finite.
        let runtime: RuntimeSettings = toml::from_str(&format!(
            r#"
            [kinematics]
            e_cm = {mass}
            [kinematics.externals]
            type = "constant"
            [kinematics.externals.data]
            momenta = [[{mass}, 0.0, 0.0, 0.0]]
            helicities = [0]
        "#
        ))?;
        let graphs = definition.generate(&model, &settings)?;
        assert_eq!(graphs.len(), 1);
        let mut process = Process::from_graph_list(
            "runtime_complex_ckm".into(),
            "default".into(),
            graphs,
            GenerationType::CrossSection,
            Some(definition),
            None,
            &model,
        )?;
        process.preprocess(&model, &settings, &(&runtime).into(), &pool)?;
        process.generate_integrands(&model, &settings, (&runtime).into(), &pool)?;
        let integrand = process.get_integrand_mut("default")?;
        let ProcessIntegrand::CrossSection(cross_section) = &*integrand else {
            unreachable!();
        };
        assert_eq!(cross_section.data.symmetrize_left_right_states, symmetrize);
        integrand.warm_up(&model)?;
        model.apply_param_card(&InputParamCard::<F<f64>>::from_str(
            "etaWS = [0.341, 0.0]".into(),
            "toml",
        )?)?;
        assert!(model.get_parameter("CKM1x3").value.unwrap().im < 0.0);
        assert_ne!(model.get_coupling("GC_43").value.unwrap().re, 0.0);
        // Updating the coupling point must not turn the user's CP assumption
        // into an automatic rejection. No CP-violating rate is certified here.
        integrand.warm_up(&model)?;
    }
    Ok(())
}

#[test]
fn nlo_dis_cross_section_counts_match_main() -> Result<()> {
    // Main's `sm_dis`: the SM restricted to massless light quarks, whose
    // electron-photon vertex alone carries the extra `LT` coupling order.
    let mut model = load_generic_model("sm");
    let mut restriction = InputParamCard::<F<f64>>::from_str(
        include_str!("../../../../assets/models/json/sm/restrict_light_quarks_massless.json")
            .to_owned(),
        "json",
    )?;
    model.simplify(&mut restriction)?;
    let mut definition = serde_json::to_value(&model)?;
    definition["name"] = serde_json::json!("sm_dis");
    definition["orders"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({"name": "LT", "expansion_order": 99, "hierarchy": 3}));
    let couplings = definition["couplings"].as_array_mut().unwrap();
    let mut lepton_tag = couplings
        .iter()
        .find(|coupling| coupling["name"] == "GC_3")
        .unwrap()
        .clone();
    lepton_tag["name"] = serde_json::json!("DIS");
    lepton_tag["orders"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!(["LT", 1]));
    couplings.push(lepton_tag);
    definition["vertex_rules"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|vertex| vertex["name"] == "V_98")
        .unwrap()["couplings"] = serde_json::json!([["DIS"]]);
    let mut model: Model = serde_json::from_value(definition)?;
    model.recompute_dependents()?;

    // Ported from main's `nlo_fs_dis`: e- p > e- X at NLO in QCD for every
    // one- and two-parton initial state, with `X` a single parton.
    let vetoed = [
        23_i64, 24, 9000001, 9000002, 9000003, 9000004, 9000005, -9000005, 12, 14, 16, 2, 4, 6, 3,
        5, 25, 250, 251, 13, 15,
    ]
    .map(ParticleSelector::Pdg)
    .to_vec();
    let partons = [-1_i64, 1, 21];
    let mut initial_states = partons.map(|parton| vec![parton]).to_vec();
    for (position, first) in partons.iter().enumerate() {
        for second in &partons[position..] {
            initial_states.push(vec![*first, *second]);
        }
    }
    let settings = GlobalSettings {
        n_cores: Parallelisation {
            feyngen: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut counts = BTreeMap::new();
    for initial in initial_states {
        let loop_count = 3 - initial.len();
        // Main grouped with `OnlyDetectZeroes`; FeynKit removes the color
        // zeros through its explicit filter and does not group.
        let mut options = GenerationOptions::default()
            .with_loop_count(loop_count, loop_count)?
            .symmetrize_initial(true)
            .symmetrize_final(true)
            .filter_zero_color(true);
        for filter in [
            GenerationFilter::Sewn(SewnFilterOptions {
                filter_tadpoles: true,
            }),
            GenerationFilter::ParticleVeto(vetoed.clone()),
            GenerationFilter::Tadpoles(TadpoleFilterOptions::default()),
            GenerationFilter::ZeroSnails(SnailFilterOptions::default()),
            GenerationFilter::PerturbativeOrders(BTreeMap::from([("QCD".to_owned(), 1)])),
            GenerationFilter::CouplingOrders(BTreeMap::from([("QCD".to_owned(), (2, Some(2)))])),
            GenerationFilter::BlobRange(1..=1),
            GenerationFilter::SpectatorRange(0..=1),
        ] {
            options = options.with_graph_filter(filter);
        }
        // Main's cut checks read only the amplitude coupling orders.
        options =
            options.with_cut_amplitude_filter(GenerationFilter::CouplingOrders(BTreeMap::from([
                ("QED".to_owned(), (2, Some(2))),
                ("LT".to_owned(), (1, Some(1))),
            ])));
        let process = ProcessDefinition {
            generation_type: GenerationType::CrossSection,
            process: GenerationProcess::new(initial.iter().copied().chain([11]), [11_i64, -1])
                .with_final_state_alternatives(partons.map(|parton| [11, parton]))?,
            generation_options: options,
            ..Default::default()
        };
        let name = initial
            .iter()
            .chain(&[11])
            .map(|pdg| Ok(model.particle_by_pdg(*pdg)?.name.clone()))
            .collect::<Result<Vec<_>>>()?
            .join(",");
        counts.insert(name, process.generate(&model, &settings)?.len());
    }
    assert_eq!(
        counts,
        BTreeMap::from(
            [
                ("d,e-", 6),
                ("d~,e-", 6),
                ("g,e-", 4),
                ("d,d,e-", 12),
                ("d~,d,e-", 12),
                ("d~,d~,e-", 12),
                ("d,g,e-", 11),
                ("d~,g,e-", 11),
                ("g,g,e-", 0),
            ]
            .map(|(name, count)| (name.to_owned(), count))
        )
    );
    Ok(())
}

/// A massive vector at rest has no momentum direction for its polarization to
/// follow; stability rotations must rotate its quantization axis with the event.
#[test]
fn rest_frame_polarizations_rotate_with_stability_rotations() -> color_eyre::Result<()> {
    test_initialise()?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
    let mut model = load_generic_model("sm");
    model.apply_param_card(&InputParamCard::<F<f64>>::from_str(
        "lamWS = [0.2253, 0.0]\nAWS = [0.808, 0.0]\nrhoWS = [0.132, 0.0]\netaWS = [0.341, 0.0]"
            .into(),
        "toml",
    )?)?;
    let definition = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        process: GenerationProcess::new([24_i64], [2_i64, -5]),
        generation_options: GenerationOptions::default()
            .with_loop_count(1, 1)?
            .symmetrize_left_right(false)
            .with_graph_filter(GenerationFilter::VertexAllow(vec![
                "V_95".into(),
                "V_125".into(),
            ])),
        ..Default::default()
    };
    let mut settings = GlobalSettings::default();
    settings.n_cores.feyngen = 1;
    settings.generation.evaluator.compile = false;
    settings.generation.uv.subtract_uv = false;
    settings.generation.uv.generate_integrated = false;
    settings.generation.threshold_subtraction.enable_thresholds = false;
    let mass = model.get_parameter("MW").value.unwrap().re;
    let runtime: RuntimeSettings = toml::from_str(&format!(
        r#"
        [kinematics]
        e_cm = {mass}
        [kinematics.externals]
        type = "constant"
        [kinematics.externals.data]
        momenta = [[{mass}, 0.0, 0.0, 0.0]]
        helicities = [0]
        [subtraction]
        disable_threshold_subtraction = true
        [stability]
        rotation_axis = []
    "#
    ))?;
    let graphs = definition.generate(&model, &settings)?;
    assert_eq!(graphs.len(), 1);
    let mut process = Process::from_graph_list(
        "rest_rotation".into(),
        "default".into(),
        graphs,
        GenerationType::CrossSection,
        Some(definition),
        None,
        &model,
    )?;
    process.preprocess(&model, &settings, &(&runtime).into(), &pool)?;
    process.generate_integrands(&model, &settings, (&runtime).into(), &pool)?;
    let pristine = process.get_integrand_mut("default")?.clone();
    for (name, axis) in [
        ("x", RotationSetting::Pi2X {}),
        ("y", RotationSetting::Pi2Y {}),
        ("z", RotationSetting::Pi2Z {}),
        (
            "euler",
            RotationSetting::EulerAngles {
                alpha: 0.3,
                beta: -0.2,
                gamma: 0.7,
            },
        ),
    ] {
        let mut values = Vec::new();
        for cache in [false, true] {
            // A fresh integrand per case, as the stack owner's probe does.
            let mut integrand = pristine.clone();
            let runtime = integrand.get_mut_settings();
            runtime.general.enable_cache = cache;
            runtime.stability.rotation_axis = vec![axis];
            runtime.stability.levels = vec![StabilityLevelSetting::default_double()];
            runtime.stability.check_on_norm = true;
            integrand.warm_up(&model)?;
            let input = MomentumSpaceEvaluationInput {
                loop_momenta: vec![ThreeMomentum {
                    px: F(11.0),
                    py: F(23.0),
                    pz: F(37.0),
                }],
                integrator_weight: F(1.0),
                graph_id: Some(0),
                group_id: None,
                orientation: None,
                channel_id: None,
            };
            let evaluation = integrand.evaluate_momentum_configuration(&model, &input, false)?;
            let stability = &evaluation.evaluation_metadata.stability_results;
            assert!(
                !stability.is_empty()
                    && stability
                        .iter()
                        .all(|result| matches!(result.status, StabilityStatus::Stable(_))),
                "{name} rotation with cache {cache}: {stability:?}"
            );
            let value = evaluation.integrand_result;
            assert!(value.re.0.is_finite() && value.re.0 != 0.0, "{value:?}");
            values.push(value);
        }
        let [off, on] = [&values[0], &values[1]];
        assert!(
            (off.re.0 - on.re.0).abs() <= 1e-12 * off.re.0.abs()
                && (off.im.0 - on.im.0).abs() <= 1e-12 * off.re.0.abs(),
            "{name} rotation: cache off {off:?}, cache on {on:?}"
        );
    }
    Ok(())
}
