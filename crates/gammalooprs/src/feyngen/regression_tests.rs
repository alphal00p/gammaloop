//! Regressions retained when the generator moved to FeynKit.
use crate::{
    initialisation::test_initialise,
    model::ModelGammaLoopExt,
    processes::ProcessDefinition,
    settings::{GlobalSettings, global::Parallelisation},
    utils::load_generic_model,
};
use feynkit_generator::{
    GenerationFilter, GenerationOptions, GenerationType, Process as GenerationProcess,
};
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
