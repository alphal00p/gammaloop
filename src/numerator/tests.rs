use spenso::{complex::Complex, iterators::IteratableTensor, structure::HasStructure};
use std::path::{Path, PathBuf};
use symbolica::{atom::Atom, domains::rational::Rational};

use crate::{
    cross_section::Amplitude,
    gammaloop_integrand::DefaultSample,
    graph::Graph,
    model::Model,
    momentum::{Dep, ExternalMomenta, Helicity},
    numerator::{ContractionSettings, ExtraInfo, GlobalPrefactor},
    tests_from_pytest::{load_amplitude_output, sample_generator, test_export_settings},
    utils::{ApproxEq, F},
    Externals, RotationSetting, Settings,
};

use super::{
    Evaluate, EvaluatorOptions, Numerator, NumeratorCompileOptions, NumeratorEvaluatorOptions,
    UnInit,
};

#[ignore]
#[test]
fn hhgghh() {
    let _ = env_logger::builder().is_test(true).try_init();
    let (model, amplitude, path) = load_amplitude_output(
        &("TEST_AMPLITUDE_".to_string() + "physical_1L_2A_final_4H_top_internal" + "/GL_OUTPUT"),
        true,
    );

    validate_gamma(amplitude.amplitude_graphs[0].graph.clone(), &model, path);
}

#[ignore]
#[test]
fn hairy_glue_box() {
    let _ = env_logger::builder().is_test(true).try_init();
    let (_, amplitude, _) = load_amplitude_output(
        &("TEST_AMPLITUDE_".to_string() + "hairy_glue_box" + "/GL_OUTPUT"),
        true,
    );

    let graph = amplitude.amplitude_graphs[0].graph.clone();
    for (i, s) in graph.bare_graph.external_slots().iter().enumerate() {
        println!("{i}:{}", s);
    }

    let color = graph.derived_data.unwrap().numerator.from_graph(&graph.bare_graph,Some(&GlobalPrefactor{color:Atom::parse("f(aind(coad(8,1),coad(8,2),coad(8,3)))*f(aind(coad(8,4),coad(8,5),coad(8,6)))*id(aind(coad(8,7),coad(8,0)))").unwrap(),colorless:Atom::new_num(1)})).color_simplify().state.color.to_dense().map_data(|a|a.to_string());

    insta::assert_ron_snapshot!(color);
}

#[test]
fn trees() {
    let _ = env_logger::builder().is_test(true).try_init();
    let tree_name = "th_th";
    let amp_name = "tree_amplitude_1_th_th";
    let file_path = PathBuf::new()
        .join("./src/test_resources/trees")
        .join(tree_name)
        .join("GL_OUTPUT");

    let model = Model::from_file(String::from(
        Path::new("./src/test_resources")
            .join("gammaloop_models/sm.yaml")
            .to_str()
            .unwrap(),
    ))
    .unwrap();

    let amplitude = Amplitude::from_file(
        &model,
        String::from(
            file_path
                .join(format!("sources/amplitudes/{}/amplitude.yaml", amp_name))
                .to_str()
                .unwrap(),
        ),
    )
    .unwrap();

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    let lmb = &graph.bare_graph.loop_momentum_basis;
    let sample: crate::gammaloop_integrand::DefaultSample<f64> = sample_generator(
        3,
        &graph.bare_graph,
        Some(vec![
            Helicity::Plus,
            Helicity::Zero,
            Helicity::Plus,
            Helicity::Zero,
        ]),
    );
    let emr = graph.bare_graph.cff_emr_from_lmb(&sample.sample, lmb);

    let three_emr = lmb.spatial_emr(&sample.sample);

    let onshell_energies = graph
        .bare_graph
        .compute_onshell_energies(sample.loop_moms(), sample.external_moms());

    for ((e, q), p) in onshell_energies.iter().zip(emr).zip(three_emr) {
        F::approx_eq_res(e, &q.temporal.value, &F(1e-10)).unwrap();

        for (a, b) in q.spatial.into_iter().zip(p) {
            F::approx_eq_res(&a, &b, &F(1e-10)).unwrap();
        }
    }

    let export_path = PathBuf::new()
        .join("./src/test_resources/trees")
        .join(tree_name);

    let contraction_settings = ContractionSettings::<Rational>::Normal;

    let mut export_settings = test_export_settings();
    for (i, s) in graph.bare_graph.external_slots().iter().enumerate() {
        println!("{i}:{}", s);
    }
    export_settings.numerator_settings.global_prefactor = Some(GlobalPrefactor {
        color: Atom::parse("id(aind(cof(3,2),coaf(3,0)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    });

    graph.generate_cff();
    let mut graph =
        graph.process_numerator(&model, contraction_settings, export_path, &export_settings);

    let external_mom = Externals::Constant {
        momenta: vec![
            ExternalMomenta::Independent([F(592.625), F(0.), F(0.), F(566.8116)]), // 1
            ExternalMomenta::Independent([F(579.375), F(0.), F(0.), F(-566.8116)]), // 2
            ExternalMomenta::Independent([F(592.625), F(125.7463), F(504.2705), F(-226.2178)]), // 3
            ExternalMomenta::Dependent(Dep::Dep),                                  // 4
        ],
        helicities: vec![
            Helicity::Minus,
            Helicity::Zero,
            Helicity::Plus,
            Helicity::Zero,
        ],
    };

    let external_signature = graph.bare_graph.external_in_or_out_signature();

    let sample: DefaultSample<f64> = DefaultSample::new(
        vec![],
        &external_mom,
        F(1.),
        &external_mom
            .generate_polarizations(&graph.bare_graph.external_particles(), &external_signature),
        &external_signature,
    );
    let settings = Settings::default();

    let val = graph
        .evaluate_fourd_expr(
            &[],
            sample.external_moms(),
            sample.polarizations(),
            &settings,
        )
        .scalar()
        .unwrap();

    let energy_product = graph
        .bare_graph
        .compute_energy_product(sample.loop_moms(), sample.external_moms());

    let valcff = graph.evaluate_cff_expression(&sample, &settings) / energy_product;
    println!("4d: {}", val);
    println!("CFF: {}", valcff);
}

#[test]
fn tree_ta_ta_1() {
    let (model, amplitude, path) = load_tree("ta_ta", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    graph.bare_graph.verify_external_edge_order().unwrap();

    for (i, s) in graph.bare_graph.external_slots().iter().enumerate() {
        println!("{i}:{}", s);
    }

    let mut test_export_settings = test_export_settings();
    test_export_settings.numerator_settings.global_prefactor = Some(GlobalPrefactor {
        color: Atom::parse("id(aind(coaf(3,0),cof(3,2)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    });

    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::<Rational>::Normal,
        path,
        &test_export_settings,
    );
    let sample: DefaultSample<f64> = sample_generator(3, &graph.bare_graph, None);

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

    let val = graph
        .evaluate_fourd_expr(
            &[],
            sample.external_moms(),
            sample.polarizations(),
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    let expected_val = Complex::new(F(0.0002488992442101011), F(0.00003005190084792535));
    let expected_cff = Complex::new(F(-0.000248899244210101), F(-0.000030051900847925438));

    val.approx_eq_res(&expected_val, &F(1e-10)).unwrap();
    cff_val.approx_eq_res(&expected_cff, &F(1e-10)).unwrap();

    let mut externals = Externals::Constant {
        momenta: vec![
            ExternalMomenta::Independent([
                F(0.5149645E+03),
                F(0.0000000E+00),
                F(0.0000000E+00),
                F(0.4850355E+03),
            ]),
            ExternalMomenta::Dependent(crate::momentum::Dep::Dep),
            ExternalMomenta::Independent([
                F(0.5149645E+03),
                F(0.1076044E+03),
                F(0.4315174E+03),
                F(-0.1935805E+03),
            ]),
            ExternalMomenta::Independent([
                F(0.4850355E+03),
                F(-0.1076044E+03),
                F(-0.4315174E+03),
                F(0.1935805E+03),
            ]),
        ],
        helicities: vec![
            Helicity::Plus,
            Helicity::Plus,
            Helicity::Plus,
            Helicity::Plus,
        ],
    };

    externals
        .set_dependent_at_end(&graph.bare_graph.external_in_or_out_signature())
        .unwrap();

    let external_signature = graph.bare_graph.external_in_or_out_signature();

    let sample: DefaultSample<f64> = DefaultSample::new(
        vec![],
        &externals,
        F(1.),
        &externals
            .generate_polarizations(&graph.bare_graph.external_particles(), &external_signature),
        &external_signature,
    );

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

    let val = graph
        .evaluate_fourd_expr(
            &[],
            sample.external_moms(),
            sample.polarizations(),
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    let rotated_sample = sample.get_rotated_sample(&RotationSetting::Pi2X.rotation_method().into());

    // println!("rotated sample: {}", rotated_sample);
    // println!("sample: {}", sample);

    let cff_val_rot = graph.evaluate_cff_expression(&rotated_sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(rotated_sample.loop_moms(), rotated_sample.external_moms());

    let expected_val = Complex::new(F(0.), F(0.));
    let expected_cff = Complex::new(F(-0.0), F(-0.0));
    let expected_cff_rot = Complex::new(F(-0.0), F(-0.0));
    val.approx_eq_res(&expected_val, &F(1e-10)).unwrap();
    cff_val.approx_eq_res(&expected_cff, &F(1e-10)).unwrap();
    cff_val_rot
        .approx_eq_res(&expected_cff_rot, &F(1e-10))
        .unwrap();

    println!(
        "{}",
        Numerator::default()
            .from_graph(
                &graph.bare_graph,
                test_export_settings
                    .numerator_settings
                    .global_prefactor
                    .as_ref()
            )
            .color_simplify()
            .gamma_simplify()
            .export()
    );
}

pub fn validate_gamma(g: Graph<UnInit>, model: &Model, path: PathBuf) {
    let num = g.derived_data.as_ref().unwrap().numerator.clone();

    let num = num.from_graph(&g.bare_graph, None);

    let path_gamma = path.join("gamma");

    let mut export_settings = test_export_settings();
    export_settings.numerator_settings.eval_settings =
        NumeratorEvaluatorOptions::Joint(EvaluatorOptions {
            cpe_rounds: Some(2),
            compile_options: NumeratorCompileOptions::NotCompiled,
        });

    let mut num_nogamma = num
        .clone()
        .color_simplify()
        // .gamma_symplify()
        .parse()
        .contract(ContractionSettings::<Rational>::Normal)
        .unwrap()
        .generate_evaluators(
            model,
            &g.bare_graph,
            &ExtraInfo {
                path,
                orientations: vec![vec![true; g.bare_graph.edges.len()]],
            },
            &export_settings,
        );
    let mut num_gamma = num
        .clone()
        .color_simplify()
        .gamma_simplify()
        .parse()
        .contract(ContractionSettings::<Rational>::Normal)
        .unwrap()
        .generate_evaluators(
            model,
            &g.bare_graph,
            &ExtraInfo {
                path: path_gamma,
                orientations: vec![vec![true; g.bare_graph.edges.len()]],
            },
            &export_settings,
        );

    for i in 0..10 {
        let s: DefaultSample<f64> = sample_generator(i, &g.bare_graph, None);

        let emr = g
            .bare_graph
            .cff_emr_from_lmb(&s.sample, &g.bare_graph.loop_momentum_basis);

        let polarizations = s.polarizations();

        let val = num_nogamma.evaluate_single(&emr, polarizations, None, &Settings::default());

        let valg = num_gamma.evaluate_single(&emr, polarizations, None, &Settings::default());

        assert!(
            Complex::approx_eq_iterator(
                valg.iter_flat().map(|(_, o)| o),
                val.iter_flat().map(|(_, o)| o),
                &F(0.0001),
            ),
            "{}: {}!={}",
            i,
            val,
            valg,
        );
    }
}

#[test]
fn tree_h_ttxaah_1() {
    let (model, amplitude, path) = load_tree("h_ttxaah", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    // for (i, s) in graph.bare_graph.external_slots().iter().enumerate() {
    //     println!("{i}:{}", s);
    // }
    let mut test_export_settings = test_export_settings();
    test_export_settings.numerator_settings.global_prefactor = Some(GlobalPrefactor {
        color: Atom::parse("id(aind(cof(3,1),coaf(3,2)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    });

    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::<Rational>::Normal,
        path,
        &test_export_settings,
    );
    let sample: DefaultSample<f64> = sample_generator(3, &graph.bare_graph, None);

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

    let val = graph
        .evaluate_fourd_expr(
            &[],
            sample.external_moms(),
            sample.polarizations(),
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    let expected_val = Complex::new(
        F(0.0000000018402069973971576),
        F(0.000000006517984974022169),
    );
    let expected_cff = Complex::new(
        F(-0.0000000018402069973971528),
        F(-0.000000006517984974022173),
    );
    val.approx_eq_res(&expected_val, &F(1e-10)).unwrap();
    cff_val.approx_eq_res(&expected_cff, &F(1e-10)).unwrap();
}

#[test]
fn tree_hh_ttxaa_1() {
    let (model, amplitude, path) = load_tree("hh_ttxaa", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    for (i, s) in graph.bare_graph.external_slots().iter().enumerate() {
        println!("{i}:{}", s);
    }
    let mut test_export_settings = test_export_settings();
    test_export_settings.numerator_settings.global_prefactor = Some(GlobalPrefactor {
        color: Atom::parse("id(aind(cof(3,2),coaf(3,3)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    });

    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::<Rational>::Normal,
        path,
        &test_export_settings,
    );
    let sample: DefaultSample<f64> = sample_generator(3, &graph.bare_graph, None);

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

    let val = graph
        .evaluate_fourd_expr(
            &[],
            sample.external_moms(),
            sample.polarizations(),
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    let expected_val = Complex::new(
        F(-0.000000007061511384861008),
        F(0.0000000024246582518362876),
    );
    let expected_cff = Complex::new(
        F(0.000000007061511384861011),
        F(-0.0000000024246582518362864),
    );
    val.approx_eq_res(&expected_val, &F(1e-10)).unwrap();
    cff_val.approx_eq_res(&expected_cff, &F(1e-10)).unwrap();

    let externals = Externals::Constant {
        momenta: vec![
            ExternalMomenta::Independent(
                [0.586E+03, 0.000E+00, 0.000E+00, 5.735_817_291_371_824E2].map(F),
            ),
            ExternalMomenta::Independent(
                [
                    5.86E2,
                    0.00000000000000000E+00,
                    0.00000000000000000E+00,
                    -5.735_817_291_371_824E2,
                ]
                .map(F),
            ),
            ExternalMomenta::Independent(
                [
                    1.953_887_163_589_759_6E2,
                    -2.266_618_827_913_845_3E1,
                    4.110_590_302_435_584E1,
                    -7.774_509_068_652_03E1,
                ]
                .map(F),
            ),
            ExternalMomenta::Independent(
                [
                    3.785_715_624_411_541E2,
                    -1.065_068_477_511_299_8E2,
                    -3.096_594_386_148_455_6E2,
                    7.845_222_334_639_679E1,
                ]
                .map(F),
            ),
            ExternalMomenta::Independent(
                [
                    1.562_565_490_076_973E2,
                    -1.085_901_723_312_049_5E2,
                    -1.002_097_685_717_689E2,
                    5.081_619_686_346_704E1,
                ]
                .map(F),
            ),
            ExternalMomenta::Dependent(crate::momentum::Dep::Dep),
        ],
        helicities: vec![
            Helicity::Zero,
            Helicity::Zero,
            Helicity::Plus,
            Helicity::Plus,
            Helicity::Plus,
            Helicity::Plus,
        ],
    };

    let external_signature = graph.bare_graph.external_in_or_out_signature();

    let sample: DefaultSample<f64> = DefaultSample::new(
        vec![],
        &externals,
        F(1.),
        &externals
            .generate_polarizations(&graph.bare_graph.external_particles(), &external_signature),
        &external_signature,
    );

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

    let val = graph
        .evaluate_fourd_expr(
            &[],
            sample.external_moms(),
            sample.polarizations(),
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    let expected_val = Complex::new(
        F(0.000000029250807246513418),
        F(-0.00000000017560186291318082),
    );
    let expected_cff = Complex::new(
        F(-0.00000002925080724651342),
        F(0.00000000017560186291318444),
    );
    val.approx_eq_res(&expected_val, &F(1e-10)).unwrap();
    cff_val.approx_eq_res(&expected_cff, &F(1e-10)).unwrap();
}

fn load_tree(tree_name: &str, amp_num: usize) -> (Model, Amplitude<UnInit>, PathBuf) {
    let _ = env_logger::builder().is_test(true).try_init();
    let file_path = PathBuf::new()
        .join("./src/test_resources/trees")
        .join(tree_name)
        .join("GL_OUTPUT");

    let model = Model::from_file(String::from(
        Path::new("./src/test_resources")
            .join("gammaloop_models/sm.yaml")
            .to_str()
            .unwrap(),
    ))
    .unwrap();

    let amplitude = Amplitude::from_file(
        &model,
        String::from(
            file_path
                .join(format!(
                    "sources/amplitudes/tree_amplitude_{}_{}/amplitude.yaml",
                    amp_num, tree_name
                ))
                .to_str()
                .unwrap(),
        ),
    )
    .unwrap();

    let export_path = file_path.join(format!(
        "sources/amplitudes/tree_amplitude_{}_{}",
        amp_num, tree_name
    ));

    (model, amplitude, export_path)
}

#[test]
fn tree_h_ttxaah_0() {
    let _ = env_logger::builder().is_test(true).try_init();
    let expr = Atom::parse("-8/3*𝑖*ee^2*vev*lam*yt*(MT*id(aind(bis(4,3),bis(4,4)))+Q(6,aind(lord(4,5)))*γ(aind(loru(4,5),bis(4,3),bis(4,4))))*(MT*id(aind(bis(4,5),bis(4,6)))+Q(7,aind(lord(4,11)))*γ(aind(loru(4,11),bis(4,5),bis(4,6))))*(ProjM(aind(bis(4,7),bis(4,6)))+ProjP(aind(bis(4,7),bis(4,6))))*sqrt(2)^-1*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,6),cof(3,7)))*id(aind(coaf(3,7),cof(3,8)))*id(aind(coaf(3,8),cof(3,9)))*id(aind(coaf(3,9),cof(3,10)))*γ(aind(lord(4,2),bis(4,3),bis(4,2)))*γ(aind(lord(4,3),bis(4,5),bis(4,4)))*ubar(1,aind(bis(4,2)))*v(2,aind(bis(4,7)))*ϵbar(3,aind(loru(4,2)))*ϵbar(4,aind(loru(4,3)))").unwrap();

    let prefactor = GlobalPrefactor {
        color: Atom::parse("id(aind(cof(3,5),coaf(3,10)))").unwrap(),
        colorless: Atom::new_num(1),
    };

    Numerator::default()
        .from_global(expr, Some(&prefactor))
        .color_simplify()
        // .color_project()
        // .gamma_symplify()
        .parse()
        .contract(ContractionSettings::<Rational>::Normal)
        .unwrap();
    // println!("{}", new.export());
}

#[test]
fn color() {
    insta::assert_snapshot!("Single color string",Numerator::default().from_global(Atom::parse("f(aind(coad(8,1),coad(8,11),coad(8,21)))*f(aind(coad(8,21),coad(8,2),coad(8,12)))*f(aind(coad(8,3),coad(8,12),coad(8,22)))*f(aind(coad(8,22),coad(8,4),coad(8,13)))*f(aind(coad(8,5),coad(8,13),coad(8,23)))*f(aind(coad(8,23),coad(8,6),coad(8,14)))*f(aind(coad(8,7),coad(8,14),coad(8,24)))*f(aind(coad(8,24),coad(8,8),coad(8,11)))*f(aind(coad(8,1),coad(8,2),coad(8,3)))*f(aind(coad(8,4),coad(8,5),coad(8,6)))*id(aind(coad(8,7),coad(8,8)))").unwrap(), None).color_simplify().export());
}
