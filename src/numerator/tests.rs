use brotli::CompressorWriter;
use insta::assert_snapshot;
use spenso::{
    complex::Complex,
    data::{DenseTensor, StorageTensor},
    iterators::IteratableTensor,
    parametric::{atomcore::TensorAtomMaps, ParamTensor},
    shadowing::ETS,
    structure::{
        representation::{BaseRepName, Bispinor, Minkowski},
        slot::IsAbstractSlot,
        HasStructure,
    },
    upgrading_arithmetic::FallibleSub,
};
use std::{
    fs::File,
    io::BufWriter,
    path::{Path, PathBuf},
};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
    function,
};

use crate::{
    cross_section::Amplitude,
    feyngen::diagram_generator::{EdgeColor, NodeColorWithVertexRule},
    gammaloop_integrand::DefaultSample,
    graph::{BareGraph, Graph},
    model::Model,
    momentum::{Dep, ExternalMomenta, Helicity},
    numerator::{ufo::UFO, ContractionSettings, ExtraInfo, GlobalPrefactor, Network},
    tests_from_pytest::{
        load_amplitude_output, load_generic_model, sample_generator, test_export_settings,
    },
    utils::{ApproxEq, F},
    Externals, RotationSetting, Settings,
};

use super::{
    Evaluate, EvaluatorOptions, GammaSimplified, Numerator, NumeratorCompileOptions,
    NumeratorEvaluatorOptions, UnInit,
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

    let color = graph.derived_data.unwrap().numerator.from_graph(&graph.bare_graph,&GlobalPrefactor{color:Atom::parse("f(aind(coad(8,1),coad(8,2),coad(8,3)))*f(aind(coad(8,4),coad(8,5),coad(8,6)))*id(aind(coad(8,7),coad(8,0)))").unwrap(),colorless:Atom::new_num(1)}).color_simplify().state.color.to_bare_dense().map_data(|a|a.to_string());

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
    export_settings.numerator_settings.global_prefactor = GlobalPrefactor {
        color: Atom::parse("id(cof(3,2),dind(cof(3,0)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    };

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

fn compare_poly_to_direct(graph: &BareGraph, prefactor: &GlobalPrefactor) -> bool {
    let color_simplified = Numerator::default()
        .from_graph(graph, prefactor)
        .color_simplify();

    let poly = color_simplified
        .clone()
        .parse_poly(graph)
        .contract()
        .unwrap()
        .to_contracted()
        .state
        .tensor;

    let direct = color_simplified
        .parse()
        .contract(ContractionSettings::<Rational>::Normal)
        .unwrap()
        .state
        .tensor;

    let zero = ParamTensor::from(DenseTensor::fill(
        poly.structure().clone(),
        Atom::new_num(0),
    ));

    zero == poly.sub_fallible(&direct).unwrap().expand()
}

#[allow(dead_code)]
pub fn save_expr(graph: &BareGraph, prefactor: &GlobalPrefactor, name: &str) {
    let color_simplified = Numerator::default()
        .from_graph(graph, prefactor)
        .color_simplify();
    let direct = color_simplified
        .parse()
        .contract::<Rational>(ContractionSettings::<Rational>::Normal)
        .unwrap()
        .state
        .tensor
        .scalar()
        .unwrap();
    let f = File::create(name).unwrap();
    let mut writer = CompressorWriter::new(BufWriter::new(f), 4096, 3, 22);

    direct.as_view().export(&mut writer).unwrap();
}

#[test]

fn t_ta() {
    let (_model, amplitude, _path) = load_tree("t_ta", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    graph.bare_graph.verify_external_edge_order().unwrap();

    let global_prefactor = GlobalPrefactor {
        color: Atom::parse("id(dind(cof(3,0)),cof(3,1))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    };
    assert!(
        compare_poly_to_direct(&graph.bare_graph, &global_prefactor),
        "Poly and direct are not equal"
    );
    // save_expr(&graph.bare_graph, Some(&global_prefactor), "t_ta.dat");
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

    let global_prefactor = GlobalPrefactor {
        color: Atom::parse("id(dind(cof(3,0)),cof(3,2))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    };
    assert!(
        compare_poly_to_direct(&graph.bare_graph, &global_prefactor),
        "Poly and direct are not equal"
    );

    // save_expr(&graph.bare_graph, &global_prefactor), "ta_ta.dat");

    let mut test_export_settings = test_export_settings();
    test_export_settings.numerator_settings.global_prefactor = global_prefactor;

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
                &test_export_settings.numerator_settings.global_prefactor
            )
            .color_simplify()
            .gamma_simplify()
            .export()
    );
}

pub fn validate_gamma(g: Graph<UnInit>, model: &Model, path: PathBuf) {
    let num = g.derived_data.as_ref().unwrap().numerator.clone();

    let num = num.from_graph(&g.bare_graph, &GlobalPrefactor::default());

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
        .contract::<Rational>(ContractionSettings::<Rational>::Normal)
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
        .contract::<Rational>(ContractionSettings::<Rational>::Normal)
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
    test_export_settings.numerator_settings.global_prefactor = GlobalPrefactor {
        color: Atom::parse("id(cof(3,1),dind(cof(3,2)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    };

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
    test_export_settings.numerator_settings.global_prefactor = GlobalPrefactor {
        color: Atom::parse("id(cof(3,2),dind(cof(3,3)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    };

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

    let num = Numerator::default();
    let expr = Atom::parse("-8/3*𝑖*ee^2*vev*lam*yt*(MT*id(bis(4,3),bis(4,4))+Q(6,mink(4,5))*γ(mink(4,5),bis(4,3),bis(4,4)))*(MT*id(bis(4,5),bis(4,6))+Q(7,mink(4,11))*γ(mink(4,11),bis(4,5),bis(4,6)))*(ProjM(bis(4,7),bis(4,6))+ProjP(bis(4,7),bis(4,6)))*sqrt(2)^-1*id(dind(cof(3,5)),cof(3,6))*id(dind(cof(3,6)),cof(3,7))*id(dind(cof(3,7)),cof(3,8))*id(dind(cof(3,8)),cof(3,9))*id(dind(cof(3,9)),cof(3,10))*γ(mink(4,2),bis(4,3),bis(4,2))*γ(mink(4,3),bis(4,5),bis(4,4))*ubar(1,bis(4,2))*v(2,bis(4,7))*ϵbar(3,mink(4,2))*ϵbar(4,mink(4,3))").unwrap();

    let prefactor = GlobalPrefactor {
        color: Atom::parse("id(cof(3,5),dind(cof(3,10)))").unwrap(),
        colorless: Atom::new_num(1),
    };

    num.from_global(expr, &prefactor)
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
    insta::assert_snapshot!("Single color string",Numerator::default().from_global(Atom::parse("f(coad(8,1),coad(8,11),coad(8,21))*f(coad(8,21),coad(8,2),coad(8,12))*f(coad(8,3),coad(8,12),coad(8,22))*f(coad(8,22),coad(8,4),coad(8,13))*f(coad(8,5),coad(8,13),coad(8,23))*f(coad(8,23),coad(8,6),coad(8,14))*f(coad(8,7),coad(8,14),coad(8,24))*f(coad(8,24),coad(8,8),coad(8,11))*f(coad(8,1),coad(8,2),coad(8,3))*f(coad(8,4),coad(8,5),coad(8,6))*id(coad(8,7),coad(8,8))").unwrap(), &GlobalPrefactor::default()).color_simplify().export());
}

#[test]
fn prefactor() {
    let mut test_export_settings = test_export_settings();
    test_export_settings.numerator_settings.global_prefactor = GlobalPrefactor {
        color: Atom::parse("id(cof(3,2),dind(cof(3,3)))/Nc").unwrap(),
        colorless: Atom::new_num(1),
    };

    println!(
        "{}",
        serde_yaml::to_string(&test_export_settings.numerator_settings).unwrap()
    );
}

#[test]
fn gamma_simplify_one() {
    fn gamma(mu: usize, i: usize, j: usize) -> Atom {
        let mink = Minkowski::rep(4);
        let bis = Bispinor::rep(4);

        function!(
            ETS.gamma,
            mink.new_slot(mu).to_atom(),
            bis.new_slot(i).to_atom(),
            bis.new_slot(j).to_atom()
        )
    }

    fn metric(mu: usize, nu: usize) -> Atom {
        let mink = Minkowski::rep(4);

        function!(
            ETS.metric,
            mink.new_slot(mu).to_atom(),
            mink.new_slot(nu).to_atom()
        )
    }

    fn test_and_assert(atom: Atom) {
        // let g_t = Network::parse_impl(atom.as_view())
        //     .contract::<Rational>(ContractionSettings::Normal)
        //     .unwrap()
        //     .tensor;
        let g_simp = GammaSimplified::gamma_symplify_impl(atom.into()).0;

        let _g_simp_t = Network::parse_impl(g_simp.as_view())
            .contract::<Rational>(ContractionSettings::Normal)
            .unwrap()
            .tensor;

        // assert_eq!(g_simp_t, g_t);

        insta::assert_snapshot!(g_simp.to_canonical_string());
    }
    let g = metric(2, 4)
        * metric(20, 40)
        * gamma(4, 1, 2)
        * gamma(40, 2, 3)
        * gamma(3, 3, 4)
        * gamma(30, 4, 5)
        * gamma(2, 5, 6)
        * gamma(20, 6, 7)
        * gamma(1, 7, 8)
        * gamma(5, 8, 1);

    test_and_assert(g);

    let g = gamma(1, 1, 2) * gamma(2, 2, 3) * gamma(3, 3, 4) * gamma(1, 4, 5);

    test_and_assert(g);

    let g = gamma(1, 1, 1);

    assert!(GammaSimplified::gamma_symplify_impl(g.into()).0.is_zero());

    let g = gamma(1, 1, 2) * gamma(2, 2, 1);

    test_and_assert(g);

    let g = gamma(1, 1, 2)
        * gamma(2, 2, 3)
        * gamma(10, 3, 4)
        * gamma(20, 4, 1)
        * gamma(1, 10, 20)
        * gamma(2, 20, 30)
        * gamma(30, 30, 40)
        * gamma(40, 40, 10);

    test_and_assert(g);

    let g = Atom::parse(
        "(Metric(mink(4,2),mink(4,3))-phat^-2*p(mink(4,2))*p(mink(4,3)))
        *Metric(mink(4,4),mink(4,5))
        *γ(mink(4,4),bis(4,9),bis(4,8))
        *γ(mink(4,28),bis(4,8),bis(4,11))
        *γ(mink(4,5),bis(4,11),bis(4,10))
        *γ(mink(4,29),bis(4,10),bis(4,7))
        *γ(mink(4,3),bis(4,7),bis(4,6))
        *γ(mink(4,30),bis(4,6),bis(4,5))
        *γ(mink(4,2),bis(4,5),bis(4,4))
        *γ(mink(4,27),bis(4,4),bis(4,9))
        *Q(4,mink(4,27))
        *Q(5,mink(4,28))
        *Q(6,mink(4,29))
        *Q(7,mink(4,30))",
    )
    .unwrap();
    test_and_assert(g);
}

#[test]
fn gamma_algebra() {
    let _ = ETS.gamma;
    let _ = UFO.t;
    let g = Atom::parse(
        "(Metric(mink(4,2),mink(4,3)))
    *(Metric(mink(4,7),mink(4,8))
        )
    *(-Metric(mink(4,4),mink(4,5))
        *Metric(mink(4,9),mink(4,6)))
    *Metric(mink(4,4),mink(4,7))
    *Metric(mink(4,5),mink(4,8))
    *Metric(mink(4,6),mink(4,11))
    *Metric(mink(4,9),mink(4,10))
    *γ(mink(4,2),bis(4,3),bis(4,2))
    *γ(mink(4,3),bis(4,5),bis(4,4))
    *γ(mink(4,10),bis(4,7),bis(4,6))
    *γ(mink(4,11),bis(4,9),bis(4,8))
    *γ(mink(4,31),bis(4,2),bis(4,5))
    *γ(mink(4,32),bis(4,8),bis(4,3))
    *γ(mink(4,33),bis(4,4),bis(4,7))
    *γ(mink(4,38),bis(4,6),bis(4,9))
    *Metric(mink(4,31),mink(4,32))
    *Metric(mink(4,33),mink(4,38))",
    )
    .unwrap();

    let g = GammaSimplified::gamma_symplify_impl(g.into()).0;
    println!(
        "{}",
        Network::parse_impl(g.as_view())
            .contract::<Rational>(ContractionSettings::Normal)
            .unwrap()
            .tensor
            .scalar()
            .unwrap()
    );
    println!("{g}");
}
#[test]
fn one_loop_lbl() {
    let (_model, amplitude, _path) = load_amplitude_output(
        &("TEST_AMPLITUDE_".to_string() + "physical_1L_6photons" + "/GL_OUTPUT"),
        true,
    );

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    let feyn = Numerator::default().from_graph(&graph.bare_graph, &GlobalPrefactor::default());

    println!("initial{:+}", feyn.get_single_atom().unwrap().0);
    println!(
        "canonized:{:+}",
        feyn.canonize_lorentz()
            .unwrap()
            .get_single_atom()
            .unwrap()
            .0
    );

    println!(
        "canonized with color:{:+}",
        feyn.color_simplify()
            .canonize_lorentz()
            .unwrap()
            .get_single_atom()
            .unwrap()
            .0
    );
}

#[test]

fn bug_check() {
    let a = Atom::parse("-1/9*𝑖*ee^2*G^2*(-TR+TR*Nc^2)*(P(0,mink(4,25))+K(1,mink(4,25)))*Metric(mink(4,0),mink(4,1))*Metric(mink(4,2),mink(4,3))*id(mink(4,2),mink(4,4))*id(mink(4,3),mink(4,5))*γ(mink(4,0),bis(4,9),bis(4,6))*γ(mink(4,1),bis(4,8),bis(4,7))*γ(mink(4,4),bis(4,5),bis(4,4))*γ(mink(4,5),bis(4,3),bis(4,2))*γ(mink(4,25),bis(4,4),bis(4,3))*γ(mink(4,27),bis(4,7),bis(4,9))*γ(mink(4,28),bis(4,6),bis(4,5))*γ(mink(4,29),bis(4,2),bis(4,8))*K(0,mink(4,27))*K(1,mink(4,28))*K(1,mink(4,29))").unwrap();
    //let b = Atom::parse("-1/9*𝑖*ee^2*G^2*(-TR+TR*Nc^2)*(P(0,mink(4,25))+K(1,mink(4,25)))*Metric(mink(4,0),mink(4,1))*Metric(mink(4,2),mink(4,3))*id(mink(4,2),mink(4,4))*id(mink(4,3),mink(4,5))*γ(mink(4,0),bis(4,9),bis(4,6))*γ(mink(4,1),bis(4,8),bis(4,7))*γ(mink(4,4),bis(4,5),bis(4,4))*γ(mink(4,5),bis(4,3),bis(4,2))*γ(mink(4,25),bis(4,4),bis(4,3))*γ(mink(4,27),bis(4,7),bis(4,9))*γ(mink(4,28),bis(4,6),bis(4,5))*γ(mink(4,29),bis(4,2),bis(4,8))*K(0,mink(4,27))*K(1,mink(4,28))*K(1,mink(4,29))").unwrap();
    let b = a.clone() * -1;
    println!("a/b={}  TT", a / b);
}

#[test]
fn bug_check_b() {
    let a = Atom::parse(
         "-4/9*ee^2*yt^2*G^4*(MT*id(bis(4,2),bis(4,5))+γ(mink(4,41),bis(4,2),bis(4,5))*K(0,mink(4,41)))*(MT*id(bis(4,3),bis(4,8))+γ(mink(4,43),bis(4,8),bis(4,3))*K(1,mink(4,43)))*(MT*id(bis(4,4),bis(4,7))-γ(mink(4,44),bis(4,4),bis(4,7))*K(1,mink(4,44)))*(MT*id(bis(4,6),bis(4,15))+γ(mink(4,46),bis(4,6),bis(4,15))*K(2,mink(4,46)))*(MT*id(bis(4,9),bis(4,12))-(K(2,mink(4,48))-K(3,mink(4,48)))*γ(mink(4,48),bis(4,12),bis(4,9)))*(MT*id(bis(4,10),bis(4,13))+(P(0,mink(4,50))+K(2,mink(4,50))-K(3,mink(4,50)))*γ(mink(4,50),bis(4,10),bis(4,13)))*(MT*id(bis(4,11),bis(4,14))-(P(0,mink(4,51))+K(2,mink(4,51)))*γ(mink(4,51),bis(4,14),bis(4,11)))*(ProjM(bis(4,3),bis(4,2))+ProjP(bis(4,3),bis(4,2)))*(ProjM(bis(4,5),bis(4,4))+ProjP(bis(4,5),bis(4,4)))*(-(-K(1,mink(4,6))-K(2,mink(4,6)))*Metric(mink(4,5),mink(4,7))+(-K(1,mink(4,7))-K(2,mink(4,7)))*Metric(mink(4,5),mink(4,6))+(K(1,mink(4,5))+K(2,mink(4,5))-K(3,mink(4,5)))*Metric(mink(4,6),mink(4,7))-(K(1,mink(4,7))+K(2,mink(4,7))-K(3,mink(4,7)))*Metric(mink(4,5),mink(4,6))+Metric(mink(4,5),mink(4,7))*K(3,mink(4,6))-Metric(mink(4,6),mink(4,7))*K(3,mink(4,5)))*sqrt(2)^-2*Metric(mink(4,0),mink(4,1))*Metric(mink(4,2),mink(4,5))*Metric(mink(4,3),mink(4,6))*Metric(mink(4,4),mink(4,7))*id(mink(4,0),mink(4,9))*id(mink(4,1),mink(4,8))*γ(mink(4,2),bis(4,7),bis(4,6))*γ(mink(4,3),bis(4,9),bis(4,8))*γ(mink(4,4),bis(4,11),bis(4,10))*γ(mink(4,8),bis(4,13),bis(4,12))*γ(mink(4,9),bis(4,15),bis(4,14))",
     )
     .unwrap();
    let indices = [
        ("mink(4,0)", "m"),
        ("mink(4,1)", "m"),
        ("mink(4,2)", "m"),
        ("mink(4,3)", "m"),
        ("mink(4,4)", "m"),
        ("mink(4,5)", "m"),
        ("mink(4,6)", "m"),
        ("mink(4,7)", "m"),
        ("mink(4,8)", "m"),
        ("mink(4,9)", "m"),
        ("mink(4,41)", "m"),
        ("mink(4,43)", "m"),
        ("mink(4,44)", "m"),
        ("mink(4,46)", "m"),
        ("mink(4,48)", "m"),
        ("mink(4,50)", "m"),
        ("mink(4,51)", "m"),
        ("bis(4,2)", "b"),
        ("bis(4,3)", "b"),
        ("bis(4,4)", "b"),
        ("bis(4,5)", "b"),
        ("bis(4,6)", "b"),
        ("bis(4,7)", "b"),
        ("bis(4,8)", "b"),
        ("bis(4,9)", "b"),
        ("bis(4,10)", "b"),
        ("bis(4,11)", "b"),
        ("bis(4,12)", "b"),
        ("bis(4,13)", "b"),
        ("bis(4,14)", "b"),
        ("bis(4,15)", "b"),
    ]
    .iter()
    .map(|(a, g)| (Atom::parse(a).unwrap(), g))
    .collect::<Vec<_>>();

    println!("res={}", a.canonize_tensors(&indices).unwrap());
}

#[test]
fn bug_check_a() {
    let a = Atom::parse(
        "-(-K(1,mink(4,6))-K(2,mink(4,6)))*Metric(mink(4,5),mink(4,7))
 +(-K(1,mink(4,7))-K(2,mink(4,7)))*Metric(mink(4,5),mink(4,6))
 +(K(1,mink(4,5))+K(2,mink(4,5))-K(3,mink(4,5)))*Metric(mink(4,6),mink(4,7))
 -(K(1,mink(4,7))+K(2,mink(4,7))-K(3,mink(4,7)))*Metric(mink(4,5),mink(4,6))
 +Metric(mink(4,5),mink(4,7))*K(3,mink(4,6))
 -Metric(mink(4,6),mink(4,7))*K(3,mink(4,5))",
    )
    .unwrap();
    let indices: [(Atom, &str); 3] = [
        (Atom::parse("mink(4,5)").unwrap(), "m"),
        (Atom::parse("mink(4,6)").unwrap(), "m"),
        (Atom::parse("mink(4,7)").unwrap(), "m"),
    ];

    println!("res={}", a.canonize_tensors(&indices).unwrap());
}

// Ignore this test by default as it's too memory intensive for the CI
#[ignore]
#[test]
fn one_loop_lbl_concretize() {
    let (_model, amplitude, _path) = load_amplitude_output(
        &("TEST_AMPLITUDE_".to_string() + "physical_1L_6photons" + "/GL_OUTPUT"),
        true,
    );

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    let feyn = Numerator::default()
        .from_graph(&graph.bare_graph, &GlobalPrefactor::default())
        .color_simplify()
        .parse();

    let reps = feyn.random_concretize_reps(None, true);
    let rep_views = reps
        .iter()
        .map(|(a, b)| (a.as_view(), b.as_view()))
        .collect::<Vec<_>>();

    println!(
        "{}",
        feyn.apply_reps(rep_views)
            .contract::<Rational>(ContractionSettings::Normal)
            .unwrap()
            .state
            .tensor
            .scalar()
            .unwrap()
            .expand()
    );
}

// #[test]
// fn color_simple() {
//     println!("{}",ColorSimplified::color_simplify_impl(&Atom::parse("(-1*T(coad(8,0),cof(3,j(17,18,19)),dind(cof(3,k(17,18,19))))*T(coad(8,1),cof(3,11),dind(cof(3,j(17,18,19))))*T(coad(8,2),cof(3,5),dind(cof(3,14)))*T(coad(8,3),cof(3,14),dind(cof(3,11)))*T(coad(8,4),cof(3,k(17,18,19)),dind(cof(3,5)))+T(coad(8,0),cof(3,11),dind(cof(3,j(17,18,19))))*T(coad(8,1),cof(3,j(17,18,19)),dind(cof(3,k(17,18,19))))*T(coad(8,2),cof(3,5),dind(cof(3,14)))*T(coad(8,3),cof(3,14),dind(cof(3,11)))*T(coad(8,4),cof(3,k(17,18,19)),dind(cof(3,5))))").unwrap().into()).map_or_else(|a|if let ColorError::NotFully(a)=a{
//         a
//     }else{
//         panic!("")
//     },|a|a))
// }

#[test]
fn dumb_four_gluon() {
    let model = load_generic_model("sm");

    let gggg = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_37"),
    };
    let mut four_gluon = symbolica::graph::Graph::new();
    let v = four_gluon.add_node(gggg);
    let g = EdgeColor::from_particle(model.get_particle("g"));

    four_gluon.add_edge(v, v, false, g).unwrap();
    four_gluon.add_edge(v, v, false, g).unwrap();

    let graph = BareGraph::from_symbolica_graph(
        &model,
        "gggg".into(),
        &four_gluon,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let num = Numerator::default().from_graph(&graph, &GlobalPrefactor::default());
    // println!("{}", num.state.color);
    // println!("{}", num.state.colorless);

    let colorsimp = num.color_simplify();
    // println!("{}", colorsimp.state.color);

    let gamma = colorsimp.clone().gamma_simplify();
    // println!("{}", gamma.state.colorless);

    let gammasingle = gamma.get_single_atom().unwrap();

    // println!("{}", gammasingle.0.expand());

    let expanded = colorsimp
        .parse()
        .contract::<Rational>(ContractionSettings::Normal)
        .unwrap()
        .state
        .tensor
        .scalar()
        .unwrap()
        .expand();

    assert!((&expanded - gammasingle.0).expand().is_zero());
    assert_snapshot!(expanded.to_canonical_string());
}
