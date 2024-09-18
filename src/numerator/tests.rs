use std::{
    fs::File,
    path::{Path, PathBuf},
};

use color_eyre::Help;
use eyre::Context;
use petgraph::graph;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use rayon::vec;
use spenso::{complex::Complex, parametric::PatternReplacement, structure::HasStructure};
use symbolica::{
    atom::Atom,
    domains::rational::Rational,
    id::{Pattern, Replacement},
    state::State,
};

use crate::{
    cross_section::Amplitude,
    gammaloop_integrand::DefaultSample,
    graph::{BareGraph, Graph, SerializableGraph},
    model::Model,
    momentum::{ExternalMomenta, FourMomentum, Helicity, Polarization},
    numerator::{Contracted, ContractionSettings, ExtraInfo},
    tests_from_pytest::{sample_generator, test_export_settings},
    utils::{approx_eq_res, f128, F},
    Externals, Settings,
};

use super::{Evaluate, Global, Numerator, UnInit};

#[test]
fn hhgghh() {
    let expr = Atom::parse("4/9*ee^2*yt^4*Nc*(MT*id(aind(bis(4,0),bis(4,3)))+Q(6,aind(lord(4,5)))*γ(aind(loru(4,5),bis(4,0),bis(4,3))))*(MT*id(aind(bis(4,2),bis(4,5)))+Q(7,aind(lord(4,11)))*γ(aind(loru(4,11),bis(4,2),bis(4,5))))*(MT*id(aind(bis(4,4),bis(4,7)))+Q(8,aind(lord(4,29)))*γ(aind(loru(4,29),bis(4,4),bis(4,7))))*(MT*id(aind(bis(4,6),bis(4,9)))+Q(9,aind(lord(4,77)))*γ(aind(loru(4,77),bis(4,6),bis(4,9))))*(MT*id(aind(bis(4,8),bis(4,11)))+Q(10,aind(lord(4,203)))*γ(aind(loru(4,203),bis(4,8),bis(4,11))))*(MT*id(aind(bis(4,10),bis(4,1)))+Q(11,aind(lord(4,533)))*γ(aind(loru(4,533),bis(4,10),bis(4,1))))*(ProjM(aind(bis(4,1),bis(4,0)))+ProjP(aind(bis(4,1),bis(4,0))))*(ProjM(aind(bis(4,3),bis(4,2)))+ProjP(aind(bis(4,3),bis(4,2))))*(ProjM(aind(bis(4,5),bis(4,4)))+ProjP(aind(bis(4,5),bis(4,4))))*(ProjM(aind(bis(4,7),bis(4,6)))+ProjP(aind(bis(4,7),bis(4,6))))*sqrt(2)^-4*γ(aind(lord(4,2),bis(4,9),bis(4,8)))*γ(aind(lord(4,3),bis(4,11),bis(4,10)))*ϵbar(2,aind(loru(4,3)))*ϵbar(3,aind(loru(4,2)))").unwrap();

    let replacements = vec![
        (
            Pattern::parse("ϵbar(a_,cind(0))").unwrap(),
            Pattern::parse("0").unwrap().into(),
        ),
        (
            Pattern::parse("ϵbar(a_,cind(2))").unwrap(),
            Pattern::parse("0").unwrap().into(),
        ),
        (
            Pattern::parse("ϵbar(a_,cind(3))").unwrap(),
            Pattern::parse("0").unwrap().into(),
        ),
        (
            Pattern::parse("Q(a__)").unwrap(),
            Pattern::parse("1").unwrap().into(),
        ),
    ];

    let reps: Vec<_> = replacements
        .iter()
        .map(|(a, b)| Replacement::new(a, b))
        .collect();
    let mut new = Numerator {
        state: Global::new(expr.clone().into()),
    }
    .color_symplify()
    .gamma_symplify()
    // .parse()
    // .contract(ContractionSettings::<Rational>::Normal)
    // .unwrap()
    // .state
    // .tensor;
    ;

    // new.replace_all_multiple_repeat_mut(&reps);

    // let new = new
    //     .tensor
    //     .try_as_dense()
    //     .unwrap()
    //     .map_data_ref(|a| a.expand().factor());

    println!("{}", new.export());
    // println!("{}", State::I);

    let n_edges = 12;
    let pol_data = vec![("ϵbar".into(), 2, 4), ("ϵbar".into(), 3, 4)];
    let mut params = Contracted::generate_kinematic_params_impl(n_edges, pol_data);

    let seed = 9;
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut param_values: Vec<Complex<F<f64>>> = (0..params.len())
        .map(|_| Complex::new_re(F(rng.gen())))
        .collect::<Vec<_>>();

    // params.push(Atom::new_var(State::I));
    // param_values.push(Complex::new_i());

    let model_params_start = params.len();
    let model_params = vec![
        ("MT", 173.0),
        ("ee", 0.3079537672443688),
        ("yt", 0.9936661458150063),
    ];

    for (name, value) in model_params {
        params.push(Atom::new_var(State::get_symbol(name)));
        param_values.push(Complex::new(F(value), F(0.)));
    }

    let quad_values: Vec<Complex<F<f128>>> =
        param_values.iter().map(|c| c.map(|f| f.higher())).collect();

    let extra_info = ExtraInfo {
        path: PathBuf::new().join("ignore"),
        orientations: vec![vec![true; n_edges]],
    };

    let mut new = Numerator {
        state: Global::new(expr.clone().into()),
    }
    .color_symplify()
    .parse()
    .contract(ContractionSettings::<Rational>::Normal)
    .unwrap()
    .generate_evaluators_from_params(
        n_edges,
        model_params_start,
        &params,
        param_values.clone(),
        quad_values.clone(),
        &extra_info,
        &test_export_settings(),
    );

    let mut newg = Numerator {
        state: Global::new(expr.into()),
    }
    .color_symplify()
    .gamma_symplify()
    .parse()
    .contract(ContractionSettings::<Rational>::Normal)
    .unwrap()
    .generate_evaluators_from_params(
        n_edges,
        model_params_start,
        &params,
        param_values,
        quad_values,
        &extra_info,
        &test_export_settings(),
    );

    for i in 0..100 {
        let emr: Vec<FourMomentum<F<f64>>> = (0..n_edges)
            .map(|_| {
                FourMomentum::from_args(F(rng.gen()), F(rng.gen()), F(rng.gen()), F(rng.gen()))
            })
            .collect::<Vec<_>>();

        let polarizations: Vec<Polarization<Complex<F<f64>>>> = vec![
            Polarization::lorentz([F(rng.gen()), F(rng.gen()), F(rng.gen()), F(rng.gen())]).cast(),
            Polarization::lorentz([F(rng.gen()), F(rng.gen()), F(rng.gen()), F(rng.gen())]).cast(),
        ];

        let val = new.evaluate_single(&emr, &polarizations, &Settings::default());

        let valg = newg.evaluate_single(&emr, &polarizations, &Settings::default());

        assert_eq!(val, valg);

        print!("{}", val);
        println!("{}", valg);
    }
}

#[test]
fn trees() {
    let _ = env_logger::builder().is_test(true).try_init();
    let tree_name = "th_th";
    let amp_name = "tree_amplitude_1_th_th";
    let file_path = PathBuf::new()
        .join("./src/test_resources/Trees")
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
    let emr = graph.bare_graph.cff_emr_from_lmb(&sample, lmb);

    let three_emr = lmb.spatial_emr(&sample);

    let onshell_energies = graph
        .bare_graph
        .compute_onshell_energies(&sample.loop_moms, &sample.external_moms);

    for ((e, q), p) in onshell_energies.iter().zip(emr).zip(three_emr) {
        approx_eq_res(e, &q.temporal.value, &F(1e-10)).unwrap();

        for (a, b) in q.spatial.into_iter().zip(p) {
            approx_eq_res(&a, &b, &F(1e-10)).unwrap();
        }
    }

    let export_path = PathBuf::new()
        .join("./src/test_resources/Trees")
        .join(tree_name);

    let contraction_settings = ContractionSettings::<Rational>::Normal;

    let export_settings = test_export_settings();

    graph.generate_cff();
    let mut graph =
        graph.process_numerator(&model, contraction_settings, export_path, &export_settings);

    let external_mom = vec![
        FourMomentum::from_args(F(592.625), F(0.), F(0.), F(566.8116)), // 1
        FourMomentum::from_args(F(579.375), F(0.), F(0.), F(-566.8116)), // 2
        FourMomentum::from_args(F(592.625), F(125.7463), F(504.2705), F(-226.2178)), // 3
    ];

    let sample = DefaultSample::new(
        vec![],
        external_mom,
        F(1.),
        &[
            Helicity::Minus,
            Helicity::Zero,
            Helicity::Plus,
            Helicity::Zero,
        ],
        &graph.bare_graph,
    );
    let settings = Settings::default();

    let val = graph
        .evaluate_fourd_expr(&[], &sample.external_moms, &sample.polarizations, &settings)
        .scalar()
        .unwrap();

    let energy_product = graph
        .bare_graph
        .compute_energy_product(&sample.loop_moms, &sample.external_moms);

    let valcff = graph.evaluate_cff_expression(&sample, &settings) / energy_product;
    println!("4d: {}", val);
    println!("CFF: {}", valcff);
}

#[test]
fn tree_ta_ta_1() {
    let (model, amplitude, path) = load_tree("ta_ta", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    let export_settings = test_export_settings();

    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::<Rational>::Normal,
        path,
        &test_export_settings(),
    );
    let sample: DefaultSample<f64> = sample_generator(3, &graph.bare_graph, None);

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);

    let val = graph
        .evaluate_fourd_expr(
            &[],
            &sample.external_moms,
            &sample.polarizations,
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    println!("4d: {}", val);
    println!("CFF: {}", cff_val);

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

    let sample = DefaultSample::new(
        vec![],
        externals.get_indep_externals(&[]).0,
        F(1.),
        &externals.get_helicities(),
        &graph.bare_graph,
    );

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);

    let val = graph
        .evaluate_fourd_expr(
            &[],
            &sample.external_moms,
            &sample.polarizations,
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    println!("4d: {}", val);
    println!("CFF: {}", cff_val);
}

#[test]
fn tree_h_ttxaah_1() {
    let (model, amplitude, path) = load_tree("h_ttxaah", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::<Rational>::Normal,
        path,
        &test_export_settings(),
    );
    let sample: DefaultSample<f64> = sample_generator(3, &graph.bare_graph, None);

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);

    let val = graph
        .evaluate_fourd_expr(
            &[],
            &sample.external_moms,
            &sample.polarizations,
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    println!("4d: {}", val);
    println!("CFF: {}", cff_val);
}


#[test]
fn tree_hh_ttxaa_1() {
    let (model, amplitude, path) = load_tree("hh_ttxaa", 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::<Rational>::Normal,
        path,
        &test_export_settings(),
    );
    let sample: DefaultSample<f64> = sample_generator(3, &graph.bare_graph, None);

    let cff_val = graph.evaluate_cff_expression(&sample, &Settings::default())
        / graph
            .bare_graph
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);

    let val = graph
        .evaluate_fourd_expr(
            &[],
            &sample.external_moms,
            &sample.polarizations,
            &Settings::default(),
        )
        .scalar()
        .unwrap();

    println!("4d: {}", val);
    println!("CFF: {}", cff_val);
}

fn load_tree(tree_name: &str, amp_num: usize) -> (Model, Amplitude<UnInit>, PathBuf) {
    let _ = env_logger::builder().is_test(true).try_init();
    let file_path = PathBuf::new()
        .join("./src/test_resources/Trees")
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

    let export_path = PathBuf::new()
        .join("./src/test_resources/Trees")
        .join(tree_name);

    (model, amplitude, export_path)
}

#[test]
fn tree_h_ttxaah_0() {
    let expr = Atom::parse("-8/3*𝑖*ee^2*vev*lam*yt*(MT*id(aind(bis(4,3),bis(4,4)))+Q(6,aind(lord(4,5)))*γ(aind(loru(4,5),bis(4,3),bis(4,4))))*(MT*id(aind(bis(4,5),bis(4,6)))+Q(7,aind(lord(4,11)))*γ(aind(loru(4,11),bis(4,5),bis(4,6))))*(ProjM(aind(bis(4,7),bis(4,6)))+ProjP(aind(bis(4,7),bis(4,6))))*sqrt(2)^-1*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,6),cof(3,7)))*id(aind(coaf(3,7),cof(3,8)))*id(aind(coaf(3,8),cof(3,9)))*id(aind(coaf(3,9),cof(3,10)))*γ(aind(lord(4,2),bis(4,3),bis(4,2)))*γ(aind(lord(4,3),bis(4,5),bis(4,4)))*ubar(1,aind(bis(4,2)))*v(2,aind(bis(4,7)))*ϵbar(3,aind(loru(4,2)))*ϵbar(4,aind(loru(4,3)))").unwrap();

    let new = Numerator {
        state: Global::new(expr.clone().into()),
    }
    .color_symplify()
    .color_project()
    // .gamma_symplify()
    .parse()
    .contract(ContractionSettings::<Rational>::Normal)
    .unwrap();
    println!("{}", new.export());
}

#[test]
fn color_three_loop_physical_photon() {
    //     let color_expr= Atom::parse(concat!(
    //         // "T(aind(coad(8,9),cof(3,8),coaf(3,7)))*T(aind(coad(8,14),cof(3,13),coaf(3,12)))*T(aind(coad(8,21),cof(3,20),coaf(3,19)))*T(aind(coad(8,26),cof(3,25),coaf(3,24)))*",
    //     // "id(aind(coaf(3,3),cof(3,4)))*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,10),cof(3,11)))*id(aind(coaf(3,15),cof(3,16)))*id(aind(coaf(3,17),cof(3,18)))*id(aind(coaf(3,22),cof(3,23)))*id(aind(cof(3,4),coaf(3,24)))*id(aind(cof(3,6),coaf(3,3)))*id(aind(cof(3,8),coaf(3,5)))*id(aind(cof(3,11),coaf(3,7)))*id(aind(cof(3,13),coaf(3,10)))*id(aind(cof(3,16),coaf(3,12)))*id(aind(cof(3,18),coaf(3,15)))*id(aind(cof(3,20),coaf(3,17)))*id(aind(cof(3,23),coaf(3,19)))*id(aind(cof(3,25),coaf(3,22)))*id(aind(coad(8,21),coad(8,9)))*id(aind(coad(8,26),coad(8,14)))",
    // "id(aind(coaf(3,3),cof(3,4)))*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,10),cof(3,11)))*id(aind(coaf(3,15),cof(3,16)))*id(aind(coaf(3,17),cof(3,18)))*id(aind(coaf(3,22),cof(3,23)))*id(aind(cof(3,4),coaf(3,24)))*id(aind(cof(3,6),coaf(3,3)))*id(aind(cof(3,8),coaf(3,5)))*id(aind(cof(3,11),coaf(3,7)))*id(aind(cof(3,13),coaf(3,10)))*id(aind(cof(3,16),coaf(3,12)))*id(aind(cof(3,18),coaf(3,15)))*id(aind(cof(3,20),coaf(3,17)))*id(aind(cof(3,23),coaf(3,19)))*id(aind(cof(3,25),coaf(3,22)))*id(aind(coad(8,21),coad(8,9)))*id(aind(coad(8,26),coad(8,14)))")).unwrap();

    //     let mut color_num = Numerator {
    //         expression: color_expr,
    //         network: None,
    //         const_map: AHashMap::new(),
    //         eval: crate::numerator::EvaluatorNumerator::UnInit,
    //         extra_info: ExtraInfo {
    //             emr_size: 17,
    //             graph_name: "srt".into(),
    //         },
    //     };

    //     color_num.process_color_simple();

    //     println!("{}", color_num.expression)
}
