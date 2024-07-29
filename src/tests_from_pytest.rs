#![allow(unused_imports)]
use crate::cff::esurface::{
    get_existing_esurfaces, Esurface, EsurfaceID, ExistingEsurfaceId, ExistingEsurfaces,
};
use crate::cff::generation::generate_cff_expression;
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::gammaloop_integrand::DefaultSample;
use crate::graph::{
    Edge, EdgeType, HasVertexInfo, InteractionVertexInfo, SerializableGraph, VertexInfo,
};
use crate::model::Model;
use crate::momentum::{FourMomentum, ThreeMomentum};
use crate::numerator::Numerator;
use crate::subtraction::overlap::{self, find_center, find_maximal_overlap};
use crate::subtraction::static_counterterm;
use crate::utils::{
    assert_approx_eq, compute_momentum, compute_three_momentum_from_four, PrecisionUpgradable,
};
use crate::utils::{f128, F};
use ahash::AHashMap;
use colored::Colorize;
use itertools::{FormatWith, Itertools};
//use libc::__c_anonymous_ptrace_syscall_info_exit;
use lorentz_vector::LorentzVector;
use petgraph::algo::greedy_matching;
use petgraph::graph;
use rayon::prelude::IndexedParallelIterator;
use serde;
use spenso::complex::Complex;
use statrs::function::evaluate;
use std::fs::File;
use std::path::Path;
use std::{clone, env};
use symbolica;
use symbolica::atom::Atom;
use symbolica::domains::float::Complex as SymComplex;

#[allow(unused)]
const LTD_COMPARISON_TOLERANCE: F<f64> = F(1.0e-12);

pub fn load_amplitude_output(output_path: &str, load_generic_model: bool) -> (Model, Amplitude) {
    let output_dir = if let Ok(pytest_output_path) = env::var("PYTEST_OUTPUT_PATH_FOR_RUST") {
        pytest_output_path
    } else {
        String::from("./src/test_resources")
    };
    let path = Path::new(&output_dir).join(output_path);
    let output_meta_data: OutputMetaData = serde_yaml::from_reader(
        File::open(path.join("output_metadata.yaml")).unwrap_or_else(|e| {
            panic!(
                "Could not open file '{:?}'. Error: {:?}",
                path.join("output_metadata.yaml"),
                e
            )
        }),
    )
    .unwrap();
    assert_eq!(output_meta_data.output_type, OutputType::Amplitudes);
    assert_eq!(output_meta_data.contents.len(), 1);

    let model = if load_generic_model {
        Model::from_file(String::from(
            Path::new(&output_dir)
                .join(format!(
                    "gammaloop_models/{}.yaml",
                    output_meta_data.model_name.as_str()
                ))
                .to_str()
                .unwrap(),
        ))
        .unwrap()
    } else {
        Model::from_file(String::from(
            path.join(format!(
                "sources/model/{}.yaml",
                output_meta_data.model_name
            ))
            .to_str()
            .unwrap(),
        ))
        .unwrap()
    };

    let amplitude = Amplitude::from_file(
        &model,
        String::from(
            path.join(format!(
                "sources/amplitudes/{}/amplitude.yaml",
                output_meta_data.contents[0]
            ))
            .to_str()
            .unwrap(),
        ),
    )
    .unwrap();
    (model, amplitude)
}

#[cfg(test)]
mod tests_scalar_massless_triangle {
    use lorentz_vector::LorentzVector;
    use rayon::prelude::IndexedParallelIterator;
    use smartstring::SmartString;
    use spenso::complex::Complex;
    use symbolica::domains::float::Complex as SymComplex;

    use crate::{
        gammaloop_integrand::DefaultSample,
        graph::EdgeType,
        momentum::{FourMomentum, ThreeMomentum},
        observables::AFBSettings,
        utils::F,
    };

    use super::*;

    #[test]
    fn pytest_massless_scalar_triangle() {
        let (model, amplitude) =
            load_amplitude_output("TEST_AMPLITUDE_massless_scalar_triangle/GL_OUTPUT", true);

        assert_eq!(model.name, "scalars");
        assert!(amplitude.amplitude_graphs.len() == 1);
        assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 6);
        assert!(
            amplitude.amplitude_graphs[0]
                .graph
                .external_connections
                .len()
                == 3
        );

        let mut graph = amplitude.amplitude_graphs[0].graph.clone();
        graph.generate_loop_momentum_bases();

        graph.generate_cff();
        graph.process_numerator(&model);
        graph.numerator_substitute_model_params(&model);
        // graph.evaluate_model_params(&model);
        graph.process_numerator(&model);
        assert_eq!(
            graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .get_num_trees(),
            6
        );
        assert_eq!(
            graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .esurfaces
                .len(),
            6
        );

        let p1 = FourMomentum::from_args(F(1.), F(3.), F(4.0), F(5.0));
        let p2 = FourMomentum::from_args(F(1.), F(6.0), F(7.), F(8.));

        let k = ThreeMomentum::new(F(1.), F(2.), F(3.));

        let energy_product = graph.compute_energy_product(&[k], &[p1, p2]);

        let sample = DefaultSample {
            loop_moms: vec![k],
            external_moms: vec![p1, p2],
            jacobian: F(1.0),
        };

        let cff_res = graph.evaluate_cff_expression(&sample, 0) / energy_product;

        // println!("res = {:+e}", res);

        // test the lmb generation;;
        assert_eq!(
            graph
                .derived_data
                .loop_momentum_bases
                .as_ref()
                .unwrap()
                .len(),
            3
        );

        let absolute_truth = Complex::new(F(0.0), F(4.531238289663331e-6));

        graph.generate_ltd();

        let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1, p2]) / energy_product;

        println!("cff_res = {:+e}", cff_res);
        println!("ltd_res = {:+e}", ltd_res);
        println!("cltd_m = {:+e}", absolute_truth);

        assert_approx_eq(&cff_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
        assert_approx_eq(&cff_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
        assert_approx_eq(&ltd_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
        assert_approx_eq(&ltd_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);

        let propagator_groups = graph.group_edges_by_signature();
        assert_eq!(propagator_groups.len(), 3);

        let generate_data = graph.generate_esurface_data();

        if let Err(e) = generate_data {
            panic!("Error: {}", e);
        }

        let existing = get_existing_esurfaces(
            &graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .esurfaces,
            graph.derived_data.esurface_derived_data.as_ref().unwrap(),
            &[p1, p2],
            &graph.loop_momentum_basis,
            0,
            F(2.0),
        );

        assert_eq!(existing.len(), 0);

        let cff = graph.get_cff();
        let unfolded = cff.expand_terms();

        assert_eq!(unfolded.len(), 6);

        for term in unfolded.iter() {
            assert_eq!(term.len(), 2);
        }
    }
}

#[test]
fn pytest_scalar_fishnet_2x2() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_fishnet_2x2/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 16);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 13);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 4
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    let cff = generate_cff_expression(&graph).unwrap();
    assert!(!cff.esurfaces.is_empty());

    // println!("basis size: {:?}", graph.loop_momentum_basis.basis);
    graph.generate_loop_momentum_bases();

    let k1 = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.));
    let k2 = ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.));
    let k3 = ThreeMomentum::new(F(17. / 19.), F(19. / 23.), F(23. / 29.));
    let k4: ThreeMomentum<F<f64>> = ThreeMomentum::new(29. / 31., 31. / 37., 37. / 41.).into();
    let p1 = FourMomentum::from_args(79. / 83., 41. / 43., 43. / 47., 47. / 53.).into();
    let p2 = FourMomentum::from_args(83. / 89., 53. / 59., 59. / 61., 61. / 67.).into();
    let p3 = FourMomentum::from_args(89. / 97., 67. / 71., 71. / 73., 73. / 79.).into();

    let emr = graph.compute_emr(&[k1, k2, k3, k4], &[p1, p2, p3]);
    let n_lmb = graph
        .clone()
        .derived_data
        .loop_momentum_bases
        .unwrap()
        .len();
    assert_eq!(n_lmb, 192);
    //println!("number of lmbs: {}", n_lmb);

    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);
    graph.generate_ltd();

    for basis in graph
        .derived_data
        .loop_momentum_bases
        .as_ref()
        .unwrap()
        .iter()
    {
        let momenta_in_basis = basis.basis.iter().map(|index| emr[*index]).collect_vec();
        let new_emr = basis
            .edge_signatures
            .iter()
            .map(|s| compute_three_momentum_from_four(s, &momenta_in_basis, &[p1, p2, p3]))
            .collect_vec();
        assert_eq!(emr.len(), new_emr.len());

        for (e1, e2) in emr.iter().zip(new_emr.iter()) {
            assert_approx_eq(&e1.px, &e2.px, &LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(&e1.py, &e2.py, &LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(&e1.pz, &e2.pz, &LTD_COMPARISON_TOLERANCE);
        }
    }

    //println!("lmb consistency check passed");

    let k1_f128: ThreeMomentum<F<f128>> = k1.cast().higher();
    let k2_f128: ThreeMomentum<F<f128>> = k2.cast().higher();
    let k3_f128: ThreeMomentum<F<f128>> = k3.cast().higher();
    let k4_f128: ThreeMomentum<F<f128>> = k4.cast().higher();
    let p1_f128: FourMomentum<F<f128>> = p1.cast().higher();
    let p2_f128: FourMomentum<F<f128>> = p2.cast().higher();
    let p3_f128: FourMomentum<F<f128>> = p3.cast().higher();

    let loop_moms_f128 = vec![k1_f128, k2_f128, k3_f128, k4_f128];
    let externals_f128 = vec![p1_f128, p2_f128, p3_f128];

    let energy_product = graph.compute_energy_product(&loop_moms_f128, &externals_f128);

    let ltd_res = graph.evaluate_ltd_expression(&loop_moms_f128, &externals_f128) / &energy_product;
    let sample = DefaultSample {
        loop_moms: loop_moms_f128,
        external_moms: externals_f128,
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / &energy_product;

    let absolute_truth: Complex<F<f128>> = Complex::new(
        F::<f128>::from_f64(0.000019991301832169422),
        F::<f128>::from_f64(0.0),
    );

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    assert_approx_eq(
        &cff_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &cff_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &ltd_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 12);
    // TODO: @Mathijs, you can put your own checks there
}

#[test]
fn pytest_scalar_sunrise() {
    let (model, amplitude) = load_amplitude_output("TEST_AMPLITUDE_scalar_sunrise/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 5);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 4);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 2
    );

    let k1 = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.));
    let k2 = ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.));

    let p1 = FourMomentum::from_args(F(0.), F(0.), F(0.), F(0.));

    let absolute_truth = Complex::new(F(0.24380172488169907), F(0.));

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_loop_momentum_bases();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);
    graph.generate_ltd();

    let energy_product = graph.compute_energy_product(&[k1, k2], &[p1]);

    let sample = DefaultSample {
        loop_moms: vec![k1, k2],
        external_moms: vec![p1],
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k1, k2], &[p1]) / energy_product;

    println!("cff_res = {:+e}", cff_res);
    println!("ltd_res = {:+e}", ltd_res);

    assert_approx_eq(&cff_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&cff_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
    println!("cff correct");
    assert_approx_eq(&ltd_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&ltd_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
    println!("ltd correct");

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 3);
}

#[test]
fn pytest_scalar_fishnet_2x3() {
    let (model, mut amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_fishnet_2x3/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 21);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 16);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 4
    );

    // generate cff
    amplitude.amplitude_graphs[0].graph.generate_cff();
    amplitude.amplitude_graphs[0]
        .graph
        .process_numerator(&model);
    amplitude.amplitude_graphs[0]
        .graph
        .numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    amplitude.amplitude_graphs[0]
        .graph
        .process_numerator(&model);
    amplitude.amplitude_graphs[0].graph.generate_ltd();

    let externals: Vec<FourMomentum<F<f128>>> = (0..3)
        .map(|i| {
            FourMomentum::from_args(
                F(1.),
                F(i as f64) + F(2.),
                F(i as f64) + F(3.),
                F(i as f64) + F(4.0),
            )
            .cast()
            .higher()
        })
        .collect_vec();

    let loop_moms: Vec<ThreeMomentum<F<f128>>> = (0..6)
        .map(|i| {
            ThreeMomentum::new(
                F(i as f64) - F(2.),
                F(i as f64) + F(3.5),
                F(i as f64) + F(4.5001),
            )
            .cast()
            .higher()
        })
        .collect_vec();

    // let before_cff = std::time::Instant::now();

    let sample = DefaultSample {
        loop_moms: loop_moms.clone(),
        external_moms: externals.clone(),
        jacobian: F(1.0),
    };

    let cff_res = amplitude.amplitude_graphs[0]
        .graph
        .evaluate_cff_expression(&sample, 0);

    // let cff_duration = before_cff.elapsed();
    // println!("cff_duration: {}", cff_duration.as_micros());

    // let before_ltd = std::time::Instant::now();
    let ltd_res = amplitude.amplitude_graphs[0]
        .graph
        .evaluate_ltd_expression(&loop_moms, &externals);

    // let ltd_duration = before_ltd.elapsed();
    // println!("ltd_duration: {}", ltd_duration.as_micros());

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    assert_approx_eq(&cff_res.re, &ltd_res.re, &ltd_comparison_tolerance128);

    assert_approx_eq(&cff_res.im, &ltd_res.im, &ltd_comparison_tolerance128);

    let propagator_groups = amplitude.amplitude_graphs[0]
        .graph
        .group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 17);
    // TODO: @Mathijs, you can put your own checks there
}

#[test]
fn pytest_scalar_cube() {
    let (model, amplitude) = load_amplitude_output("TEST_AMPLITUDE_scalar_cube/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 20);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 16);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 8
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_loop_momentum_bases();
    graph.generate_cff();
    graph.generate_ltd();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let ext = graph
        .edges
        .iter()
        .filter(|e| e.edge_type != EdgeType::Virtual)
        .count();

    assert_eq!(ext, 8);

    let mut external_momenta: Vec<FourMomentum<F<f128>>> = Vec::with_capacity(7);
    for i in 0..7 {
        external_momenta.push(
            FourMomentum::from_args(
                F(1.),
                F(3.) + F(i as f64),
                F(5.0) + F(i as f64),
                F(4.0) + F(i as f64),
            )
            .cast()
            .higher(),
        );
    }

    assert_eq!(graph.loop_momentum_basis.edge_signatures[0].1.len(), 8);

    assert_eq!(graph.loop_momentum_basis.edge_signatures[0].0.len(), 5);

    let mut loop_momenta: Vec<ThreeMomentum<F<f128>>> = Vec::with_capacity(5);
    for i in 0..5 {
        loop_momenta.push(
            ThreeMomentum::new(
                F(1.) + F(i as f64),
                F(2.) + F(i as f64),
                F(3.) + F(i as f64),
            )
            .cast()
            .higher(),
        );
    }

    let sample = DefaultSample {
        loop_moms: loop_momenta.clone(),
        external_moms: external_momenta.clone(),
        jacobian: F(1.0),
    };

    let ltd_res = graph.evaluate_ltd_expression(&loop_momenta, &external_momenta);
    let cff_res = graph.evaluate_cff_expression(&sample, 0);
    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&cff_res.re, &ltd_res.re, &ltd_comparison_tolerance128);

    assert_approx_eq(&cff_res.im, &ltd_res.im, &ltd_comparison_tolerance128);

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 12);
}

#[test]
fn pytest_scalar_bubble() {
    let (model, amplitude) = load_amplitude_output("TEST_AMPLITUDE_scalar_bubble/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 4);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 4);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 2
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    let p1 = FourMomentum::from_args(F(17. / 19.), F(7. / 11.), F(11. / 13.), F(13. / 17.));
    let k = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.));

    let onshell_energies = graph.compute_onshell_energies(&[k], &[p1, p1]);

    assert_eq!(onshell_energies.len(), 4);

    graph.generate_ltd();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let energy_product = graph.compute_energy_product(&[k], &[p1, p1]);

    let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1]) / energy_product;

    let sample = DefaultSample {
        loop_moms: vec![k],
        external_moms: vec![p1],
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / energy_product;

    let absolute_truth = Complex::new(F(0.), -F(0.052955801144924944));

    assert_approx_eq(&cff_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&cff_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&ltd_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&ltd_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 2);
}

#[test]
fn pytest_scalar_massless_box() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_massless_box/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 8);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 8);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 4
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_ltd();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let p1: FourMomentum<F<f128>> =
        FourMomentum::from_args(79. / 83., 41. / 43., 43. / 47., 47. / 53.)
            .cast()
            .higher();
    let p2: FourMomentum<F<f128>> =
        FourMomentum::from_args(83. / 89., 53. / 59., 59. / 61., 61. / 67.)
            .cast()
            .higher();
    let p3: FourMomentum<F<f128>> =
        FourMomentum::from_args(89. / 97., 67. / 71., 71. / 73., 73. / 79.)
            .cast()
            .higher();

    let externals = vec![p1, p2, p3];

    let absolute_truth = Complex::new(
        F::<f128>::from_f64(0.0),
        F::<f128>::from_f64(-1.5735382832053006e-6),
    );

    let k: ThreeMomentum<F<f128>> = ThreeMomentum::new(
        F::<f128>::from_f64(1.),
        F::<f128>::from_f64(2.),
        F::<f128>::from_f64(3.),
    )
    .cast();

    let loop_moms = vec![k];

    let energy_product = graph.compute_energy_product(&loop_moms, &externals);

    let ltd_res = graph.evaluate_ltd_expression(&loop_moms, &externals) / &energy_product;
    let sample = DefaultSample {
        loop_moms,
        external_moms: externals,
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / &energy_product;

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    assert_approx_eq(
        &cff_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &cff_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &ltd_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &ltd_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 4);
    graph.generate_esurface_data().unwrap();

    let box4_e = [
        FourMomentum::from_args(F(14.0), F(-6.6), F(-40.0), F(0.0)),
        FourMomentum::from_args(F(43.0), F(-15.2), F(-33.0), F(0.0)),
        FourMomentum::from_args(F(17.9), F(50.0), F(-11.8), F(0.0)),
    ];

    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let existing = get_existing_esurfaces(
        esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &box4_e,
        &graph.loop_momentum_basis,
        0,
        F(57.0),
    );

    let edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .collect_vec();

    find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing,
        esurfaces,
        &edge_masses,
        &box4_e,
        0,
    );

    assert_eq!(existing.len(), 4);

    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing,
        esurfaces,
        &edge_masses,
        &box4_e,
        0,
    );

    assert_eq!(maximal_overlap.overlap_groups.len(), 4);
    for overlap in maximal_overlap.overlap_groups.iter() {
        assert_eq!(overlap.existing_esurfaces.len(), 2);
    }
}

#[test]
fn pytest_scalar_double_triangle() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_double_triangle/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 7);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 6);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 2
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_ltd();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let absolute_truth = Complex::new(F(0.00009115369712210525), F(0.)).higher();

    let p1 = FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.)
        .cast()
        .higher();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.))
        .cast()
        .higher();
    let k1 = ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.))
        .cast()
        .higher();

    let loop_moms = vec![k0, k1];
    let externals = vec![p1];

    let energy_product = graph.compute_energy_product(&loop_moms, &externals);

    let ltd_res = graph.evaluate_ltd_expression(&loop_moms, &externals) / &energy_product;
    let sample = DefaultSample {
        loop_moms,
        external_moms: externals,
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / &energy_product;

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    assert_approx_eq(
        &cff_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &cff_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );
    assert_approx_eq(
        &ltd_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 5);
}

#[test]
fn pytest_scalar_mercedes() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_mercedes/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 9);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 7);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 3
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_ltd();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let absolute_truth = Complex::new(F(0.0), F(2.3081733247975594e-13)).higher();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(3., 4., 5.).cast().higher();
    let k1: ThreeMomentum<F<f128>> = ThreeMomentum::new(7., 7., 9.).cast().higher();
    let k2: ThreeMomentum<F<f128>> = ThreeMomentum::new(9., 3., 1.).cast().higher();

    let p1: FourMomentum<F<f128>> = FourMomentum::from_args(1., 12., 13., 14.).cast().higher();
    let p2: FourMomentum<F<f128>> = FourMomentum::from_args(2., 15., 17., 19.).cast().higher();

    let loop_moms = vec![k0, k1, k2];
    let externals = vec![p1, p2];

    let energy_product = graph.compute_energy_product(&loop_moms, &externals);
    let ltd_res = graph.evaluate_ltd_expression(&loop_moms, &externals) / &energy_product;

    let sample = DefaultSample {
        loop_moms,
        external_moms: externals,
        jacobian: F(1.0),
    };
    let cff_res = graph.evaluate_cff_expression(&sample, 0) / &energy_product;

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    assert_approx_eq(
        &cff_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &cff_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 6);
}

#[test]
fn pytest_scalar_triangle_box() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_triangle_box/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 9);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 8);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 3
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_ltd();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let absolute_truth = Complex::new(
        F::<f128>::from_f64(-1.264_354_742_167_213_3e-7),
        F::<f128>::from_f64(0.),
    );

    let p1: FourMomentum<F<f128>> =
        FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.)
            .cast()
            .higher();

    let p2: FourMomentum<F<f128>> = FourMomentum::from_args(2., 15., 17., 19.).cast().higher();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.))
        .cast()
        .higher();
    let k1: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.))
        .cast()
        .higher();

    let loop_moms = vec![k0, k1];
    let externals = vec![p1, p2];

    let energy_product = graph.compute_energy_product(&loop_moms, &externals);
    let ltd_res = graph.evaluate_ltd_expression(&loop_moms, &externals) / &energy_product;
    let sample = DefaultSample {
        loop_moms,
        external_moms: externals,
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / &energy_product;

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    assert_approx_eq(
        &cff_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &cff_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 6);
}

#[test]
fn pytest_scalar_isopod() {
    let (model, amplitude) = load_amplitude_output("TEST_AMPLITUDE_scalar_isopod/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);
    assert!(amplitude.amplitude_graphs[0].graph.edges.len() == 12);
    assert!(amplitude.amplitude_graphs[0].graph.vertices.len() == 10);
    assert!(
        amplitude.amplitude_graphs[0]
            .graph
            .external_connections
            .len()
            == 3
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_ltd();
    graph.generate_cff();
    graph.process_numerator(&model);
    graph.numerator_substitute_model_params(&model);
    // graph.evaluate_model_params(&model);
    graph.process_numerator(&model);

    let absolute_truth = Complex::new(F(0.0), F(-2.9299520787585056e-23)).higher();

    let p1: FourMomentum<F<f128>> =
        FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.)
            .cast()
            .higher();

    let p2: FourMomentum<F<f128>> = FourMomentum::from_args(2., 15., 17., 19.).cast().higher();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.))
        .cast()
        .higher();
    let k1: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.))
        .cast()
        .higher();

    let k2: ThreeMomentum<F<f128>> = ThreeMomentum::new(8., 9., 10.).cast().higher();

    let loop_moms = [k0, k1, k2];
    let externals = [p1, p2];

    let energy_product = graph.compute_energy_product(&loop_moms, &externals);

    let ltd_res = graph.evaluate_ltd_expression(&loop_moms, &externals) / &energy_product;
    let sample = DefaultSample {
        loop_moms: loop_moms.to_vec(),
        external_moms: externals.to_vec(),
        jacobian: F(1.0),
    };

    let cff_res = graph.evaluate_cff_expression(&sample, 0) / &energy_product;

    println!("cff_res = {:+e}", cff_res);
    println!("ltd_res = {:+e}", ltd_res);

    let ltd_comparison_tolerance128 = F::<f128>::from_ff64(LTD_COMPARISON_TOLERANCE);

    println!("cff_res.re = {:+e}", cff_res.re);
    assert_approx_eq(
        &cff_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &cff_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.re,
        &absolute_truth.re,
        &ltd_comparison_tolerance128,
    );

    assert_approx_eq(
        &ltd_res.im,
        &absolute_truth.im,
        &ltd_comparison_tolerance128,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 9);
}

#[test]
fn pytest_scalar_raised_triangle() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_raised_triangle/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 5);
}

#[test]
fn pytest_scalar_hexagon() {
    let (model, amplitude) = load_amplitude_output("TEST_AMPLITUDE_scalar_hexagon/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let kinematics = [
        FourMomentum::from_args(F(24.), F(-21.2), F(71.), F(0.)),
        FourMomentum::from_args(F(50.4), F(15.8), F(-18.8), F(0.)),
        FourMomentum::from_args(F(-0.2), F(46.2), F(8.6), F(0.)),
        -FourMomentum::from_args(F(-33.2), F(2.6), F(-70.8), F(0.)),
        -FourMomentum::from_args(F(-80.), F(-5.6), F(-40.0), F(0.0)),
    ];

    let existing_esurface = get_existing_esurfaces(
        esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &kinematics,
        &graph.loop_momentum_basis,
        0,
        F(75.),
    );

    assert_eq!(existing_esurface.len(), 6);

    let edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .collect_vec();

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurface,
        esurfaces,
        &edge_masses,
        &kinematics,
        0,
    );
    let duration = now.elapsed();
    println!("duration: {}", duration.as_micros());

    assert_eq!(maximal_overlap.overlap_groups.len(), 4);
    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        3
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        3
    );
    assert_eq!(
        maximal_overlap.overlap_groups[2].existing_esurfaces.len(),
        3
    );
    assert_eq!(
        maximal_overlap.overlap_groups[3].existing_esurfaces.len(),
        2
    );

    let hexagon_10_e = [
        FourMomentum::from_args(F(-80.), F(29.), F(-70.), F(0.)),
        FourMomentum::from_args(F(83.5), F(14.0), F(70.0), F(0.0)),
        FourMomentum::from_args(F(88.5), F(6.5), F(-6.), F(0.)),
        -FourMomentum::from_args(F(36.5), F(-71.), F(97.5), F(0.)),
        -FourMomentum::from_args(F(12.5), F(-83.5), F(-57.5), F(0.)),
    ];

    let existing_esurfaces = get_existing_esurfaces(
        esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &hexagon_10_e,
        &graph.loop_momentum_basis,
        0,
        F(88.),
    );

    assert_eq!(existing_esurfaces.len(), 10);

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurfaces,
        esurfaces,
        &edge_masses,
        &hexagon_10_e,
        0,
    );
    let duration = now.elapsed();
    println!("duration: {}", duration.as_micros());

    assert_eq!(maximal_overlap.overlap_groups.len(), 4);
    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        8
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        8
    );
    assert_eq!(
        maximal_overlap.overlap_groups[2].existing_esurfaces.len(),
        7
    );
    assert_eq!(
        maximal_overlap.overlap_groups[3].existing_esurfaces.len(),
        7
    );
}

#[test]
fn pytest_scalar_ltd_topology_c() {
    let (model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_ltd_topology_c/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let kinematics = [
        FourMomentum::from_args(F(9.0), F(0.0), F(0.0), F(8.94427190999916)),
        FourMomentum::from_args(F(9.0), F(0.0), F(0.0), F(-8.94427190999916)),
        -FourMomentum::from_args(
            F(1.83442509122858),
            F(-0.383828222192743),
            F(0.69085529916260),
            F(-1.31653190094982),
        ),
        -FourMomentum::from_args(
            F(5.78920098524940),
            F(-1.80358221330469),
            F(-5.24375913342836),
            F(1.328506453),
        ),
        -FourMomentum::from_args(F(2.82869), F(-1.83886), F(-1.6969477), F(0.8605192)),
    ];

    let _edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re.0 == 0.0 {
                    None
                } else {
                    Some(value)
                }
            } else {
                mass
            }
        })
        .collect_vec();

    let existing_esurfaces = get_existing_esurfaces(
        esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &kinematics,
        &graph.loop_momentum_basis,
        2,
        F(18.),
    );

    let overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurfaces,
        esurfaces,
        &_edge_masses,
        &kinematics,
        0,
    );

    assert_eq!(overlap.overlap_groups.len(), 9);
    assert_eq!(overlap.overlap_groups[0].existing_esurfaces.len(), 22);
    assert_eq!(overlap.overlap_groups[1].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[2].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[3].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[4].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[5].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[6].existing_esurfaces.len(), 20);
    assert_eq!(overlap.overlap_groups[7].existing_esurfaces.len(), 17);
    assert_eq!(overlap.overlap_groups[8].existing_esurfaces.len(), 17);
}

#[test]
fn pytest_scalar_massless_pentabox() {
    let (_model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_massless_pentabox/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let rescaling = F(1.0e-3);
    let kinematics = [
        &FourMomentum::from_args(
            F(5.980_260_048_915_123e2),
            F(0.0),
            F(0.0),
            F(5.724_562_014_045_295e2),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(5.980_260_048_915_123e2),
            F(0.0),
            F(0.0),
            F(-5.724_562_014_045_295e2),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-5.394_473_213_122_507e2),
            F(-1.971_081_698_462_961e2),
            F(-4.416_135_519_343_869e2),
            F(2.250_822_886_064_787e2),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-2.255_538_754_188_549e2),
            F(1.757_868_459_829_899e2),
            F(3.716_353_112_335_996e1),
            F(-1.013_763_093_935_658e2),
        ) * &rescaling,
    ];

    let edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re.0 == 0.0 {
                    None
                } else {
                    Some(value)
                }
            } else {
                mass
            }
        })
        .collect_vec();

    let existing_esurfaces = get_existing_esurfaces(
        &graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &kinematics,
        &graph.loop_momentum_basis,
        0,
        F(1.0),
    );

    assert_eq!(existing_esurfaces.len(), 17);

    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurfaces,
        &graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        &edge_masses,
        &kinematics,
        0,
    );

    assert_eq!(maximal_overlap.overlap_groups.len(), 3);

    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        16
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        13
    );
    assert_eq!(
        maximal_overlap.overlap_groups[2].existing_esurfaces.len(),
        13
    );
}

#[test]
fn pytest_scalar_massless_3l_pentabox() {
    let (_model, amplitude) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_massless_3l_pentabox/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let rescaling = F(1.0e0);
    let kinematics = [
        &FourMomentum::from_args(
            F(0.149500000000000E+01),
            F(0.000000000000000E+00),
            F(0.000000000000000E+00),
            F(0.149165176901313E+01),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(0.150500000000000E+01),
            F(0.000000000000000E+00),
            F(0.000000000000000E+00),
            F(-0.149165176901313E+01),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-0.126041949101381e+01),
            F(-0.452362952912639e+00),
            F(-0.101350243653045e+01),
            F(0.516563513332600e+00),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-0.105098730574850e+01),
            F(0.489324061520790e-01),
            F(0.928212188578101e+00),
            F(-0.283905035967510e+00),
        ) * &rescaling,
    ];

    let edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re.0 == 0.0 {
                    None
                } else {
                    Some(value)
                }
            } else {
                mass
            }
        })
        .collect_vec();

    let existing_esurfaces = get_existing_esurfaces(
        &graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &kinematics,
        &graph.loop_momentum_basis,
        0,
        F(1.0),
    );

    assert_eq!(graph.loop_momentum_basis.basis.len(), 3);
    assert_eq!(existing_esurfaces.len(), 28);

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurfaces,
        &graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        &edge_masses,
        &kinematics,
        0,
    );
    let _elapsed = now.elapsed();

    assert_eq!(maximal_overlap.overlap_groups.len(), 2);
    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        27
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        26
    );
}

#[test]
fn pytest_lbl_box() {
    let (model, amplitude) = load_amplitude_output("TEST_AMPLITUDE_lbl_box/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.process_numerator(&model);
    // println!();

    // for v in graph
    //     .vertices
    //     .iter()
    //     .filter(|v| v.vertex_info.get_type() == "interacton_vertex_info")
    // {
    //     println!("vertex: {}", v.name);

    //     println!("From edges: ");
    //     for (i, e) in v.edges.clone().iter().enumerate() {
    //         println!("{} : {:?}", i, graph.edges[*e].particle.name)
    //     }
    //     println!("From vertex info: ");
    //     if let VertexInfo::InteractonVertexInfo(s) = &v.vertex_info {
    //         s.vertex_rule
    //             .particles
    //             .iter()
    //             .enumerate()
    //             .for_each(|(i, p)| println!("{} : {:?}", i, p.name));
    //     }
    // }

    // for e in graph.edges.iter() {
    //     println!("edge: {}", e.name);
    //     for v in e.vertices {
    //         if e.is_incoming_to(v) {
    //             println!("incoming to vertex: {}", graph.vertices[v].name);
    //         } else {
    //             println!("outgoing to vertex: {}", graph.vertices[v].name);
    //         }
    //         let i = graph.vertices[v]
    //             .edges
    //             .iter()
    //             .enumerate()
    //             .filter(|(_, &es)| es == graph.get_edge_position(&e.name).unwrap())
    //             .map(|(i, _)| i)
    //             .collect::<Vec<usize>>();

    //         if let VertexInfo::InteractonVertexInfo(s) = &graph.vertices[v].vertex_info {
    //             let p = &s.vertex_rule.particles[i[0]];
    //             println!("{:?}", p.name);
    //         }
    //     }
    // }

    // println!("{}", graph.derived_data.numerator.unwrap().expression);
}

#[test]
#[allow(non_snake_case)]
fn pytest_physical_3L_6photons_topology_A_inspect() {
    env_logger::init();
    let (model, amplitude) = load_amplitude_output(
        "TEST_AMPLITUDE_physical_3L_6photons_topology_A/GL_OUTPUT",
        true,
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_numerator();

    graph.process_numerator(&model);

    println!(
        "{}",
        graph.derived_data.numerator.as_ref().unwrap().expression
    );

    println!(
        "{}",
        graph.derived_data.numerator.unwrap().network.unwrap().dot()
    );
    // let mut onlycolor = Numerator {
    //     expression: Atom::parse("T(aind(coad(8,9),cof(3,8),coaf(3,7)))*T(aind(coad(8,14),cof(3,13),coaf(3,12)))*T(aind(coad(8,21),cof(3,20),coaf(3,19)))*T(aind(coad(8,26),cof(3,25),coaf(3,24)))*id(aind(coaf(3,3),cof(3,4)))*id(aind(coaf(3,4),cof(3,24)))*id(aind(coaf(3,5),cof(3,6)))*id(aind(coaf(3,6),cof(3,3)))*id(aind(coaf(3,8),cof(3,5)))*id(aind(coaf(3,10),cof(3,11)))*id(aind(coaf(3,11),cof(3,7)))*id(aind(coaf(3,13),cof(3,10)))*id(aind(coaf(3,15),cof(3,16)))*id(aind(coaf(3,16),cof(3,12)))*id(aind(coaf(3,17),cof(3,18)))*id(aind(coaf(3,18),cof(3,15)))*id(aind(coaf(3,20),cof(3,17)))*id(aind(coaf(3,22),cof(3,23)))*id(aind(coaf(3,23),cof(3,19)))*id(aind(coaf(3,25),cof(3,22)))*id(aind(coad(8,21),coad(8,9)))*id(aind(coad(8,26),coad(8,14)))").unwrap(),
    //     network: None,
    //     const_map: AHashMap::new(),
    // };

    // onlycolor.fill_network();
    // println!("{}", onlycolor.network.as_ref().unwrap().dot());

    // onlycolor.process_color_simple();
    // println!("{}", onlycolor.expression);
    // let a = graph
    //     .derived_data
    //     .numerator
    //     .as_ref()
    //     .unwrap()
    //     .network
    //     .as_ref()
    //     .unwrap();

    // println!("{}", a.dot_nodes());
}
