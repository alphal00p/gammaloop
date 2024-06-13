#![allow(unused_imports)]
use crate::cff::generate_cff_expression;
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::graph::{Edge, EdgeType, HasVertexInfo, InteractionVertexInfo, VertexInfo};
use crate::model::Model;
use crate::momentum::{FourMomentum, ThreeMomentum};
use crate::utils::{assert_approx_eq, compute_momentum, compute_three_momentum_from_four};
use crate::utils::{f128, F};
use colored::Colorize;
use itertools::{FormatWith, Itertools};
use lorentz_vector::LorentzVector;
use rayon::prelude::IndexedParallelIterator;
use serde;
use std::fs::File;
use std::path::Path;
use std::{clone, env};
use symbolica;
use symbolica::domains::float::Complex;

#[allow(unused)]
const LTD_COMPARISON_TOLERANCE: F<f64> = F(1.0e-12);

pub fn load_amplitude_output(output_path: &str) -> (Model, Amplitude) {
    let path = Path::new(output_path);
    let output_meta_data: OutputMetaData =
        serde_yaml::from_reader(File::open(path.join("output_metadata.yaml")).unwrap()).unwrap();
    assert_eq!(output_meta_data.output_type, OutputType::Amplitudes);
    assert_eq!(output_meta_data.contents.len(), 1);

    let model = Model::from_file(String::from(
        path.join(format!(
            "sources/model/{}.yaml",
            output_meta_data.model_name
        ))
        .to_str()
        .unwrap(),
    ))
    .unwrap();
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
    use symbolica::domains::float::Complex;

    use crate::{
        graph::EdgeType,
        momentum::{FourMomentum, ThreeMomentum},
        observables::AFBSettings,
        utils::F,
    };

    use super::*;

    #[test]
    #[ignore] // Important since this test will only run successfully when called from with pytest where the massless_triangle_generation fixture will be run
    fn pytest_massless_scalar_triangle() {
        assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

        let (model, amplitude) =
            load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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
        assert_eq!(
            graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .terms
                .len(),
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

        let cff_res = graph.evaluate_cff_expression(&[k], &[p1, p2]) / energy_product;

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
        // TODO: @Mathijs, you can put your own checks there
    }
}

#[test]
#[ignore] // Important since this test will only run successfully when called from with pytest where the massless_triangle_generation fixture will be run
fn pytest_scalar_fishnet_2x2() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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
    assert!(!cff.terms.is_empty());

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

    let k1_f128: ThreeMomentum<F<f128>> = k1.cast();
    let k2_f128: ThreeMomentum<F<f128>> = k2.cast();
    let k3_f128: ThreeMomentum<F<f128>> = k3.cast();
    let k4_f128: ThreeMomentum<F<f128>> = k4.cast();
    let p1_f128: FourMomentum<F<f128>> = p1.cast();
    let p2_f128: FourMomentum<F<f128>> = p2.cast();
    let p3_f128: FourMomentum<F<f128>> = p3.cast();

    let loop_moms_f128 = [k1_f128, k2_f128, k3_f128, k4_f128];
    let externals_f128 = [p1_f128, p2_f128, p3_f128];

    let energy_product = graph.compute_energy_product(&loop_moms_f128, &externals_f128);

    let cff_res = graph.evaluate_cff_expression(&loop_moms_f128, &externals_f128) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&loop_moms_f128, &externals_f128) / energy_product;

    let absolute_truth: Complex<F<f64>> = Complex::new(F(0.000019991301832169422), F(0.0));

    assert_approx_eq(
        &cff_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &cff_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &ltd_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
    // TODO: @Mathijs, you can put your own checks there
}

#[test]
#[ignore]
fn pytest_scalar_sunrise() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let p1 = FourMomentum::from_args(F(0.), F(0.), F(0.), F(0.)).into();

    let absolute_truth = Complex::new(F(0.24380172488169907), F(0.));

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_loop_momentum_bases();
    graph.generate_cff();
    graph.generate_ltd();

    let energy_product = graph.compute_energy_product(&[k1, k2], &[p1]);

    let cff_res = graph.evaluate_cff_expression(&[k1, k2], &[p1]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k1, k2], &[p1]) / energy_product;

    println!("cff_res = {:+e}", cff_res);
    println!("ltd_res = {:+e}", ltd_res);

    assert_approx_eq(&cff_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&cff_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
    println!("cff correct");
    assert_approx_eq(&ltd_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&ltd_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
    println!("ltd correct");
}

#[test]
#[ignore] // Important since this test will only run successfully when called from with pytest where the massless_triangle_generation fixture will be run
fn pytest_scalar_fishnet_2x3() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, mut amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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
        })
        .collect_vec();

    // let before_cff = std::time::Instant::now();

    let cff_res = amplitude.amplitude_graphs[0]
        .graph
        .evaluate_cff_expression(&loop_moms, &externals);

    // let cff_duration = before_cff.elapsed();
    // println!("cff_duration: {}", cff_duration.as_micros());

    // let before_ltd = std::time::Instant::now();
    let ltd_res = amplitude.amplitude_graphs[0]
        .graph
        .evaluate_ltd_expression(&loop_moms, &externals);

    // let ltd_duration = before_ltd.elapsed();
    // println!("ltd_duration: {}", ltd_duration.as_micros());

    assert_approx_eq(
        &cff_res.re.into(),
        &ltd_res.re.into(),
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &cff_res.im.into(),
        &ltd_res.im.into(),
        &LTD_COMPARISON_TOLERANCE,
    );
    // TODO: @Mathijs, you can put your own checks there
}

#[test]
#[ignore]
fn pytest_scalar_cube() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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
            .cast(),
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
            .cast(),
        );
    }

    let ltd_res = graph.evaluate_ltd_expression(&loop_momenta, &external_momenta);
    let cff_res = graph.evaluate_cff_expression(&loop_momenta, &external_momenta);
    assert_approx_eq(
        &cff_res.re.into(),
        &ltd_res.re.into(),
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &cff_res.im.into(),
        &ltd_res.im.into(),
        &LTD_COMPARISON_TOLERANCE,
    );
}

#[test]
#[ignore]
fn pytest_scalar_bubble() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let energy_product = graph.compute_energy_product(&[k], &[p1, p1]);

    let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1]) / energy_product;
    let cff_res = graph.evaluate_cff_expression(&[k], &[p1]) / energy_product;

    let absolute_truth = Complex::new(F(0.), -F(0.052955801144924944));

    assert_approx_eq(&cff_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&cff_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&ltd_res.re, &absolute_truth.re, &LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(&ltd_res.im, &absolute_truth.im, &LTD_COMPARISON_TOLERANCE);
}

#[test]
#[ignore]
fn pytest_massless_scalar_box() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let p1: FourMomentum<F<f128>> =
        FourMomentum::from_args(79. / 83., 41. / 43., 43. / 47., 47. / 53.).cast();
    let p2: FourMomentum<F<f128>> =
        FourMomentum::from_args(83. / 89., 53. / 59., 59. / 61., 61. / 67.).cast();
    let p3: FourMomentum<F<f128>> =
        FourMomentum::from_args(89. / 97., 67. / 71., 71. / 73., 73. / 79.).cast();

    let externals = [p1, p2, p3];

    let absolute_truth = Complex::new(F(0.0), F(-1.5735382832053006e-6));

    let k: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(1.), F(2.), F(3.)).cast();

    let energy_product = graph.compute_energy_product(&[k], &externals);

    let cff_res = graph.evaluate_cff_expression(&[k], &externals) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k], &externals) / energy_product;

    assert_approx_eq(
        &cff_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &cff_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &ltd_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &ltd_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
}

#[test]
#[ignore]
fn pytest_scalar_double_triangle() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let absolute_truth = Complex::new(F(0.00009115369712210525), F(0.));

    let p1 = FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.).cast();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.)).cast();
    let k1 = ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.)).cast();

    let energy_product = graph.compute_energy_product(&[k0, k1], &[p1]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1], &[p1]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1], &[p1]) / energy_product;

    assert_approx_eq(
        &cff_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &cff_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        &ltd_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
}

#[test]
#[ignore]
fn pytest_scalar_mercedes() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let absolute_truth = Complex::new(F(0.0), F(2.3081733247975594e-13));

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(3., 4., 5.).cast();
    let k1: ThreeMomentum<F<f128>> = ThreeMomentum::new(7., 7., 9.).cast();
    let k2: ThreeMomentum<F<f128>> = ThreeMomentum::new(9., 3., 1.).cast();

    let p1: FourMomentum<F<f128>> = FourMomentum::from_args(1., 12., 13., 14.).cast();
    let p2: FourMomentum<F<f128>> = FourMomentum::from_args(2., 15., 17., 19.).cast();

    let energy_product = graph.compute_energy_product(&[k0, k1, k2], &[p1, p2]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;

    assert_approx_eq(
        &cff_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &cff_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
}

#[test]
#[ignore]
fn pytest_scalar_triangle_box() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let absolute_truth = Complex::new(F(-1.264_354_742_167_213_3e-7), F(0.));

    let p1: FourMomentum<F<f128>> =
        FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.).cast();

    let p2: FourMomentum<F<f128>> = FourMomentum::from_args(2., 15., 17., 19.).cast();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.)).cast();
    let k1: ThreeMomentum<F<f128>> =
        ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.)).cast();

    let energy_product = graph.compute_energy_product(&[k0, k1], &[p1, p2]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1], &[p1, p2]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1], &[p1, p2]) / energy_product;

    assert_approx_eq(
        &cff_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &cff_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
}

#[test]
#[ignore]
fn pytest_scalar_isopod() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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

    let absolute_truth = Complex::new(F(0.0), F(-2.9299520787585056e-23));

    let p1: FourMomentum<F<f128>> =
        FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.).cast();

    let p2: FourMomentum<F<f128>> = FourMomentum::from_args(2., 15., 17., 19.).cast();

    let k0: ThreeMomentum<F<f128>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.)).cast();
    let k1: ThreeMomentum<F<f128>> =
        ThreeMomentum::new(F(2. / 11.), F(11. / 13.), F(13. / 17.)).cast();

    let k2: ThreeMomentum<F<f128>> = ThreeMomentum::new(8., 9., 10.).cast();

    let energy_product = graph.compute_energy_product(&[k0, k1, k2], &[p1, p2]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;

    println!("cff_res = {:+e}", cff_res);
    println!("ltd_res = {:+e}", ltd_res);

    assert_approx_eq(
        &cff_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &cff_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.re.into(),
        &absolute_truth.re,
        &LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        &ltd_res.im.into(),
        &absolute_truth.im,
        &LTD_COMPARISON_TOLERANCE,
    );
}

#[test]
#[ignore]

fn pytest_lbl_box() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_numerator(&model);
    println!();

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

    println!("{}", graph.derived_data.numerator.unwrap());
}
