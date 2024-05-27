#![allow(unused_imports)]
use crate::cff::esurface::{
    get_existing_esurfaces, Esurface, EsurfaceID, ExistingEsurfaceId, ExistingEsurfaces,
};
use crate::cff::generation::generate_cff_expression;
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::graph::{Edge, EdgeType, HasVertexInfo, InteractionVertexInfo, VertexInfo};
use crate::model::Model;
use crate::subtraction::overlap::{self, find_center, find_maximal_overlap};
use crate::subtraction::static_counterterm;
use crate::utils::{assert_approx_eq, compute_momentum, upgrade_lorentz_vector};
use colored::Colorize;
use itertools::Itertools;
use libc::__c_anonymous_ptrace_syscall_info_exit;
use lorentz_vector::LorentzVector;
use num::Complex;
use petgraph::algo::greedy_matching;
use petgraph::graph;
use rayon::prelude::IndexedParallelIterator;
use serde;
use statrs::function::evaluate;
use std::fs::File;
use std::path::Path;
use std::{clone, env};
use symbolica;

#[allow(unused)]
const LTD_COMPARISON_TOLERANCE: f64 = 1.0e-12;

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
    use num::complex::Complex64;
    use rayon::prelude::IndexedParallelIterator;
    use smartstring::SmartString;

    use crate::{cff::esurface::get_existing_esurfaces, graph::EdgeType, observables::AFBSettings};

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

        let p1 = LorentzVector::from_args(1.0, 3.0, 4.0, 5.0);
        let p2 = LorentzVector::from_args(1.0, 6.0, 7.0, 8.0);

        let k = LorentzVector::from_args(0.0, 1.0, 2.0, 3.0);

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

        let absolute_truth = Complex::new(0.0, 4.531238289663331e-6);

        graph.generate_ltd();

        let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1, p2]) / energy_product;

        println!("cff_res = {:+e}", cff_res);
        println!("ltd_res = {:+e}", ltd_res);
        println!("cltd_m = {:+e}", absolute_truth);

        assert_approx_eq(cff_res.re, absolute_truth.re, LTD_COMPARISON_TOLERANCE);
        assert_approx_eq(cff_res.im, absolute_truth.im, LTD_COMPARISON_TOLERANCE);
        assert_approx_eq(ltd_res.re, absolute_truth.re, LTD_COMPARISON_TOLERANCE);
        assert_approx_eq(ltd_res.im, absolute_truth.im, LTD_COMPARISON_TOLERANCE);

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
            2.0,
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
    assert!(!cff.esurfaces.is_empty());

    // println!("basis size: {:?}", graph.loop_momentum_basis.basis);
    graph.generate_loop_momentum_bases();

    let k1 = LorentzVector::from_args(0., 2. / 3., 3. / 5., 5. / 7.);
    let k2 = LorentzVector::from_args(0., 7. / 11., 11. / 13., 13. / 17.);
    let k3 = LorentzVector::from_args(0., 17. / 19., 19. / 23., 23. / 29.);
    let k4 = LorentzVector::from_args(0., 29. / 31., 31. / 37., 37. / 41.);
    let p1 = LorentzVector::from_args(79. / 83., 41. / 43., 43. / 47., 47. / 53.);
    let p2 = LorentzVector::from_args(83. / 89., 53. / 59., 59. / 61., 61. / 67.);
    let p3 = LorentzVector::from_args(89. / 97., 67. / 71., 71. / 73., 73. / 79.);

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
            .map(|s| compute_momentum(s, &momenta_in_basis, &[p1, p2, p3]))
            .collect_vec();
        assert_eq!(emr.len(), new_emr.len());

        for (e1, e2) in emr.iter().zip(new_emr.iter()) {
            assert_approx_eq(e1.t, e2.t, LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(e1.x, e2.x, LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(e1.y, e2.y, LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(e1.z, e2.z, LTD_COMPARISON_TOLERANCE);
        }
    }

    //println!("lmb consistency check passed");

    let k1_f128 = upgrade_lorentz_vector(&k1);
    let k2_f128 = upgrade_lorentz_vector(&k2);
    let k3_f128 = upgrade_lorentz_vector(&k3);
    let k4_f128 = upgrade_lorentz_vector(&k4);
    let p1_f128 = upgrade_lorentz_vector(&p1);
    let p2_f128 = upgrade_lorentz_vector(&p2);
    let p3_f128 = upgrade_lorentz_vector(&p3);

    let loop_moms_f128 = [k1_f128, k2_f128, k3_f128, k4_f128];
    let externals_f128 = [p1_f128, p2_f128, p3_f128];

    let energy_product = graph.compute_energy_product(&loop_moms_f128, &externals_f128);

    let cff_res = graph.evaluate_cff_expression(&loop_moms_f128, &externals_f128) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&loop_moms_f128, &externals_f128) / energy_product;

    let absolute_truth = Complex::new(0.000019991301832169422, 0.0);

    assert_approx_eq(
        cff_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        cff_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        ltd_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 12);
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

    let k1 = LorentzVector::from_args(0., 2. / 3., 3. / 5., 5. / 7.);
    let k2 = LorentzVector::from_args(0., 7. / 11., 11. / 13., 13. / 17.);

    let p1 = LorentzVector::from_args(0., 0., 0., 0.);

    let absolute_truth = Complex::new(0.24380172488169907, 0.);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_loop_momentum_bases();
    graph.generate_cff();
    graph.generate_ltd();

    let energy_product = graph.compute_energy_product(&[k1, k2], &[p1]);

    let cff_res = graph.evaluate_cff_expression(&[k1, k2], &[p1]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k1, k2], &[p1]) / energy_product;

    println!("cff_res = {:+e}", cff_res);
    println!("ltd_res = {:+e}", ltd_res);

    assert_approx_eq(cff_res.re, absolute_truth.re, LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(cff_res.im, absolute_truth.im, LTD_COMPARISON_TOLERANCE);
    println!("cff correct");
    assert_approx_eq(ltd_res.re, absolute_truth.re, LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(ltd_res.im, absolute_truth.im, LTD_COMPARISON_TOLERANCE);
    println!("ltd correct");

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 3);
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

    let externals = (0..3)
        .map(|i| {
            upgrade_lorentz_vector(&LorentzVector::from_args(
                1.0,
                i as f64 + 2.0,
                i as f64 + 3.0,
                i as f64 + 4.0,
            ))
        })
        .collect_vec();

    let loop_moms = (0..6)
        .map(|i| {
            upgrade_lorentz_vector(&LorentzVector::from_args(
                1.0,
                i as f64 - 2.0,
                i as f64 + 3.5,
                i as f64 + 4.5001,
            ))
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
        cff_res.re.into(),
        ltd_res.re.into(),
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        cff_res.im.into(),
        ltd_res.im.into(),
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = amplitude.amplitude_graphs[0]
        .graph
        .group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 17);
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

    let mut external_momenta = Vec::with_capacity(7);
    for i in 0..7 {
        external_momenta.push(upgrade_lorentz_vector(&LorentzVector::from_args(
            1.0,
            3.0 + i as f64,
            5.0 + i as f64,
            4.0 + i as f64,
        )));
    }

    assert_eq!(graph.loop_momentum_basis.edge_signatures[0].1.len(), 8);

    assert_eq!(graph.loop_momentum_basis.edge_signatures[0].0.len(), 5);

    let mut loop_momenta = Vec::with_capacity(5);
    for i in 0..5 {
        loop_momenta.push(upgrade_lorentz_vector(&LorentzVector::from_args(
            0.0,
            1.5 + i as f64,
            2.5 + i as f64,
            4.5 + i as f64,
        )));
    }

    let ltd_res = graph.evaluate_ltd_expression(&loop_momenta, &external_momenta);
    let cff_res = graph.evaluate_cff_expression(&loop_momenta, &external_momenta);
    assert_approx_eq(
        cff_res.re.into(),
        ltd_res.re.into(),
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        cff_res.im.into(),
        ltd_res.im.into(),
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 12);
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

    let p1 = LorentzVector::from_args(17. / 19., 7. / 11., 11. / 13., 13. / 17.);
    let k = LorentzVector::from_args(0.0, 2. / 3., 3. / 5., 5. / 7.);

    let onshell_energies = graph.compute_onshell_energies(&[k], &[p1, p1]);

    assert_eq!(onshell_energies.len(), 4);

    graph.generate_ltd();
    graph.generate_cff();

    let energy_product = graph.compute_energy_product(&[k], &[p1, p1]);

    let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1]) / energy_product;
    let cff_res = graph.evaluate_cff_expression(&[k], &[p1]) / energy_product;

    let absolute_truth = Complex::new(0., -0.052955801144924944);

    assert_approx_eq(cff_res.re, absolute_truth.re, LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(cff_res.im, absolute_truth.im, LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(ltd_res.re, absolute_truth.re, LTD_COMPARISON_TOLERANCE);
    assert_approx_eq(ltd_res.im, absolute_truth.im, LTD_COMPARISON_TOLERANCE);

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 2);
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

    let p1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        79. / 83.,
        41. / 43.,
        43. / 47.,
        47. / 53.,
    ));
    let p2 = upgrade_lorentz_vector(&LorentzVector::from_args(
        83. / 89.,
        53. / 59.,
        59. / 61.,
        61. / 67.,
    ));
    let p3 = upgrade_lorentz_vector(&LorentzVector::from_args(
        89. / 97.,
        67. / 71.,
        71. / 73.,
        73. / 79.,
    ));

    let externals = [p1, p2, p3];

    let absolute_truth = Complex::new(0.0, -1.5735382832053006e-6);

    let k = upgrade_lorentz_vector(&LorentzVector::from_args(0.0, 1.0, 2.0, 3.0));

    let energy_product = graph.compute_energy_product(&[k], &externals);

    let cff_res = graph.evaluate_cff_expression(&[k], &externals) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k], &externals) / energy_product;

    assert_approx_eq(
        cff_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        cff_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        ltd_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        ltd_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 4);
    graph.generate_esurface_data().unwrap();

    let box4_e = [
        LorentzVector::from_args(14.0, -6.6, -40.0, 0.0),
        LorentzVector::from_args(43.0, -15.2, -33.0, 0.0),
        LorentzVector::from_args(17.9, 50.0, -11.8, 0.0),
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
        57.0,
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
    );

    assert_eq!(existing.len(), 4);

    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing,
        esurfaces,
        &edge_masses,
        &box4_e,
    );

    assert_eq!(maximal_overlap.len(), 4);
    for overlap in maximal_overlap.iter() {
        assert_eq!(overlap.0.len(), 2);
    }
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

    let absolute_truth = Complex::new(0.00009115369712210525, 0.0);

    let p1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        53. / 59.,
        41. / 43.,
        43. / 47.,
        47. / 53.,
    ));

    let k0 = upgrade_lorentz_vector(&LorentzVector::from_args(0., 2. / 3., 3. / 5., 5. / 7.));
    let k1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        0.,
        7. / 11.,
        11. / 13.,
        13. / 17.,
    ));

    let energy_product = graph.compute_energy_product(&[k0, k1], &[p1]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1], &[p1]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1], &[p1]) / energy_product;

    assert_approx_eq(
        cff_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        cff_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );
    assert_approx_eq(
        ltd_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 5);
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

    let absolute_truth = Complex::new(0.0, 2.3081733247975594e-13);

    let k0 = upgrade_lorentz_vector(&LorentzVector::from_args(0.0, 3., 4., 5.));
    let k1 = upgrade_lorentz_vector(&LorentzVector::from_args(0.0, 7., 7., 9.));
    let k2 = upgrade_lorentz_vector(&LorentzVector::from_args(0.0, 9., 3., 1.));

    let p1 = upgrade_lorentz_vector(&LorentzVector::from_args(1.0, 12., 13., 14.));
    let p2 = upgrade_lorentz_vector(&LorentzVector::from_args(2.0, 15., 17., 19.));

    let energy_product = graph.compute_energy_product(&[k0, k1, k2], &[p1, p2]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;

    assert_approx_eq(
        cff_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        cff_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 6);
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

    let absolute_truth = Complex::new(-1.264_354_742_167_213_3e-7, 0.0);

    let p1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        53. / 59.,
        41. / 43.,
        43. / 47.,
        47. / 53.,
    ));

    let p2 = upgrade_lorentz_vector(&LorentzVector::from_args(2., 15., 17., 19.));

    let k0 = upgrade_lorentz_vector(&LorentzVector::from_args(0., 2. / 3., 3. / 5., 5. / 7.));
    let k1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        0.,
        7. / 11.,
        11. / 13.,
        13. / 17.,
    ));

    let energy_product = graph.compute_energy_product(&[k0, k1], &[p1, p2]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1], &[p1, p2]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1], &[p1, p2]) / energy_product;

    assert_approx_eq(
        cff_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        cff_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 6);
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

    let absolute_truth = Complex::new(0.0, -2.9299520787585056e-23);

    let p1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        53. / 59.,
        41. / 43.,
        43. / 47.,
        47. / 53.,
    ));

    let p2 = upgrade_lorentz_vector(&LorentzVector::from_args(2., 15., 17., 19.));

    let k0 = upgrade_lorentz_vector(&LorentzVector::from_args(0., 2. / 3., 3. / 5., 5. / 7.));
    let k1 = upgrade_lorentz_vector(&LorentzVector::from_args(
        0.,
        7. / 11.,
        11. / 13.,
        13. / 17.,
    ));

    let k2 = upgrade_lorentz_vector(&LorentzVector::from_args(0., 8., 9., 10.));

    let energy_product = graph.compute_energy_product(&[k0, k1, k2], &[p1, p2]);

    let cff_res = graph.evaluate_cff_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;
    let ltd_res = graph.evaluate_ltd_expression(&[k0, k1, k2], &[p1, p2]) / energy_product;

    println!("cff_res = {:+e}", cff_res);
    println!("ltd_res = {:+e}", ltd_res);

    assert_approx_eq(
        cff_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        cff_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.re.into(),
        absolute_truth.re,
        LTD_COMPARISON_TOLERANCE,
    );

    assert_approx_eq(
        ltd_res.im.into(),
        absolute_truth.im,
        LTD_COMPARISON_TOLERANCE,
    );

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 9);
}

#[test]
#[ignore]
fn pytest_raised_triangle() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 5);
}

#[test]
#[ignore]
fn pytest_hexagon() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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
        LorentzVector::from_args(24., -21.2, 71., 0.),
        LorentzVector::from_args(50.4, 15.8, -18.8, 0.),
        LorentzVector::from_args(-0.2, 46.2, 8.6, 0.),
        -LorentzVector::from_args(-33.2, 2.6, -70.8, 0.),
        -LorentzVector::from_args(-80., -5.6, -40.0, 0.0),
    ];

    let existing_esurface = get_existing_esurfaces(
        esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &kinematics,
        &graph.loop_momentum_basis,
        0,
        75.,
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
    );
    let duration = now.elapsed();
    println!("duration: {}", duration.as_micros());

    assert_eq!(maximal_overlap.len(), 4);
    assert_eq!(maximal_overlap[0].0.len(), 3);
    assert_eq!(maximal_overlap[1].0.len(), 3);
    assert_eq!(maximal_overlap[2].0.len(), 3);
    assert_eq!(maximal_overlap[3].0.len(), 2);

    let hexagon_10_e = [
        LorentzVector::from_args(-80., 29., -70., 0.),
        LorentzVector::from_args(83.5, 14.0, 70.0, 0.0),
        LorentzVector::from_args(88.5, 6.5, -6., 0.),
        -LorentzVector::from_args(36.5, -71., 97.5, 0.),
        -LorentzVector::from_args(12.5, -83.5, -57.5, 0.),
    ];

    let existing_esurfaces = get_existing_esurfaces(
        esurfaces,
        graph.derived_data.esurface_derived_data.as_ref().unwrap(),
        &hexagon_10_e,
        &graph.loop_momentum_basis,
        0,
        88.,
    );

    assert_eq!(existing_esurfaces.len(), 10);

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurfaces,
        esurfaces,
        &edge_masses,
        &hexagon_10_e,
    );
    let duration = now.elapsed();
    println!("duration: {}", duration.as_micros());

    assert_eq!(maximal_overlap.len(), 4);
    assert_eq!(maximal_overlap[0].0.len(), 8);
    assert_eq!(maximal_overlap[1].0.len(), 8);
    assert_eq!(maximal_overlap[2].0.len(), 7);
    assert_eq!(maximal_overlap[3].0.len(), 7);
}

#[test]
#[ignore]
fn pytest_topology_c() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

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
        LorentzVector::from_args(9.0, 0.0, 0.0, 8.94427190999916),
        LorentzVector::from_args(9.0, 0.0, 0.0, -8.94427190999916),
        -LorentzVector::from_args(
            1.83442509122858,
            -0.383828222192743,
            0.69085529916260,
            -1.31653190094982,
        ),
        -LorentzVector::from_args(
            5.78920098524940,
            -1.80358221330469,
            -5.24375913342836,
            1.328506453,
        ),
        -LorentzVector::from_args(2.82869, -1.83886, -1.6969477, 0.8605192),
    ];

    let _edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re == 0.0 {
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
        18.,
    );

    let overlap = find_maximal_overlap(
        &graph.loop_momentum_basis,
        &existing_esurfaces,
        esurfaces,
        &_edge_masses,
        &kinematics,
    );

    assert_eq!(overlap.len(), 9);
    assert_eq!(overlap[0].0.len(), 22);
    assert_eq!(overlap[1].0.len(), 21);
    assert_eq!(overlap[2].0.len(), 21);
    assert_eq!(overlap[3].0.len(), 21);
    assert_eq!(overlap[4].0.len(), 21);
    assert_eq!(overlap[5].0.len(), 21);
    assert_eq!(overlap[6].0.len(), 20);
    assert_eq!(overlap[7].0.len(), 17);
    assert_eq!(overlap[8].0.len(), 17);
}

#[test]
#[ignore]
fn pytest_massless_pentabox() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (_model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let rescaling = 1.0e-3;
    let kinematics = [
        LorentzVector::from_args(5.980_260_048_915_123e2, 0.0, 0.0, 5.724_562_014_045_295e2)
            * rescaling,
        LorentzVector::from_args(5.980_260_048_915_123e2, 0.0, 0.0, -5.724_562_014_045_295e2)
            * rescaling,
        LorentzVector::from_args(
            -5.394_473_213_122_507e2,
            -1.971_081_698_462_961e2,
            -4.416_135_519_343_869e2,
            2.250_822_886_064_787e2,
        ) * rescaling,
        LorentzVector::from_args(
            -2.255_538_754_188_549e2,
            1.757_868_459_829_899e2,
            3.716_353_112_335_996e1,
            -1.013_763_093_935_658e2,
        ) * rescaling,
    ];

    let edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re == 0.0 {
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
        1.0,
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
    );

    assert_eq!(maximal_overlap.len(), 3);

    assert_eq!(maximal_overlap[0].0.len(), 16);
    assert_eq!(maximal_overlap[1].0.len(), 13);
    assert_eq!(maximal_overlap[2].0.len(), 13);
}

#[test]
#[ignore]
fn pytest_massless_3l_pentabox() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let (_model, amplitude) =
        load_amplitude_output(&env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap());

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let rescaling = 1.0e0;
    let kinematics = [
        LorentzVector::from_args(
            0.149500000000000E+01,
            0.000000000000000E+00,
            0.000000000000000E+00,
            0.149165176901313E+01,
        ) * rescaling,
        LorentzVector::from_args(
            0.150500000000000E+01,
            0.000000000000000E+00,
            0.000000000000000E+00,
            -0.149165176901313E+01,
        ) * rescaling,
        LorentzVector::from_args(
            -0.126041949101381e+01,
            -0.452362952912639e+00,
            -0.101350243653045e+01,
            0.516563513332600e+00,
        ) * rescaling,
        LorentzVector::from_args(
            -0.105098730574850e+01,
            0.489324061520790e-01,
            0.928212188578101e+00,
            -0.283905035967510e+00,
        ) * rescaling,
    ];

    let edge_masses = graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re == 0.0 {
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
        1.0,
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
    );
    let _elapsed = now.elapsed();

    assert_eq!(maximal_overlap.len(), 2);
    assert_eq!(maximal_overlap[0].0.len(), 27);
    assert_eq!(maximal_overlap[1].0.len(), 26);
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
