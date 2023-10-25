#![allow(unused_imports)]
use crate::cff::generate_cff_expression;
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::graph::EdgeType;
use crate::model::Model;
use crate::utils::{assert_approx_eq, compute_momentum};
use colored::Colorize;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use serde;
use std::fs::File;
use std::path::Path;
use std::{clone, env};
use symbolica;

#[allow(unused)]
const LTD_COMPARISON_TOLERANCE: f64 = 1.0e-15;

pub fn load_amplitude_output(
    output_path: String,
    sb_state: &mut symbolica::state::State,
    sb_workspace: &symbolica::state::Workspace,
) -> (Model, Amplitude) {
    let path = Path::new(&output_path);
    let output_meta_data: OutputMetaData =
        serde_yaml::from_reader(File::open(path.join("output_metadata.yaml")).unwrap()).unwrap();
    assert_eq!(output_meta_data.output_type, OutputType::Amplitudes);
    assert_eq!(output_meta_data.contents.len(), 1);

    let model = Model::from_file(
        String::from(
            path.join(format!(
                "sources/model/{}.yaml",
                output_meta_data.model_name
            ))
            .to_str()
            .unwrap(),
        ),
        sb_state,
        sb_workspace,
    )
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

    use crate::graph::EdgeType;

    use super::*;

    #[test]
    #[ignore] // Important since this test will only run successfully when called from with pytest where the massless_triangle_generation fixture will be run
    fn pytest_massless_scalar_triangle() {
        assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());
        let mut sb_state = symbolica::state::State::new();
        let sb_workspace = symbolica::state::Workspace::new();
        let (model, amplitude) = load_amplitude_output(
            env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap(),
            &mut sb_state,
            &sb_workspace,
        );

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

        let p2 = LorentzVector::from_args(1.0, 3.0, 4.0, 5.0);
        let p1 = LorentzVector::from_args(1.0, 6.0, 7.0, 8.0);

        let k = LorentzVector::from_args(0.0, 1.0, 2.0, 3.0);

        let onshell_energies = graph.compute_onshell_energies(&[k], &[p1, p2, p1 - p2]);

        assert_eq!(onshell_energies.len(), 6);

        let edge_types = graph.get_edge_type_list();

        let energy_product = onshell_energies
            .iter()
            .zip(edge_types.iter())
            .filter(|(_e, t)| **t == EdgeType::Virtual)
            .map(|(e, _)| 2.0 * e)
            .product::<f64>()
            .recip();

        let res = energy_product * graph.evaluate_cff_expression(&[k], &[p1, p2]);

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

        graph.generate_ltd();

        let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1, p2, LorentzVector::new()]);
        for edge in graph.edges.iter() {
            let _position = graph.get_edge_position(&edge.name).unwrap();
            // println!(
            //     "edge name: {}\n edge position: {}\n vertices: {:?}\n signature: {:?}",
            //     edge.name,
            //     position,
            //     edge.vertices,
            //     graph.loop_momentum_basis.edge_signatures[position],
            // );
        }
        // println!("ltd_res = {:+e}", ltd_res);

        assert_approx_eq(res, ltd_res, 1.0e-5);
        // TODO: @Mathijs, you can put your own checks there
    }
}

#[test]
#[ignore] // Important since this test will only run successfully when called from with pytest where the massless_triangle_generation fixture will be run
fn pytest_scalar_fishnet_2x2() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());
    let mut sb_state = symbolica::state::State::new();
    let sb_workspace = symbolica::state::Workspace::new();
    let (model, amplitude) = load_amplitude_output(
        env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap(),
        &mut sb_state,
        &sb_workspace,
    );

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

    graph.generate_loop_momentum_bases();

    let k1 = LorentzVector::from_args(4.0, 1.0, 2.0, 3.0);
    let k2 = LorentzVector::from_args(5.0, 1.0, 6.0, 3.0);
    let k3 = LorentzVector::from_args(7.0, 1.0, 8.0, 0.0);
    let k4 = LorentzVector::from_args(9.0, 1.0, 10.0, 89.0);

    let zero_vector = LorentzVector::from_args(0.0, 0.0, 0.0, 0.0);

    let emr = graph.compute_emr(&[k1, k2, k3, k4], &[zero_vector; 4]);
    let n_lmb = graph
        .clone()
        .derived_data
        .loop_momentum_bases
        .unwrap()
        .len();
    assert_eq!(n_lmb, 192);
    //println!("number of lmbs: {}", n_lmb);

    for basis in graph.derived_data.loop_momentum_bases.unwrap() {
        let momenta_in_basis = basis.basis.iter().map(|index| emr[*index]).collect_vec();
        let new_emr = basis
            .edge_signatures
            .iter()
            .map(|s| compute_momentum(s, &momenta_in_basis, &[zero_vector; 4]))
            .collect_vec();
        assert_eq!(emr.len(), new_emr.len());
        for (e1, e2) in emr.iter().zip(new_emr.iter()) {
            assert_approx_eq(e1.t, e2.t, LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(e1.x, e2.x, LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(e1.y, e2.y, LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(e1.z, e2.z, LTD_COMPARISON_TOLERANCE);
        }
    }

    // TODO: @Mathijs, you can put your own checks there
}

#[test]
#[ignore] // Important since this test will only run successfully when called from with pytest where the massless_triangle_generation fixture will be run
fn pytest_scalar_fishnet_2x3() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());
    let mut sb_state = symbolica::state::State::new();
    let sb_workspace = symbolica::state::Workspace::new();
    let (model, amplitude) = load_amplitude_output(
        env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap(),
        &mut sb_state,
        &sb_workspace,
    );

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

    // TODO: @Mathijs, you can put your own checks there
}

#[test]
#[ignore]
fn pytest_scalar_cube() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let mut sb_state = symbolica::state::State::new();
    let sb_workspace = symbolica::state::Workspace::new();
    let (model, amplitude) = load_amplitude_output(
        env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap(),
        &mut sb_state,
        &sb_workspace,
    );

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

    let graph = &amplitude.amplitude_graphs[0].graph;
    let cff = generate_cff_expression(graph).unwrap();
    assert!(!cff.terms.is_empty());
}

#[test]
#[ignore]
fn pytest_scalar_bubble() {
    assert!(env::var("PYTEST_OUTPUT_PATH_FOR_RUST").is_ok());

    let mut sb_state = symbolica::state::State::new();
    let sb_workspace = symbolica::state::Workspace::new();
    let (model, amplitude) = load_amplitude_output(
        env::var("PYTEST_OUTPUT_PATH_FOR_RUST").unwrap(),
        &mut sb_state,
        &sb_workspace,
    );

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
    let cff = generate_cff_expression(&graph).unwrap();

    let p1 = LorentzVector::from_args(1.0, 3.0, 4.0, 5.0);
    let k = LorentzVector::from_args(0.0, 1.0, 2.0, 3.0);

    let onshell_energies = graph.compute_onshell_energies(&[k], &[p1, p1]);

    assert_eq!(onshell_energies.len(), 4);

    //  let inverse_energy_product = graph
    //      .edges
    //      .iter()
    //      .filter(|edge| edge.edge_type == EdgeType::Virtual)
    //      .map(|edge| graph.get_edge_position(&edge.name).unwrap())
    //      .map(|index| 2.0 * onshell_energies[index])
    //      .product::<f64>()
    //      .recip();

    // let cff_res = inverse_energy_product * cff.evaluate(&[k], &[p1, p1], &graph);

    graph.generate_ltd();

    //    let ltd_res = graph.evaluate_ltd_expression(&[k], &[p1, p1]);

    //assert_eq!(ltd_res, cff_res);

    assert_eq!(cff.terms.len(), 2);
}
