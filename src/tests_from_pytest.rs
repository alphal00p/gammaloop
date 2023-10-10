#![allow(unused_imports)]
use crate::cff::generate_cff_expression;
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::model::Model;
use crate::tests::approx_eq;
use colored::Colorize;
use serde;
use std::env;
use std::fs::File;
use std::path::Path;
use symbolica;

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

        let graph = &amplitude.amplitude_graphs[0].graph;
        let cff = generate_cff_expression(graph).unwrap();
        assert_eq!(cff.terms.len(), 6);
        assert_eq!(cff.esurfaces.len(), 6);

        println!("what is this {:?}", graph.edge_signatures);

        let p1 = LorentzVector::from_args(1.0, 3.0, 4.0, 5.0);
        let p2 = -LorentzVector::from_args(1.0, 6.0, 7.0, 8.0);

        let k = LorentzVector::from_args(0.0, 1.0, 2.0, 3.0);

        let onshell_energies = graph.compute_onshell_energies(&[k], &[p1, p2, p1 - p2]);

        assert_eq!(onshell_energies.len(), 6);
        assert_eq!(
            onshell_energies[graph.get_edge_position(&SmartString::from("p1")).unwrap()],
            p1.t
        );

        assert_eq!(
            onshell_energies[graph.get_edge_position(&SmartString::from("p2")).unwrap()],
            p2.t
        );
        assert_eq!(
            onshell_energies[graph.get_edge_position(&SmartString::from("p3")).unwrap()],
            (p1 - p2).t
        );
        assert_eq!(
            onshell_energies[graph.get_edge_position(&SmartString::from("q1")).unwrap()],
            (k - p2).spatial_distance()
        );
        assert_eq!(
            onshell_energies[graph.get_edge_position(&SmartString::from("q2")).unwrap()],
            k.spatial_distance()
        );
        assert_eq!(
            onshell_energies[graph.get_edge_position(&SmartString::from("q3")).unwrap()],
            (k - p1 - p2).spatial_distance()
        );

        println!("onshell energies = {:?}", onshell_energies);
        assert_eq!(onshell_energies.len(), 6);

        let edge_types = graph.get_edge_type_list();

        let energy_product = onshell_energies
            .iter()
            .zip(edge_types.iter())
            .filter(|(_e, t)| **t == EdgeType::Virtual)
            .map(|(e, _)| 2.0 * e)
            .product::<f64>()
            .recip();

        let res = energy_product * cff.evaluate(&onshell_energies, &edge_types);
        let benchmark_val = 4.646388885763160e-6;

        println!("res = {:+e}", res);
        assert!(approx_eq(res, benchmark_val, 1e-15));

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

    let graph = &amplitude.amplitude_graphs[0].graph;
    let cff = generate_cff_expression(graph).unwrap();
    assert!(cff.terms.len() > 0);

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
    assert!(cff.terms.len() > 0);
}
