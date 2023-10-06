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

        // this is a bit handcrafted at the moment
        // In the future this is generated automatically when the full evaluation stack is there
        let energy_cache = [
            1.0,
            1.0,
            -2.0,
            10.770329614269007,
            22.9128784747792,
            3.7416573867739413,
        ];

        let benchmark_val = 0.011605113880815013;
        let res = cff.evaluate(&energy_cache);

        println!("res = {}", res);
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
