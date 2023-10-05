#![allow(unused_imports)]
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::model::Model;
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
