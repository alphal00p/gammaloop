use std::{env, path::Path};

use crate::model::Model;

pub(crate) fn load_generic_model(name: &str) -> Model {
    Model::from_file(String::from(
        Path::new(&output_dir())
            .join(format!("gammaloop_models/{}.yaml", name))
            .to_str()
            .unwrap(),
    ))
    .unwrap()
}

fn output_dir() -> String {
    if let Ok(pytest_output_path) = env::var("PYTEST_OUTPUT_PATH_FOR_RUST") {
        pytest_output_path
    } else {
        String::from("./src/test_resources")
    }
}
