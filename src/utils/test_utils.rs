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

#[cfg(test)]
pub mod test {
    use insta::assert_snapshot;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::atom::AtomCore;

    use crate::{symbolica_ext::CallSymbol, utils::GS};

    #[test]
    fn normalization() {
        let a = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.f(&[0]));
        let b = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.f(&[1]));

        assert_snapshot!(a.to_canonical_string(),@"0");
        assert_snapshot!(b.to_canonical_string(),@"_gammaloop::Q(1,spenso::cind(1))");
    }
}
