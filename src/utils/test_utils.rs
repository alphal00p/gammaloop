use std::{env, path::PathBuf};

use crate::model::Model;

use include_dir::{include_dir, Dir};

static BUILTIN_MODELS: Dir = include_dir!("$CARGO_MANIFEST_DIR/models/json");

pub fn load_generic_model(name: &str) -> Model {
    if let Some(file) = BUILTIN_MODELS.get_file(format!("{}/{}.json", name, name)) {
        Model::from_str(file.contents_utf8().unwrap().into(), "json").unwrap()
    } else {
        panic!("Model {} not found in built-in models.", name);
    }
}

pub fn output_dir() -> PathBuf {
    if let Ok(pytest_output_path) = env::var("PYTEST_OUTPUT_PATH_FOR_RUST") {
        pytest_output_path.into()
    } else {
        PathBuf::from("./src/test_resources")
    }
}

#[cfg(test)]
pub mod test {
    use insta::assert_snapshot;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::{atom::AtomCore, parse_lit};

    use crate::utils::{symbolica_ext::CallSymbol, GS};

    #[test]
    fn normalization() {
        let a = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.f(&[0]));
        let b = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.f(&[1]));

        assert_snapshot!(a.to_canonical_string(),@"0");
        assert_snapshot!(b.to_canonical_string(),@"_gammaloop::Q(1,spenso::cind(1))");

        let c = GS.delta_vec(0, GS.cind(1));
        assert_snapshot!(c.to_canonical_string(),@"0");
        let c = GS.delta_vec(0, GS.cind(0));
        assert_snapshot!(c.to_canonical_string(),@"1");

        let expr = parse_lit!(f(a + p + r));
        assert_snapshot!(GS.linearize.f(&[expr]).to_canonical_string(),@"_gammaloop::f(_gammaloop::a)+_gammaloop::f(_gammaloop::p)+_gammaloop::f(_gammaloop::r)");
    }
}
