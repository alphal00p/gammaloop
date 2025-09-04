use std::{env, path::PathBuf};

use crate::model::Model;

pub fn load_generic_model(name: &str) -> Model {
    Model::from_file(output_dir().join(format!("gammaloop_models/{}.json", name))).unwrap()
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
    use idenso::color::ColorSimplifier;
    use insta::assert_snapshot;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::{atom::AtomCore, parse_lit};

    use crate::{symbolica_ext::CallSymbol, utils::GS};

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
        c.simplify_color();

        let expr = parse_lit!(f(a + p + r));
        assert_snapshot!(GS.linearize.f(&[expr]).to_canonical_string(),@"_gammaloop::f(_gammaloop::a)+_gammaloop::f(_gammaloop::p)+_gammaloop::f(_gammaloop::r)");
    }
}
