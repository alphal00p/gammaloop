use std::{env, path::PathBuf};

use crate::model::Model;

use include_dir::{Dir, include_dir};

static BUILTIN_MODELS: Dir = include_dir!("$CARGO_MANIFEST_DIR/../../assets/models/json");

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
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/test_resources")
    }
}

#[cfg(test)]
pub mod test {
    use insta::assert_snapshot;
    use itertools::Itertools;
    use linnet::half_edge::involution::EdgeIndex;
    use momtrop::assert_approx_eq;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::{atom::AtomCore, parse_lit};

    use crate::{
        momentum::ThreeMomentum,
        settings::runtime::ParameterizationSettings,
        utils::{F, GS, global_inv_parameterize, global_parameterize, symbolica_ext::CallSymbol},
    };

    #[test]
    fn normalization() {
        let a = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.f(&[0]));
        let b = GS.emr_vec_index(EdgeIndex(1), AIND_SYMBOLS.cind.f(&[1]));

        assert_snapshot!(a.to_canonical_string(),@"0");
        assert_snapshot!(b.to_canonical_string(),@"gammalooprs::{}::Q(1,spenso::{}::cind(1))");

        let c = GS.delta_vec(0, GS.cind(1));
        assert_snapshot!(c.to_canonical_string(),@"0");
        let c = GS.delta_vec(0, GS.cind(0));
        assert_snapshot!(c.to_canonical_string(),@"1");

        let expr = parse_lit!(f(a + p + r));
        assert_snapshot!(GS.linearize.f(&[expr]).to_canonical_string(),@"gammalooprs::{}::f(gammalooprs::{}::a)+gammalooprs::{}::f(gammalooprs::{}::p)+gammalooprs::{}::f(gammalooprs::{}::r)");
    }

    #[test]
    fn test_inv_param() {
        let x = [
            F(0.1),
            F(0.2),
            F(0.3),
            F(0.4),
            F(0.5),
            F(0.6),
            F(0.7),
            F(0.8),
            F(0.9),
        ];
        let e_cm = F(42.2);
        let param_settings = ParameterizationSettings::default();

        let (momenta, jac_1) = global_parameterize(&x, e_cm.clone(), &param_settings);
        let actual_momenta = momenta
            .iter()
            .map(|p| ThreeMomentum {
                px: p[0],
                py: p[1],
                pz: p[2],
            })
            .collect_vec();

        let (xs, jac_2) = global_inv_parameterize(&actual_momenta, e_cm, &param_settings);

        for i in 0..9 {
            assert_approx_eq(&x[i], &xs[i], &F(1.0e-14));
        }

        let prod = jac_1 * jac_2;
        assert_approx_eq(&prod, &F(1.0), &F(1.0e-14));
    }
}
