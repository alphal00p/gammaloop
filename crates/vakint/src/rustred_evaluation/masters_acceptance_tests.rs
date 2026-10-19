//! Exercise coefficients beyond the historical six-term export through the
//! actual package registry, symbolic expansion, and numerical evaluator.

use super::super::data::RustRedData;
use super::*;
use symbolica::domains::float::Real;

fn check_added_coefficient(data: &RustRedData, loops: usize, name: &str, power: i64) {
    crate::Vakint::initialize_vakint_symbols();
    let metadata = data.loops[&loops]
        .manifest
        .masters
        .iter()
        .find(|master| master["name"].as_str() == Some(name))
        .expect("acceptance master must be declared by its package");
    let declared = metadata["coefficients"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["epsilon_power"].as_i64() == Some(power))
        .expect("new epsilon order must survive generated package metadata");
    assert_eq!(declared["symbolic"], true);
    assert_eq!(declared["numerical"], true);
    let catalog = &data.masters;
    let (symbol, master) = catalog
        .masters
        .iter()
        .find(|(_, master)| Some(master.identity.as_str()) == metadata["master_id"].as_str())
        .expect("declared master must reach the runtime catalog");
    let settings = VakintSettings {
        epsilon_symbol: "ep".into(),
        run_time_decimal_precision: 80,
        number_of_terms_in_epsilon_expansion: loops as i64 + power + 1,
        ..Default::default()
    };
    let symbolic = catalog
        .finalize(
            symbol.clone(),
            &settings,
            loops,
            RustRedMasterMode::Symbolic,
        )
        .unwrap();
    assert!(!symbolic.contains_symbol(vk_symbol!("RustRedCoefficient")));
    let numerical = catalog
        .finalize(
            symbol.clone(),
            &settings,
            loops,
            RustRedMasterMode::Numerical,
        )
        .unwrap();
    assert_eq!(
        catalog.resolve_numerical(symbolic, &settings).unwrap(),
        numerical,
        "symbolic-first and direct numerical modes must use the same exact coefficients"
    );
    let epsilon = vk_parse!("ep").unwrap();
    let coefficient = numerical
        .coefficient_list::<i32>(std::slice::from_ref(&epsilon))
        .into_iter()
        .find(|(monomial, _)| *monomial == epsilon.pow(Atom::num(power)))
        .expect("requested highest epsilon coefficient must not be truncated")
        .1;
    let bits = settings.get_binary_precision();
    let actual = coefficient
        .evaluate_with_prec(&HashMap::<Atom, Complex<Float>>::new(), bits)
        .unwrap();
    let expected = Float::parse(
        master.coefficients[&power].numerical.as_ref().unwrap(),
        Some(bits),
    )
    .unwrap();
    let error = (actual.re.clone() - expected.clone()).norm() + actual.im.norm();
    let tolerance =
        Float::parse("1e-70", Some(bits)).unwrap() * (expected.norm() + Float::with_val(bits, 1));
    assert!(
        error < tolerance,
        "{loops}L {name} epsilon^{power}: error {error}"
    );

    // The next absent order must still fail instead of becoming an exact zero.
    let unavailable = VakintSettings {
        number_of_terms_in_epsilon_expansion: settings.number_of_terms_in_epsilon_expansion + 1,
        ..settings
    };
    let error = catalog
        .finalize(
            symbol.clone(),
            &unavailable,
            loops,
            RustRedMasterMode::Numerical,
        )
        .unwrap_err()
        .to_string();
    assert!(
        error.contains(&format!("unavailable epsilon coefficient {}", power + 1)),
        "{error}"
    );
}

#[test]
fn vendored_sunset_epsilon_four_is_available() {
    check_added_coefficient(&RustRedData::vendored().unwrap(), 2, "P2", 4);
}

#[test]
fn vendored_mibn_epsilon_four_is_available() {
    check_added_coefficient(&RustRedData::vendored().unwrap(), 3, "P13", 4);
}

#[test]
fn vendored_five_line_epsilon_three_is_available() {
    check_added_coefficient(&RustRedData::vendored().unwrap(), 3, "P14", 3);
}

#[test]
fn external_pr2_epsilon_three_is_available() {
    let path = super::super::cache_tests::external_path();
    let data = RustRedData::from_directory(&path, "latest", 4).unwrap();
    check_added_coefficient(&data, 4, "P2", 3);
}
