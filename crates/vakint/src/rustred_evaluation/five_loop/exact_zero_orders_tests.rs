//! Numerical-mode regression for Laurent orders that vanish identically.
//!
//! From external data v0.3.0 on, the former five-loop masters P31, P125 and
//! P131 reduce to other masters. Their 1/ε^5..1/ε orders are identically zero
//! (the Symbolic evaluation has no such term), but in Numerical mode the exact
//! master coefficients of those orders were evaluated as independent finite
//! values whose cancellation left an uncertain zero, which the precision
//! check refused. These orders must now be exact zeros, the finite part must
//! agree with the values these integrals had as published masters (v0.2.3,
//! to all printed digits within the requested precision), and orders that do
//! not vanish (the exact poles of the former master P108) must be kept.
//!
//! Needs `VAKINT_TEST_EXTERNAL_DATA` pointing to a v0.3.0 or later checkout;
//! with an older checkout (or none) the test prints why and returns.
//! `cargo nextest run -p vakint --run-ignored only -E 'test(exact_zero_orders_tests)'`.
use super::*;
use crate::get_integer_from_atom;
use std::collections::HashMap;
use symbolica::domains::float::{Complex, Float, Real};

/// The data version of an external-data checkout, as (major, minor, patch).
fn data_version(root: &std::path::Path) -> Option<(u64, u64, u64)> {
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("manifest.json")).ok()?).ok()?;
    let version = manifest["version"].as_str()?.trim_start_matches('v');
    let mut parts = version
        .split(|c: char| !c.is_ascii_digit())
        .map(|part| part.parse::<u64>().ok());
    Some((parts.next()??, parts.next()??, parts.next()??))
}

fn epsilon_orders(value: &Atom) -> BTreeMap<i64, Atom> {
    let epsilon = Atom::var(vk_symbol!("ε"));
    value
        .coefficient_list::<i32>(std::slice::from_ref(&epsilon))
        .into_iter()
        .map(|(monomial, coefficient)| {
            let power = if monomial == Atom::num(1) {
                0
            } else if monomial == epsilon {
                1
            } else {
                match monomial.as_view() {
                    AtomView::Pow(power) => get_integer_from_atom(power.get_exp())
                        .unwrap_or_else(|| panic!("non-integer epsilon power {monomial}")),
                    _ => panic!("unexpected epsilon monomial {monomial}"),
                }
            };
            (power, coefficient)
        })
        .collect()
}

/// The significant decimal digits of an exported floating coefficient.
fn exported_digits(coefficient: &Atom) -> u32 {
    let text = coefficient.to_canonical_string();
    let mantissa = text
        .trim_start_matches('-')
        .split(['e', 'E', '`'])
        .next()
        .unwrap();
    let digits = mantissa
        .chars()
        .filter(char::is_ascii_digit)
        .skip_while(|c| *c == '0');
    u32::try_from(digits.count()).unwrap()
}

fn real_value(coefficient: &Atom) -> Float {
    coefficient
        .evaluate_with_prec(&HashMap::<Atom, Complex<Float>>::default(), 256)
        .unwrap_or_else(|error| panic!("{coefficient}: {error}"))
        .re
}

#[test]
#[ignore = "loads a v0.3.0 five-loop package and evaluates four integrals (about 10 minutes, ~25 GB); needs VAKINT_TEST_EXTERNAL_DATA (v0.3.0 or later)"]
fn former_masters_with_vanishing_poles_evaluate_numerically() {
    let Some(root) = std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(PathBuf::from) else {
        eprintln!("skipped: VAKINT_TEST_EXTERNAL_DATA is unset");
        return;
    };
    match data_version(&root) {
        Some(version) if version >= (0, 3, 0) => {}
        version => {
            eprintln!(
                "skipped: {} declares data version {version:?}; this regression needs v0.3.0 or later",
                root.display()
            );
            return;
        }
    }
    let mut vakint = Vakint::new().unwrap();
    vakint.load_external_data(&root, 5).unwrap();
    let settings = VakintSettings {
        form_exe_path: "/exact-zero-orders-test-must-not-use-form".into(),
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Numerical,
            },
        )]),
        number_of_terms_in_epsilon_expansion: 6,
        ..Default::default()
    };
    let digits = settings.run_time_decimal_precision;
    let bits = crate::precision::decimal_bits(digits + 10);
    // The published v0.2.3 finite parts of these integrals as masters (all
    // their pole coefficients were "0"), and the digits the fixed evaluation
    // exports on the v0.3.0 RC tree minus two: P31 certifies all 32 requested
    // digits, P125 and P131 certify 28 (the tracked enclosure of their finite
    // part; `to_atom` keeps the certified bits less two guard bits).
    for (name, powers, published, floor) in [
        (
            "P31",
            [0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
            "-2.21096283270677178249823255964778140410637976127e+0",
            30,
        ),
        (
            "P125",
            [2, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
            "-1.522366493525372497169378588027918096329e-2",
            26,
        ),
        (
            "P131",
            [3, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
            "5.17609395952911714479391175457927012585452527727170634e-3",
            26,
        ),
    ] {
        let input = super::admission_census_tests::input_for_key(
            &IntegralKey::try_new(powers.to_vec()).unwrap(),
        )
        .unwrap();
        let value = vakint
            .evaluate_integral(&settings, input.expression().as_view())
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        let orders = epsilon_orders(&value);
        assert_eq!(
            orders.keys().copied().collect::<Vec<_>>(),
            vec![0],
            "{name}: only the finite part is nonzero: {value}"
        );
        let finite = real_value(&orders[&0]);
        let published = Float::parse(published, Some(bits)).unwrap();
        // The exported value carries only its certified digits (at most the
        // requested ones): at least the floor, and agreement to all of them
        // but the last.
        let exported = exported_digits(&orders[&0]);
        assert!(
            exported >= floor,
            "{name}: only {exported} digits: {}",
            orders[&0]
        );
        let tolerance = published.clone().norm()
            * Float::parse(
                &format!("1e{}", 2 - i64::from(exported.min(digits))),
                Some(bits),
            )
            .unwrap();
        assert!(
            (finite.clone() - &published).norm() <= tolerance,
            "{name}: finite part {finite} disagrees with the published {published}"
        );
    }
    // P108's poles are exact and nonzero (1/ε^2: -ζ(3)/2): kept, not dropped.
    let input = super::admission_census_tests::input_for_key(
        &IntegralKey::try_new(vec![1, 1, 2, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1]).unwrap(),
    )
    .unwrap();
    let value = vakint
        .evaluate_integral(&settings, input.expression().as_view())
        .unwrap();
    let orders = epsilon_orders(&value);
    assert_eq!(
        orders.keys().copied().collect::<Vec<_>>(),
        vec![-2, -1, 0],
        "P108: {value}"
    );
    let half_zeta3 = Float::parse("0.60102845157979714269986908075573", Some(bits)).unwrap();
    assert!(
        (real_value(&orders[&-2]) + half_zeta3).norm() < Float::parse("1e-30", Some(bits)).unwrap(),
        "P108 1/ε^2: {}",
        orders[&-2]
    );
}
