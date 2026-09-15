use super::*;
use crate::{LoopNormalizationFactor, VakintSettings};

fn finalizer(terms: i64, digits: u32) -> FMFT {
    Vakint::initialize_vakint_symbols();
    FMFT::with_settings(VakintSettings {
        form_exe_path: "/definitely/not/a/form/executable".into(),
        number_of_terms_in_epsilon_expansion: terms,
        run_time_decimal_precision: digits,
        use_dot_product_notation: true,
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        ..VakintSettings::default()
    })
}

#[test]
fn native_pr_finite_value_uses_the_existing_table_without_form() {
    for name in ["PR0", "PR12", "PR11d"] {
        // PR11d's shipped finite part carries 26 decimal digits, not 32.
        let fmft = if name == "PR11d" {
            finalizer(5, 25)
        } else {
            finalizer(5, 32)
        };
        let result = fmft
            .finalize_native_reduced_masters(
                vk_parse!(name).unwrap(),
                &Atom::num(1),
                &FMFTOptions::default(),
            )
            .unwrap();
        let source = vk_parse!(format!("{name}ep0").as_str()).unwrap();
        let expected = fmft.substitute_masters(source.as_view()).unwrap();
        assert_eq!(result, expected, "{name}");
    }
}

#[test]
fn native_pr_missing_laurent_order_and_spurious_poles_fail_closed() {
    for (terms, expression) in [(5, "PR9x"), (6, "PR9d"), (5, "PR11d/ep")] {
        let error = finalizer(terms, 25)
            .finalize_native_reduced_masters(
                vk_parse!(expression).unwrap(),
                &Atom::num(1),
                &FMFTOptions::default(),
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("beyond expansion depth supported"),
            "{error}"
        );
    }
}

#[test]
fn native_pr_dimension_and_custom_epsilon_poles_expose_unknown_orders() {
    let mut fmft = finalizer(5, 25);
    fmft.settings.epsilon_symbol = "eps_test".into();
    for expression in ["PR11d/(d-4)", "PR11d/eps_test"] {
        let error = fmft
            .finalize_native_reduced_masters(
                vk_parse!(expression).unwrap(),
                &Atom::num(1),
                &FMFTOptions::default(),
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("beyond expansion depth supported"),
            "{expression}: {error}",
        );
    }
}

#[test]
fn native_pr_dimension_and_custom_epsilon_coefficients_match_internal_ep() {
    let mut fmft = finalizer(6, 32);
    fmft.settings.epsilon_symbol = "eps_test".into();
    let results = ["(d-4)*PR12", "-2*eps_test*PR12", "-2*ep*PR12"].map(|expression| {
        fmft.finalize_native_reduced_masters(
            vk_parse!(expression).unwrap(),
            &Atom::num(1),
            &FMFTOptions::default(),
        )
        .unwrap()
    });
    assert_eq!(results[0], results[1]);
    assert_eq!(results[1], results[2]);
    assert!(results[0].contains_symbol(vk_symbol!("eps_test")));
    assert!(!results[0].contains_symbol(vk_symbol!("d")));
    assert!(!results[0].contains_symbol(vk_symbol!("ep")));
}

#[test]
fn native_pr_rejects_unavailable_precision_without_changing_legacy_policy() {
    let fmft = finalizer(5, 32);
    let value = vk_parse!("PR11d").unwrap();
    let error = fmft
        .finalize_native_reduced_masters(value.clone(), &Atom::num(1), &FMFTOptions::default())
        .unwrap_err();
    assert!(error.to_string().contains("PR11dep0"), "{error}");
    assert!(error.to_string().contains("stored precision"), "{error}");
    assert!(
        fmft.finalize_master_expression(value, 4, &Atom::num(1), &FMFTOptions::default(), false)
            .is_ok()
    );

    let error = finalizer(5, 20_000)
        .finalize_native_reduced_masters(
            vk_parse!("PR12").unwrap(),
            &Atom::num(1),
            &FMFTOptions::default(),
        )
        .unwrap_err();
    assert!(error.to_string().contains("stored precision"), "{error}");
}

#[test]
fn native_pr_checks_only_constants_surviving_requested_truncation() {
    // PR11d starts at epsilon^0: no numerical value is needed through epsilon^-1.
    let result = finalizer(4, 20_000)
        .finalize_native_reduced_masters(
            vk_parse!("PR11d").unwrap(),
            &Atom::num(1),
            &FMFTOptions::default(),
        )
        .unwrap();
    assert!(result.is_zero());
}

#[test]
fn native_pr_retains_symbolic_basis_when_substitution_is_disabled() {
    let options = FMFTOptions {
        expand_masters: false,
        susbstitute_masters: false,
    };
    let value = vk_parse!("PR11d").unwrap();
    let result = finalizer(5, 20_000)
        .finalize_native_reduced_masters(value.clone(), &Atom::num(1), &options)
        .unwrap();
    assert!(result.contains_symbol(vk_symbol!("PR11d")));
    assert!(!result.contains_symbol(vk_symbol!("PR11dep0")));
}

#[test]
fn native_pr_preserves_mass_and_custom_epsilon_conventions() {
    let mut fmft = finalizer(5, 32);
    fmft.settings.epsilon_symbol = "eps_test".into();
    let mass = vk_parse!("13/10").unwrap();
    let value = mass.clone().pow(Atom::num(-1)) * vk_parse!("PR12").unwrap();
    let native = fmft
        .finalize_native_reduced_masters(value.clone(), &mass, &FMFTOptions::default())
        .unwrap();
    let legacy = fmft
        .finalize_master_expression(value, 4, &mass, &FMFTOptions::default(), false)
        .unwrap();
    assert_eq!(native, legacy);
    let expected_finite = fmft
        .substitute_masters(vk_parse!("PR12ep0").unwrap().as_view())
        .unwrap()
        / mass.clone();
    assert_eq!(native, expected_finite);
    assert!(!native.contains_symbol(vk_symbol!("ep")));
    assert!(!native.contains_symbol(vk_symbol!("Oep")));

    let pole = fmft
        .finalize_native_reduced_masters(vk_parse!("PR11").unwrap(), &mass, &FMFTOptions::default())
        .unwrap();
    assert!(pole.contains_symbol(vk_symbol!("eps_test")));
    assert!(!pole.contains_symbol(vk_symbol!("ep")));
    assert!(!pole.contains_symbol(vk_symbol!("ε")));
}
