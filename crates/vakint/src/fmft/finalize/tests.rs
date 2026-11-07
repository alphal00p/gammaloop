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
        let expected = fmft
            .substitute_additional_constants(expected.as_view())
            .unwrap();
        assert_eq!(result, expected, "{name}");
    }
}

#[test]
fn native_pr_exact_cancellation_precedes_approximate_master_tables() {
    let fmft = finalizer(5, 25);
    let value = vk_parse!("(PR11d+PR12)/3-PR11d/3-PR12/3").unwrap();
    assert!(!value.is_zero(), "regression input must remain unexpanded");
    assert!(value.expand().is_zero(), "input cancels exactly");
    let pre_table = fmft
        .expand_masters(value.as_view())
        .unwrap()
        .series(vk_symbol!("ep"), Atom::Zero.as_view(), Rational::from(0))
        .unwrap()
        .to_atom();
    assert!(
        !pre_table.is_zero(),
        "the series must retain the unexpanded cancellation"
    );
    assert!(pre_table.expand().is_zero());
    // On the pinned native stack, substituting this unexpanded coefficient
    // first leaves -2.1e-25 at 25 digits. Its exact expansion must precede
    // approximation; neither a tolerance nor increased precision fixes that.
    let result = fmft
        .finalize_native_reduced_masters(value, &Atom::num(1), &FMFTOptions::default())
        .unwrap();
    assert!(
        result.is_zero(),
        "exact cancellation must not leave a float residual"
    );
}

#[test]
fn native_pr_missing_laurent_order_and_spurious_poles_fail_closed() {
    // (terms, expression, master, first missing order, last needed order)
    for (terms, expression, master, first, last) in [
        (5, "PR9x", "I_PR9x", 0, 0),
        (6, "PR9d", "I_PR9d", 1, 1),
        (5, "PR11d/ep", "I_PR11d", 1, 1),
    ] {
        let error = finalizer(terms, 25)
            .finalize_native_reduced_masters(
                vk_parse!(expression).unwrap(),
                &Atom::num(1),
                &FMFTOptions::default(),
            )
            .unwrap_err();
        let missing = crate::truncation_guard::expect_missing_orders(error);
        assert_eq!(missing.backend, "FMFT", "{expression}");
        assert_eq!(missing.requested_order, Some(terms - 5), "{expression}");
        assert_eq!(missing.masters.len(), 1, "{expression}: {missing}");
        let missing = &missing.masters[0];
        assert_eq!(
            (missing.master.as_str(), missing.first_missing_order),
            (master, first),
            "{expression}"
        );
        assert_eq!(missing.last_needed_order, Some(last), "{expression}");
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
        let missing = crate::truncation_guard::expect_missing_orders(error);
        assert_eq!(missing.masters[0].master, "I_PR11d", "{expression}");
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

#[test_log::test]
fn native_pr_warns_and_preserves_requested_precision_and_legacy_values() {
    let fmft = finalizer(5, 32);
    let value = vk_parse!("PR11d").unwrap();
    let native = fmft
        .finalize_native_reduced_masters(value.clone(), &Atom::num(1), &FMFTOptions::default())
        .unwrap();
    let legacy = fmft
        .finalize_master_expression(value, 4, &Atom::num(1), &FMFTOptions::default())
        .unwrap();
    assert_eq!(native, legacy);
    assert_eq!(fmft.settings.run_time_decimal_precision, 32);

    let high = finalizer(5, 20_000);
    let native = high
        .finalize_native_reduced_masters(
            vk_parse!("PR12").unwrap(),
            &Atom::num(1),
            &FMFTOptions::default(),
        )
        .unwrap();
    let expected = high
        .substitute_masters(vk_parse!("PR12ep0").unwrap().as_view())
        .unwrap();
    let expected = high
        .substitute_additional_constants(expected.as_view())
        .unwrap();
    assert_eq!(native, expected);
    assert_eq!(high.settings.run_time_decimal_precision, 20_000);
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
fn native_pr_symbolic_poles_retain_canonical_exact_zeta() {
    let result = finalizer(4, 20_000)
        .finalize_native_reduced_masters(
            vk_parse!("PR11").unwrap(),
            &Atom::num(1),
            &FMFTOptions {
                expand_masters: true,
                susbstitute_masters: false,
            },
        )
        .unwrap();
    assert_eq!(
        result,
        vk_parse!("5*symbolica::zeta(5)/ε").unwrap(),
        "symbolic output must cancel exact Forcer pole coefficients"
    );
}

#[test]
fn native_pr_symbolic_poles_match_exact_lower_loop_masters() {
    let fmft = finalizer(4, 32);
    let actual = fmft
        .finalize_native_reduced_masters(
            vk_parse!("PR13").unwrap(),
            &Atom::num(1),
            &FMFTOptions {
                expand_masters: true,
                susbstitute_masters: false,
            },
        )
        .unwrap();
    // PR13 factorizes into the three-loop miD6 master and a one-loop
    // tadpole, with the sign appropriate to the FMFT master convention.
    let expected = crate::matad::MATAD::with_settings(fmft.settings.clone())
        .finalize_reduced_masters(
            vk_parse!("miD6*Gam(1,1)/(ep*(1-ep))").unwrap(),
            4,
            None,
            &Atom::num(1),
            &crate::MATADOptions {
                direct_numerical_substition: false,
                substitute_hpls: false,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(
        (actual - expected).expand().is_zero(),
        "FMFT's D6 and MATAD's exact finite miD6 coefficient must cancel"
    );
}

#[test]
fn native_pr_laurent_expansion_preserves_open_metric_and_scalar_evanescences() {
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::{BaseRepName, Minkowski, initialize},
    };

    let mut fmft = finalizer(5, 32);
    fmft.settings.epsilon_symbol = "eps_metric".into();
    initialize();
    let dimension = vk_parse!("d").unwrap();
    let slot = |index| Minkowski::selfless_symbol().call((&dimension, index));
    let metric = ETS.metric(slot(vk_symbol!("open_mu")), slot(vk_symbol!("open_nu")));
    let input = &metric * vk_parse!("PR11*(1+ep)/(4-2*ep)").unwrap();
    let result = fmft
        .finalize_native_reduced_masters(
            input,
            &Atom::num(1),
            &FMFTOptions {
                expand_masters: true,
                susbstitute_masters: false,
            },
        )
        .unwrap();
    let metric = metric
        .replace(dimension.to_pattern())
        .with(vk_parse!("4-2*eps_metric").unwrap());
    let expected = metric
        * vk_parse!("5*symbolica::zeta(5)/(4*eps_metric)+15/8*symbolica::zeta(5)+PR11ep0/4")
            .unwrap();
    assert!(!result.contains_symbol(Symbol::DERIVATIVE));
    assert!(!result.contains_symbol(vk_symbol!("ep")));
    assert!((result - expected).expand().is_zero());
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
        .finalize_master_expression(value, 4, &mass, &FMFTOptions::default())
        .unwrap();
    assert_eq!(native, legacy);
    let expected_finite = fmft
        .substitute_masters(vk_parse!("PR12ep0").unwrap().as_view())
        .unwrap()
        / mass.clone();
    let expected_finite = fmft
        .substitute_additional_constants(expected_finite.as_view())
        .unwrap();
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
