use symbolica::domains::float::{Complex, Float};
use vakint::{NumericalEvaluationResult, VakintSettings};

#[test]
fn statistical_errors_preserve_distinct_laurent_orders_and_component_quadrature() {
    let settings = VakintSettings {
        run_time_decimal_precision: 40,
        ..VakintSettings::default()
    };
    let values = |input: &[(i64, &str, &str)]| {
        NumericalEvaluationResult::from_vec(
            input
                .iter()
                .map(|&(power, real, imaginary)| (power, (real.to_owned(), imaginary.to_owned())))
                .collect(),
            &settings,
        )
        .unwrap()
    };
    let first = values(&[(-2, "3", "4"), (2, "0", "7")]);
    let second = values(&[(-2, "4", "3"), (0, "5", "12")]);
    let expected = values(&[(-2, "5", "5"), (0, "5", "12"), (2, "0", "7")]);

    for combined in [
        first.aggregate_errors(&second).unwrap(),
        second.aggregate_errors(&first).unwrap(),
    ] {
        assert_eq!(
            combined.coefficients().keys().copied().collect::<Vec<_>>(),
            [-2, 0, 2]
        );
        let (matches, message) = combined.does_approx_match(&expected, None, 1e-30, 0.0);
        assert!(matches, "{message}");
        assert!(combined.coefficients()[&2].precision().real().exact_zero());
    }
}

#[test]
fn malformed_decimal_input_returns_an_error() {
    let result = NumericalEvaluationResult::from_vec(
        vec![(0, ("not-a-number".to_owned(), "0".to_owned()))],
        &VakintSettings::default(),
    );
    assert!(result.is_err());
}

#[test]
fn repeated_laurent_orders_are_combined_with_their_uncertainty() {
    let settings = VakintSettings {
        run_time_decimal_precision: 40,
        ..Default::default()
    };
    let combined = NumericalEvaluationResult::from_vec(
        vec![(2, ("3".into(), "4".into())), (2, ("4".into(), "3".into()))],
        &settings,
    )
    .unwrap();
    let expected =
        NumericalEvaluationResult::from_vec(vec![(2, ("7".into(), "7".into()))], &settings)
            .unwrap();
    let (matches, message) = combined.does_approx_match(&expected, None, 1e-30, 0.0);
    assert!(matches, "{message}");
    assert_eq!(combined.coefficients().len(), 1);

    let cancellation = NumericalEvaluationResult::from_vec(
        vec![
            (2, ("3".into(), "0".into())),
            (2, ("-3".into(), "0".into())),
        ],
        &settings,
    )
    .unwrap();
    assert!(
        !cancellation.coefficients()[&2]
            .precision()
            .real()
            .exact_zero()
    );
    assert!(cancellation.coefficients()[&2].to_atom().is_err());
    assert!(cancellation.warnings().iter().any(|warning| {
        warning.epsilon_power == Some(2)
            && warning.achieved_decimal_precision == Some(0)
            && warning.data_version == "explicit numerical inputs"
            && warning.sources.len() == 2
            && warning
                .sources
                .iter()
                .all(|source| source.epsilon_power == Some(2))
    }));
}

#[test]
fn imported_float_precision_limits_are_named_and_exact_zero_does_not_warn() {
    let result = NumericalEvaluationResult::from_values(
        vec![(-3, Complex::new(Float::with_val(53, 1.25), Float::new(53)))],
        100,
    );
    assert!(
        result.coefficients()[&-3]
            .precision()
            .real()
            .significant_digits()
            .unwrap()
            <= 15
    );
    assert!(result.warnings().iter().any(|warning| {
        warning.epsilon_power == Some(-3)
            && warning.requested_decimal_precision == 100
            && warning.achieved_decimal_precision.is_none()
            && warning.sources.iter().any(|source| {
                source.identity.contains("imported epsilon^-3 input 0")
                    && source
                        .available_decimal_precision
                        .is_some_and(|digits| digits <= 15)
            })
    }));
    assert!(
        result
            .warnings()
            .iter()
            .any(|warning| warning.achieved_decimal_precision.is_some())
    );

    let exact_zero = NumericalEvaluationResult::from_values(
        vec![(0, Complex::new(Float::new(53), Float::new(53)))],
        100,
    );
    assert!(exact_zero.is_zero());
    assert!(exact_zero.warnings().is_empty());
}

#[test]
fn explicit_decimal_precision_suffixes_do_not_gain_requested_digits() {
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        ..Default::default()
    };
    let result = NumericalEvaluationResult::from_vec(
        vec![(1, ("1.2345`5".into(), "2.34567890123456789`200".into()))],
        &settings,
    )
    .unwrap();
    let coefficient = &result.coefficients()[&1];
    let declared = Float::parse("1.2345`5", None).unwrap();
    assert!(coefficient.value().re.prec() <= declared.prec());
    assert!(coefficient.precision().real().significant_digits().unwrap() <= 5);
    assert!(
        coefficient
            .precision()
            .imaginary()
            .significant_digits()
            .unwrap()
            > 60
    );
    assert!(result.warnings().iter().any(|warning| {
        warning.requested_decimal_precision == 80
            && warning.sources.iter().any(|source| {
                source
                    .available_decimal_precision
                    .is_some_and(|digits| digits <= 5)
            })
    }));
}

#[test]
fn suffix_free_decimals_have_the_declared_input_precision() {
    let settings = VakintSettings {
        run_time_decimal_precision: 40,
        ..Default::default()
    };
    let result = NumericalEvaluationResult::from_vec(
        vec![(
            0,
            (
                "1.23456789012345678901234567890123456789".into(),
                "0".into(),
            ),
        )],
        &settings,
    )
    .unwrap();
    assert!(!result.warnings().iter().any(|warning| {
        warning.achieved_decimal_precision.is_none()
            && warning.sources.iter().any(|source| {
                source
                    .available_decimal_precision
                    .is_some_and(|digits| digits < 40)
            })
    }));
}
