use super::*;

fn settings(digits: u32) -> VakintSettings {
    VakintSettings {
        run_time_decimal_precision: digits,
        ..Default::default()
    }
}

#[test]
fn public_numerical_evaluation_collects_orders_beyond_the_old_i8_range() {
    let engine = Vakint::new().unwrap();
    let settings = VakintSettings {
        epsilon_symbol: "vakint::test_regulator".into(),
        ..settings(32)
    };
    for (input, expected) in [
        ("test_regulator^128", vec![(128, 1)]),
        ("test_regulator^-129", vec![(-129, 1)]),
        (
            "(test_regulator^100+2*test_regulator^-130)*(test_regulator^30+3)",
            vec![(-130, 6), (-100, 2), (100, 3), (130, 1)],
        ),
        (
            "(test_regulator^2-1)/(test_regulator-1)",
            vec![(0, 1), (1, 1)],
        ),
    ] {
        let expression = vk_parse!(input).unwrap();
        let (result, error) = engine
            .numerical_evaluation(
                &settings,
                expression.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                None,
            )
            .unwrap();
        assert!(error.is_none());
        assert_eq!(
            result.coefficients().keys().copied().collect::<Vec<_>>(),
            expected.iter().map(|(power, _)| *power).collect::<Vec<_>>(),
            "{input}"
        );
        for (power, value) in expected {
            let coefficient = &result.coefficients()[&power];
            assert_eq!(coefficient.value().re.to_rational(), value, "{input}");
            assert!(coefficient.precision().imaginary().exact_zero());
        }
    }
}

#[test]
fn public_numerical_evaluation_preserves_namespaced_regulator_identity() {
    let engine = Vakint::new().unwrap();
    for name in [
        "ordered_carrier_control::epsilon",
        "another_regulator::epsilon",
    ] {
        let settings = VakintSettings {
            epsilon_symbol: name.into(),
            ..settings(32)
        };
        let epsilon = Atom::var(vk_symbol!(name));
        let expression = epsilon.clone().pow(-129)
            + Atom::num(2) * epsilon.clone().pow(-2)
            + Atom::num(3)
            + Atom::num(4) * epsilon.clone().pow(2)
            + Atom::num(5) * epsilon.pow(128);
        let (result, error) = engine
            .numerical_evaluation(
                &settings,
                expression.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                None,
            )
            .unwrap();
        assert!(error.is_none());
        assert_eq!(
            result.coefficients().keys().copied().collect::<Vec<_>>(),
            [-129, -2, 0, 2, 128]
        );
        for (power, value) in [(-129, 1), (-2, 2), (0, 3), (2, 4), (128, 5)] {
            assert_eq!(
                result.coefficients()[&power].value().re.to_rational(),
                value
            );
            assert!(
                result.coefficients()[&power]
                    .precision()
                    .imaginary()
                    .exact_zero()
            );
        }
        let homonym = Atom::var(vk_symbol!("unconfigured_regulator::epsilon")).pow(2);
        let error = engine
            .numerical_evaluation(
                &settings,
                homonym.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                None,
            )
            .unwrap_err();
        assert!(!error.to_string().is_empty());
        assert!(
            !error
                .to_string()
                .contains("epsilon power is not an integer monomial")
        );
    }
}

#[test]
fn namespaced_regulators_reject_noninteger_laurent_powers_gracefully() {
    let engine = Vakint::new().unwrap();
    let settings = VakintSettings {
        epsilon_symbol: "ordered_carrier_control::epsilon".into(),
        ..settings(32)
    };
    for exponent in [vk_parse!("1/2").unwrap(), vk_parse!("-1/2").unwrap()] {
        let expression = Atom::var(vk_symbol!(&settings.epsilon_symbol)).pow(exponent);
        let error = engine
            .numerical_evaluation(
                &settings,
                expression.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                None,
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("not an integer Laurent monomial")
        );
    }
}

#[test]
fn public_numerical_evaluation_rejects_oversized_degrees_before_expansion() {
    let engine = Vakint::new().unwrap();
    let settings = VakintSettings {
        epsilon_symbol: "vakint::test_regulator".into(),
        ..settings(32)
    };
    for input in [
        "test_regulator^2147483648",
        "test_regulator^-2147483648",
        "test_regulator^9223372036854775808",
        "test_regulator^-9223372036854775808",
        "((test_regulator+1)^1073741824)^2",
        "(test_regulator^1073741824+1)*(test_regulator^1073741824+1)",
        "(test_regulator+1)^-2147483648",
        "test_regulator^2147483647*(1+test_regulator)",
    ] {
        let expression = vk_parse!(input).unwrap();
        for error in [
            engine
                .numerical_evaluation(
                    &settings,
                    expression.as_view(),
                    &HashMap::default(),
                    &HashMap::default(),
                    None,
                )
                .unwrap_err(),
            NumericalEvaluationResult::from_atom(
                expression.as_view(),
                vk_symbol!(&settings.epsilon_symbol),
                &settings,
            )
            .unwrap_err(),
        ] {
            assert!(
                error.to_string().contains("epsilon degree budget"),
                "{input}: {error}"
            );
        }
    }
}

#[test]
fn public_numerical_evaluation_rejects_noninteger_laurent_powers_gracefully() {
    let engine = Vakint::new().unwrap();
    let settings = settings(32);
    for input in ["ε^(1/2)", "ε^(4294967296/3)", "ε^(4294967296+x)"] {
        let expression = vk_parse!(input).unwrap();
        let error = engine
            .numerical_evaluation(
                &settings,
                expression.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                None,
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("not an integer Laurent monomial"),
            "{input}: {error}"
        );
    }
    let cancelled = vk_parse!("(ε^(1/2))^2-ε").unwrap();
    let (result, _) = engine
        .numerical_evaluation(
            &settings,
            cancelled.as_view(),
            &HashMap::default(),
            &HashMap::default(),
            None,
        )
        .unwrap();
    assert!(result.is_zero());
}

fn source(name: &str, value: &str, digits: u32) -> (Atom, FiniteSource) {
    (
        Atom::var(vk_symbol!(name)),
        FiniteSource {
            value: value.into(),
            significant_digits: digits,
            identity: name.into(),
            master_id: Some(format!("test::{name}")),
            epsilon_power: Some(0),
            data_version: "precision-test-v1".into(),
            source_sha256: "test-payload".into(),
            absolute_error: None,
            error_is_bounded: false,
        },
    )
}

#[test]
fn raw_center_and_original_error_survive_until_checked_symbolic_conversion() {
    let tracked = TrackedComplex::from_decimal_with_error(
        "1.23456789123456789123456789",
        "1e-10",
        256,
        ErrorQuality::Estimated,
    )
    .unwrap();
    let coefficient = NumericalCoefficient::from_tracked(tracked.clone(), 16);
    assert_eq!(coefficient.value().re, tracked.ball.re.center);
    assert_eq!(
        coefficient.precision().real().absolute_error(),
        &tracked.ball.re.radius
    );
    assert_eq!(
        coefficient.precision().real().quality(),
        ErrorQuality::Estimated
    );
    assert_eq!(
        coefficient.precision().imaginary().absolute_error(),
        &tracked.ball.im.radius
    );
    assert!(coefficient.precision().imaginary().exact_zero());
    assert!(coefficient.precision().real().significant_digits().unwrap() < 12);

    let encoded = coefficient.to_atom().unwrap();
    let AtomView::Num(number) = encoded.as_view() else {
        panic!("expected encoded Float")
    };
    let symbolica::coefficient::CoefficientView::Float(real, imaginary) = number.get_coeff_view()
    else {
        panic!("expected encoded Float")
    };
    let real = real.to_float();
    assert!(imaginary.to_float().is_zero());
    assert_ne!(real.to_rational(), coefficient.value().re.to_rational());
    assert!(real.prec() <= coefficient.precision().real().significant_bits().unwrap());
    let rounded = tracked.rounded(Some(decimal_bits(16))).unwrap();
    assert_eq!(real, rounded.ball.re.center);
    let shift = (real.to_rational() - tracked.ball.re.center.to_rational()).abs();
    assert!(shift + tracked.ball.re.radius.to_rational() <= rounded.ball.re.radius.to_rational());
}

#[test]
fn explicit_presentation_caps_do_not_reduce_the_raw_center_or_increase_export_precision() {
    let tracked = TrackedComplex::from_decimal_with_error(
        "1.23456789123456789123456789",
        "1e-70",
        256,
        ErrorQuality::Bounded,
    )
    .unwrap();
    let coefficient = NumericalCoefficient::from_tracked(tracked.clone(), 16);
    let presented = coefficient.capped_decimal_precision(50);
    assert_eq!(coefficient.value().re, tracked.ball.re.center);
    assert!(presented.value().re.prec() <= decimal_bits(16));
    assert_eq!(presented.maximum_export_bits, decimal_bits(16));
    assert!(presented.precision().real().absolute_error() >= &tracked.ball.re.radius);
    assert_eq!(
        presented.precision().real().quality(),
        ErrorQuality::Bounded
    );
    assert!(coefficient.to_real_float().unwrap().prec() <= decimal_bits(16));
}

#[test]
fn independent_master_cancellation_loses_measured_digits() {
    Vakint::initialize_vakint_symbols();
    let sources = BTreeMap::from([
        source("precision_a", "1.0000000000000000000000001", 50),
        source("precision_b", "1", 50),
    ]);
    let expression = vk_parse!("precision_a-precision_b").unwrap();
    let (coefficient, warnings) = evaluate_coefficient(
        expression.as_view(),
        &sources,
        &HashMap::new(),
        &settings(50),
    )
    .unwrap();
    let digits = coefficient.precision().real().significant_digits().unwrap();
    assert!((20..=25).contains(&digits), "achieved {digits} digits");
    // The central estimate still retains its guard digits; its uncertainty is
    // independently limited by the two validated finite sources.
    let expected_center = Float::parse("1e-25", Some(320)).unwrap().to_rational();
    let center_error = (coefficient.value().re.to_rational() - expected_center).abs();
    assert!(center_error < Float::parse("1e-70", Some(320)).unwrap().to_rational());
    assert!(coefficient.precision().imaginary().exact_zero());
    assert!(
        warnings
            .iter()
            .any(|warning| warning.achieved_decimal_precision == Some(digits))
    );
    assert!(
        warnings
            .iter()
            .flat_map(|warning| &warning.sources)
            .any(|source| source.identity == "precision_a")
    );
    assert!(
        warnings
            .iter()
            .flat_map(|warning| &warning.sources)
            .any(|source| source.identity == "precision_b")
    );
    assert!(coefficient.to_atom().is_ok());
}

#[test]
fn equal_independent_estimates_are_an_uncertain_zero() {
    Vakint::initialize_vakint_symbols();
    let sources = BTreeMap::from([
        source("precision_c", "1", 50),
        source("precision_d", "1", 50),
    ]);
    let expression = vk_parse!("precision_c-precision_d").unwrap();
    let (coefficient, warnings) = evaluate_coefficient(
        expression.as_view(),
        &sources,
        &HashMap::new(),
        &settings(50),
    )
    .unwrap();
    assert!(coefficient.value().re.is_zero());
    assert!(!coefficient.precision().real().exact_zero());
    assert!(!coefficient.precision().real().absolute_error().is_zero());
    assert!(
        coefficient
            .to_atom()
            .unwrap_err()
            .to_string()
            .contains("no reliable")
    );
    assert!(
        warnings
            .iter()
            .any(|warning| warning.achieved_decimal_precision == Some(0))
    );
}

#[test]
fn identical_source_cancellation_is_exact_and_has_no_warning() {
    Vakint::initialize_vakint_symbols();
    let sources = BTreeMap::from([source("precision_same", "1", 4)]);
    let expression = vk_parse!("precision_same-precision_same").unwrap();
    let (coefficient, warnings) = evaluate_coefficient(
        expression.as_view(),
        &sources,
        &HashMap::new(),
        &settings(200),
    )
    .unwrap();
    assert!(coefficient.precision().real().exact_zero());
    assert!(coefficient.precision().imaginary().exact_zero());
    assert!(warnings.is_empty());
    assert_eq!(coefficient.to_atom().unwrap(), Atom::Zero);
}

#[test]
fn standalone_polynomial_cancellation_is_exact_after_constant_normalization() {
    Vakint::initialize_vakint_symbols();
    let expression = vk_parse!("(z3+1)^2-Zeta(3)^2-2*Zeta(3)-1").unwrap();
    let (coefficient, warnings) = evaluate_coefficient(
        expression.as_view(),
        &BTreeMap::new(),
        &HashMap::new(),
        &settings(32),
    )
    .unwrap();
    assert!(coefficient.precision().real().exact_zero());
    assert!(coefficient.precision().imaginary().exact_zero());
    assert!(warnings.is_empty());
    assert_eq!(coefficient.to_atom().unwrap(), Atom::Zero);

    let result =
        NumericalEvaluationResult::from_atom(expression.as_view(), vk_symbol!("ep"), &settings(32))
            .unwrap();
    assert!(result.is_zero());
    assert!(result.warnings().is_empty());
    assert_eq!(result.to_atom(vk_symbol!("ep")).unwrap(), Atom::Zero);
}

#[test]
fn standalone_polynomial_preserves_independent_finite_literals() {
    Vakint::initialize_vakint_symbols();
    let expression = vk_parse!("(1.25`30+z3)^2-z3^2-1.5625`30-2.5`30*z3").unwrap();
    let result =
        NumericalEvaluationResult::from_atom(expression.as_view(), vk_symbol!("ep"), &settings(40))
            .unwrap();
    let coefficient = &result.coefficients()[&0];
    assert!(!coefficient.precision().real().exact_zero());
    assert!(!coefficient.precision().real().absolute_error().is_zero());
    assert!(coefficient.to_atom().is_err());
    assert!(result.warnings().iter().any(|warning| {
        warning.achieved_decimal_precision == Some(0) && warning.sources.len() == 3
    }));
}

#[test]
fn cancelled_snapshot_references_do_not_require_their_catalog() {
    Vakint::initialize_vakint_symbols();
    let master =
        vk_parse!("RustRedCoefficient(test_snapshot,RustRedMaster(test_family,1),-2)").unwrap();
    let expression =
        (&master + Atom::one()).pow(2) - master.pow(2) - &master * Atom::num(2) - Atom::one();
    let result =
        NumericalEvaluationResult::from_atom(expression.as_view(), vk_symbol!("ep"), &settings(32))
            .unwrap();
    assert!(result.is_zero());
    assert!(result.warnings().is_empty());
    assert_eq!(result.to_atom(vk_symbol!("ep")).unwrap(), Atom::Zero);
}

#[test]
fn standalone_evaluation_lists_snapshot_references_without_reading_identity_tags() {
    Vakint::initialize_vakint_symbols();
    let expression = vk_parse!("RustRedMaster(test_family,1,-1)+RustRedCoefficient(test_snapshot,RustRedMaster(other_family,2,1),-2)").unwrap();
    let error =
        NumericalEvaluationResult::from_atom(expression.as_view(), vk_symbol!("ep"), &settings(30))
            .unwrap_err()
            .to_string();
    assert!(error.contains("RustRedMaster"), "{error}");
    assert!(error.contains("RustRedCoefficient"), "{error}");
    assert!(error.contains("Vakint instance"), "{error}");
}

#[test]
fn pysecdec_normalization_tracks_independent_epsilon_coefficients() {
    Vakint::initialize_vakint_symbols();
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        number_of_terms_in_epsilon_expansion: 4,
        integral_normalization_factor: crate::LoopNormalizationFactor::Custom(
            "(1-eps)/(𝑖*𝜋^(2-eps))".into(),
        ),
        ..Default::default()
    };
    let expression = Vakint::process_pysecdec_output(
        &settings,
        1,
        7,
        &["1.0+1.0*ep".into(), "0.01+0.01*ep".into()],
    )
    .unwrap();
    let (central, error) =
        evaluate_with_errors(expression.as_view(), &HashMap::new(), &settings, None).unwrap();
    for result in [central, error.unwrap()] {
        let cancellation = &result.coefficients()[&1];
        assert!(cancellation.value().re.is_zero());
        assert!(!cancellation.precision().real().exact_zero());
        assert!(!cancellation.precision().real().absolute_error().is_zero());
        assert!(cancellation.to_atom().is_err());
        assert!(
            result
                .warnings()
                .iter()
                .any(|warning| warning.achieved_decimal_precision == Some(0))
        );
    }
}

#[test]
fn finite_custom_normalization_is_tracked_before_series_expansion() {
    Vakint::initialize_vakint_symbols();
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        number_of_terms_in_epsilon_expansion: 4,
        integral_normalization_factor: crate::LoopNormalizationFactor::Custom(
            "(0.1-0.1*eps)/(𝑖*𝜋^(2-eps))".into(),
        ),
        ..Default::default()
    };
    let expression =
        Vakint::process_pysecdec_output(&settings, 1, 8, &["1+ep".into(), "0".into()]).unwrap();
    let (result, _) =
        evaluate_with_errors(expression.as_view(), &HashMap::new(), &settings, None).unwrap();
    let cancellation = &result.coefficients()[&1];
    assert!(cancellation.value().re.is_zero());
    assert!(!cancellation.precision().real().exact_zero());
    assert!(cancellation.to_atom().is_err());
    assert!(
        result
            .warnings()
            .iter()
            .flat_map(|warning| &warning.sources)
            .any(|source| source.identity.contains("normalization_input"))
    );
}

#[test]
fn symbolic_rustred_keeps_finite_input_uncertainty_until_numerical_evaluation() {
    Vakint::initialize_vakint_symbols();
    let engine = Vakint::new().unwrap();
    let settings = VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Symbolic,
            },
        )]),
        number_of_terms_in_epsilon_expansion: 2,
        run_time_decimal_precision: 40,
        integral_normalization_factor: crate::LoopNormalizationFactor::FMFTandMATAD,
        ..Default::default()
    };
    let input = vk_parse!(format!(
        "(1.25-1.25*{})*topo(I1L(muvsq,1))",
        settings.epsilon_symbol,
    ))
    .unwrap();
    let symbolic = engine.evaluate(&settings, input.as_view()).unwrap();
    assert!(
        symbolic.contains_symbol(vk_symbol!("NumericalSource")),
        "{symbolic}"
    );
    let (numerical, _) = engine
        .numerical_evaluation(
            &settings,
            symbolic.as_view(),
            &ahash::HashMap::from_iter([("muvsq".into(), settings.real_to_prec("1"))]),
            &ahash::HashMap::default(),
            None,
        )
        .unwrap();
    let coefficient = &numerical.coefficients()[&0];
    assert!(coefficient.value().re.is_zero());
    assert!(!coefficient.precision().real().exact_zero());
    assert!(coefficient.to_atom().is_err());
}

#[test]
fn finite_source_shortfall_warns_and_returns_achieved_accuracy() {
    Vakint::initialize_vakint_symbols();
    let sources = BTreeMap::from([source("precision_limited", "2.34567890123456789", 16)]);
    let (coefficient, warnings) = evaluate_coefficient(
        Atom::var(vk_symbol!("precision_limited")).as_view(),
        &sources,
        &HashMap::new(),
        &settings(100),
    )
    .unwrap();
    assert!(coefficient.precision().real().significant_digits().unwrap() < 20);
    assert!(
        warnings
            .iter()
            .any(|warning| warning.requested_decimal_precision == 100
                && warning.sources[0].available_decimal_precision == Some(16))
    );
    assert!(coefficient.to_real_float().is_ok());
}

#[test]
fn a_tighter_bound_does_not_promote_the_declared_source_cap() {
    Vakint::initialize_vakint_symbols();
    let (symbol, mut finite) = source("precision_guard_digits", "1.23456789", 30);
    finite.absolute_error = Some("1e-38".into());
    finite.error_is_bounded = true;
    let sources = BTreeMap::from([(symbol.clone(), finite)]);
    let (coefficient, warnings) =
        evaluate_coefficient(symbol.as_view(), &sources, &HashMap::new(), &settings(40)).unwrap();
    assert!(coefficient.precision().real().significant_digits().unwrap() <= 30);
    assert!(warnings.iter().any(|warning| {
        warning
            .sources
            .iter()
            .any(|source| source.available_decimal_precision == Some(30))
    }));
}

#[test]
fn exact_constants_have_no_artificial_requested_digit_shortfall() {
    Vakint::initialize_vakint_symbols();
    for digits in [16, 32, 200] {
        let expression = vk_parse!("z3+sqrt3+log(2)").unwrap();
        let (coefficient, warnings) = evaluate_coefficient(
            expression.as_view(),
            &BTreeMap::new(),
            &HashMap::new(),
            &settings(digits),
        )
        .unwrap();
        assert!(coefficient.precision().real().significant_digits().unwrap() >= digits);
        assert!(coefficient.precision().imaginary().exact_zero());
        assert!(
            warnings.is_empty(),
            "unexpected precision warnings: {warnings:?}"
        );
    }
}

#[test]
fn unused_low_precision_input_does_not_limit_output() {
    Vakint::initialize_vakint_symbols();
    let parameters = HashMap::from([(
        Atom::var(vk_symbol!("precision_unused")),
        Complex::new(Float::with_val(12, 3), Float::new(12)),
    )]);
    let expression = vk_parse!("1/7").unwrap();
    let (coefficient, warnings) = evaluate_coefficient(
        expression.as_view(),
        &BTreeMap::new(),
        &parameters,
        &settings(100),
    )
    .unwrap();
    assert!(coefficient.precision().real().significant_digits().unwrap() >= 100);
    assert!(warnings.is_empty());
}

#[test]
fn external_contractions_track_cancellation_before_becoming_floats() {
    Vakint::initialize_vakint_symbols();
    let externals = HashMap::<_, _, ahash::RandomState>::from_iter([(
        1,
        crate::Momentum::Real((
            Float::with_val(150, 1),
            Float::with_val(150, 1),
            Float::new(150),
            Float::new(150),
        )),
    )]);
    let expression = vk_parse!("dot(p(1),p(1))*ep^-2").unwrap();
    let (expression, inputs) = evaluation_inputs(
        expression.as_view(),
        vk_symbol!("ep"),
        &HashMap::default(),
        &HashMap::default(),
        Some(&externals),
    )
    .unwrap();
    let result = evaluate_series(
        expression.as_view(),
        vk_symbol!("ep"),
        &inputs,
        &settings(40),
        None,
    )
    .unwrap();
    let coefficient = result.coefficients().get(&-2).unwrap();
    assert!(coefficient.value().re.is_zero());
    assert!(!coefficient.precision().real().exact_zero());
    assert!(result.to_atom(vk_symbol!("ep")).is_err());
}

#[test]
fn f64_helpers_preserve_their_input_precision() {
    let inputs = HashMap::<_, _, ahash::RandomState>::from_iter([("precision_f64".into(), 1.25)]);
    assert_eq!(
        crate::params_from_f64(&inputs, 200)["precision_f64"].prec(),
        53
    );
    let inputs = HashMap::<_, _, ahash::RandomState>::from_iter([(
        "precision_f64".into(),
        Complex::new(1.25, 2.5),
    )]);
    let result = crate::params_from_complex_f64(&inputs, 200);
    assert_eq!(result["precision_f64"].re.prec(), 53);
    assert_eq!(result["precision_f64"].im.prec(), 53);
}

#[test]
fn partial_evaluation_preserves_exact_real_and_gaussian_rationals() {
    Vakint::initialize_vakint_symbols();
    let expression =
        vk_parse!("x/3+(2/7+5*𝑖/11)*y+z/123456789012345678901234567890123456789").unwrap();
    let partial = Vakint::partial_numerical_evaluation(
        &settings(50),
        expression.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )
    .unwrap();
    assert_eq!(partial.expand(), expression.expand());
    partial.visitor(&mut |part| {
        if let AtomView::Num(number) = part {
            assert!(!matches!(
                number.get_coeff_view(),
                symbolica::coefficient::CoefficientView::Float(..)
            ));
        }
        true
    });

    let constant = partial
        .replace(vk_parse!("x").unwrap())
        .with(Atom::one())
        .replace(vk_parse!("y").unwrap())
        .with(Atom::one())
        .replace(vk_parse!("z").unwrap())
        .with(Atom::one());
    let (numeric, _) = Vakint::full_numerical_evaluation(
        &settings(50),
        constant.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )
    .unwrap();
    let exported = numeric.coefficients()[&0].to_atom().unwrap();
    let AtomView::Num(number) = exported.as_view() else {
        panic!("final numerical evaluation must export a numerical coefficient")
    };
    assert!(matches!(
        number.get_coeff_view(),
        symbolica::coefficient::CoefficientView::Float(..)
    ));
}

#[test]
fn partial_evaluation_does_not_promote_a_float_one_to_an_exact_rational() {
    Vakint::initialize_vakint_symbols();
    let expression = vk_parse!("1.0`20*x").unwrap();
    let partial = Vakint::partial_numerical_evaluation(
        &settings(80),
        expression.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )
    .unwrap();
    let value = partial.replace(vk_parse!("x").unwrap()).with(Atom::one());
    let AtomView::Num(number) = value.as_view() else {
        panic!("expected the explicitly finite input coefficient")
    };
    let symbolica::coefficient::CoefficientView::Float(real, _) = number.get_coeff_view() else {
        panic!("a Float midpoint of one must remain approximate")
    };
    assert_eq!(real.to_float().to_rational(), 1);
    assert!(real.to_float().prec() <= Float::parse("1.0`20", None).unwrap().prec());

    // Keep the literal inside an unexpanded polynomial: constructing
    // `1.0*x-x` already combines its midpoint before this API receives it.
    let cancellation = vk_parse!("(1.0`20+x)^2-x^2-2*x-1").unwrap();
    assert!(!cancellation.is_zero());
    let error = Vakint::partial_numerical_evaluation(
        &settings(80),
        cancellation.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )
    .unwrap_err();
    assert!(error.to_string().contains("no reliable"), "{error}");
}

#[test]
fn no_decimal_digits_cannot_be_exported_as_an_atom() {
    Vakint::initialize_vakint_symbols();
    let coefficient =
        NumericalCoefficient::from_value(Complex::new(Float::with_val(3, 2), Float::new(3)), 16);
    assert!(coefficient.to_atom().is_err());
}

#[test]
fn requested_precision_does_not_promote_a_weak_complex_component() {
    Vakint::initialize_vakint_symbols();
    let coefficient = NumericalCoefficient::from_value(
        Complex::new(
            Float::parse("1.234567890123456789", Some(200)).unwrap(),
            Float::parse("2.3456789", Some(60)).unwrap(),
        ),
        40,
    );
    let value = coefficient.as_tracked().to_float_with_precision().unwrap();
    assert!(value.re.prec() > value.im.prec());
    assert!(value.im.prec() <= 60);
}

#[test]
fn tracked_imaginary_alias_preserves_symbolic_input_and_evaluates_scalars_and_series() {
    let engine = Vakint::new().unwrap();
    let unit = Atom::var(crate::symbols::S.cmplx_i);
    assert_ne!(unit.clone().pow(2), Atom::num(-1));
    let epsilon = Atom::var(vk_symbol!("ε"));
    for digits in [16, 80] {
        let settings = settings(digits);
        let (scalar, error) = engine
            .numerical_evaluation(
                &settings,
                unit.as_view(),
                &Default::default(),
                &Default::default(),
                None,
            )
            .unwrap();
        assert!(error.is_none());
        let value = &scalar.coefficients()[&0];
        assert!(value.precision().real().exact_zero());
        assert_eq!(value.value().im.to_rational(), 1);
        assert!(value.precision().imaginary().absolute_error().is_zero());
        assert!(scalar.warnings().is_empty());

        let expression =
            &unit * epsilon.clone().pow(-2) + Atom::one() + unit.clone().pow(2) * &epsilon;
        let (series, error) = engine
            .numerical_evaluation(
                &settings,
                expression.as_view(),
                &Default::default(),
                &Default::default(),
                None,
            )
            .unwrap();
        assert!(error.is_none());
        assert_eq!(
            series.coefficients().keys().copied().collect::<Vec<_>>(),
            [-2, 0, 1]
        );
        assert!(series.coefficients()[&-2].precision().real().exact_zero());
        assert_eq!(series.coefficients()[&-2].value().im.to_rational(), 1);
        assert_eq!(series.coefficients()[&0].value().re.to_rational(), 1);
        assert!(
            series.coefficients()[&0]
                .precision()
                .imaginary()
                .exact_zero()
        );
        assert_eq!(series.coefficients()[&1].value().re.to_rational(), -1);
        assert!(
            series.coefficients()[&1]
                .precision()
                .imaginary()
                .exact_zero()
        );
        assert!(series.warnings().is_empty());
    }
    assert_ne!(unit.clone().pow(2), Atom::num(-1));
}

#[test]
fn tracked_imaginary_alias_default_precedes_explicit_parameter_overrides() {
    let unit = Atom::var(crate::symbols::S.cmplx_i);
    for (real, imaginary) in [(7, 0), (3, 4)] {
        let parameters = HashMap::from([(
            unit.clone(),
            Complex::new(Float::with_val(256, real), Float::with_val(256, imaginary)),
        )]);
        let (value, _) = evaluate_coefficient(
            (&unit + Atom::one()).as_view(),
            &BTreeMap::new(),
            &parameters,
            &settings(32),
        )
        .unwrap();
        assert_eq!(value.value().re.to_rational(), real + 1);
        assert_eq!(value.value().im.to_rational(), imaginary);
    }
}

#[test]
fn tracked_imaginary_alias_preserves_independent_component_accuracy() {
    let (identity, finite) = source("imaginary_alias_finite", "1.23456789012345", 12);
    let sources = BTreeMap::from([(identity.clone(), finite)]);
    let expression = Atom::num(2) + Atom::var(crate::symbols::S.cmplx_i) * identity;
    let (value, warnings) = evaluate_coefficient(
        expression.as_view(),
        &sources,
        &HashMap::new(),
        &settings(50),
    )
    .unwrap();
    assert_eq!(value.value().re.to_rational(), 2);
    assert!(value.precision().real().absolute_error().is_zero());
    assert!(!value.precision().imaginary().absolute_error().is_zero());
    assert!(value.precision().imaginary().significant_digits().unwrap() < 50);
    assert!(warnings.iter().any(|warning| {
        warning.requested_decimal_precision == 50
            && warning
                .sources
                .iter()
                .any(|source| source.identity == "imaginary_alias_finite")
    }));
}

#[test]
fn tracked_imaginary_alias_does_not_retry_a_limited_tracked_callback() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    use symbolica::atom::EvaluationInfo;

    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let limited = symbolica::symbol!(
        "vakint_test::imaginary_alias_limited_callback",
        eval = EvaluationInfo::new().register_tracked(move |tags, arguments, bits| {
            if !tags.is_empty() || !arguments.is_empty() {
                return Err("test constant accepts no tags or arguments".into());
            }
            observed.fetch_add(1, Ordering::SeqCst);
            TrackedComplex::from_source_decimal("1.25", 8, bits).map_err(|error| error.to_string())
        })
    );
    let expression = Atom::var(crate::symbols::S.cmplx_i) * Atom::var(limited);
    let (value, warnings) = evaluate_coefficient(
        expression.as_view(),
        &BTreeMap::new(),
        &HashMap::new(),
        &settings(60),
    )
    .unwrap();
    assert_eq!(calls.load(Ordering::SeqCst), 1);
    assert!(value.precision().real().exact_zero());
    assert!(value.precision().imaginary().significant_digits().unwrap() < 60);
    assert!(
        warnings
            .iter()
            .any(|warning| warning.achieved_decimal_precision.is_some())
    );
}
