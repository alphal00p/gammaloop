use super::*;

#[test]
fn exact_complex_constant_denominators_preserve_real_master_coefficients() {
    crate::Vakint::initialize_vakint_symbols();
    let source = vk_parse!(
        "(symbolica::polylog(2,(1+3^(1/2)*𝑖)/2)+symbolica::polylog(2,(1-3^(1/2)*𝑖)/2))*(3^(1/2)-𝑖)^-1*(3^(1/2)+𝑖)^-1"
    )
    .unwrap();
    let normalized = normalize_constants(source.as_view()).unwrap();
    assert_eq!(
        normalized,
        vk_parse!("PolylogReal(2,(1+3^(1/2)*𝑖)/2)/2").unwrap()
    );
    let settings = VakintSettings {
        run_time_decimal_precision: 100,
        ..Default::default()
    };
    let (number, warnings) = crate::precision::evaluate_coefficient(
        source.as_view(),
        &BTreeMap::new(),
        &std::collections::HashMap::new(),
        &settings,
    )
    .unwrap();
    assert!(number.precision().imaginary().exact_zero());
    assert!(number.real_midpoint().is_ok());
    assert!(warnings.is_empty());
}

#[test]
fn constant_denominator_normalization_retains_finite_input_accuracy() {
    crate::Vakint::initialize_vakint_symbols();
    let source = vk_parse!("1.25`40*(3^(1/2)-𝑖)^-1*(3^(1/2)+𝑖)^-1").unwrap();
    // Raw finite values must be staged by the tracked evaluator first.
    assert_eq!(
        normalize_constant_denominators(source.clone()).unwrap(),
        source
    );
    let inert = vk_parse!("inert((3^(1/2)-𝑖)^-1*(3^(1/2)+𝑖)^-1)").unwrap();
    assert_eq!(
        normalize_constant_denominators(inert.clone()).unwrap(),
        inert
    );
    let settings = VakintSettings {
        run_time_decimal_precision: 100,
        ..Default::default()
    };
    let (number, warnings) = crate::precision::evaluate_coefficient(
        source.as_view(),
        &BTreeMap::new(),
        &std::collections::HashMap::new(),
        &settings,
    )
    .unwrap();
    assert!(number.precision().imaginary().exact_zero());
    let digits = number.precision().real().significant_digits().unwrap();
    assert!((30..=40).contains(&digits));
    assert!(!warnings.is_empty());
}

#[test]
fn pr14_endpoint_one_projections_use_checked_native_evaluators() {
    crate::Vakint::initialize_vakint_symbols();
    // Nine new components in the authenticated PR14 epsilon-one formulas.
    // Their numerical values are independently checked in Symbolica's kernel
    // tests; this checks Vakint's exact letter recognition and public export.
    let projections = [
        "GPLImag(0,(1-3^(1/2)*𝑖)/2,-1,1)",
        "GPLImag(0,0,0,1,(1-3^(1/2)*𝑖)/2,1)",
        "GPLReal(0,0,(1-3^(1/2)*𝑖)/2,-1,1)",
        "GPLReal((1-3^(1/2)*𝑖)/2,-1,1)",
        "GPLReal(0,0,1,1,-1,1)",
        "GPLReal(0,0,1,(1-3^(1/2)*𝑖)/2,-1,1)",
        "GPLReal(0,0,1,(1-3^(1/2)*𝑖)/2,(-1+3^(1/2)*𝑖)/2,1)",
        "GPLReal(0,0,1,(1-3^(1/2)*𝑖)/2,(-1-3^(1/2)*𝑖)/2,1)",
        "GPLReal(0,0,(1-3^(1/2)*𝑖)/2,1,-1,1)",
    ];
    for source in projections {
        let source = vk_parse!(source).unwrap();
        let native = to_symbolica_expression(source.as_view()).unwrap();
        let value = native
            .evaluate_tracked(
                &std::collections::HashMap::<Atom, symbolica::domains::float::TrackedComplex>::new(
                ),
                256,
            )
            .unwrap();
        assert!(value.precision().real.significant_digits().unwrap() >= 60);
        assert!(value.precision().imaginary.exact_zero());
    }
    for invalid in ["GPLReal(0,2,1)", "GPLReal(1,-1,1)", "GPLImag(-1,0,1)"] {
        assert!(to_symbolica_expression(vk_parse!(invalid).unwrap().as_view()).is_err());
    }
}

#[test]
fn historical_and_native_hpl_notation_share_exact_values() {
    crate::Vakint::initialize_vakint_symbols();
    let input = vk_parse!("cl2+HI(1,0,om1)").unwrap();
    assert!(normalize_constants(input.as_view()).unwrap().is_zero());
    let input = vk_parse!("9*Sqrt(3)*S2-4*cl2").unwrap();
    assert!(
        normalize_constants(input.as_view())
            .unwrap()
            .together()
            .is_zero()
    );
    let input = vk_parse!("HR(1,0,om1)+HI(1,0,(1-3^(1/2)*𝑖)/2)").unwrap();
    let native = normalize_constants(input.as_view()).unwrap();
    assert_eq!(normalize_constants(native.as_view()).unwrap(), native);
    let foreign = vk_parse!("user_space::HI(2,x)+user_space::z3").unwrap();
    assert_eq!(normalize_constants(foreign.as_view()).unwrap(), foreign);
}

#[test]
fn numerical_replacement_preserves_laurent_and_spectator_factors() {
    crate::Vakint::initialize_vakint_symbols();
    let settings = VakintSettings {
        run_time_decimal_precision: 100,
        ..Default::default()
    };
    let input = vk_parse!("user_space::spectator/ep^3*HI(1,0,om1)").unwrap();
    let value = substitute_numerically(input.as_view(), &settings).unwrap();
    let constant =
        substitute_numerically(vk_parse!("HI(1,0,om1)").unwrap().as_view(), &settings).unwrap();
    assert_eq!(
        value,
        vk_parse!("user_space::spectator/ep^3").unwrap() * constant
    );
}

#[test]
fn nonlinear_real_hpl_projections_have_exactly_zero_imaginary_part() {
    crate::Vakint::initialize_vakint_symbols();
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        ..Default::default()
    };
    let input =
        vk_parse!("(HR(1,-1,0,om1)+3*HI(0,1,1,om1))^3 + 7*HR(1,0,om1)*HI(1,-1,0,om1)^2").unwrap();
    let symbolic = normalize_constants(input.as_view()).unwrap().expand();
    assert!(symbolic.contains_symbol(vk_symbol!("HR")));
    assert!(symbolic.contains_symbol(vk_symbol!("HI")));
    assert!(!symbolic.contains_symbol(hpl_at_sixth_root()));
    let expanded = substitute_numerically(symbolic.as_view(), &settings).unwrap();
    let unexpanded = substitute_numerically(input.as_view(), &settings).unwrap();
    let parameters = ahash::HashMap::<Atom, Complex<Float>>::default();
    for expression in [expanded, unexpanded] {
        let value = expression
            .evaluate_with_prec(&parameters, settings.get_binary_precision())
            .unwrap();
        assert!(value.im.as_raw().is_zero(), "{value}");
    }
    let conjugate = vk_parse!("HR(1,0,(1-3^(1/2)*𝑖)/2)+HI(1,0,(1-3^(1/2)*𝑖)/2)").unwrap();
    assert_eq!(
        normalize_constants(conjugate.as_view()).unwrap(),
        vk_parse!("HR(1,0,om1)-HI(1,0,om1)").unwrap()
    );
}

#[test]
fn exact_tagged_constants_are_evaluated_without_visiting_their_tags() {
    crate::Vakint::initialize_vakint_symbols();
    let native = sixth_root_hpl(&[1, 0], 1).unwrap();
    // A tagged evaluator is more capable than the conservative structural predicate.
    assert!(!native.is_constant());
    let settings = VakintSettings {
        run_time_decimal_precision: 100,
        ..Default::default()
    };
    let value = substitute_numerically(native.as_view(), &settings).unwrap();
    assert!(matches!(value.as_view(), AtomView::Num(_)));
    let conjugate =
        substitute_numerically(sixth_root_hpl(&[1, 0], -1).unwrap().as_view(), &settings).unwrap();
    assert_eq!(value.conj(), conjugate);
}

#[test]
fn malformed_hpls_and_nonfinite_constants_fail_gracefully() {
    crate::Vakint::initialize_vakint_symbols();
    for source in ["HI(2,om1)", "HR(1,0,1/2)", "HI(1,1,1,1,1,1,1,om1)"] {
        assert!(
            normalize_constants(vk_parse!(source).unwrap().as_view()).is_err(),
            "{source}"
        );
    }
    assert!(
        substitute_numerically(
            vk_parse!("Zeta(1)").unwrap().as_view(),
            &VakintSettings::default()
        )
        .is_err()
    );
}

#[test]
fn elementary_aliases_use_native_arbitrary_precision_functions() {
    crate::Vakint::initialize_vakint_symbols();
    let input = vk_parse!("z2+z3+z4+z5+z6+sqrt3+Sqrt(5)").unwrap();
    let expected=vk_parse!("symbolica::zeta(2)+symbolica::zeta(3)+symbolica::zeta(4)+symbolica::zeta(5)+symbolica::zeta(6)+3^(1/2)+5^(1/2)").unwrap();
    assert!(
        (normalize_constants(input.as_view()).unwrap() - expected)
            .expand()
            .is_zero()
    );
}

#[test]
fn conjugate_polylogarithm_master_formulas_remain_exactly_real() {
    crate::Vakint::initialize_vakint_symbols();
    let positive = vk_parse!("symbolica::polylog(2,(1+3^(1/2)*𝑖)/2)").unwrap();
    let negative = vk_parse!("symbolica::polylog(2,(1-3^(1/2)*𝑖)/2)").unwrap();
    let difference = &positive - &negative;
    let formula = Atom::i() * &difference + difference.pow(2) + positive * negative;
    let normalized =
        crate::rustred_evaluation::combine_exact(normalize_constants(formula.as_view()).unwrap());
    assert_eq!(
        normalize_constants(normalized.as_view()).unwrap(),
        normalized
    );
    let settings = VakintSettings {
        run_time_decimal_precision: 40,
        ..Default::default()
    };
    let (number, _) = crate::precision::evaluate_coefficient(
        normalized.as_view(),
        &BTreeMap::new(),
        &std::collections::HashMap::new(),
        &settings,
    )
    .unwrap();
    assert!(number.precision().imaginary().exact_zero());
    assert!(number.to_real_float().is_ok());
    assert!(projection_kernel(vk_parse!("PolylogReal(2,2)").unwrap().as_view()).is_err());
}

#[test]
fn exported_period_formulas_use_standalone_symbolica_evaluators() {
    crate::Vakint::initialize_vakint_symbols();
    let input =
        vk_parse!("HR(1,0,om1)+GPLReal(0,0,(-1-3^(1/2)*𝑖)/2,1,1)+PolylogImag(2,(1+3^(1/2)*𝑖)/2)")
            .unwrap();
    let exported = to_symbolica_expression(input.as_view()).unwrap();
    for name in [
        "HR",
        "HI",
        "GPLReal",
        "GPLImag",
        "PolylogReal",
        "PolylogImag",
    ] {
        assert!(!exported.contains_symbol(vk_symbol!(name)));
    }
    let value = exported
        .evaluate_with_prec(
            &std::collections::HashMap::<Atom, Complex<Float>>::new(),
            256,
        )
        .unwrap();
    assert!(value.re.as_raw().is_finite());
    assert!(value.im.as_raw().is_zero());
    let tracked = exported
        .evaluate_tracked(
            &std::collections::HashMap::<Atom, symbolica::domains::float::TrackedComplex>::new(),
            256,
        )
        .unwrap();
    assert!(tracked.precision().imaginary().exact_zero());
}

#[test]
fn inverse_sqrt_three_trilogarithms_use_the_exact_hpl_transform() {
    crate::Vakint::initialize_vakint_symbols();
    let point = vk_parse!("3^(1/2)*𝑖/3").unwrap();
    let projection = vk_symbol!("PolylogReal").call((3, point.clone()));
    let (kernel, imaginary) = projection_kernel(projection.as_view()).unwrap().unwrap();
    assert!(!imaginary);
    assert!(!kernel.contains_symbol(polylog()));
    let parameters = std::collections::HashMap::<Atom, Complex<Float>>::new();
    let actual = kernel.evaluate_with_prec(&parameters, 720).unwrap();
    let independent = polylog()
        .call((3, point))
        .evaluate_with_prec(&parameters, 840)
        .unwrap();
    let tolerance = Float::parse("1e-200", Some(720)).unwrap();
    assert!((actual.re.clone() - &independent.re).as_raw().clone().abs() < *tolerance.as_raw());
    assert!((actual.im.clone() - &independent.im).as_raw().clone().abs() < *tolerance.as_raw());
}

#[test]
fn scalar_projection_products_are_not_interpreted_as_indexed_vectors() {
    crate::Vakint::initialize_vakint_symbols();
    let formula = vk_parse!("PolylogImag(2,(1+3^(1/2)*𝑖)/2)^2 + HR(1,om1)*HI(1,om1)").unwrap();
    let converted = crate::Vakint::convert_to_dot_notation(formula.as_view());
    assert_eq!(converted, formula);
    assert!(!converted.contains_symbol(crate::symbols::S.dot));
}
