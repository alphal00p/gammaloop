use super::*;

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
