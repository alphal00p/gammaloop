use vakint::{
    RustRedOptions, TensorReductionError, TensorReductionMode, Vakint, VakintSettings, vakint_parse,
};

#[test]
fn tensor_reducer_defaults_to_the_form_backend() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings::default();
    let reducer = vakint.tensor_reducer(&settings);

    assert_eq!(reducer.selected_mode(), &TensorReductionMode::Form);
    assert_eq!(
        reducer
            .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
            .selected_mode(),
        &TensorReductionMode::RustRed(RustRedOptions::new())
    );
}

#[test]
fn rustred_rank_two_projection_uses_the_matcher_and_never_touches_form() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(99,mu)*k(99,nu)*topo(prop(7,edge(3,3),k(99),muvsq,1))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(1))*g(mu,nu)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}

#[test]
fn rustred_rank_two_projection_respects_vakints_indexed_output_setting() {
    let vakint = Vakint::new().unwrap();
    let rustred_settings = VakintSettings {
        form_exe_path: "/this/path/must/not_be_invoked_by_rustred".into(),
        use_dot_product_notation: false,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*k(1,nu)*topo(I1L(muvsq,1))").unwrap();

    let reduced = vakint
        .tensor_reducer(&rustred_settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    // Frozen after an independent FORM 5.0.0 compatibility run. This
    // RustRed-mode test deliberately has no executable FORM dependency.
    let expected = vakint_parse!(
        "-g(mu,nu)*k(1,dot_dummy_ind(2))^2*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )/(-4+2*ε)"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}

#[test]
fn rustred_unsupported_family_is_typed_and_never_falls_back_to_form() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*k(1,nu)*topo(I2L(muvsq,1,1,1))").unwrap();

    let error = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap_err();

    assert!(matches!(
        error,
        TensorReductionError::RustRedUnsupportedFamily { loop_count: 2, .. }
    ));
}

#[test]
fn rustred_rejects_non_bare_explicit_input_routings_before_canonicalization() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    let inputs = [
        vakint_parse!("k(99,mu)^2*topo(prop(7,edge(3,3),0,muvsq,1))").unwrap(),
        vakint_parse!("k(99,mu)^2*topo(prop(7,edge(3,3),2*k(99),muvsq,1))").unwrap(),
        vakint_parse!("k(99,mu)^2*topo(prop(7,edge(3,3),k(99)+k(98),muvsq,1))").unwrap(),
        vakint_parse!("k(99,mu)^2*topo(prop(7,edge(3,3),k(99)+p(1),muvsq,1))").unwrap(),
    ];

    for input in inputs {
        let error = vakint
            .tensor_reducer(&settings)
            .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
            .reduce(input.as_view())
            .unwrap_err();
        assert!(matches!(
            error,
            TensorReductionError::RustRedUnsupportedMomentum { .. }
        ));
    }
}

#[test]
fn rustred_rejects_a_zero_scale_with_a_precise_adapter_error() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*k(1,nu)*topo(prop(1,edge(1,1),k(1),0,1))").unwrap();

    let error = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap_err();

    assert!(matches!(
        error,
        TensorReductionError::RustRedUnsupportedMass { .. }
    ));
}

#[test]
fn rustred_rejects_a_nonsymbolic_epsilon_before_projection() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        epsilon_symbol: "2".to_owned(),
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)^2*topo(I1L(muvsq,1))").unwrap();

    let error = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap_err();

    assert!(matches!(
        error,
        TensorReductionError::RustRedInvalidDimension { .. }
    ));
}

#[test]
fn rustred_adapter_preserves_scalar_rank_and_exposes_the_even_rank_frontier() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };

    let scalar = vakint_parse!("spectator(x)*topo(I1L(muvsq,1))").unwrap();
    let scalar_reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(scalar.as_view())
        .unwrap();
    let scalar_expected =
        vakint_parse!("spectator(x)*topo(prop(1,edge(1,1),k(1),muvsq,1))").unwrap();
    assert_eq!(scalar_reduced, scalar_expected);

    for odd in [
        vakint_parse!("k(1,mu)*topo(I1L(muvsq,1))").unwrap(),
        vakint_parse!("k(1,mu)*k(1,nu)*k(1,rho)*topo(I1L(muvsq,1))").unwrap(),
    ] {
        assert!(
            vakint
                .tensor_reducer(&settings)
                .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
                .reduce(odd.as_view())
                .unwrap()
                .is_zero()
        );
    }

    let rank_four =
        vakint_parse!("k(1,mu)*k(1,nu)*k(1,rho)*k(1,sigma)*topo(I1L(muvsq,1))").unwrap();
    assert!(matches!(
        vakint
            .tensor_reducer(&settings)
            .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
            .reduce(rank_four.as_view()),
        Err(TensorReductionError::RustRedTensor { .. })
    ));
}

#[test]
fn rustred_adapter_retains_exact_numeric_mass_and_nonunit_power() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(7,mu)*k(7,nu)*topo(prop(5,edge(4,4),k(7),3/2,-2))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(1))*g(mu,nu)*topo(\
            prop(1,edge(1,1),k(1),3/2,-2)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}
