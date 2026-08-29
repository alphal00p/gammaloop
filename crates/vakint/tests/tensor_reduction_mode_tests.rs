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
fn rustred_disjoins_overlapping_vector_ids_and_normalizes_indexed_contractions() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    // p(2) occurs only in an already scalar dot. The repeated rho index must
    // become dot(k(1),p(1)) before RustRed sees the numerator, while the two
    // different vectors sharing numeric ID 1 remain distinct.
    let input =
        vakint_parse!("dot(p(1),p(2))*p(1,rho)*k(1,rho)*k(1,mu)*topo(I1L(muvsq,1))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(1))*dot(p(1),p(2))*p(1,mu)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
    let rendered = reduced.to_string();
    assert!(!rendered.contains("rustred_loop_label"));
    assert!(!rendered.contains("rustred_external_label"));
}

#[test]
fn rustred_normalized_contraction_respects_indexed_output_without_label_leaks() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: false,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("p(1,rho)*k(1,rho)*k(1,mu)*topo(I1L(muvsq,1))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-p(1,mu)*k(1,dot_dummy_ind(2))^2*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)\
        )/(-4+2*ε)"
    )
    .unwrap();

    assert_eq!(reduced, expected);
    let rendered = reduced.to_string();
    assert!(!rendered.contains("rustred_loop_label"));
    assert!(!rendered.contains("rustred_external_label"));
}

#[test]
fn rustred_disjoins_two_loop_overlapping_ids_in_indexed_contractions() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    // The repeated indices normalize to dot(k(1),p(2))*dot(k(2),p(3)).
    // In particular, k(2) and p(2) must remain different vectors in the
    // single RustRed label namespace.
    let input =
        vakint_parse!("k(1,alpha)*k(2,beta)*p(2,alpha)*p(3,beta)*topo(I2L(muvsq,1,2,1))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(2))*dot(p(2),p(3))*topo(\
            prop(1,edge(1,2),k(1),muvsq,1)*\
            prop(2,edge(1,2),k(2),muvsq,2)*\
            prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
    let rendered = reduced.to_string();
    assert!(!rendered.contains("rustred_loop_label"));
    assert!(!rendered.contains("rustred_external_label"));
}

#[test]
fn rustred_projects_two_loop_mixed_momenta_without_form() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*k(2,nu)*topo(I2L(muvsq,1,2,1))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(2))*g(mu,nu)*topo(\
            prop(1,edge(1,2),k(1),muvsq,1)*\
            prop(2,edge(1,2),k(2),muvsq,2)*\
            prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}

#[test]
fn rustred_projects_an_isp_completed_two_loop_pinch_without_form() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*k(2,nu)*topo(I2L_pinch_3(muvsq,1,1,0))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(2))*g(mu,nu)*topo(\
            prop(1,edge(1,1),k(1),muvsq,1)*\
            prop(2,edge(1,1),k(2),muvsq,1)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}

#[test]
fn rustred_projects_an_isp_completed_three_loop_pinch_without_form() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let input = vakint_parse!("k(1,mu)*k(3,nu)*topo(I3L_pinch_1_6(muvsq,0,1,1,1,2,0))").unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(3))*g(mu,nu)*topo(\
            prop(2,edge(1,3),k(1),muvsq,1)*\
            prop(3,edge(3,1),k(2),muvsq,1)*\
            prop(4,edge(1,3),k(3),muvsq,1)*\
            prop(5,edge(1,3),-k(1)+k(2)-k(3),muvsq,2)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}

#[test]
fn rustred_rejects_caller_spelled_private_momentum_labels() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    let attacks = [
        (
            vakint_parse!("dot(rustred_loop_label(1),p(1))*k(1,mu)*topo(I1L(muvsq,1))").unwrap(),
            "rustred_loop_label",
        ),
        (
            vakint_parse!("dot(k(1),rustred_external_label(1))*k(1,mu)*topo(I1L(muvsq,1))")
                .unwrap(),
            "rustred_external_label",
        ),
        (
            vakint_parse!("k(rustred_loop_label(1),mu)*k(1,nu)*topo(I1L(muvsq,1))").unwrap(),
            "rustred_loop_label",
        ),
        (
            vakint_parse!("p(rustred_external_label(1),mu)*k(1,nu)*topo(I1L(muvsq,1))").unwrap(),
            "rustred_external_label",
        ),
    ];

    for (input, expected_head) in attacks {
        let original = input.clone();
        let error = vakint
            .tensor_reducer(&settings)
            .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
            .reduce(input.as_view())
            .unwrap_err();

        match error {
            TensorReductionError::RustRedReservedMomentumLabel { term, head } => {
                assert_eq!(term, 0);
                assert_eq!(head, expected_head);
            }
            other => panic!("expected a reserved-label error, got {other}"),
        }
        assert_eq!(input, original);
    }
}

#[test]
fn rustred_accepts_a_replayable_alternate_two_loop_routing() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    // Two edges expose the input basis and the third is their difference.
    // Vakint's chosen graph match supplies the signed canonical basis map;
    // replaying all three squared momenta proves the routing equivalent.
    let input = vakint_parse!(
        "k(11,mu)*k(11,nu)*topo(\
            prop(9,edge(7,10),k(11),muvsq,1)*\
            prop(33,edge(7,10),-k(11)+k(22),muvsq,2)*\
            prop(55,edge(10,7),k(22),muvsq,1)\
        )"
    )
    .unwrap();

    let reduced = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap();
    let expected = vakint_parse!(
        "-(2*ε-4)^-1*dot(k(1),k(1))*g(mu,nu)*topo(\
            prop(1,edge(1,2),k(1),muvsq,1)*\
            prop(2,edge(1,2),k(2),muvsq,1)*\
            prop(3,edge(2,1),k(1)+k(2),muvsq,2)\
        )"
    )
    .unwrap();

    assert_eq!(reduced, expected);
}

#[test]
fn rustred_rejects_a_non_equivalent_two_loop_routing_before_canonicalization() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred".into(),
        ..VakintSettings::default()
    };
    // With the same basis witness as the valid alternate routing, the middle
    // edge would replay to 2*k(1)+k(2), not the matched canonical k(2).
    let input = vakint_parse!(
        "k(11,mu)^2*topo(\
            prop(9,edge(7,10),k(11),muvsq,1)*\
            prop(33,edge(7,10),k(11)+k(22),muvsq,2)*\
            prop(55,edge(10,7),k(22),muvsq,1)\
        )"
    )
    .unwrap();

    let error = vakint
        .tensor_reducer(&settings)
        .mode(TensorReductionMode::RustRed(RustRedOptions::new()))
        .reduce(input.as_view())
        .unwrap_err();

    assert!(matches!(
        error,
        TensorReductionError::RustRedUnsupportedMomentum { propagator: 33, .. }
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
