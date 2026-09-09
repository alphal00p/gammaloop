mod test_utils;
use test_utils::{compare_output, get_vakint};
use vakint::VakintSettings;
use vakint::vakint_parse;

#[test_log::test]
fn test_reduction_1l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });

    let integral = vakint
        .to_canonical(
            vakint_parse!(
                "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        vakint_parse!(
            "(\
                -(2*ε-4)^-1*dot(k(1),k(1))*g(1,2)\
            )*topo(I1L(muvsq,1))"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_reduction_1l_b() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });

    let integral = vakint
        .to_canonical(
            vakint_parse!(
                "((k(1,1)*k(1,2))^2*g(1,2)+k(1,3)*p(1,3)+k(1,1)*k(1,2)*p(2,1)*p(3,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        vakint_parse!(
            "(\
                dot(k(1),k(1))^2*g(1,2)-(2*ε-4)^-1*dot(p(2),p(3))*dot(k(1),k(1))\
            )*topo(I1L(muvsq,1))"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_reduction_2l_a() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });

    let integral = vakint
        .to_canonical(
            vakint_parse!(
                "(\
                    (k(1,1)*k(2,2))^2*g(1,2)+k(2,3)*p(1,3)+k(1,1)*k(2,2)*p(2,1)*p(3,2)\
                )*topo(I2L(mUVsq,1,2,1))"
            )
            .unwrap()
            .as_view(),
            true,
        )
        .unwrap();

    _ = compare_output(
        vakint
            .tensor_reduce(integral.as_view())
            .as_ref()
            .map(|a| a.as_view()),
        vakint_parse!(
            "(\
                -(2*ε-4)^-1*dot(p(2),p(3))*dot(k(1),k(2))+dot(k(1),k(1))*dot(k(2),k(2))*g(1,2)\
            )*topo(I2L(mUVsq,1,2,1))"
        )
        .unwrap(),
    );
}

#[allow(dead_code)]
fn run_tensor_reduction_tests() {
    test_reduction_1l_a();
    test_reduction_1l_b();
    test_reduction_2l_a();
}

#[test_log::test]
fn dot_conversion_preserves_factorized_scalar_powers() {
    use symbolica::atom::AtomCore;
    use vakint::Vakint;

    let _ = vakint::symbols::S.dot;
    for expression in [
        "(dot(p(5),k(1))+dot(k(1),k(2)))^2",
        "(dot(p(5),k(1))+dot(k(1),k(2)))^3",
        "(dot(p(5),k(1))+dot(k(1),k(2)))^4",
        "(dot(p(5),k(1))+dot(k(1),k(2)))^2*dot(p(5),k(1))",
        "((dot(p(5),k(1))+dot(k(1),k(2)))^2+dot(k(1),k(1)))^3",
        "(dot(p(5),k(1))+dot(k(1),k(2)))*(dot(p(5),k(1))+dot(k(1),k(1)))",
    ] {
        let numerator = vakint_parse!(expression).unwrap();
        let indexed = Vakint::convert_from_dot_notation(numerator.as_view());
        let round_trip = Vakint::convert_to_dot_notation(indexed.as_view());
        // Every copy of a scalar contraction, including both copies in a
        // square, needs independent summed indices. Expansion is confined to
        // this diagnostic comparison.
        assert!(
            (&round_trip - &numerator).expand().is_zero(),
            "conversion changed the scalar contractions of {expression}"
        );
    }
}

#[test_log::test]
fn tensor_reduction_of_powered_scalar_sum_matches_angular_average() {
    use symbolica::atom::AtomCore;

    let vakint = get_vakint(VakintSettings {
        use_dot_product_notation: true,
        epsilon_symbol: "eps".into(),
        ..VakintSettings::default()
    });
    let input = vakint_parse!("(dot(k(1),p(1))+dot(k(1),p(2)))^2*topo(I1L(muvsq,1))").unwrap();
    // Vacuum isotropy gives <k_mu k_nu> = g_mu_nu k^2/D.
    // In particular, both diagonal terms need the same 1/D as the cross term.
    let expected = vakint_parse!(
        "dot(k(1),k(1))*(dot(p(1),p(1))+2*dot(p(1),p(2))+dot(p(2),p(2)))/(4-2*eps)*topo(I1L(muvsq,1))"
    )
    .unwrap();
    let input = vakint.to_canonical(input.as_view(), false).unwrap();
    let expected = vakint.to_canonical(expected.as_view(), false).unwrap();
    let actual = vakint.tensor_reduce(input.as_view()).unwrap();
    assert!((actual - expected).together().cancel().is_zero());
}

#[test_log::test]
fn tensor_reduction_preserves_factorized_scalar_cancellation() {
    use symbolica::{
        atom::{Atom, AtomCore},
        id::Replacement,
    };

    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });
    let a = vakint_parse!("dot(p(5),k(1))").unwrap();
    let b = vakint_parse!("dot(k(1),k(2))").unwrap();
    let c = vakint_parse!("dot(k(1),k(1))").unwrap();
    let numerator: symbolica::atom::Atom = 8 * (&a + &b) * &a * c.clone().pow(2)
        - 4 * (&a + &b).pow(2) * c.clone().pow(2)
        - 4 * a.pow(2) * c.clone().pow(2);
    let scalar: symbolica::atom::Atom = -4 * b.pow(2) * c.pow(2);
    // Both displayed expressions have degree at most two in each of a, b, c.
    // Their exact values on three distinct points per variable therefore
    // certify the complete identity, with no polynomial conversion or expansion.
    for a_value in [-1, 0, 1] {
        for b_value in [-1, 0, 1] {
            for c_value in [-1, 0, 1] {
                let replacements = [
                    Replacement::new(a.clone(), Atom::num(a_value)),
                    Replacement::new(b.clone(), Atom::num(b_value)),
                    Replacement::new(c.clone(), Atom::num(c_value)),
                ];
                assert_eq!(
                    numerator.replace_multiple(&replacements),
                    scalar.replace_multiple(&replacements),
                );
            }
        }
    }
    let topology =
        vakint_parse!("topo(prop(1,edge(1,1),k(1),muvsq,5)*prop(2,edge(1,1),k(2),muvsq,3))")
            .unwrap();
    let expected = vakint
        .to_canonical((&scalar * &topology).as_view(), false)
        .unwrap();

    // A numerator containing only internal scalar products is already tensor
    // reduced. The complete GL06 source is identically this monomial, with no
    // external p(5) dependence, before any angular or radial integration.
    for (label, input) in [
        ("scalar control", scalar),
        ("complete factorized source", numerator),
    ] {
        let canonical = vakint
            .to_canonical((input * &topology).as_view(), false)
            .unwrap();
        let actual = vakint.tensor_reduce(canonical.as_view()).unwrap();
        assert!(
            (&actual - &expected).together().is_zero(),
            "{label} tensor reduction differs from the exact scalar source: {}",
            (&actual - &expected).together()
        );
    }
}

#[test_log::test]
fn dot_conversion_preserves_existing_indices() {
    use symbolica::atom::AtomCore;
    use vakint::Vakint;

    let _ = vakint::symbols::S.dot;
    let indexed = vakint_parse!("(p(5,dot_dummy_ind(9))+k(1,dot_dummy_ind(9)))^2").unwrap();
    let scalar = vakint_parse!("(dot(p(5),k(1))+dot(k(1),k(2)))^2").unwrap();
    for (label, input) in [
        ("already indexed", indexed.clone()),
        ("mixed old and new indices", indexed * scalar),
    ] {
        let converted = Vakint::convert_from_dot_notation(input.as_view());
        let twice = Vakint::convert_from_dot_notation(converted.as_view());
        let expected = Vakint::convert_to_dot_notation(input.as_view());
        let actual = Vakint::convert_to_dot_notation(converted.as_view());
        assert_eq!(converted, twice, "dot conversion must be idempotent");
        if label == "already indexed" {
            assert_eq!(
                converted, input,
                "no dot notation must require no conversion"
            );
        }
        assert!(
            (actual - expected).expand().is_zero(),
            "{label} changed a pre-existing contraction"
        );
    }
}

#[test_log::test]
fn dot_conversion_preserves_reciprocal_powers() {
    use symbolica::atom::AtomCore;
    use vakint::Vakint;

    let _ = vakint::symbols::S.dot;
    for expression in [
        "dot(p(5),k(1))^-2",
        "(dot(p(5),k(1))+dot(k(1),k(2)))^-1",
        "(dot(p(5),k(1))+dot(k(1),k(2)))^-2",
        "(dot(p(5),k(1))+dot(k(1),k(2)))^-3",
    ] {
        let input = vakint_parse!(expression).unwrap();
        let converted = Vakint::convert_from_dot_notation(input.as_view());
        assert_eq!(
            converted,
            Vakint::convert_from_dot_notation(converted.as_view()),
            "dot conversion must be idempotent"
        );
        // Inverting the scalar coefficient exposes every repeated contraction
        // to the ordinary polynomial round-trip check, without integrating an
        // input with a non-polynomial energy numerator.
        let reciprocal = Vakint::convert_to_dot_notation(converted.pow(-1).as_view());
        let expected = input.pow(-1);
        assert!(
            (&reciprocal - &expected).together().is_zero(),
            "{expression} changed its reciprocal scalar contractions: {reciprocal}"
        );
    }
}

#[test_log::test]
fn loop_normalization_uses_numeric_imaginary_coefficients() {
    use symbolica::atom::{Atom, AtomCore};
    use vakint::LoopNormalizationFactor;

    let settings = VakintSettings::default();
    // Canonical user symbols can be imported from a persisted Symbolica state.
    // Their registration must not change a built-in normalization's algebra.
    let symbolic_i = vakint_parse!("symbolica::{}::𝑖").unwrap();
    assert_ne!(symbolic_i, Atom::i());
    assert_eq!(vakint_parse!("1𝑖").unwrap(), Atom::i());
    for normalization in [
        LoopNormalizationFactor::pySecDec,
        LoopNormalizationFactor::FMFTandMATAD,
    ] {
        let expression = normalization.to_atom(&settings).unwrap();
        assert!(
            !expression
                .get_all_symbols(true)
                .contains(&symbolic_i.get_symbol().unwrap())
        );
        for loops in [1, 2, 3] {
            let actual = expression
                .replace(vakint_parse!("n_loops").unwrap())
                .with(Atom::num(loops))
                .replace(vakint_parse!(&settings.epsilon_symbol).unwrap())
                .with(Atom::Zero);
            let expected = vakint_parse!("1𝑖*𝜋^2").unwrap().pow(-loops);
            assert!((actual - expected).together().is_zero());
        }
    }
}

#[test_log::test]
fn tensor_reduction_preserves_symbolica_user_namespaces() {
    use symbolica::atom::AtomCore;

    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    });
    // A same-spelled user variable/function must retain its full identity when
    // escaped for FORM; only Vakint's namespace is stripped before escaping.
    for coefficient in [
        "symbolica::{}::x",
        "symbolica::{}::𝑖",
        "symbolica::{}::f(vakint::{}::x)",
    ] {
        let coefficient_atom = vakint_parse!(coefficient).unwrap();
        let (header, sanitized, indices) = vakint
            .vakint
            .sanitize_user_expressions(&vakint.settings, coefficient_atom.as_view(), false, &[])
            .unwrap();
        let round_trip = vakint
            .vakint
            .process_form_output(
                &vakint.settings,
                sanitized.clone(),
                indices,
                Default::default(),
            )
            .unwrap();
        println!(
            "FORM_NAMESPACE_ROUNDTRIP coefficient={} header={header:?} sanitized={sanitized:?} actual={} expected={}",
            coefficient,
            round_trip.to_canonical_string(),
            coefficient_atom.to_canonical_string(),
        );
        let input = vakint_parse!(format!(
            "({coefficient})*k(1,1)*k(1,2)*topo(prop(1,edge(1,1),k(1),muvsq,1))"
        ))
        .unwrap();
        let canonical = vakint.to_canonical(input.as_view(), true).unwrap();
        let reduced = vakint.tensor_reduce(canonical.as_view()).unwrap();
        let expected = vakint_parse!(format!(
            "-({coefficient})*(2*ε-4)^-1*dot(k(1),k(1))*g(1,2)*topo(I1L(muvsq,1))"
        ))
        .unwrap();
        assert!(
            (&round_trip - &coefficient_atom).together().is_zero(),
            "no-FORM namespace round trip changed {coefficient}: actual={}, expected={}; tensor actual={}, expected={}",
            round_trip.to_canonical_string(),
            coefficient_atom.to_canonical_string(),
            reduced.to_canonical_string(),
            expected.to_canonical_string(),
        );
        assert!(
            (&reduced - &expected).together().is_zero(),
            "namespace round trip changed {coefficient}: actual={}, expected={}",
            reduced.to_canonical_string(),
            expected.to_canonical_string(),
        );
    }
}
