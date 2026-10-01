mod test_utils;
use log::debug;
use symbolica::try_parse;
use test_utils::{compare_output, get_vakint};
use vakint::{VakintError, VakintSettings};

#[test_log::test]
fn test_1l_matching() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    debug!("Topologies:\n{}", vakint.topologies);
    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(lm(8,2)*lm(8,2)+lm(8,33)*pext(12,33))*vakint::topo(\
                    vakint::prop(77,vakint::edge(42,42),lm(8),muvsq,1)\
                )"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(1,2)^2+vakint::k(1,33)*pext(12,33))*vakint::topo(\
                vakint::prop(1,vakint::edge(1,1),vakint::k(1),muvsq,1)\
            )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(lm(8,2)*lm(8,2)+lm(8,33)*pext(12,33))*vakint::topo(\
                    vakint::prop(77,vakint::edge(42,42),lm(8),muvsq,1)\
                )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(1,2)^2+vakint::k(1,33)*pext(12,33))*vakint::topo(vakint::I1L(muvsq,1))"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!("(lm(1,2)^2+lm(1,33)*pext(12,33))*vakint::topo(vakint::I1L(muvsq,-3))")
                    .unwrap()
                    .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(lm(1,2)^2+lm(1,33)*pext(12,33))*vakint::topo(vakint::prop(1,vakint::edge(1,1),vakint::k(1),muvsq,-3))")
            .unwrap(),
    );
}

#[test_log::test]
fn test_2l_matching_3prop() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    //println!("Topologies:\n{}", vakint.topologies);

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*vakint::topo(\
                            vakint::prop(9,vakint::edge(7,10),k(11),mUVsq,1)*\
                            vakint::prop(33,vakint::edge(7,10),k(22),mUVsq,2)*\
                            vakint::prop(55,vakint::edge(7,10),k(11)+k(22),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(1,2)*vakint::k(1,2)+vakint::k(1,77)*vakint::k(2,77)+vakint::k(2,33)*p(42,33))*vakint::topo(\
                        vakint::prop(1,vakint::edge(1,2),vakint::k(1),mUVsq,1)*\
                        vakint::prop(2,vakint::edge(1,2),vakint::k(2),mUVsq,2)*\
                        vakint::prop(3,vakint::edge(2,1),vakint::k(1)+vakint::k(2),mUVsq,1)\
                    )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(vakint::k(11,2)*vakint::k(11,2)+vakint::k(11,77)*vakint::k(22,77)+vakint::k(22,33)*p(42,33))*vakint::topo(\
                            vakint::prop(9,vakint::edge(7,10),vakint::k(11),mUVsq,1)*\
                            vakint::prop(33,vakint::edge(7,10),vakint::k(22),mUVsq,2)*\
                            vakint::prop(55,vakint::edge(7,10),vakint::k(11)+vakint::k(22),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(vakint::k(1,2)*vakint::k(1,2)+vakint::k(1,77)*vakint::k(2,77)+vakint::k(2,33)*p(42,33))*vakint::topo(vakint::I2L(mUVsq,1,2,1))")
            .unwrap(),
    );

    _ = compare_output(
        vakint.to_canonical(
            try_parse!("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*vakint::topo(vakint::I2L(mUVsq,1,2,1))")
                .unwrap()
                .as_view(),
            false,
        ).as_ref().map(|a| a.as_view()),
        try_parse!("(k(1,2)*k(1,2)+k(1,77)*k(2,77)+k(2,33)*p(42,33))*vakint::topo(\
                                                vakint::prop(1,vakint::edge(1,2),vakint::k(1),mUVsq,1)*vakint::prop(2,vakint::edge(1,2),vakint::k(2),mUVsq,2)*\
                                                vakint::prop(3,vakint::edge(2,1),vakint::k(1)+vakint::k(2),mUVsq,1)\
                                            )").unwrap(),
    );
}

#[test_log::test]
fn partly_massless_sunsets_preserve_mass_labels_and_exclude_alphaloop() {
    use symbolica::atom::AtomCore;
    use vakint::{EvaluationOrder, vakint_parse};

    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        evaluation_order: EvaluationOrder::alphaloop_only(),
        ..VakintSettings::default()
    });
    for (family, second_mass, massless_count) in [("I2L_MM0", "muvsq", 1), ("I2L_M00", "0", 2)] {
        // Put the massless line first and use arbitrary labels. Canonicalizing
        // must permute its momentum and incidence with its physical mass.
        let input = vakint_parse!(format!(
            "(a+b)*(c+d)*topo(\
             prop(9,edge(7,10),k(11),0,1)*\
             prop(33,edge(7,10),k(22),muvsq,2)*\
             prop(55,edge(10,7),k(11)+k(22),{second_mass},1))"
        ))
        .unwrap();
        let canonical = vakint.to_canonical(input.as_view(), false).unwrap();
        assert_eq!(
            canonical
                .pattern_match(
                    &vakint_parse!("prop(id_,edge(a_,b_),q_,0,power_)")
                        .unwrap()
                        .to_pattern(),
                    None,
                    None,
                )
                .count(),
            massless_count
        );
        for factor in ["a+b", "c+d"] {
            assert!(
                canonical
                    .pattern_match(&vakint_parse!(factor).unwrap().to_pattern(), None, None,)
                    .next()
                    .is_some()
            );
        }
        let short = vakint.to_canonical(input.as_view(), true).unwrap();
        assert!(
            short
                .pattern_match(
                    &vakint_parse!(format!("topo({family}(muvsq,powers__))"))
                        .unwrap()
                        .to_pattern(),
                    None,
                    None,
                )
                .next()
                .is_some()
        );
        assert!(matches!(
            vakint.evaluate_integral(short.as_view()),
            Err(VakintError::NoEvaluationMethodFound(_, _))
        ));
    }
}

#[test_log::test]
fn test_2l_matching_pinched() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*vakint::topo(\
                            vakint::prop(33,vakint::edge(10,10),k(22),mUVsq,2)*\
                            vakint::prop(55,vakint::edge(10,10),k(11),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(2,2)^2+vakint::k(1,33)*p(42,33)+vakint::k(1,77)*vakint::k(2,77))*vakint::topo(\
                        vakint::prop(1,vakint::edge(1,1),vakint::k(1),mUVsq,2)*\
                        vakint::prop(2,vakint::edge(1,1),vakint::k(2),mUVsq,1)
                    )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(k(11,2)*k(11,2)+k(11,77)*k(22,77)+k(22,33)*p(42,33))*vakint::topo(\
                            vakint::prop(33,vakint::edge(10,10),k(22),mUVsq,2)*\
                            vakint::prop(55,vakint::edge(10,10),k(11),mUVsq,1)\
                        )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("(vakint::k(2,2)^2+vakint::k(1,33)*p(42,33)+vakint::k(1,77)*vakint::k(2,77))*vakint::topo(vakint::I2L_pinch_3(mUVsq,2,1,0))")
            .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(vakint::k(2,2)^2+vakint::k(1,33)*p(42,33)+vakint::k(1,77)*vakint::k(2,77))*vakint::topo(vakint::I2L(mUVsq,2,1,0))"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(2,2)^2+vakint::k(1,33)*p(42,33)+vakint::k(1,77)*vakint::k(2,77))*vakint::topo(\
                        vakint::prop(1,vakint::edge(1,2),vakint::k(1),mUVsq,2)*\
                        vakint::prop(2,vakint::edge(1,2),vakint::k(2),mUVsq,1)*\
                        vakint::prop(3,vakint::edge(2,1),vakint::k(1)+vakint::k(2),mUVsq,0)
                )"
        )
        .unwrap(),
    );

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "(vakint::k(2,2)^2+vakint::k(1,33)*p(42,33)+vakint::k(1,77)*vakint::k(2,77))*vakint::topo(vakint::I2L_pinch_3(mUVsq,2,1,0))"
                )
                .unwrap()
                .as_view(),
                false,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(2,2)^2+vakint::k(1,33)*p(42,33)+vakint::k(1,77)*vakint::k(2,77))*vakint::topo(\
                        vakint::prop(1,vakint::edge(1,1),vakint::k(1),mUVsq,2)*\
                        vakint::prop(2,vakint::edge(1,1),vakint::k(2),mUVsq,1)
                )"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_3l_matching_with_zero_powers_in_short_form() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    debug!("Topologies:\n{}", vakint.topologies);

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!(
                    "( 1 )*vakint::topo(\
                vakint::prop(1,vakint::edge(1,2),k(1),muvsq,1)\
              * vakint::prop(2,vakint::edge(1,2),k(2),muvsq,1)\
              * vakint::prop(3,vakint::edge(1,2),k(3),muvsq,1)\
              * vakint::prop(4,vakint::edge(2,1),k(1)+k(2)+k(3),muvsq,2)\
            )"
                )
                .unwrap()
                .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("vakint::topo(vakint::I3L_pinch_1_6(muvsq,0,1,1,1,2,0))").unwrap(),
    );
}

#[test_log::test]
fn test_unknown_integrals() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });
    let vakint_with_unknown_integrals = get_vakint(VakintSettings {
        allow_unknown_integrals: true,
        ..VakintSettings::default()
    });

    let unknown_integral = try_parse!(
        "(vakint::k(1,2)*vakint::k(1,2)+vakint::k(1,77)*vakint::k(2,77)+vakint::k(2,33)*p(42,33))*vakint::topo(\
                    vakint::prop(1,vakint::edge(7,10),vakint::k(2),mA,2)*\
                    vakint::prop(2,vakint::edge(7,10),vakint::k(1),mB,1)\
                )"
    )
    .unwrap();

    assert!(matches!(
        vakint.to_canonical(unknown_integral.as_view(), false),
        Err(VakintError::UnreckognizedIntegral(_))
    ));

    // println!("Topologies:\n{}", vakint_with_unknown_integrals.topologies);

    _ = compare_output(
        vakint_with_unknown_integrals
            .to_canonical(unknown_integral.as_view(), false)
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!(
            "(vakint::k(1,2)^2+vakint::k(1,77)*vakint::k(2,77)+vakint::k(2,33)*p(42,33))*vakint::topo(vakint::UNKNOWN(\
                        vakint::prop(1,vakint::edge(7,10),vakint::k(2),mA,2)*\
                        vakint::prop(2,vakint::edge(7,10),vakint::k(1),mB,1))\
                    )"
        )
        .unwrap(),
    );
}

#[test_log::test]
fn test_2l_pinched_matching() {
    let vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    });

    // println!("Topologies:\n{}", vakint_with_unknown_integrals.topologies);

    _ = compare_output(
        vakint
            .to_canonical(
                try_parse!("vakint::topo(vakint::prop(1,vakint::edge(1,1),k(1),muvsq,1)*vakint::prop(2,vakint::edge(1,1),k(2),muvsq,1))")
                    .unwrap()
                    .as_view(),
                true,
            )
            .as_ref()
            .map(|a| a.as_view()),
        try_parse!("vakint::topo(vakint::I2L_pinch_3(muvsq,1,1,0))").unwrap(),
    );
}

#[allow(dead_code)]
fn run_input_matching_tests() {
    test_1l_matching();
    test_2l_matching_3prop();
    test_2l_matching_pinched();
    test_3l_matching_with_zero_powers_in_short_form();
    test_unknown_integrals();
    test_2l_pinched_matching();
}
