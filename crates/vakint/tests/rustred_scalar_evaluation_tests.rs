use std::collections::HashMap;
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    EvaluationMethod, EvaluationOrder, MATADOptions, RustRedEvaluationOptions, Vakint,
    VakintSettings, vakint_parse,
};

fn rustred_settings(substitute_masters: bool) -> VakintSettings {
    VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by/rustred-scalar".to_owned(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { substitute_masters },
        )]),
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn matad_raw_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder::matad_only(Some(MATADOptions {
            expand_masters: false,
            ..MATADOptions::default()
        })),
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn matad_substituted_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder::matad_only(None),
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

fn assert_exactly_equal(left: Atom, right: Atom) {
    let difference = (left - right).together().expand();
    assert!(
        difference.is_zero(),
        "nonzero exact difference: {difference}"
    );
}

#[test]
fn scalar_rustred_reduces_without_a_form_executable() {
    let vakint = Vakint::new().unwrap();
    let result = vakint
        .evaluate_integral(
            &rustred_settings(false),
            vakint_parse!("topo(I1L(muvsq,4))").unwrap().as_view(),
        )
        .unwrap();

    assert!(result.to_canonical_string().contains("Gam(1,1)"));

    let settings = rustred_settings(false);
    let scalar_identities = [
        (
            vakint_parse!("dot(k(1),k(1))*topo(I1L(muvsq,2))").unwrap(),
            vakint_parse!("topo(I1L(muvsq,1))+muvsq*topo(I1L(muvsq,2))").unwrap(),
        ),
        (
            vakint_parse!("dot(k(1),k(2))*topo(I2L(muvsq,1,1,1))").unwrap(),
            vakint_parse!(
                "1/2*(\
                    topo(I2L(muvsq,1,1,0))\
                   -topo(I2L(muvsq,0,1,1))\
                   -topo(I2L(muvsq,1,0,1))\
                   -muvsq*topo(I2L(muvsq,1,1,1))\
                )"
            )
            .unwrap(),
        ),
    ];
    for (scalar_numerator, affine_expansion) in scalar_identities {
        let lowered = vakint
            .evaluate_integral(&settings, scalar_numerator.as_view())
            .unwrap();
        let expanded = vakint
            .evaluate_integral(&settings, affine_expansion.as_view())
            .unwrap();
        assert_exactly_equal(lowered, expanded);
    }
}

#[test]
fn scalar_rustred_matches_raw_matad_through_two_loops() {
    let vakint = Vakint::new().unwrap();
    let rustred = rustred_settings(false);
    let matad = matad_raw_settings();
    let mut inputs = (1..=6)
        .map(|power| vakint_parse!(format!("topo(I1L(muvsq,{power}))").as_str()).unwrap())
        .collect::<Vec<_>>();
    inputs.extend([
        vakint_parse!("topo(I2L(muvsq,1,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,1,2,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,1,1,2))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,2,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,3,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,-1,2,2))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,0,1,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,0,2,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,0,1))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,1,0))").unwrap(),
        vakint_parse!("topo(I2L_pinch_3(muvsq,1,1,0))").unwrap(),
        vakint_parse!("topo(I2L_pinch_3(muvsq,2,1,0))").unwrap(),
    ]);

    for input in inputs {
        let rustred_result = vakint.evaluate_integral(&rustred, input.as_view()).unwrap();
        let matad_result = vakint.evaluate_integral(&matad, input.as_view()).unwrap();
        assert_exactly_equal(rustred_result, matad_result);
    }
}

#[test]
fn scalar_rustred_reuses_matcher_routing_and_preserves_spectators() {
    let vakint = Vakint::new().unwrap();
    let settings = rustred_settings(false);
    let equivalent_inputs = [
        (
            vakint_parse!(
                "7*user_space::c*topo(\
                    prop(77,edge(42,42),k(8),user_space::mass_squared,4)\
                )"
            )
            .unwrap(),
            vakint_parse!("7*user_space::c*topo(I1L(user_space::mass_squared,4))").unwrap(),
        ),
        (
            vakint_parse!(
                "user_space::c*topo(\
                    prop(55,edge(7,10),k(22),user_space::mass_squared,1)*\
                    prop(9,edge(7,10),k(11),user_space::mass_squared,2)*\
                    prop(33,edge(10,7),k(11)+k(22),user_space::mass_squared,1)\
                )"
            )
            .unwrap(),
            vakint_parse!("user_space::c*topo(I2L(user_space::mass_squared,2,1,1))").unwrap(),
        ),
        (
            vakint_parse!(
                "user_space::c*topo(\
                    prop(33,edge(10,10),k(22),user_space::mass_squared,2)*\
                    prop(55,edge(10,10),k(11),user_space::mass_squared,1)\
                )"
            )
            .unwrap(),
            vakint_parse!("user_space::c*topo(I2L_pinch_3(user_space::mass_squared,2,1,0))")
                .unwrap(),
        ),
    ];

    for (graph_input, short_input) in equivalent_inputs {
        let graph_result = vakint
            .evaluate_integral(&settings, graph_input.as_view())
            .unwrap();
        let short_result = vakint
            .evaluate_integral(&settings, short_input.as_view())
            .unwrap();
        assert_exactly_equal(graph_result, short_result);
    }
}

#[test]
fn scalar_rustred_accepts_literal_unit_mass_without_a_mass_symbol() {
    let vakint = Vakint::new().unwrap();
    let settings = rustred_settings(false);
    let cases = [
        (
            vakint_parse!("topo(I1L(1,4))").unwrap(),
            vakint_parse!("topo(I1L(user_space::mass_squared,4))").unwrap(),
        ),
        (
            vakint_parse!("topo(I2L(1,2,1,1))").unwrap(),
            vakint_parse!("topo(I2L(user_space::mass_squared,2,1,1))").unwrap(),
        ),
        (
            vakint_parse!("topo(I2L_pinch_3(1,2,1,0))").unwrap(),
            vakint_parse!("topo(I2L_pinch_3(user_space::mass_squared,2,1,0))").unwrap(),
        ),
    ];

    for (unit_mass, symbolic_mass) in cases {
        let unit_result = vakint
            .evaluate_integral(&settings, unit_mass.as_view())
            .unwrap();
        let symbolic_result = vakint
            .evaluate_integral(&settings, symbolic_mass.as_view())
            .unwrap()
            .replace(
                vakint_parse!("user_space::mass_squared")
                    .unwrap()
                    .to_pattern(),
            )
            .with(Atom::num(1).to_pattern());
        assert_exactly_equal(unit_result, symbolic_result);
    }

    let caller_m_result = vakint
        .evaluate_integral(
            &settings,
            vakint_parse!("topo(I1L(M,4))").unwrap().as_view(),
        )
        .unwrap();
    let reference = vakint
        .evaluate_integral(
            &settings,
            vakint_parse!("topo(I1L(user_space::mass_squared,4))")
                .unwrap()
                .as_view(),
        )
        .unwrap()
        .replace(
            vakint_parse!("user_space::mass_squared")
                .unwrap()
                .to_pattern(),
        )
        .with(vakint_parse!("M").unwrap().to_pattern());
    assert_exactly_equal(caller_m_result, reference);
}

#[test]
fn scalar_rustred_matches_substituted_matad_through_two_loops() {
    let vakint = Vakint::new().unwrap();
    let rustred = rustred_settings(true);
    let matad = matad_substituted_settings();
    let parameter_values = HashMap::from([("muvsq".to_owned(), 1.0), ("mursq".to_owned(), 1.0)]);
    let rustred_parameters = vakint.params_from_f64(&rustred, &parameter_values);
    let matad_parameters = vakint.params_from_f64(&matad, &parameter_values);
    let inputs = [
        vakint_parse!("topo(I1L(muvsq,4))").unwrap(),
        vakint_parse!("topo(I2L(muvsq,2,2,1))").unwrap(),
        vakint_parse!("topo(I2L_pinch_3(muvsq,2,1,0))").unwrap(),
    ];

    for input in inputs {
        let rustred_result = vakint.evaluate_integral(&rustred, input.as_view()).unwrap();
        let matad_result = vakint.evaluate_integral(&matad, input.as_view()).unwrap();
        let (rustred_numerical, _) = Vakint::full_numerical_evaluation(
            &rustred,
            rustred_result.as_view(),
            &rustred_parameters,
            &HashMap::default(),
            None,
        )
        .unwrap();
        let (matad_numerical, _) = Vakint::full_numerical_evaluation(
            &matad,
            matad_result.as_view(),
            &matad_parameters,
            &HashMap::default(),
            None,
        )
        .unwrap();
        let (matches, detail) =
            rustred_numerical.does_approx_match(&matad_numerical, None, 1.0e-25, 0.0);
        assert!(matches, "{detail}");
    }
}

#[test]
fn form_tensor_prepass_then_rustred_scalar_matches_matad() {
    let vakint = Vakint::new().unwrap();
    let mut rustred = rustred_settings(false);
    rustred.form_exe_path = "form".to_owned();
    let matad = matad_raw_settings();
    let inputs = [
        vakint_parse!(
            "(k(1,1)*k(1,2)+k(1,3)*p(1,3))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(k(1,11)*p(1,11)*k(1,12)*p(1,12))*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(user_space::A*k(1,11)*p(1,11)*k(1,12)*p(1,12)+user_space::B)*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(k(1,1)*p(1,1)*k(1,2)*p(2,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "(\
                k(1,11)*k(2,22)*k(1,11)*k(2,22)\
                +p(1,11)*k(1,11)*k(1,22)*p(1,22)\
                +p(1,11)*p(2,11)*k(2,22)*k(2,22)\
            )*topo(\
                prop(1,edge(1,2),k(1),muvsq,1)*\
                prop(2,edge(1,2),k(2),muvsq,1)*\
                prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )"
        )
        .unwrap(),
        vakint_parse!(
            "k(11,mu)*k(22,nu)*topo(\
                prop(33,edge(10,10),k(22),muvsq,1)*\
                prop(55,edge(10,10),k(11),muvsq,1)\
            )"
        )
        .unwrap(),
    ];

    for input in inputs {
        let rustred_result = vakint.evaluate(&rustred, input.as_view()).unwrap();
        let matad_result = vakint.evaluate(&matad, input.as_view()).unwrap();
        assert_exactly_equal(rustred_result, matad_result);
    }
}
