mod test_utils;
use vakint::{
    AlphaLoopOptions, EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor,
    MATADOptions,
};

use std::vec;

use log::debug;
use std::collections::HashMap;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
};
use test_utils::{compare_numerical_output, compare_output, get_vakint};
use vakint::{Vakint, VakintSettings};

use crate::test_utils::compare_vakint_evaluation_vs_reference;
use vakint::{externals_from_f64, params_from_f64, vakint_parse};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS: u32 = 32;

#[test_log::test]
fn test_integrate_1l_a() {
    let mut vakint = get_vakint(VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: true,
        mu_r_sq_symbol: format!("{}::mursq", vakint::NAMESPACE),
        integral_normalization_factor: vakint::LoopNormalizationFactor::MSbar,
        run_time_decimal_precision: 16,
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::AlphaLoop(AlphaLoopOptions {
            susbstitute_masters: false,
        })]),
        ..VakintSettings::default()
    });

    let mut integral = vakint
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

    integral = vakint.tensor_reduce(integral.as_view()).unwrap();

    let evaluated_integral_res = vakint.evaluate_integral(integral.as_view());
    let evaluated_integral = evaluated_integral_res.unwrap();

    //let evaluated_integral = evaluated_integral_res.unwrap().expand();

    let mut targets: HashMap<Atom, Atom> = HashMap::default();

    for (eps_term, trgt) in [
        ("ε^-1", "1/64*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2)/𝜋^2"),
        (
            "1",
            "-1/64*vakint::𝑖*vakint::muvsq^2*log(vakint::muvsq)*vakint::g(1,2)/𝜋^2+1/64*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)*vakint::g(1,2)/𝜋^2+3/128*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2)/𝜋^2",
        ),
        (
            "ε",
            "-(1/64*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)*vakint::g(1,2)/𝜋^2+3/128*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2)/𝜋^2)*log(vakint::muvsq)+1/16*vakint::𝑖*(7/16*vakint::muvsq^2*vakint::g(1,2)+1/48*𝜋^2*vakint::muvsq^2*vakint::g(1,2))/𝜋^2+3/128*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)*vakint::g(1,2)/𝜋^2+1/128*vakint::𝑖*vakint::muvsq^2*log(vakint::muvsq)^2*vakint::g(1,2)/𝜋^2+1/128*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)^2*vakint::g(1,2)/𝜋^2",
        ),
        (
            "ε^2",
            "1/2*(1/64*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)*vakint::g(1,2)/𝜋^2+3/128*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2)/𝜋^2)*log(vakint::muvsq)^2-(1/16*vakint::𝑖*(7/16*vakint::muvsq^2*vakint::g(1,2)+1/48*𝜋^2*vakint::muvsq^2*vakint::g(1,2))/𝜋^2+3/128*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)*vakint::g(1,2)/𝜋^2+1/128*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)^2*vakint::g(1,2)/𝜋^2)*log(vakint::muvsq)+1/16*vakint::𝑖*(7/16*vakint::muvsq^2*vakint::g(1,2)+1/48*𝜋^2*vakint::muvsq^2*vakint::g(1,2))*log(vakint::mursq)/𝜋^2+1/16*vakint::𝑖*(-1/12*vakint::muvsq^2*vakint::z3*vakint::g(1,2)+15/32*vakint::muvsq^2*vakint::g(1,2)+1/32*𝜋^2*vakint::muvsq^2*vakint::g(1,2))/𝜋^2-1/384*vakint::𝑖*vakint::muvsq^2*log(vakint::muvsq)^3*vakint::g(1,2)/𝜋^2+3/256*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)^2*vakint::g(1,2)/𝜋^2+1/384*vakint::𝑖*vakint::muvsq^2*log(vakint::mursq)^3*vakint::g(1,2)/𝜋^2",
        ),
    ] {
        targets.insert(
            vakint_parse!(eps_term).unwrap(),
            vakint_parse!(trgt).unwrap(),
        );
    }
    for (v, c) in evaluated_integral
        .coefficient_list::<i8>(&[vakint_parse!("ε").unwrap()])
        .iter()
    {
        println!("{} -> {}", v, c.to_canonical_string());
        _ = compare_output(
            Ok(c.as_view()),
            targets.get(&v.to_owned()).unwrap().to_owned(),
        );
    }

    // Now let's re-evaluate but now substituting the masters with their numerical expressions so as to proceed with numerical comparisons.
    vakint.settings.evaluation_order =
        EvaluationOrder(vec![EvaluationMethod::AlphaLoop(AlphaLoopOptions {
            susbstitute_masters: true,
        })]);
    let evaluated_integral_res = vakint.evaluate_integral(integral.as_view());
    let evaluated_integral = evaluated_integral_res.unwrap();

    //debug!("Evaluated integral: {}", evaluated_integral);
    let mut params = HashMap::default();
    params.insert("muvsq".into(), vakint.settings.real_to_prec("1"));
    params.insert("mursq".into(), vakint.settings.real_to_prec("1"));

    //evaluated_integral = evaluated_integral.expand();

    let numerical_partial_eval = Vakint::partial_numerical_evaluation(
        &vakint.settings,
        evaluated_integral.as_view(),
        &params,
        &HashMap::default(),
        None,
    );
    // numerical_partial_eval = numerical_partial_eval.expand();

    // This test is too unstable as the printout at fixed precision is not accurate enough
    // let numerical_partial_eval_canonical_str = numerical_partial_eval.to_canonical_string();
    // assert_eq!(
    //     numerical_partial_eval_canonical_str,
    //     "-12.0696723514860*g(1,2)*ε^2*𝑖+-5.36845814123893*g(1,2)*𝑖+2.46740110027234*g(1,2)*ε^-1*𝑖+9.41170424471097*g(1,2)*ε*𝑖"
    // );
    let fractional_precision: Rational = (
        1,
        10_u64.pow((vakint.settings.run_time_decimal_precision - 4) as u32),
    )
        .into();
    // println!(
    //     "TARGET {}",
    //     numerical_partial_eval
    //         .rationalize_coefficients(&Fraction::from(
    //             0.1_f64.powi((vakint.settings.run_time_decimal_precision - 4) as i32)
    //         ))
    //         .to_canonical_string()
    // );
    assert_eq!(
        numerical_partial_eval.rationalize(&fractional_precision),
       vakint_parse!("62865𝑖/15436144*vakint::ε*vakint::g(1,2)+44525𝑖/10385624*vakint::ε^2*vakint::g(1,2)+34571𝑖/21836934*vakint::g(1,2)/vakint::ε+34571𝑖/14557956*vakint::g(1,2)").unwrap()
    );

    let prec: Rational = (
        1,
        10_u64.pow((vakint.settings.run_time_decimal_precision - 4) as u32),
    )
        .into();
    //println!("numerical_partial_eval=\n{}", numerical_partial_eval);
    _ = compare_output(
        Ok(numerical_partial_eval.rationalize(&prec).as_view()),
        vakint_parse!(
            format!(
                "2.374715241617292e-3`{prec}𝑖*g(1,2)\
             +4.072584448553507e-3`{prec}𝑖*ε*g(1,2)\
             +1.583143494411528e-3`{prec}𝑖*ε^-1*g(1,2)\
             +4.287176196638421e-3`{prec}𝑖*ε^2*g(1,2)",
                prec = vakint.settings.run_time_decimal_precision - 1
            )
            .as_str()
        )
        .unwrap()
        .rationalize(&prec),
    );

    params.insert("g(1,2)".into(), vakint.settings.real_to_prec("1"));
    let numerical_full_eval = Vakint::full_numerical_evaluation_without_error(
        &vakint.settings,
        evaluated_integral.as_view(),
        &params,
        &HashMap::default(),
        None,
    );
    let numerical_full_eval_ref = numerical_full_eval.as_ref();
    debug!(
        "Full eval (metric substituted with 1):\n{}",
        numerical_full_eval_ref.unwrap()
    );
    compare_numerical_output(
        numerical_full_eval_ref,
        vec![
            (-1, ("0.0".into(), "1.583143494411528e-3".into())),
            (0, ("0.0".into(), "2.374715241617292e-3".into())),
            (1, ("0.0".into(), "4.072584448553507e-3".into())),
            (2, ("0.0".into(), "4.287176196638421e-3".into())),
        ],
        vakint.settings.run_time_decimal_precision,
    );
}

#[test_log::test]
fn test_integrate_1l_simple() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 2, integral_normalization_factor: LoopNormalizationFactor::pySecDec,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "( 1 )*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        HashMap::default(),
        vec![
            (-1, ("1.0".into(), "0.0".into()),),
            (0,  ("4.227843350984671393934879099176e-1".into(),  "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_simple_squared_mass() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 2, integral_normalization_factor: LoopNormalizationFactor::pySecDec,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(  2*user_space::muv - user_space::muv^2 )*topo(\
                prop(1,edge(1,1),k(1),user_space::muv^2,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("user_space::muv".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        HashMap::default(),
        vec![
            (-1, ("1.0".into(), "0.0".into()),),
            (0,  ("4.227843350984671393934879099176e-1".into(),  "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_cross_product() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(k(1,11)*p(1,11)*k(1,12)*p(1,12))*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-9.892747067878755034798287150883e-3".into()),),
            (0,  ("0.0".into(),  "-9.892747067878755034798287150883e-3".into()),),
            (1,  ("0.0".into(), "-1.802920540121208908815751124632e-2".into()),),
            (2,  ("0.0".into(),  "-1.406532376649359143630292751833e-2".into()),),
            (3,  ("0.0".into(),  "-2.008809563968960604519605888938e-2".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_cross_product_with_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(user_space::A*k(1,11)*p(1,11)*k(1,12)*p(1,12)+user_space::B)*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0), ("user_space::A".into(), 3.0), ("user_space::B".into(), 4.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-4.347945293051822243424995650216e-3".into()),),
            (0,  ("0.0".into(), "-2.967824120363626510439486145265e-2".into()),),
            (1,  ("0.0".into(), "-3.325428287030293393113920040564e-2".into()),),
            (2,  ("0.0".into(), "-5.234545698561186456899262604538e-2".into()),),
            (3,  ("0.0".into(), "-4.484303004236669529365928448084e-2".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_dot_product_external() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(k(1,1)*p(1,1)*k(1,2)*p(2,2))*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-4.682938456469298873921803940225".into()),),
            (0,  ("0.0".into(),  "-7.024407684703948310882705910337".into()),),
            (1,  ("0.0".into(), "-12.04670479882127302936315689539".into()),),
            (2,  ("0.0".into(),  "-12.68146718965645007677038182264".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_2l() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::analytic_only(),
        vakint_parse!(
            "(1)*topo(\
                prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(1,2),k(2),muvsq,1)\
                *prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-2, ("-6.015223977354102649894208945578e-5".into(), "0.0".into()),),
            (-1, ("-1.804567193206230794968262683673e-4".into(),  "0.0".into()),),
            (0,  ("-3.790208765869208760159732351246e-4".into(), "0.0".into()),),
            (1,  ("-9.708141619945564042866174113010e-4".into(),  "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-6.105142965702799027412986390958e-7".into()),),
            ( 0, ("0.0".into(),  "2.548415539172476714996527540285e-6".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        vakint_parse!(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),  "1.314137110381316809259740852289e-3".into()),),
            (-2, ("0.0".into(),  "1.027363035879144648390245670512e-2".into()),),
            (-1, ("0.0".into(),  "4.561360375915191500552582268384e-2".into()),),
            ( 0, ("0.0".into(),  "1.409480043924165037345830826758e-1".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        vakint_parse!(
            "(
                  user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[
            ("muvsq".into(), 1.0), ("mursq".into(), 2.0),
            ("user_space::A".into(), 5.0), ("user_space::B".into(), 7.0), ("user_space::C".into(), 11.0)
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "1.370681288615627774073730957852e-2".into()),),
            (-2, ("0.0".into(),   "1.068331591381801491615911929647e-1".into()),),
            (-1, ("0.0".into(),   "4.730917320481948967945496964352e-1".into()),),
            ( 0, ("0.0".into(),   "1.457780409165155417799457987928".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_matad() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        EvaluationOrder::matad_only(Some(MATADOptions {direct_numerical_substition: true,..MATADOptions::default()})),
        vakint_parse!(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 2.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),   "1.314137110381316809259740852289e-3".into()),),
            (-2, ("0.0".into(),   "1.027363035879144648390245670512e-2".into()),),
            (-1, ("0.0".into(),   "4.561360375915191500552582268384e-2".into()),),
            ( 0, ("0.0".into(),   "1.409480043924165037345830826758e-1".into()),),
            ( 1, ("0.0".into(),   "5.102719154490984941123728761288e-1".into()),),
            ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_matad_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        EvaluationOrder::matad_only(Some(MATADOptions {direct_numerical_substition: true,..MATADOptions::default()})),
        vakint_parse!(
            "(
                  user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[
            ("muvsq".into(), 1.0), ("mursq".into(), 2.0),
            ("user_space::A".into(), 5.0), ("user_space::B".into(), 7.0), ("user_space::C".into(), 11.0)
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-3, ("0.0".into(),  "1.370681288615627774073730957852e-2".into()),),
            (-2, ("0.0".into(),  "1.068331591381801491615911929647e-1".into()),),
            (-1, ("0.0".into(),  "4.730917320481948967945496964352e-1".into()),),
            ( 0, ("0.0".into(),  "1.457780409165155417799457987928   ".into()),),
            ( 1, ( "0.0".into(), "5.282245494843259168589935087942   ".into()),),
            ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_3l_matad() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5,..VakintSettings::default()},
        EvaluationOrder::matad_only(None),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(),  "-6.105142965702799027412986390958e-7".into()),),
            ( 0, ("0.0".into(),  "2.548415539172476714996527540285e-6".into()),),
            ( 1, ("0.0".into(),  "-1.063440725963320582103097545956e-5".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            // This does not have an analytical expression yet (FMFT not implemented yet)
            (0,  ("-2.169283452273432986058475530569e-9".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h_squared_mass() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 2*user_space::muv - user_space::muv^2 )*topo(\
                 prop(1,edge(5,1),k(1),user_space::muv^2,1)\
                *prop(2,edge(2,6),k(2),user_space::muv^2,1)\
                *prop(3,edge(6,5),k(3),user_space::muv^2,1)\
                *prop(4,edge(3,4),k(4),user_space::muv^2,1)\
                *prop(5,edge(4,5),k(1)-k(3),user_space::muv^2,1)\
                *prop(6,edge(6,3),k(2)-k(3),user_space::muv^2,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),user_space::muv^2,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),user_space::muv^2,1)\
                *prop(9,edge(1,2),k(3)+k(4),user_space::muv^2,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("user_space::muv".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            // This does not have an analytical expression yet (FMFT not implemented yet)
            (0,  ("-2.169283452273432986058475530569e-9".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h_rank_4() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(
                  k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 5.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("1.809145974886785501452557650622e-9".into(), "0.0".into()),),
            (-3,  ("1.862208677723707446525921998109e-8".into(), "0.0".into()),),
            (-2,  ("8.865577059648962609058045604434e-8".into(), "0.0".into()),),
            (-1,  ("4.064224364096375531168559391726e-7".into(), "0.0".into()),),
            ( 0,  ("7.705260630861442737917312763495e-6".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_4l_h_rank_4_additional_symbols_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(
                  user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 5.0), ("user_space::A".into(), 5.0), ("user_space::B".into(), 7.0), ("user_space::C".into(), 11.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("9.045729874433927507262788253108e-9".into(), "0.0".into()),),
            (-3,  ("9.311043388618537232629609990544e-8".into(), "0.0".into()),),
            (-2,  ("3.740608759535330337057785114301e-7".into(), "0.0".into()),),
            (-1,  ("2.152059306128503126262830041472e-6".into(), "0.0".into()),),
            ( 0,  ("6.989489067178944692297775010846e-5".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR9d_from_H() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,0)\
                *prop(4,edge(3,4),k(4),muvsq,2)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,0)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR9d_from_X() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,0)\
                *prop(4,edge(4,3),k(4),muvsq,2)\
                *prop(5,edge(3,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,4),k(2)-k(3),muvsq,1)\
                *prop(7,edge(3,2),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(1,4),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(2,1),k(3)-k(1)-k(2)+k(4),muvsq,0)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR9d_from_H_pinch() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,2),k(1),muvsq,1)\
                *prop(2,edge(2,5),k(2),muvsq,1)\
                *prop(4,edge(3,4),k(3),muvsq,2)\
                *prop(5,edge(4,5),k(4),muvsq,1)\
                *prop(6,edge(5,3),k(4)+k(2)-k(1),muvsq,1)\
                *prop(7,edge(4,2),k(3)-k(4),muvsq,1)\
                *prop(8,edge(2,3),k(1)-k(2)+k(3)-k(4),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR9d_from_FG() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,3),k(1),muvsq,1)\
                *prop(2,edge(3,4),k(2),muvsq,2)\
                *prop(3,edge(4,5),k(3),muvsq,1)\
                *prop(4,edge(2,1),k(1)-k(3),muvsq,0)\
                *prop(5,edge(5,1),k(4),muvsq,1)\
                *prop(6,edge(4,2),k(2)-k(3),muvsq,1)\
                *prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)\
                *prop(8,edge(3,2),k(1)-k(2),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR9d_from_FG_pinch() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD, number_of_terms_in_epsilon_expansion: 5,
                         run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(5,3),k(1),muvsq,1)\
                *prop(2,edge(3,4),k(2),muvsq,2)\
                *prop(3,edge(4,5),k(3),muvsq,1)\
                *prop(5,edge(5,1),k(4),muvsq,1)\
                *prop(6,edge(4,1),k(2)-k(3),muvsq,1)\
                *prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)\
                *prop(8,edge(3,1),k(1)-k(2),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR11d() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(1,2),k(1),muvsq,2)\
                *prop(2,edge(2,5),k(2),muvsq,1)\
                *prop(3,edge(3,4),k(3),muvsq,1)\
                *prop(4,edge(4,5),k(4),muvsq,1)\
                *prop(5,edge(2,3),k(1)-k(2),muvsq,1)\
                *prop(6,edge(4,1),k(3)-k(4),muvsq,1)\
                *prop(7,edge(5,3),k(2)+k(3)-k(1),muvsq,1)\
                *prop(8,edge(1,5),k(3)-k(4)-k(1),muvsq,1)\
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (0,  ("-2.906486288643112641819206002127".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_clover() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings {
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 1)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("1.000000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("4.000000000000000000000000000000".into(), "0.0".into()),),
            (-2,  ("13.28986813369645287294483033329".into(), "0.0".into()),),
            (-1,  ("31.55672999723968577791300378449".into(), "0.0".into()),),
            (0,   ("67.98165058904685502307905531744".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_clover_with_non_unit_scales() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 1)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 7.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("81.00000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("-218.9682569565641868485693739496".into(), "0.0".into()),),
            (-2,  ("724.4490568207455247307151297051".into(), "0.0".into()),),
            (-1,  ("-1446.834846767122729283863130264".into(), "0.0".into()),),
            (0,   ("3106.843653546628093141699453745".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_dotted_clover() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { number_of_terms_in_epsilon_expansion: 5, run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 2)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("27.00000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("-99.98941898552139561618979131653".into(), "0.0".into()),),
            (-2,  ("314.4724379257699038597615012182".into(), "0.0".into()),),
            (-1,  ("-723.7613011959560846715260866565".into(), "0.0".into()),),
            (0,   ("1517.892833437916940808520861336".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_clover_with_numerator() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings {
            number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS,
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            ..VakintSettings::default() },
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {..FMFTOptions::default()} )]),
        vakint_parse!(
            "(
                  user_space::A * k(1,11)*k(2,11)*k(1,22)*k(2,22)
                + user_space::B * p(1,11)*k(3,11)*k(3,22)*p(2,22)
                + user_space::C * p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
             )*topo(
                prop(1, edge(1, 1), k(1), muvsq, 2)*\
                prop(2, edge(1, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 1), k(3), muvsq, 1)*\
                prop(4, edge(1, 1), k(4), muvsq, 1)
            )"
        ).unwrap().as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 0.3), ("mursq".into(), 0.7), ("user_space::A".into(), 3.0), ("user_space::B".into(), 4.0), ("user_space::C".into(), 5.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("-1.897149599999999855007182247846e-1".into(), "0.0".into()),),
            (-3,  ("-1.495259819655131009380566817668".into(), "0.0".into()),),
            (-2,  ("-6.805240907875078933713325181389".into(), "0.0".into()),),
            (-1,  ("-22.56027900679456203938234477552".into(), "0.0".into()),),
            (0,   ("-60.49337040949871593265194938449".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[allow(dead_code)]
fn run_integral_evaluation_analytic_tests() {
    test_integrate_1l_a();
    test_integrate_1l_simple();
    test_integrate_1l_simple_squared_mass();
    test_integrate_1l_cross_product();
    test_integrate_1l_cross_product_with_additional_symbols_numerator();
    test_integrate_1l_dot_product_external();
    test_integrate_2l();
    test_integrate_3l();
    test_integrate_3l_rank_4();
    test_integrate_3l_rank_4_additional_symbols_numerator();
    test_integrate_3l_rank_4_matad();
    test_integrate_3l_rank_4_matad_additional_symbols_numerator();
    test_integrate_3l_matad();
    test_integrate_4l_h();
    test_integrate_4l_h_squared_mass();
    test_integrate_4l_h_rank_4();
    test_integrate_4l_h_rank_4_additional_symbols_numerator();
    test_integrate_4l_PR9d_from_H();
    test_integrate_4l_PR9d_from_X();
    test_integrate_4l_PR9d_from_H_pinch();
    test_integrate_4l_PR9d_from_FG();
    test_integrate_4l_PR9d_from_FG_pinch();
    test_integrate_4l_PR11d();
    test_integrate_4l_clover();
    test_integrate_4l_clover_with_non_unit_scales();
    test_integrate_4l_dotted_clover();
    test_integrate_4l_clover_with_numerator();
}
