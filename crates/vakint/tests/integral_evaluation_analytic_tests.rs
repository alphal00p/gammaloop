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
        ("ε^-1", "1/4*𝜋^2*𝑖*muvsq^2*g(1,2)"),
        (
            "1",
            "1/4*vakint::muvsq^2*(-𝜋^2*vakint::𝑖*log(𝜋)+𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+3/8*𝜋^2*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2)-1/4*𝜋^2*vakint::𝑖*vakint::muvsq^2*log(vakint::muvsq)*vakint::g(1,2)",
        ),
        (
            "ε",
            "-(1/4*vakint::muvsq^2*(-𝜋^2*vakint::𝑖*log(𝜋)+𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+3/8*𝜋^2*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2))*log(vakint::muvsq)+𝜋^2*vakint::𝑖*(7/16*vakint::muvsq^2*vakint::g(1,2)+1/48*𝜋^2*vakint::muvsq^2*vakint::g(1,2))+3/8*vakint::muvsq^2*(-𝜋^2*vakint::𝑖*log(𝜋)+𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+1/4*vakint::muvsq^2*(1/2*𝜋^2*vakint::𝑖*log(𝜋)^2+1/2*𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq)^2-𝜋^2*vakint::𝑖*log(𝜋)*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+1/8*𝜋^2*vakint::𝑖*vakint::muvsq^2*log(vakint::muvsq)^2*vakint::g(1,2)",
        ),
        (
            "ε^2",
            "(7/16*vakint::muvsq^2*vakint::g(1,2)+1/48*𝜋^2*vakint::muvsq^2*vakint::g(1,2))*(-𝜋^2*vakint::𝑖*log(𝜋)+𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq))+1/2*(1/4*vakint::muvsq^2*(-𝜋^2*vakint::𝑖*log(𝜋)+𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+3/8*𝜋^2*vakint::𝑖*vakint::muvsq^2*vakint::g(1,2))*log(vakint::muvsq)^2-(𝜋^2*vakint::𝑖*(7/16*vakint::muvsq^2*vakint::g(1,2)+1/48*𝜋^2*vakint::muvsq^2*vakint::g(1,2))+3/8*vakint::muvsq^2*(-𝜋^2*vakint::𝑖*log(𝜋)+𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+1/4*vakint::muvsq^2*(1/2*𝜋^2*vakint::𝑖*log(𝜋)^2+1/2*𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq)^2-𝜋^2*vakint::𝑖*log(𝜋)*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2))*log(vakint::muvsq)+𝜋^2*vakint::𝑖*(-1/12*vakint::muvsq^2*vakint::z3*vakint::g(1,2)+15/32*vakint::muvsq^2*vakint::g(1,2)+1/32*𝜋^2*vakint::muvsq^2*vakint::g(1,2))+3/8*vakint::muvsq^2*(1/2*𝜋^2*vakint::𝑖*log(𝜋)^2+1/2*𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq)^2-𝜋^2*vakint::𝑖*log(𝜋)*log(1/4*𝜋^-1*vakint::mursq))*vakint::g(1,2)+1/4*vakint::muvsq^2*(-1/6*𝜋^2*vakint::𝑖*log(𝜋)^3+1/6*𝜋^2*vakint::𝑖*log(1/4*𝜋^-1*vakint::mursq)^3+1/2*𝜋^2*vakint::𝑖*log(𝜋)^2*log(1/4*𝜋^-1*vakint::mursq)-1/2*𝜋^2*vakint::𝑖*log(𝜋)*log(1/4*𝜋^-1*vakint::mursq)^2)*vakint::g(1,2)-1/24*𝜋^2*vakint::𝑖*vakint::muvsq^2*log(vakint::muvsq)^3*vakint::g(1,2)",
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
        // println!("{} -> {}", v, c.to_canonical_string());
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
       vakint_parse!("-2879700𝑖/536411*g(1,2)+3726809𝑖/395976*ε*g(1,2)+1075967𝑖/436073*ε^-1*g(1,2)-4041047𝑖/334810*ε^2*g(1,2)").unwrap()
    );

    let prec: Rational = (
        1,
        10_u64.pow((vakint.settings.run_time_decimal_precision - 4) as u32),
    )
        .into();

    _ = compare_output(
        Ok(numerical_partial_eval.rationalize(&prec).as_view()),
        vakint_parse!(
            format!(
                "-5.36845814123893`{prec}𝑖*g(1,2)\
             +9.4117042447109682`{prec}𝑖*ε*g(1,2)\
             +2.46740110027234`{prec}𝑖*ε^-1*g(1,2)\
             -12.0696723514860`{prec}𝑖*ε^2*g(1,2)",
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
            (-1, ("0.0".into(), "2.46740110027234".into())),
            (0, ("0.0".into(), "-5.36845814123893".into())),
            (1, ("0.0".into(), "9.41170424471097".into())),
            (2, ("0.0".into(), "-12.0696723514860".into())),
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
            (-1, ("0.0".into(), "-15.41829599538179739876641630989".into()),),
            (0,  ("0.0".into(),  "41.25556923066471696715797280407".into()),),
            (1,  ("0.0".into(), "-75.58506810083682014451198579381".into()),),
            (2,  ("0.0".into(),  "104.8268973321405776943695037314".into()),),
            (3,  ("0.0".into(),  "-130.2125934660588794479583559489".into()),),
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
            (-1, ("0.0".into(), "-6.776470381787957720961284930165".into()),),
            (0,  ("0.0".into(), "-21.34624897436485396509954629830".into()),),
            (1,  ("0.0".into(), "72.41426780477804749630966479154".into()),),
            (2,  ("0.0".into(), "-147.4626326303287005964359960851".into()),),
            (3,  ("0.0".into(), "211.1788648410998597813244933220".into()),),
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
            (-1, ("0.0".into(), "-7298.572454605580698628106094408".into()),),
            (0,  ("0.0".into(),  "15879.89918178474997676561014674".into()),),
            (1,  ("0.0".into(), "-27839.82115585504758906774113191".into()),),
            (2,  ("0.0".into(),  "35702.09081569554005452137401203".into()),),
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
            (-2, ("-146.1136365510036558546604990331".into(), "0.0".into()),),
            (-1, ("635.8146971740286947808753047759".into(),  "0.0".into()),),
            (0,  ("-1646.531034471454483109201793220".into(), "0.0".into()),),
            (1,  ("2240.516116133454318298096346441".into(),  "0.0".into()),),
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
            (-1, ("0.0".into(), "-2311.289033520460340396770711738".into()),),
            ( 0, ("0.0".into(),  "35134.99893627257345553503414002".into()),),
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
            (-3, ("0.0".into(),   "4975068.903103003548576756843470".into()),),
            (-2, ("0.0".into(),  "-15967412.96033300288621485195252".into()),),
            (-1, ("0.0".into(),   "46275660.33034806160550888444548".into()),),
            ( 0, ("0.0".into(),  "-117731367.8844665539198405383934".into()),),
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
            (-3, ("0.0".into(),   "51891342.24417464772786998080343".into()),),
            (-2, ("0.0".into(),  "-167769771.6403982559590402103780".into()),),
            (-1, ("0.0".into(),   "486068084.1087102468110182518293".into()),),
            ( 0, ("0.0".into(),  "-1237710265.809503746583431557798".into()),),
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
            (-3, ("0.0".into(),   "4975068.903103003548576756843470".into()),),
            (-2, ("0.0".into(),  "-15967412.96033300288621485195252".into()),),
            (-1, ("0.0".into(),   "46275660.33034806160550888444548".into()),),
            ( 0, ("0.0".into(),  "-117731367.8844665539198405383934".into()),),
            ( 1, ("0.0".into(),  "919780256.2106401071849409513709".into()),),
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
            (-3, ("0.0".into(),   "51891342.24417464772786998080343".into()),),
            (-2, ("0.0".into(),  "-167769771.6403982559590402103780".into()),),
            (-1, ("0.0".into(),   "486068084.1087102468110182518293".into()),),
            ( 0, ("0.0".into(),  "-1237710265.809503746583431557798".into()),),
            ( 1, ( "0.0".into(),  "9616506746.507901578686816427718".into()),),
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
            (-1, ("0.0".into(), "-2311.289033520460340396770711738".into()),),
            ( 0, ("0.0".into(),  "35134.99893627257345553503414002".into()),),
            ( 1, ("0.0".into(),  "-287175.6919292485174272232526581".into()),),
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
            (0,  ("-12799.53514305961548130719263292".into(), "0.0".into()),),
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
            (0,  ("-12799.53514305961548130719263292".into(), "0.0".into()),),
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
            (-4,  ("+10674.59739307939575801964744189".into(), "0.0".into()),),
            (-3,  ("-47071.92195516143798974687670839".into(), "0.0".into()),),
            (-2,  ("+61389.19910667747809242280299166".into(), "0.0".into()),),
            (-1,  ("+928571.7264940044122518812620547".into(), "0.0".into()),),
            ( 0,  ("+29325841.99455641744844958694905".into(), "0.0".into()),),
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
            (-4,  ("+53372.98696539697879009823720947".into(), "0.0".into()),),
            (-3,  ("-235359.6097758071899487343835420".into(), "0.0".into()),),
            (-2,  ("-101464.4077196306917939828600505".into(), "0.0".into()),),
            (-1,  ("+11355453.75154459334966095599239".into(), "0.0".into()),),
            ( 0,  ("+277164323.3312205706615743239968".into(), "0.0".into()),),
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
