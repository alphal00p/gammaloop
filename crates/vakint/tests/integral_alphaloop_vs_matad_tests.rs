mod test_utils;

use ahash::HashMap;
use log::debug;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::rational::Rational,
    function,
};
use test_utils::{TestVakint, compare_two_evaluations, get_vakint};
use vakint::{
    EvaluationOrder, NumericalEvaluationResult, Vakint, VakintSettings, matad::MATAD, vakint_parse,
    vakint_symbol,
};

use vakint::{externals_from_f64, params_from_f64};

const N_DIGITS_ANLYTICAL_EVALUATION: u32 = 32;
const COMPARISON_REL_THRESHOLD: f64 = 1.0e-25;
// No error on anylic expressions
const MAX_PULL: f64 = 0.0e0;

#[test_log::test]
fn test_integrate_1l_no_numerator() {
    for pow in 1..=6 {
        #[rustfmt::skip]
        compare_two_evaluations(
            VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 4, ..VakintSettings::default()},
            ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
            vakint_parse!(
                format!("( 1 )
                *topo(\
                    prop(1,edge(1,1),k(1),muvsq,{})
                )",pow).as_str()
            ).unwrap().as_view(),
            params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            externals_from_f64(
            &(1..=2)
                .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
                //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
                .collect(),
                N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
            true,
        );
    }
}

#[test_log::test]
fn test_integrate_1l_no_numerator_squared_mass() {
    for pow in 1..=2 {
        #[rustfmt::skip]
        compare_two_evaluations(
            VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 4, ..VakintSettings::default()},
            ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
            vakint_parse!(
                format!("( 2*user_space::muv - user_space::muv^2 )
                *topo(\
                    prop(1,edge(1,1),k(1),user_space::muv^2,{})
                )",pow).as_str()
            ).unwrap().as_view(),
            params_from_f64(&[("user_space::muv".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            externals_from_f64(
            &(1..=2)
                .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
                //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
                .collect(),
                N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
            true,
        );
    }
}

#[test_log::test]
fn test_integrate_2l_no_numerator() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 4, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
        vakint_parse!(
            "(1)*topo(\
            prop(1,edge(1,2),k(1),muvsq,1)\
            * prop(2,edge(1,2),k(2),muvsq,1)\
            * prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
        )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
        COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_basketball_a() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 3, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
        vakint_parse!(
            "( k(1, 1)^2 )*topo(
                prop(1, edge(2, 1), k(1), muvsq, 1)*\
                prop(2, edge(2, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 2), k(3), muvsq, 1)*\
                prop(4, edge(1, 2), k(1) + k(2) - k(3), muvsq, 1)
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_basketball_b() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 3, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
        vakint_parse!(
            "( k(3, 1)^2 )*topo(
                prop(1, edge(2, 1), k(1), muvsq, 1)*\
                prop(2, edge(2, 1), k(2), muvsq, 1)*\
                prop(3, edge(1, 2), k(3), muvsq, 1)*\
                prop(4, edge(1, 2), k(1) + k(2) - k(3), muvsq, 1)
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_no_numerator() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 4, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
        vakint_parse!(
            "( 1 )
            *topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        ).unwrap().as_view(),
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 4, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
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
        N_DIGITS_ANLYTICAL_EVALUATION),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

#[test_log::test]
fn test_integrate_3l_rank_4_different_scales() {
    #[rustfmt::skip]
    compare_two_evaluations(
        VakintSettings { run_time_decimal_precision: N_DIGITS_ANLYTICAL_EVALUATION, number_of_terms_in_epsilon_expansion: 4, ..VakintSettings::default()},
        ((&EvaluationOrder::alphaloop_only() ,true),(&EvaluationOrder::matad_only(None) ,true)),
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
        params_from_f64(&[("muvsq".into(), 3.0), ("mursq".into(), 5.0)].iter().cloned().collect(),
        N_DIGITS_ANLYTICAL_EVALUATION),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            //.map(|i| (i, (17.0*((0) as f64), 4.0*((0) as f64), 3.0*((0) as f64), 12.0*((0) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION),
            COMPARISON_REL_THRESHOLD, MAX_PULL,
        true,
    );
}

pub fn evaluate_expression_with_matad(
    vakint: &TestVakint,
    input: AtomView,
    direct_masters_substitution: bool,
) -> NumericalEvaluationResult {
    let mut integral = input.to_owned();

    let matad = MATAD::with_settings(vakint.settings.clone());
    if direct_masters_substitution {
        integral = matad
            .substitute_masters_directly(integral.as_view())
            .unwrap();
    } else {
        integral = matad.expand_matad_masters(integral.as_view()).unwrap();
    }
    // println!(
    //     "NEW vakint.settings.number_of_terms_in_epsilon_expansion:\n{}",
    //     vakint.settings.number_of_terms_in_epsilon_expansion
    // );
    // println!("NEW Integral:\n{}", integral);
    // integral = integral.expand()
    /*
    println!(
        "OLD vakint.settings.number_of_terms_in_epsilon_expansion:\n{}",
        vakint.settings.number_of_terms_in_epsilon_expansion
    );
    println!("OLD Integral:\n{}", integral);
    */
    integral = integral
        .series(
            vakint_symbol!("ep"),
            Atom::Zero.as_view(),
            Rational::from(vakint.settings.number_of_terms_in_epsilon_expansion - 3),
            true,
        )
        .unwrap()
        .to_atom();

    if direct_masters_substitution {
        // Expanding here is important to improve efficiency and avoid symbolica bugs with floating point coefficients
        integral = integral.expand();
        integral = matad.substitute_poly_gamma(integral.as_view()).unwrap();
        integral = matad
            .substitute_additional_constants(integral.as_view())
            .unwrap();
    } else {
        integral = matad.substitute_masters(integral.as_view()).unwrap();
        // Expanding here is important to improve efficiency and avoid symbolica bugs with floating point coefficients
        integral = integral.expand();
        integral = matad.substitute_hpls(integral.as_view()).unwrap();
        integral = matad.substitute_poly_gamma(integral.as_view()).unwrap();
        integral = matad
            .substitute_additional_constants(integral.as_view())
            .unwrap();
    }
    let muv_sq_symbol = vakint_symbol!("muvsq");
    let log_muv_mu_sq = function!(
        Symbol::LOG,
        Atom::var(muv_sq_symbol)
            / Atom::var(vakint_symbol!(vakint.settings.mu_r_sq_symbol.as_str()))
    );

    let log_mu_sq = function!(
        Symbol::LOG,
        Atom::var(vakint_symbol!(vakint.settings.mu_r_sq_symbol.as_str()))
    );

    integral = integral
        .replace(vakint_parse!("logmUVmu").unwrap().to_pattern())
        .with((log_muv_mu_sq).to_pattern());
    integral = integral
        .replace(vakint_parse!("log_mu_sq").unwrap().to_pattern())
        .with((log_mu_sq).to_pattern());

    let test_params = params_from_f64(
        &[("muvsq".into(), 1.0), ("mursq".into(), 1.0)]
            .iter()
            .cloned()
            .collect(),
        vakint.settings.run_time_decimal_precision,
    );

    integral = integral
        .replace(vakint_parse!("ep").unwrap().to_pattern())
        .with(Atom::var(vakint_symbol!(vakint.settings.epsilon_symbol.as_str())).to_pattern());

    let (res, _error) = Vakint::full_numerical_evaluation(
        &vakint.settings,
        integral.as_view(),
        &test_params,
        &HashMap::default(),
        None,
    )
    .unwrap();

    res
}

#[test_log::test]
pub fn test_eval_matad_masters() {
    let mut vakint = get_vakint(VakintSettings {
        evaluation_order: EvaluationOrder::empty(),
        run_time_decimal_precision: 40,
        number_of_terms_in_epsilon_expansion: 4,
        ..VakintSettings::default()
    });

    #[rustfmt::skip]
    let expression_to_tests = vec![
        (
            vakint_parse!("miD6").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -1, ("  2.40411380631918857079947632302289998153".into(), "0.0".into()) ),
                (  0, ("-10.03527847976878917191470068515890023865".into(), "0.0".into()) ),
                (  1, (" 41.87670208303157617433490267097099146643".into(), "0.0".into()) ),
                (  2, ("-146.80128953959941603962375965123680914875".into(), "0.0".into()) ),
                ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miT111").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -2, ("-1.5".into(), "0.0".into()) ),
                ( -1, ("-4.5".into(), "0.0".into()) ),
                (  0, ("-9.451540242238151318806279316660622180185".into(), "0.0".into()) ),
                (  1, ("-24.20892802120359267872133821957094892580".into(), "0.0".into()) ),
                (  2, ("-38.71759974491583887231661364177794370942".into(), "0.0".into()) ),
                (  3, ("-101.953991009597114422662478576975778730".into(), "0.0".into()) ),
            ], &vakint.settings),
            6
        ),
        (
            vakint_parse!("miD5").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -1, ("2.404113806319188570799476323022899981530".into(), "0.0".into()) ),
                (  0, ("-8.216859817508738062913398338601085824970".into(), "0.0".into()) ),
                (  1, ("36.47368421155096825994471856976365848559".into(), "0.0".into()) ),
                (  2, ("-122.5028439280736143862645277852881354128".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miBN").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -3, ("2.0".into(), "0.0".into()) ),
                ( -2, ("7.666666666666666666666666666666666666667".into(), "0.0".into()) ),
                ( -1, ("22.43480220054467930941724549993807556766".into(), "0.0".into()) ),
                (  0, ("39.42929462910208211529996476007305636115".into(), "0.0".into()) ),
                (  1, ("62.92709375535910070547798992048699891696".into(), "0.0".into()) ),
                (  2, ("-126.3366653990120700722498217033370728331".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("Gam(1,1)").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                (  0, ("1.0".into(), "0.0".into()) ),
                (  1, ("0.0".into(), "0.0".into()) ),
                (  2, (" 0.8224670334241132182362075833230125946".into(), "0.0".into()) ),
                (  3, ("-0.40068563438653142846657938717048333026".into(), "0.0".into()) ),
                (  4, (" 0.60880681896251523272775207930440694543".into(), "0.0".into()) ),
            ], &vakint.settings),
            7
        ),
        (
            vakint_parse!("iGam(1,1)").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                (  0, ("1.0".into(), "0.0".into()) ),
                (  1, ("0.0".into(), "0.0".into()) ),
                (  2, ("-0.8224670334241132182362075833230125946".into(), "0.0".into()) ),
                (  3, (" 0.40068563438653142846657938717048333026".into(), "0.0".into()) ),
                (  4, (" 0.067645202106946136969750231033822993923423".into(), "0.0".into()) ),
            ], &vakint.settings),
            7
        ),
        (
            vakint_parse!("Gam(1,2)").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                (  0, ("1.0".into(), "0.0".into()) ),
                (  1, ("0.0".into(), "0.0".into()) ),
                (  2, (" 3.289868133696452872944830333292050378438".into(), "0.0".into()) ),
                (  3, ("-3.205485075092251427732635097363866642040".into(), "0.0".into()) ),
                (  4, (" 9.740909103400243723644033268870511124973".into(), "0.0".into()) ),
            ], &vakint.settings),
            7
        ),
        (
            vakint_parse!("iGam(1,2)").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                (  0, ("1.0".into(), "0.0".into()) ),
                (  1, ("0.0".into(), "0.0".into()) ),
                (  2, ("-3.289868133696452872944830333292050378438".into(), "0.0".into()) ),
                (  3, (" 3.205485075092251427732635097363866642040".into(), "0.0".into()) ),
                (  4, (" 1.082323233711138191516003696541167902775".into(), "0.0".into()) ),
            ], &vakint.settings),
            7
        ),
        (
            vakint_parse!("miD4").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -1, ("2.404113806319188570799476323022899981530".into(), "0.0".into()) ),
                (  0, ("-5.913204783884020530495717892535405026883".into(), "0.0".into()) ),
                (  1, ("31.79387586520335091502703130598293231890".into(), "0.0".into()) ),
                (  2, ("-95.53186858506048154153099646099149596351".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miD3").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -1, ("2.404113806319188570799476323022899981530".into(), "0.0".into()) ),
                (  0, ("-3.027009493987652019786374701758957286151".into(), "0.0".into()) ),
                (  1, ("28.73643511952363680901046995862299631519".into(), "0.0".into()) ),
                (  2, ("-63.46100358831692185768876871928863693891".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miDM").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -1, ("2.404113806319188570799476323022899981530".into(), "0.0".into()) ),
                (  0, ("-2.860862224139327350272784567773241917561".into(), "0.0".into()) ),
                (  1, ("29.00667443783775908331902631581722417518".into(), "0.0".into()) ),
                (  2, ("-62.36139634229625160648107039345983094078".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miDN").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -1, ("2.404113806319188570799476323022899981530".into(), "0.0".into()) ),
                (  0, ("1.120248397039242082272516548224209526276".into(), "0.0".into()) ),
                (  1, ("30.68103527534589055078588298235627556530".into(), "0.0".into()) ),
                (  2, ("-13.30346064085824836188894234098890655258".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miE3").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -3, ("-0.6666666666666666666666666666666666666667".into(), "0.0".into()) ),
                ( -2, ("-3.666666666666666666666666666666666666667".into(), "0.0".into()) ),
                ( -1, ("-13.77400727566226453704248689998363477479".into(), "0.0".into()) ),
                (  0, ("-55.65962246120633017136139508642412198254".into(), "0.0".into()) ),
                (  1, ("-151.9352362017674553145984030189617080164".into(), "0.0".into()) ),
                (  2, ("-574.6540576129672528661511219788586831297".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
        (
            vakint_parse!("miBN1").unwrap(),
            NumericalEvaluationResult::from_vec(vec![
                ( -3, ("1.0".into(), "0.0".into()) ),
                ( -2, ("3.75".into(), "0.0".into()) ),
                ( -1, ("10.59240110027233965470862274996903778383".into(), "0.0".into()) ),
                (  0, ("21.76198850991296192361111230083506510406".into(), "0.0".into()) ),
                (  1, ("7.851742836425531175775595591525685030810".into(), "0.0".into()) ),
                (  2, ("-71.05207017591209503800257696579714623490".into(), "0.0".into()) ),
            ], &vakint.settings),
            5
        ),
    ];
    for (input, target, number_of_expansion_terms) in &expression_to_tests {
        vakint.settings.number_of_terms_in_epsilon_expansion = *number_of_expansion_terms;
        let res = evaluate_expression_with_matad(&vakint, input.as_view(), true);
        let (matches, msg) = res.does_approx_match(
            target,
            None,
            0.1_f64.powi((vakint.settings.run_time_decimal_precision - 5) as i32),
            0.,
        );
        println!(
            "MATAD evaluation with direct master substitution of expression: {}\nyields:\n{}",
            input, res
        );
        assert!(matches, "mismatch in evaluation of {}: {}", input, msg);
    }

    for (input, target, number_of_expansion_terms) in &expression_to_tests {
        vakint.settings.number_of_terms_in_epsilon_expansion = *number_of_expansion_terms;
        let res = evaluate_expression_with_matad(&vakint, input.as_view(), false);
        let (matches, msg) = res.does_approx_match(
            target,
            None,
            0.1_f64.powi((vakint.settings.run_time_decimal_precision - 5) as i32),
            0.,
        );
        println!(
            "MATAD evaluation without direct master substitution of expression: {}\nyields:\n{}",
            input, res
        );
        assert!(matches, "mismatch in evaluation of {}: {}", input, msg);
    }
}

#[test_log::test]
pub fn test_eval_matad_one_master_combination() {
    let vakint = get_vakint(VakintSettings {
        evaluation_order: EvaluationOrder::empty(),
        run_time_decimal_precision: 40,
        number_of_terms_in_epsilon_expansion: 4,
        ..VakintSettings::default()
    });

    let mut input = vakint_parse!("M^2*miD6+16*M^2*(1736*(-2*ep+4)^2-718*(-2*ep+4)^3+165*(-2*ep+4)^4-20*(-2*ep+4)^5+(-2*ep+4)^6-2208*(-2*ep+4)+1152)^-1*Gam(1,1)^3+4*M^2*miT111*((-2*ep+4)^2-7*(-2*ep+4)+12)^-1*Gam(1,1)+M^2*miD5*(2*(-2*ep+4)-6)^-1*(3*(-2*ep+4)-12)+M^2*miBN*(-3*(-2*ep+4)+8)*(8*(-2*ep+4)-24)^-1").unwrap();
    input = input
        .replace(vakint_parse!("M").unwrap().to_pattern())
        .with(vakint_parse!("muvsq^(1/2)").unwrap().to_pattern());
    let res = evaluate_expression_with_matad(&vakint, input.as_view(), true);
    #[rustfmt::skip]
    let (matches, msg) = res.does_approx_match(
        &NumericalEvaluationResult::from_vec(vec![
            ( -3, ("1.0".into(), "0.0".into()) ),
            ( -2, ("5.66666666666666666666666666666666666667".into(),   "0.0".into()) ),
            ( -1, ("20.17312652385648488703674553970843989141".into(),  "0.0".into()) ),
            (  0, ("48.81971591461289153001703696320240598328".into(),  "0.0".into()) ),
            (  1, ("225.02689067182955254127148271373064571626".into(), "0.0".into()) ),
        ], &vakint.settings),
        None,
        0.1_f64.powi((vakint.settings.run_time_decimal_precision -2) as i32),
        0.,
    );
    debug!(
        "MATAD evaluation with direct master substitution of expression: {}\nyields:\n{}",
        input, res
    );
    assert!(matches, "{}", msg);
    let res = evaluate_expression_with_matad(&vakint, input.as_view(), false);
    #[rustfmt::skip]
    let (matches, msg) = res.does_approx_match(
        &NumericalEvaluationResult::from_vec(vec![
            ( -3, ("1.0".into(), "0.0".into()) ),
            ( -2, ("5.66666666666666666666666666666666666667".into(),   "0.0".into()) ),
            ( -1, ("20.17312652385648488703674553970843989141".into(),  "0.0".into()) ),
            (  0, ("48.81971591461289153001703696320240598328".into(),  "0.0".into()) ),
            (  1, ("225.02689067182955254127148271373064571626".into(), "0.0".into()) ),
        ], &vakint.settings),
        None,
        0.1_f64.powi((vakint.settings.run_time_decimal_precision -2) as i32),
        0.,
    );
    debug!(
        "MATAD evaluation without direct master substitution of expression: {}\nyields:\n{}",
        input, res
    );
    assert!(matches, "{}", msg);
}

#[allow(dead_code)]
fn run_integral_alphaloop_vs_matad_tests() {
    // Convenience runner to execute all tests in this module.
    test_integrate_1l_no_numerator();
    test_integrate_1l_no_numerator_squared_mass();
    test_integrate_2l_no_numerator();
    test_integrate_3l_basketball_a();
    test_integrate_3l_basketball_b();
    test_integrate_3l_no_numerator();
    test_integrate_3l_rank_4();
    test_integrate_3l_rank_4_different_scales();
    test_eval_matad_masters();
    test_eval_matad_one_master_combination();
}
