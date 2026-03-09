use ahash::HashMap;
use vakint::{
    EvaluationOrder, LoopNormalizationFactor, NumericalEvaluationResult, Vakint, VakintExpression,
    VakintSettings, vakint_parse, vakint_symbol,
};

fn main() {
    // Set vakint parameters
    let settings = VakintSettings {
        evaluation_order: EvaluationOrder::alphaloop_only(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        run_time_decimal_precision: 32,
        ..VakintSettings::default()
    };
    let vakint = Vakint::new().unwrap();
    vakint.validate_settings(&settings).unwrap();

    let mut integral = vakint_parse!(
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
    )
    .unwrap();
    let mut vakint_expr = VakintExpression::try_from(integral.clone()).unwrap();
    println!("\nInput integral:\n{}\n", vakint_expr);
    // Convert the numerator of the first integral to a dot notation
    vakint_expr.0[0].numerator =
        Vakint::convert_to_dot_notation(vakint_expr.0[0].numerator.as_view());
    println!("\nInput integral in dot notation:\n{}\n", vakint_expr);

    integral = vakint.evaluate(&settings, integral.as_view()).unwrap();
    println!("Evaluated integral:\n{}\n", integral.clone());

    // Set some value for the mass parameters
    let params = vakint.params_from_f64(
        &settings,
        &[("muvsq".into(), 3.0), ("mursq".into(), 5.0)]
            .iter()
            .cloned()
            .collect(),
    );

    // And for the external momenta part of the numerator
    let externals = vakint.externals_from_f64(
        &settings,
        &(1..=2)
            .map(|i| {
                (
                    i,
                    (
                        0.17 * ((i + 1) as f64),
                        0.4 * ((i + 2) as f64),
                        0.3 * ((i + 3) as f64),
                        0.12 * ((i + 4) as f64),
                    ),
                )
            })
            .collect(),
    );

    let (eval, error) = vakint
        .numerical_evaluation(
            &settings,
            integral.as_view(),
            &params,
            &HashMap::default(),
            Some(&externals),
        )
        .unwrap();
    println!("Numerical evaluation:\n{}\n", eval);
    let eval_atom = eval.to_atom(vakint_symbol!(settings.epsilon_symbol.clone()));
    println!("Numerical evaluation as atom:\n{}\n", eval_atom);
    #[rustfmt::skip]
    let target_eval =  NumericalEvaluationResult::from_vec(
    vec![
            (-3, ( "0.0".into(),  "-10202.59860843888064555902993586".into()),),
            (-2, ( "0.0".into(),  "62122.38565651740465978420334366".into()),),
            (-1, ( "0.0".into(),  "-188670.2193437045050954664088623".into()),),
            ( 0, ( "0.0".into(),  "148095.4883501202267659938351786".into()),),
        ],
        &settings);
    let (matches, match_msg) = target_eval.does_approx_match(
        &eval,
        error.as_ref(),
        10.0_f64.powi(-((settings.run_time_decimal_precision - 4) as i32)),
        1.0,
    );
    if matches {
        println!("Numerical evaluation matches target result.");
    } else {
        println!(
            "Numerical evaluation does not match target result:\n{}",
            match_msg
        );
    }
}
