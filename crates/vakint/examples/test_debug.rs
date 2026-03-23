use test_log::env_logger;
use vakint::{InputFloatRationalizationPrecision, Vakint, VakintExpression, VakintSettings};

fn main() {
    env_logger::init();
    let vakint = Vakint::new().unwrap();

    let settings = VakintSettings {
        allow_unknown_integrals: false,
        use_dot_product_notation: false,
        integral_normalization_factor: vakint::LoopNormalizationFactor::MSbar,
        run_time_decimal_precision: 16,
        number_of_terms_in_epsilon_expansion: 2,
        mu_r_sq_symbol: "gammalooprs::μᵣ²".to_string(),
        epsilon_symbol: "gammalooprs::ε".to_string(),
        verify_numerator_identification: false,
        precision_for_input_float_rationalization:
            InputFloatRationalizationPrecision::FullPrecision,
        evaluation_order: vakint::EvaluationOrder::matad_only(None),
        ..VakintSettings::default()
    };

    println!("{:#?}", settings);

    // let input = fs::read_to_string("./examples/test_debug_input_lucien.txt")
    //     .expect("failed to read ./examples/test_debug_input_lucien.txt");

    // let mut integral = //symbolica::atom::Atom::Zero;
    // symbolica::parse!(input.trim());

    let mut integral = symbolica::parse!(
        "(
             1.1243423413786486276348761294716297346178236`100
        )*vakint::topo(\
        vakint::prop(9,vakint::edge(66,66),vakint::k(3),MUVsq,1)\
    )"
    );

    let full = vakint.evaluate(&settings, integral.as_view()).unwrap();

    println!(
        "\nInput integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint
        .to_canonical(&settings, integral.as_view(), true)
        .unwrap();
    println!(
        "Matched integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint.tensor_reduce(&settings, integral.as_view()).unwrap();
    println!(
        "Tensor reduced integral:\n{}\n",
        VakintExpression::try_from(integral.clone()).unwrap()
    );

    integral = vakint
        .evaluate_integral(&settings, integral.as_view())
        .unwrap();
    println!("Evaluated integral:\n{}\n", integral);

    // let mut params = HashMap::default();
    // params.insert(
    //     "oneloop_evaluation::MUVsq".into(),
    //     vakint.settings.real_to_prec("1.0"),
    // );
    // params.insert("test::myMUV".into(), vakint.settings.real_to_prec("1.0"));
    // params.insert("test::myMUVI".into(), vakint.settings.real_to_prec("1.0"));
    // params.insert("test::mu_r_sq".into(), vakint.settings.real_to_prec("1.0"));
    // let numerical_partial_eval = Vakint::partial_numerical_evaluation(
    //     &vakint.settings,
    //     integral.as_view(),
    //     &params,
    //     &HashMap::default(),
    //     None,
    // );
    // println!("Partial eval:\n{}\n", numerical_partial_eval);

    // let externals = vakint.externals_from_f64(
    //     &(1..=2)
    //         .map(|i| {
    //             (
    //                 i,
    //                 (
    //                     0.17 * ((i + 1) as f64),
    //                     0.4 * ((i + 2) as f64),
    //                     0.3 * ((i + 3) as f64),
    //                     0.12 * ((i + 4) as f64),
    //                 ),
    //             )
    //         })
    //         .collect(),
    // );

    // params.insert(
    //     "vakint::g(oneloop_evaluation::mink4(4,22),oneloop_evaluation::mink4(4,11))".into(),
    //     vakint.settings.real_to_prec("1.0"),
    // );
    // params.insert(
    //     "vakint::p(1,oneloop_evaluation::mink4(4,99))".into(),
    //     vakint.settings.real_to_prec("1.0"),
    // );
    // params.insert(
    //     "vakint::p(2,oneloop_evaluation::mink4(4,101))".into(),
    //     vakint.settings.real_to_prec("1.0"),
    // );
    // let numerical_full_eval = Vakint::full_numerical_evaluation_without_error(
    //     &vakint.settings,
    //     integral.as_view(),
    //     &params,
    //     &HashMap::default(),
    //     Some(&externals),
    // )
    // .unwrap();
    // println!(
    //     "Full eval (tensor structures left all substituted with 1):\n{}\n",
    //     numerical_full_eval
    // );
}
