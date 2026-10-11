use super::{ArbPrec, Complex, F, QuadFloat};
use symbolica::{
    atom::{Atom, AtomCore, EvaluationInfo},
    domains::{
        float::{DoubleFloat, Float, Real},
        rational::Rational,
    },
    evaluate::{EvaluationDomain, FunctionMap, OptimizationSettings},
    function, symbol,
};

#[test]
fn direct_wrapper_callbacks_take_precedence() {
    let info = EvaluationInfo::new()
        .register(|args: &[f64]| args[0])
        .register(|args: &[F<f64>]| args[1])
        .register(|args: &[DoubleFloat]| args[0])
        .register(|args: &[QuadFloat]| args[1])
        .register(|args: &[Float]| args[0].clone())
        .register(|args: &[ArbPrec]| args[1].clone());
    assert_eq!(
        F::<f64>::resolve_function(&[], &info).unwrap()(&[F(2.0), F(3.0)]),
        F(3.0)
    );
    let inputs = [Rational::from(2), Rational::from(3)];
    let quad = inputs.each_ref().map(QuadFloat::from);
    assert_eq!(
        QuadFloat::resolve_function(&[], &info).unwrap()(&quad),
        quad[1]
    );
    let precise = inputs.each_ref().map(ArbPrec::from);
    assert_eq!(
        ArbPrec::resolve_function(&[], &info).unwrap()(&precise),
        precise[1]
    );
}

#[test]
fn underlying_callbacks_receive_tags_and_all_arguments() {
    let info = EvaluationInfo::new()
        .with_tags(1)
        .register_tagged(|tags| {
            assert_eq!(tags, &[Atom::num(7).as_view()]);
            Box::new(|args: &[f64]| args[0] - args[1])
        })
        .register_tagged(|tags| {
            assert_eq!(tags, &[Atom::num(7).as_view()]);
            Box::new(|args: &[DoubleFloat]| args[0] - args[1])
        })
        .register_tagged(|tags| {
            assert_eq!(tags, &[Atom::num(7).as_view()]);
            Box::new(|args: &[Float]| args[0].clone() - &args[1])
        });
    let tag = Atom::num(7);
    let tags = [tag.as_view()];
    assert_eq!(
        F::<f64>::resolve_function(&tags, &info).unwrap()(&[F(7.0), F(2.0)]),
        F(5.0)
    );
    let inputs = [Rational::from(7), Rational::from(2)];
    let quad = inputs.each_ref().map(QuadFloat::from);
    assert_eq!(
        QuadFloat::resolve_function(&tags, &info).unwrap()(&quad),
        QuadFloat::from(&Rational::from(5))
    );
    let precise = inputs.each_ref().map(ArbPrec::from);
    assert_eq!(
        ArbPrec::resolve_function(&tags, &info).unwrap()(&precise),
        ArbPrec::from(&Rational::from(5))
    );
}

#[test]
fn unregistered_callbacks_remain_unresolved() {
    let info = EvaluationInfo::new();
    assert!(F::<f64>::resolve_function(&[], &info).is_none());
    assert!(QuadFloat::resolve_function(&[], &info).is_none());
    assert!(ArbPrec::resolve_function(&[], &info).is_none());
}

#[test]
fn native_hyperbolic_evaluators_preserve_wrapper_precision() {
    use symbolica::transcendental::{coth, csch, sech, tanh};

    let x = symbol!("evaluation_domain_x").to_atom();
    let input = Rational::from((1_i128, 1_i128 << 80)) + Rational::from(1);
    let precise_input = F::<ArbPrec>::from(&input);
    let cases = [
        (tanh(), precise_input.tanh()),
        (coth(), precise_input.tanh().inv()),
        (sech(), precise_input.cosh().inv()),
        (csch(), precise_input.sinh().inv()),
    ];
    let tolerance = F::<QuadFloat>::from(&Rational::from((1_i128, 1_i128 << 98)));
    for (symbol, expected) in cases {
        let evaluator = function!(symbol, x.clone())
            .as_view()
            .to_evaluation_tree(&FunctionMap::new(), std::slice::from_ref(&x))
            .unwrap()
            .linearize(&OptimizationSettings::default());
        let result = evaluator
            .clone()
            .map_coeff(&|c| F::<f64>::from(&c.re))
            .evaluate_single(&[F::<f64>::from(&input)]);
        assert!((result.0 - expected.into_ff64().0).abs() < 1e-15);

        let result = evaluator
            .clone()
            .map_coeff(&|c| F::<QuadFloat>::from(&c.re))
            .evaluate_single(&[F::<QuadFloat>::from(&input)]);
        let expected_quad = F(QuadFloat::from(Float::from(expected.clone())));
        assert!((result - expected_quad).abs() < tolerance);

        let result = evaluator
            .map_coeff(&|c| F::<ArbPrec>::from(&c.re))
            .evaluate_single(std::slice::from_ref(&precise_input));
        // Equivalent hyperbolic formulas need not round identically.
        let tolerance = expected.epsilon() * expected.abs() * expected.from_usize(4);
        assert!((result - expected).abs() <= tolerance);
    }
}

#[test]
fn native_hyperbolic_evaluators_support_complex_wrappers() {
    use symbolica::transcendental::{coth, csch, sech, tanh};

    let x = symbol!("evaluation_domain_complex_x").to_atom();
    let input = Complex::new(F(0.75), F(0.5));
    let cases = [
        (tanh(), input.tanh()),
        (coth(), input.tanh().inv()),
        (sech(), input.cosh().inv()),
        (csch(), input.sinh().inv()),
    ];
    for (symbol, expected) in cases {
        let evaluator = function!(symbol, x.clone())
            .as_view()
            .to_evaluation_tree(&FunctionMap::new(), std::slice::from_ref(&x))
            .unwrap()
            .linearize(&OptimizationSettings::default());
        let double = evaluator
            .clone()
            .map_coeff(&|c| Complex::new(F::<f64>::from(&c.re), F::<f64>::from(&c.im)))
            .evaluate_single(&[input]);
        let quad = evaluator
            .clone()
            .map_coeff(&|c| Complex::new(F::<QuadFloat>::from(&c.re), F::<QuadFloat>::from(&c.im)))
            .evaluate_single(&[Complex::new(
                F::<QuadFloat>::from_ff64(input.re),
                F::<QuadFloat>::from_ff64(input.im),
            )]);
        let precise = evaluator
            .map_coeff(&|c| Complex::new(F::<ArbPrec>::from(&c.re), F::<ArbPrec>::from(&c.im)))
            .evaluate_single(&[Complex::new(
                F::<ArbPrec>::from_ff64(input.re),
                F::<ArbPrec>::from_ff64(input.im),
            )]);
        for result in [
            double,
            Complex::new(quad.re.into_ff64(), quad.im.into_ff64()),
            Complex::new(precise.re.into_ff64(), precise.im.into_ff64()),
        ] {
            assert!((result.re - expected.re).abs().0 < 1e-14);
            assert!((result.im - expected.im).abs().0 < 1e-14);
        }
    }
}

#[test]
fn native_tanh_remains_finite_at_uv_scales_and_preserves_derivatives() {
    use symbolica::{
        domains::dual::HyperDual,
        prelude::{Complex as SymComplex, Dualizer},
    };

    crate::initialisation::test_initialise().unwrap();
    let argument = Atom::var(symbol!("native_tanh_argument"));
    // Native tanh must preserve real and complex values at UV scales, with
    // analytic derivatives at x=0.
    let evaluator = symbolica::transcendental::tanh()
        .call_args([argument.clone()])
        .evaluator(std::slice::from_ref(&argument))
        .build()
        .unwrap();
    let mut real = evaluator
        .clone()
        .map_coeff(&|coefficient| coefficient.re.to_f64());
    let dualizer = Dualizer::new(
        HyperDual::<SymComplex<Rational>>::new(
            crate::utils::hyperdual_utils::simple_n_deriv_shape(2),
        ),
        vec![],
    );
    let mut dual = evaluator
        .clone()
        .vectorize(&dualizer)
        .unwrap()
        .map_coeff(&|coefficient| coefficient.re.to_f64());
    for value in [-1e12_f64, -1000., -1., 0., 1., 1000., 1e12] {
        let tanh = value.tanh();
        let first_derivative = 1. - tanh * tanh;
        let expected = [tanh, first_derivative, -tanh * first_derivative];
        assert!((real.evaluate_single(&[value]) - tanh).abs() < 1e-14);
        let mut actual = [0.; 3];
        dual.evaluate(&[value, 1., 0.], &mut actual);
        for (actual, expected) in actual.into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-14, "tanh at {value}");
        }
    }

    let mut complex = evaluator.map_coeff(&|coefficient| {
        SymComplex::new(coefficient.re.to_f64(), coefficient.im.to_f64())
    });
    for value in [-1e12_f64, 1e12] {
        let actual = complex.evaluate_single(&[SymComplex::new(value, 1.)]);
        assert_eq!(actual.re, value.signum());
        assert_eq!(actual.im, 0.);
    }
}
