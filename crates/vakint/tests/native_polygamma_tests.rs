//! Native special-function integration and the removal of fixed-precision tables.

use ahash::HashMap;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::float::{Complex, Float, Real},
    domains::rational::Rational,
    function,
};
use vakint::{Vakint, VakintSettings, fmft::FMFT, matad::MATAD, vakint_parse, vakint_symbol};

// Independent regression data from the former MATAD table at revision
// 24cb1ceb6e906452e29816a11737ef855b7150a3, stored at approximately 20,100
// significant digits. These are test-only values, not evaluator substitutions.
const PSI2: &str = include_str!("resources/native_polygamma/polygamma_2_at_1.txt");
const PSI4: &str = include_str!("resources/native_polygamma/polygamma_4_at_1.txt");

fn gamma_expansions(input: &Atom, settings: &VakintSettings) -> Vec<Atom> {
    let matad = MATAD::with_settings(settings.clone());
    let fmft = FMFT::with_settings(settings.clone());
    vec![
        matad.substitute_masters_directly(input.as_view()).unwrap(),
        matad.expand_matad_masters(input.as_view()).unwrap(),
        fmft.expand_masters(input.as_view()).unwrap(),
    ]
}

fn evaluate(expression: &Atom, settings: &VakintSettings) -> Complex<Float> {
    let result = Vakint::full_numerical_evaluation_without_error(
        settings,
        expression.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )
    .unwrap();
    assert_eq!(result.0.len(), 1);
    assert_eq!(result.0[0].0, 0);
    let value = result.0[0].1.clone();
    assert!(value.re.as_raw().is_finite());
    assert!(value.im.as_raw().is_zero());
    value
}

fn assert_relative_digits(actual: &Complex<Float>, expected: &Complex<Float>, digits: u32) {
    let error = (actual - expected).norm().re;
    let threshold = Float::parse(&format!("1e-{digits}"), Some(actual.re.prec())).unwrap();
    assert!(
        error < expected.norm().re * threshold,
        "native result misses {digits} digits"
    );
}

#[test]
fn gamma_expansions_use_native_constants_and_keep_unknown_orders() {
    Vakint::initialize_vakint_symbols();
    let settings = VakintSettings::default();
    let input = vakint_parse!("Gam(1,1)").unwrap();
    let expected = [
        "1",
        "0",
        "𝜋^2/12",
        "symbolica::polygamma(2,1)/6",
        "𝜋^4/160",
        "symbolica::polygamma(2,1)*𝜋^2/72+symbolica::polygamma(4,1)/120",
    ];
    for expanded in gamma_expansions(&input, &settings) {
        assert!(
            !expanded
                .get_all_symbols(true)
                .iter()
                .any(|symbol| symbol.get_namespace() == "vakint"
                    && ["PolyGamma", "EulerGamma", "Gamma"].contains(&symbol.get_stripped_name()))
        );
        assert!(expanded.contains(vakint_parse!("Oep(6,I_Gamma)").unwrap().as_view()));
        let series = expanded
            .series(
                vakint_symbol!("ep"),
                Atom::Zero.as_view(),
                Rational::from(5),
            )
            .unwrap();
        for (power, expected) in expected.iter().enumerate() {
            let coefficient = series.coefficient(Rational::from(power)).unwrap();
            assert!(
                (coefficient - vakint_parse!(*expected).unwrap())
                    .expand()
                    .is_zero()
            );
        }
        assert!(series.coefficient(series.absolute_order()).is_none());
    }
}

#[test]
fn scalar_special_functions_are_not_lorentz_vectors() {
    Vakint::initialize_vakint_symbols();
    let scalar = vakint_parse!(
        "symbolica::polygamma(2,z)^2+symbolica::polygamma(1,z)*symbolica::polygamma(2,z)\
         +g(mu,nu)*symbolica::polygamma(2,mu)*symbolica::polygamma(4,nu)"
    )
    .unwrap();
    assert_eq!(
        Vakint::convert_to_dot_notation(scalar.as_view()),
        scalar.expand()
    );
    let mixed = &scalar * vakint_parse!("k(1,101)*p(2,101)").unwrap();
    let expected = scalar * vakint_parse!("dot(k(1),p(2))").unwrap();
    assert_eq!(
        Vakint::convert_to_dot_notation(mixed.as_view()),
        expected.expand()
    );
    for (input, expected) in [
        ("k(1,mu)^2", "dot(k(1),k(1))"),
        ("k(1,mu)^4", "dot(k(1),k(1))^2"),
        ("g(mu,nu)*k(1,mu)*p(2,nu)", "dot(k(1),p(2))"),
    ] {
        assert_eq!(
            Vakint::convert_to_dot_notation(vakint_parse!(input).unwrap().as_view()),
            vakint_parse!(expected).unwrap(),
            "{input}"
        );
    }
}

#[test]
fn positive_integer_gamma_arguments_use_native_numerical_evaluation() {
    let settings = VakintSettings {
        run_time_decimal_precision: 140,
        ..Default::default()
    };
    let input = vakint_parse!("Gam(3,1)").unwrap();
    let base = gamma_expansions(&vakint_parse!("Gam(1,1)").unwrap(), &settings).remove(0);
    let expected = (base * vakint_parse!("(1+ep)*(2+ep)").unwrap())
        .series(
            vakint_symbol!("ep"),
            Atom::Zero.as_view(),
            Rational::from(5),
        )
        .unwrap();
    for expanded in gamma_expansions(&input, &settings) {
        let series = expanded
            .series(
                vakint_symbol!("ep"),
                Atom::Zero.as_view(),
                Rational::from(5),
            )
            .unwrap();
        for power in 0..=5 {
            let actual = evaluate(
                &series.coefficient(Rational::from(power)).unwrap(),
                &settings,
            );
            let expected = evaluate(
                &expected.coefficient(Rational::from(power)).unwrap(),
                &settings,
            );
            assert_relative_digits(&actual, &expected, 100);
        }
    }
}

#[test]
fn native_polygamma_matches_20000_digit_reference_values() {
    let settings = VakintSettings {
        run_time_decimal_precision: 20032,
        ..Default::default()
    };
    for (order, reference) in [(2, PSI2), (4, PSI4)] {
        let expression = function!(symbolica::transcendental::polygamma(), order, 1);
        let expected = settings.complex_to_prec(reference.trim(), "0");
        assert_relative_digits(&evaluate(&expression, &settings), &expected, 20000);
    }
}

#[test]
fn gamma_expansion_coefficients_evaluate_at_20000_digits() {
    let settings = VakintSettings {
        run_time_decimal_precision: 20032,
        ..Default::default()
    };
    let psi2 = Atom::num(settings.real_to_prec(PSI2.trim()));
    let psi4 = Atom::num(settings.real_to_prec(PSI4.trim()));
    let expected = [
        (3, psi2.clone() / 6),
        (5, psi2 * Atom::var(Symbol::PI).pow(2) / 72 + psi4 / 120),
    ];
    for expanded in gamma_expansions(&vakint_parse!("Gam(1,1)").unwrap(), &settings) {
        let series = expanded
            .series(
                vakint_symbol!("ep"),
                Atom::Zero.as_view(),
                Rational::from(5),
            )
            .unwrap();
        for (power, expected) in &expected {
            let actual = evaluate(
                &series.coefficient(Rational::from(*power)).unwrap(),
                &settings,
            );
            assert_relative_digits(&actual, &evaluate(expected, &settings), 20000);
        }
    }
}
