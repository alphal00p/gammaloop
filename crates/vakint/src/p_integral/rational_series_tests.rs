use super::{Laurent, rational_series::RationalSeries};
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::domains::rational::Rational;

#[test]
fn exact_recurrence_admits_captured_four_loop_background_coefficients() {
    compare_captured_coefficients(
        include_bytes!("../../tests/resources/p_integral/qcd-l4-graph17-series.json.gz"),
        6,
    );
}

#[test]
fn exact_recurrence_admits_large_repeated_four_loop_coefficients() {
    compare_captured_coefficients(
        include_bytes!("../../tests/resources/p_integral/qcd-l4-graph41-series.json.gz"),
        1,
    );
}

fn compare_captured_coefficients(compressed: &[u8], expected_cases: usize) {
    use std::io::Read;
    crate::Vakint::initialize_vakint_symbols();
    let mut fixture = String::new();
    flate2::read::GzDecoder::new(compressed)
        .read_to_string(&mut fixture)
        .unwrap();
    let fixture: serde_json::Value = serde_json::from_str(&fixture).unwrap();
    let cases = fixture["cases"].as_array().unwrap();
    assert_eq!(cases.len(), expected_cases);
    for case in cases {
        let input = Atom::parse(
            case["input"].as_str().unwrap(),
            "vakint",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let epsilon = Atom::parse(
            case["epsilon"].as_str().unwrap(),
            "vakint",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
        .get_symbol()
        .unwrap();
        // These are complete captured coefficients, including production
        // rational prefactors, complex phases and open Lorentz tensors.
        // Compare the full requested Laurent polynomial with the unchanged
        // generic expansion; do not replace or specialize spectator tensors.
        compare(&input, epsilon, case["order"].as_i64().unwrap());
    }
}

fn compare(input: &Atom, epsilon: Symbol, order: i64) {
    compare_expected_valuation(
        input,
        epsilon,
        order,
        super::valuation(input, epsilon).unwrap(),
    );
}

fn compare_expected_valuation(input: &Atom, epsilon: Symbol, order: i64, expected: Option<i64>) {
    let prepared = RationalSeries::new(input, epsilon).expect("fixture must use exact recurrence");
    assert_eq!(
        prepared.valuation(),
        expected,
        "incorrect exact valuation for {input}"
    );
    let actual = prepared
        .series(order)
        .expect("fixture must remain inside exact work bounds")
        .to_atom(epsilon);
    let reference = input
        .series(epsilon, Atom::Zero.as_view(), Rational::from(order))
        .unwrap()
        .to_atom();
    assert!(
        (&actual - &reference).expand().together().is_zero(),
        "{input}: {actual} != {reference}"
    );
}

#[test]
fn exact_recurrence_admits_frozen_four_loop_forcer_coefficients_and_complex_tensor_factors() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let battery: serde_json::Value =
        serde_json::from_str(include_str!("../../tests/fixtures/forcer/p4l_battery.json")).unwrap();
    let tensor = Atom::parse(
        "UFO::{}::gs_minus^4*vakint::{symmetric}::g(gammalooprs::{spenso::index}::hedge(11,1),gammalooprs::{spenso::index}::hedge(20,1))",
        "vakint", Default::default()).unwrap();
    let mut coefficients = 0;
    for case in battery["cases"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|case| matches!(case["name"].as_str(), Some("la4_top" | "la4_dot")))
    {
        for output in case["outputs"].as_array().unwrap() {
            let reduced =
                crate::p_integral::forcer::parse_forcer_output(output.as_str().unwrap(), epsilon)
                    .unwrap();
            let mut masters = Vec::new();
            reduced.visitor(&mut |view| {
                if let symbolica::atom::AtomView::Fun(function) = view
                    && [
                        crate::p_integral::pg_symbol(),
                        crate::p_integral::p_master_symbol(),
                        crate::p_integral::p_gamma_symbol(),
                    ]
                    .contains(&function.get_symbol())
                {
                    masters.push(view.to_owned());
                    return false;
                }
                true
            });
            masters.sort();
            masters.dedup();
            for (_, coefficient) in reduced.coefficient_list::<i16>(&masters) {
                for prefactor in [Atom::one(), Atom::i(), &tensor * Atom::i()] {
                    let input = &coefficient * &prefactor;
                    for order in [-1, 2] {
                        compare(&input, epsilon, order);
                        compare(&input.expand(), epsilon, order);
                    }
                }
                coefficients += 1;
            }
        }
    }
    assert!(
        coefficients >= 4,
        "real master coefficients must be exercised"
    );
}

#[test]
fn exact_recurrence_respects_shift_cancellation_and_scalar_denominator() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    // Explicit mathematical orders are independent of either implementation.
    // The legacy relative-depth series may retain an algebraically zero first
    // coefficient and therefore report only a conservative structural bound:
    // the first input reduces to epsilon^-2/(1-2*epsilon)^3, not order -5.
    // Complete coefficients are still compared with the original, unchanged
    // generic series of each uncanonicalized input in compare_expected_valuation.
    for (input, expected) in [
        ("((a+1)^2-a^2-2*a-1+ε^3)/(ε^5*(1-2*ε)^3)", Some(-2)),
        ("(a+ε*b+ε^2*c)/((a+b)*(1-ε+2*ε^2))", Some(0)),
        ("1/(a-1)*(a^2-1-(a-1)*(a+1))/ε^4", None),
        ("(ε^7+2*ε^8)/(ε^2*(1-ε))", Some(5)),
    ] {
        let input = vk_parse!(input).unwrap();
        for order in [-3, 0, 3, 6] {
            compare_expected_valuation(&input, epsilon, order, expected);
        }
    }
}

#[test]
fn exact_recurrence_rejects_unsupported_dependencies_and_work() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    for input in [
        "1/(a+ε)",
        "sin(ε)",
        "f(ε)/(1-ε)",
        "1/(1-ε)^65",
        "ε^(1/2)",
        "a^(4294967296/3)/(1-ε)",
        "(a+b)^(100000*x)/(1-ε)",
    ] {
        assert!(
            RationalSeries::new(&vk_parse!(input).unwrap(), epsilon).is_none(),
            "{input}"
        );
    }
    let floating = Atom::num(0.25_f64) / (Atom::one() - Atom::var(epsilon));
    assert!(RationalSeries::new(&floating, epsilon).is_none());
    let input = vk_parse!("1/(1-ε)").unwrap();
    let prepared = RationalSeries::new(&input, epsilon).unwrap();
    assert!(prepared.series(100).is_none());
}

#[test]
fn exact_recurrence_preserves_existing_gamma_split_and_truncated_products() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let input = crate::Vakint::normalize_gamma_functions(
        vk_parse!("Gamma(1+ε)^2*Gamma(1-2*ε)/Gamma(2-ε)/(ε^3*(1-2*ε))")
            .unwrap()
            .as_view(),
    );
    let split = super::split(&input, epsilon);
    let prepared = RationalSeries::new(&split.rational, epsilon).unwrap();
    let valuation = prepared.valuation().unwrap();
    let order = 1;
    let exponential = Laurent::exp(&split.exponent(order - valuation), order - valuation);
    let actual = prepared
        .series(order)
        .unwrap()
        .mul(&exponential, order)
        .to_atom(epsilon);
    let reference = Laurent::series(&split.rational, epsilon, order)
        .unwrap()
        .mul(&exponential, order)
        .to_atom(epsilon);
    assert!((actual - reference).expand().together().is_zero());
}

#[test]
fn exact_recurrence_admits_pi_without_hiding_epsilon_dependence() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    for (input, expected) in [
        ("𝜋*ε/(1-ε)", Some(1)),
        ("𝑖/(16*𝜋^2*ε*(1-2*ε))", Some(-1)),
        ("(1+𝜋*ε+𝜋^2*ε^2)/(𝜋^4*ε^2*(1-3*ε+ε^2))", Some(-2)),
        ("((𝜋+1)^2-𝜋^2-2*𝜋-1+ε^3)/(𝜋^2*ε^5*(1-2*ε))", Some(-2)),
    ] {
        let input = vk_parse!(input).unwrap();
        for order in [-1, 0, 2, 5] {
            compare_expected_valuation(&input, epsilon, order, expected);
        }
    }
    // Independent expected coefficients ensure PI*epsilon is not detached
    // as a constant factor or promoted to an opaque variable hiding epsilon.
    let input = vk_parse!("𝜋*ε/(1-ε)").unwrap();
    let actual = RationalSeries::new(&input, epsilon)
        .unwrap()
        .series(3)
        .unwrap()
        .to_atom(epsilon);
    let expected = vk_parse!("𝜋*ε+𝜋*ε^2+𝜋*ε^3").unwrap();
    assert!((actual - expected).expand().is_zero());
}

#[test]
fn exact_recurrence_admits_the_production_msbar_prefactor() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let settings = crate::VakintSettings {
        integral_normalization_factor: crate::LoopNormalizationFactor::MSbar,
        ..Default::default()
    };
    // Construct the same normalization and measure as finalize(), rather
    // than testing only the earlier pure rational FORCER table coefficients.
    let normalization = settings
        .get_integral_normalization_factor_atom()
        .unwrap()
        .replace(crate::symbols::S.n_loops.to_pattern())
        .with(Atom::num(1).to_pattern());
    let measure =
        vk_parse!("(𝑖*𝜋^(2-ε)*Gamma(1+ε)*Gamma(1-ε)^2/Gamma(2-2*ε))*exp(-logQ2mu-log_mu_sq)^ε")
            .unwrap();
    let prefactor = crate::Vakint::normalize_gamma_functions((normalization * measure).as_view());
    let split = super::split(&prefactor, epsilon);
    assert!(split.rational.contains_symbol(Symbol::PI));
    for order in [0, 2, 4] {
        compare(&split.rational, epsilon, order);
    }
}

#[test]
fn analytic_spectators_preserve_epsilon_and_exactness_guards() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    for input in [
        "1/(𝜋+ε)",
        "1/(symbolica::zeta(3)+ε)",
        "1/(sin(x)+ε)",
        "sin(ε)/(1-ε)",
        "symbolica::zeta(3+ε)/(1-ε)",
        "f(sin(ε))/(1-ε)",
        "f(𝜋*ε)/(1-ε)",
        "1/((sin(x)^2+cos(x)^2-1)*(1-ε))",
    ] {
        assert!(
            RationalSeries::new(&vk_parse!(input).unwrap(), epsilon).is_none(),
            "{input}"
        );
    }
    let denominator = Atom::one() - Atom::var(epsilon);
    let floating = Atom::num(0.25_f64) / denominator;
    assert!(RationalSeries::new(&floating, epsilon).is_none());
}

#[test]
fn exact_recurrence_admits_epsilon_free_analytic_spectators() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let e = Atom::var(epsilon);
    for spectator in [
        symbolica::transcendental::zeta().call(3),
        Atom::var(symbolica::transcendental::euler_gamma()),
        Atom::var(Symbol::E),
        vk_parse!("sin(x)+cos(y)").unwrap(),
        vk_parse!("log(x)+symbolica::gamma(y)").unwrap(),
    ] {
        let input = &spectator * (Atom::one() + &e) / (&e * (Atom::one() - &e));
        for order in [-2, -1, 0, 3] {
            compare_expected_valuation(&input, epsilon, order, Some(-1));
        }
    }
}

#[test]
fn analytic_spectator_relations_preserve_semantic_valuation() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let identity = vk_parse!("sin(x)^2+cos(x)^2-1").unwrap();
    assert!(!identity.is_zero());
    assert!(identity.check_zero().is_true());
    let e = Atom::var(epsilon);
    let denominator = e.pow(5) * (Atom::one() - &e);
    for (numerator, expected_valuation) in
        [(identity.clone(), None), (&identity + e.pow(3), Some(-2))]
    {
        let input = numerator / &denominator;
        let prepared = RationalSeries::new(&input, epsilon).expect("analytic relation is admitted");
        assert_eq!(prepared.valuation(), expected_valuation);
        for order in [-3, -1, 0, 2] {
            let actual = prepared.series(order).unwrap().to_atom(epsilon);
            let reference = input
                .series(epsilon, Atom::Zero.as_view(), Rational::from(order))
                .unwrap()
                .to_atom();
            assert!((&actual - reference).check_zero().is_true());
            let expected = if expected_valuation.is_some() && order >= -2 {
                Atom::add_many((-2..=order).map(|power| e.pow(power)))
            } else {
                Atom::Zero
            };
            assert!((actual - expected).check_zero().is_true());
        }
    }
}

#[test]
fn semantic_zero_denominators_are_checked_before_formal_factor_cancellation() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let zero = vk_parse!("sin(x)^2+cos(x)^2-1").unwrap();
    assert!(!zero.is_zero() && zero.check_zero().is_true());
    // In independent s=sin(x), c=cos(x), the expanded numerator contains
    // (s^2+c^2-1)(s^2+c^2+1). A polynomial GCD must not cancel the
    // identically-zero denominator and turn this division into a finite value.
    let input =
        vk_parse!("(sin(x)^4+2*sin(x)^2*cos(x)^2+cos(x)^4-1)/((sin(x)^2+cos(x)^2-1)*(1-ε))")
            .unwrap();
    assert!(RationalSeries::new(&input, epsilon).is_none());
}

#[test]
fn nonzero_analytic_denominator_factors_may_still_cancel() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    // Here the common factor sin(x)^2+cos(x)^2+1 equals two, so cancellation
    // is valid and the result is exactly 2/(1-epsilon).
    let input = vk_parse!(
        "(sin(x)^4+2*sin(x)^2*cos(x)^2+cos(x)^4+2*sin(x)^2+2*cos(x)^2+1)/((sin(x)^2+cos(x)^2+1)*(1-ε))"
    ).unwrap();
    let actual = RationalSeries::new(&input, epsilon)
        .unwrap()
        .series(3)
        .unwrap()
        .to_atom(epsilon);
    let expected = vk_parse!("2+2*ε+2*ε^2+2*ε^3").unwrap();
    assert!((actual - expected).check_zero().is_true());
    compare(&input, epsilon, 3);
}
