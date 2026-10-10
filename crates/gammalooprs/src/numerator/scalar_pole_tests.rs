use spenso::algebra::complex::Complex;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    evaluate::{FunctionMap, OptimizationSettings},
    function, parse_lit,
    poly::series::SeriesDepth,
    symbol,
};

use crate::{
    cff::expression::OrientationID,
    initialisation::test_initialise,
    integrands::process::{GenericEvaluator, GenericEvaluatorFloat},
    processes::EvaluatorSettings,
    utils::{F, GS},
};

use super::NumeratorAtomExt;

#[test]
fn scalar_poles_preserve_floating_coefficients_in_arithmetic() {
    let x = parse_lit!(pole_domain::x);
    let removable = (x.pow(2) - 1) / (&x - 1);
    let family = function!(symbol!("pole_domain::N"; Scalar), &x);
    for coefficient in [Atom::num(0.25f64), Atom::i() * Atom::num(0.25f64)] {
        for scalar in [
            &coefficient * &removable,
            (&x + &coefficient) * &removable,
            (&x + &coefficient).pow(-1) * &removable,
        ] {
            let source = &family * scalar;
            assert_eq!(
                source.cancel_scalar_poles(std::slice::from_ref(&family)),
                source,
            );
        }
    }
    // Function arguments are opaque to rational conversion. A float inside
    // one must not disable cancellation of the surrounding exact coefficient.
    let opaque = function!(symbol!("pole_domain::f"; Scalar), Atom::num(0.25f64));
    let source = &opaque * &removable;
    assert_eq!(source.cancel_scalar_poles(&[]), opaque * (&x + 1));
    for coefficient in [Atom::one() / 4, Atom::i() / 4] {
        let source = &coefficient * &removable;
        assert_eq!(source.cancel_scalar_poles(&[]), coefficient * (&x + 1));
    }
}

// These oracles normalize only scalar coefficients of literal retained calls;
// no graph numerator definition is introduced or expanded.
fn assert_scalar_coefficients_equal(actual: &Atom, expected: &Atom, keys: &[Atom]) {
    for (_, coefficient) in (actual - expected)
        .coefficient_list_exact(keys)
        .expect("the scalar oracle must have supported factorized coefficients")
    {
        assert!(keys.iter().all(|key| !coefficient.contains(key)));
        assert!(coefficient.together().cancel().is_zero());
    }
}

#[test]
fn scalar_poles_cancel_retained_family_laurent_coefficients() {
    let t_symbol = symbol!("pole_jet::t");
    let t = Atom::var(t_symbol);
    let s = parse_lit!(pole_jet::s);
    let m = parse_lit!(pole_jet::m);
    let a = parse_lit!(pole_jet::a);
    let q = parse_lit!(pole_jet::q);
    let keys = [0, 1].map(|order| function!(symbol!("pole_jet::N"; Scalar), order, &q));
    let [n0, n1] = &keys;
    let source =
        (n0 + n1 * &t) / (&s * &t) * ((&m + &s + &a * &t).pow(-1) - (&m + &a * &t).pow(-1));
    let series = source
        .series(t_symbol, Atom::Zero, SeriesDepth::absolute(0))
        .unwrap();
    let mut result = Atom::Zero;
    for (power, coefficient) in series.terms() {
        assert!(coefficient.contains(s.pow(-1)));
        let reduced = coefficient.cancel_scalar_poles(&keys);
        assert!(!reduced.contains(s.pow(-1)));
        assert_eq!(reduced.cancel_scalar_poles(&keys), reduced);
        assert_scalar_coefficients_equal(&reduced, coefficient, &keys);
        result += reduced * t.pow(Atom::num(power));
    }
    let expected = -n0 / (&t * &m * (&m + &s)) - n1 / (&m * (&m + &s))
        + n0 * &a * (&s + 2 * &m) / (m.pow(2) * (&m + &s).pow(2));
    assert_scalar_coefficients_equal(&result, &expected, &keys);
    assert!(result.contains(t.pow(-1)));
    assert!(
        result.contains(n1),
        "the numerator's Taylor derivative was lost"
    );
}

#[test]
fn scalar_poles_preserve_retained_products_and_distinct_arguments() {
    let s = parse_lit!(pole_product::s);
    let m = parse_lit!(pole_product::m);
    let q = parse_lit!(pole_product::q);
    let p = parse_lit!(pole_product::p);
    let family = symbol!("pole_product::N"; Scalar);
    let n = function!(family, &q);
    let shifted = function!(family, &q + &p);
    let other = function!(symbol!("pole_product::M"; Scalar), &p);
    let keys = [n.clone(), shifted.clone(), other.clone()];
    let scalar = (m.pow(-1) - (&m + &s).pow(-1)) / &s;
    let source = &n * other.pow(2) * &scalar;
    let result = source.cancel_scalar_poles(&keys);
    let expected = &n * other.pow(2) / (&m * (&m + &s));
    assert!(!result.contains(s.pow(-1)));
    assert!(result.contains(other.pow(2)));
    assert_scalar_coefficients_equal(&result, &expected, &keys);

    let distinct = (&n / &m - &shifted / (&m + &s)) / &s;
    let result = distinct.cancel_scalar_poles(&keys);
    assert!(result.contains(s.pow(-1)));
    assert!(result.contains(&n));
    assert!(result.contains(&shifted));
    assert_scalar_coefficients_equal(&result, &distinct, &keys);
}

#[test]
fn scalar_poles_preserve_momentum_dependence_for_an_outer_taylor_operation() {
    let t_symbol = symbol!("pole_momentum::t");
    let t = Atom::var(t_symbol);
    let q = parse_lit!(pole_momentum::q);
    let s = parse_lit!(pole_momentum::s);
    let m = parse_lit!(pole_momentum::m);
    let momentum = &q + &t;
    let numerator = function!(symbol!("pole_momentum::N"; Scalar), &momentum);
    let keys = std::slice::from_ref(&numerator);
    let source = &numerator / &s * (m.pow(-1) - (&m + &s).pow(-1));
    let result = source.cancel_scalar_poles(keys);
    assert!(result.contains(&numerator));
    assert!(!result.contains(s.pow(-1)));

    // A local scalar toy definition supplies a nonzero outer derivative. The
    // cancellation pass itself never sees or substitutes this definition.
    let specialized = result
        .replace(numerator.to_pattern())
        .with(momentum.pow(2).to_pattern());
    let actual = specialized
        .series(t_symbol, Atom::Zero, SeriesDepth::absolute(1))
        .unwrap()
        .to_atom();
    let expected = (q.pow(2) + 2 * &q * &t) / (&m * (&m + &s));
    assert_scalar_coefficients_equal(&actual, &expected, &[]);
}

#[test]
fn scalar_poles_preserve_owned_energies_and_their_outer_derivatives() {
    test_initialise().unwrap();
    let t_symbol = symbol!("pole_energy::t");
    let t = Atom::var(t_symbol);
    let q = parse_lit!(pole_energy::q);
    let s = parse_lit!(pole_energy::s);
    let mass = parse_lit!(pole_energy::mass);
    let invariant = mass.pow(2) + (&q + &t).pow(2);
    let energy = function!(GS.on_shell_energy, 3, &invariant);
    let other_owner = function!(GS.on_shell_energy, 7, &invariant);
    let source = (energy.pow(-1) - (&energy + &s).pow(-1)) / &s;
    let result = source.cancel_scalar_poles(&[]);
    assert!(!result.contains(s.pow(-1)));
    assert!(result.contains(&energy));
    let actual_jet = result
        .series(t_symbol, Atom::Zero, SeriesDepth::absolute(1))
        .unwrap()
        .to_atom();
    let original_jet = source
        .series(t_symbol, Atom::Zero, SeriesDepth::absolute(1))
        .unwrap()
        .to_atom();
    assert_scalar_coefficients_equal(&actual_jet, &original_jet, &[]);
    assert!(!actual_jet.derivative(t_symbol).is_zero());

    for owned_difference in [
        energy.pow(2) - &invariant,
        energy.pow(-1) - other_owner.pow(-1),
    ] {
        let source = owned_difference / &s;
        let result = source.cancel_scalar_poles(&[]);
        assert!(!result.is_zero(), "an outer-U energy owner was erased");
        assert!(result.contains(&energy));
        assert!(result.contains(s.pow(-1)));
        assert_scalar_coefficients_equal(&result, &source, &[]);
    }
}

#[test]
fn scalar_poles_do_not_cross_guards_or_evaluate_inactive_singularities() {
    test_initialise().unwrap();
    let q = parse_lit!(pole_guard::q);
    let active = parse_lit!(pole_guard::active);
    let other = parse_lit!(pole_guard::other);
    let source = Symbol::IF.call_args([active.clone(), q.pow(-1), Atom::Zero])
        + Symbol::IF.call_args([other.clone(), q.pow(-2), Atom::Zero]);
    let result = source.cancel_scalar_poles(&[]);
    assert_eq!(result, source);
    for guard in [
        OrientationID(4).atom(),
        function!(GS.theta, &active),
        function!(GS.orientation_delta, &active, &other),
    ] {
        let guarded = &guard / &q - guard / (&q + 1);
        assert_eq!(guarded.cancel_scalar_poles(&[]), guarded);
    }
    let mut evaluator = GenericEvaluator::new_from_raw_params(
        [result],
        &[active, other, q],
        &FunctionMap::default(),
        vec![],
        OptimizationSettings::default(),
        None,
        &EvaluatorSettings::default(),
    )
    .unwrap();
    for (inputs, expected) in [
        ([0.0, 0.0, 0.0], 0.0),
        ([1.0, 0.0, 2.0], 0.5),
        ([0.0, 1.0, 2.0], 0.25),
    ] {
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(
            &inputs.map(|value| Complex::new_re(F(value))),
        );
        assert_eq!(actual, Complex::new_re(F(expected)));
    }
}

#[test]
fn scalar_poles_bound_common_denominator_expansion() {
    let x = parse_lit!(pole_budget::x);
    let y = parse_lit!(pole_budget::y);
    let s = parse_lit!(pole_budget::s);
    let family = function!(symbol!("pole_budget::N"; Scalar), &x);
    let source = &family * ((&x + &y).pow(100_000) - 1) / &s;
    assert_eq!(
        source.cancel_scalar_poles(std::slice::from_ref(&family)),
        source
    );
}

#[test]
fn scalar_poles_keep_products_of_numerator_sums_and_hidden_families_intact() {
    let q = parse_lit!(pole_opaque::q);
    let s = parse_lit!(pole_opaque::s);
    let m = parse_lit!(pole_opaque::m);
    let family = symbol!("pole_opaque::N"; Scalar);
    let keys = [0, 1, 2, 3].map(|index| function!(family, index, &q));
    let [a, b, c, d] = &keys;
    let removable = (m.pow(-1) - (&m + &s).pow(-1)) / &s;
    for numerator in [
        (a + b) * (c + d),
        (a + b).pow(2),
        a * (b + c),
        function!(symbol!("pole_opaque::outer"; Scalar), a),
        a.pow(-1),
    ] {
        let source = &numerator * &removable;
        assert_eq!(
            source.cancel_scalar_poles(&keys),
            source,
            "unsupported numerator structure was redistributed: {numerator}"
        );
    }
}

#[test]
fn scalar_poles_preserve_guards_inside_retained_arguments() {
    let q = parse_lit!(pole_nested_guard::q);
    let s = parse_lit!(pole_nested_guard::s);
    let m = parse_lit!(pole_nested_guard::m);
    let guarded_argument = Symbol::IF.call_args([q.clone(), q.pow(-1), Atom::Zero]);
    let family = function!(symbol!("pole_nested_guard::N"; Scalar), guarded_argument);
    let source = &family / &s * (m.pow(-1) - (&m + &s).pow(-1));
    assert_eq!(
        source.cancel_scalar_poles(std::slice::from_ref(&family)),
        source
    );
}

#[test]
fn exact_coefficients_preserve_nonzero_rational_differences() {
    let x = parse_lit!(exact_coefficient_small::x);
    let q = parse_lit!(exact_coefficient_small::q);
    let delta = Atom::num(10).pow(-30);
    let small = (&x + 1).pow(-1) - (&x + 1 + &delta).pow(-1);
    let n = function!(symbol!("exact_coefficient_small::N"; Scalar), &q);
    let m = function!(symbol!("exact_coefficient_small::M"; Scalar), &q);
    let keys = [n.clone(), m.clone()];
    for numerator in [n.clone(), &n * m.pow(2)] {
        let source = &numerator * &small;
        let coefficients = source.coefficient_list_exact(&keys).unwrap();
        assert_eq!(coefficients.len(), 1, "a nonzero coefficient was discarded");
        let (key, coefficient) = &coefficients[0];
        assert_eq!(key, &numerator);
        assert!(!coefficient.is_zero());
        assert!(keys.iter().all(|key| !coefficient.contains(key)));
        // An exact rational sample distinguishes this from the zero which a
        // statistical coefficient test can assign to the uncombined input.
        let at_one = coefficient.replace(x.to_pattern()).with(Atom::num(1));
        assert_eq!(at_one, &delta / (2 * (2 + &delta)));
    }
}

#[test]
fn exact_coefficients_retain_factorized_zero_until_scalar_cancellation() {
    let t = parse_lit!(exact_coefficient_zero::t);
    let x = parse_lit!(exact_coefficient_zero::x);
    let y = parse_lit!(exact_coefficient_zero::y);
    let q = parse_lit!(exact_coefficient_zero::q);
    let n = function!(symbol!("exact_coefficient_zero::N"; Scalar), &q);
    let m = function!(symbol!("exact_coefficient_zero::M"; Scalar), &q);
    let zero: Atom = (&x + &y).pow(2) - x.pow(2) - 2 * &x * &y - y.pow(2);
    assert!(
        !zero.is_zero(),
        "the fixture must retain its factored syntax"
    );
    let source = t.pow(-2) * (&n + &m * &t) * &zero;
    let coefficients = source
        .coefficient_list_exact(std::slice::from_ref(&t))
        .unwrap()
        .into_iter()
        .collect::<std::collections::BTreeMap<_, _>>();
    assert_eq!(coefficients.len(), 2);
    assert_eq!(coefficients.get(&t.pow(-2)), Some(&(&n * &zero)));
    assert_eq!(coefficients.get(&t.pow(-1)), Some(&(&m * &zero)));
    assert!(
        coefficients
            .values()
            .all(|coefficient| !coefficient.is_zero())
    );
    for coefficient in coefficients.values() {
        // Add a genuine scalar denominator so the cancellation pass, rather
        // than statistical coefficient extraction, owns removal of the zero.
        let rational = coefficient / (&x + 1);
        assert!(
            rational
                .cancel_scalar_poles(&[n.clone(), m.clone()])
                .is_zero()
        );
    }
}

#[test]
fn exact_coefficients_preserve_signed_powers_and_family_arguments() {
    let t = parse_lit!(exact_coefficient_laurent::t);
    let q = parse_lit!(exact_coefficient_laurent::q);
    let p = parse_lit!(exact_coefficient_laurent::p);
    let family = symbol!("exact_coefficient_laurent::N"; Scalar);
    let n = function!(family, &q);
    let shifted = function!(family, &q + &p);
    let spectator = (parse_lit!(exact_coefficient_laurent::x)
        + parse_lit!(exact_coefficient_laurent::y))
    .pow(3);
    let source = &spectator * (n.pow(2) * &shifted / t.pow(3) + &n * &shifted * t.pow(2));
    let coefficients = source
        .coefficient_list_exact(std::slice::from_ref(&t))
        .unwrap()
        .into_iter()
        .collect::<std::collections::BTreeMap<_, _>>();
    assert_eq!(coefficients.len(), 2);
    assert_eq!(
        coefficients.get(&t.pow(-3)),
        Some(&(&spectator * n.pow(2) * &shifted))
    );
    assert_eq!(
        coefficients.get(&t.pow(2)),
        Some(&(&spectator * &n * &shifted))
    );
    assert!(
        coefficients
            .values()
            .all(|coefficient| coefficient.contains(&spectator))
    );
    for coefficient in coefficients.values() {
        let grouped = coefficient
            .coefficient_list_exact(&[n.clone(), shifted.clone()])
            .unwrap();
        assert_eq!(grouped.len(), 1);
        assert_eq!(grouped[0].1, spectator);
    }
    assert!(
        Atom::Zero
            .coefficient_list_exact(std::slice::from_ref(&t))
            .unwrap()
            .is_empty()
    );
}

#[test]
fn exact_coefficients_refuse_non_laurent_key_dependence() {
    let t = parse_lit!(exact_coefficient_unsupported::t);
    let q = parse_lit!(exact_coefficient_unsupported::q);
    let x = parse_lit!(exact_coefficient_unsupported::x);
    for source in [
        function!(symbol!("exact_coefficient_unsupported::f"; Scalar), &t),
        (&t + &q).pow(-1),
        (&t + &q).pow(2),
        t.pow(Atom::num((1, 2))),
        t.pow(x),
    ] {
        assert!(
            source
                .coefficient_list_exact(std::slice::from_ref(&t))
                .is_none(),
            "unsupported dependence was hidden in a coefficient: {source}"
        );
    }
}
