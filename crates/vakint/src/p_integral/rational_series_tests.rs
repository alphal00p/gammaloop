use super::{
    Laurent,
    rational_series::{RationalPart, RationalSeries},
};
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

/// The split multiplies its input out, which cancels a denominator that also
/// divides every numerator term. A 0/0 input (the single path refuses its
/// zero denominator) must stay refused there, not come back as the series of
/// the cancelled form; the same shape with a nonzero factor is still split.
#[test]
fn a_semantic_zero_denominator_is_refused_before_the_split_expands_it_away() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let quotient = |factor: &str| {
        vk_parse!(&format!("({factor}*a+{factor}*b*ε)/({factor})")).unwrap()
    };
    for zero in ["sin(𝜋)", "(sin(x)^2+cos(x)^2-1)"] {
        let factor = vk_parse!(zero).unwrap();
        assert!(!factor.is_zero() && factor.check_zero().is_true(), "{zero}");
        let input = quotient(zero);
        assert!(RationalPart::new(&input, epsilon).is_none(), "{zero}");
        assert!(RationalSeries::new(&input, epsilon).is_none(), "{zero}");
    }
    let input = quotient("sin(x)");
    let actual = RationalSeries::new(&input, epsilon)
        .expect("a nonzero common factor is admitted")
        .series(2)
        .unwrap()
        .to_atom(epsilon);
    assert!((actual - vk_parse!("a+b*ε").unwrap()).check_zero().is_true());
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

/// More spectator monomials × numerator degree than one exact part admits
/// (Q L4 entry 1621, two-point row at depth 7: 1927 monomials × degree 17):
/// the input is split per spectator monomial and the parts' series summed.
/// 1100 monomials `a^i b^j` over 30 distinct denominators `1 + kε`
/// (support 1100 × 31 > 32768), against the closed form
/// `Σ mono/(1+kε) = Σ mono Σ_n (−kε)^n`; and the phase and an ε-free
/// spectator denominator (`𝜋^-2`) carried as in the declined input.
#[test]
fn many_spectator_monomials_are_expanded_per_monomial() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let e = Atom::var(epsilon);
    let (a, b) = (vk_parse!("a").unwrap(), vk_parse!("b").unwrap());
    let pi2 = vk_parse!("𝜋^(-2)").unwrap();
    let terms = (0..1100i64)
        .map(|i| {
            let k = i % 30 + 1;
            let monomial = a.pow(Atom::num(i % 40)) * b.pow(Atom::num(i / 40));
            (monomial, k)
        })
        .collect::<Vec<_>>();
    let sum = Atom::add_many(
        terms
            .iter()
            .map(|(monomial, k)| monomial * &pi2 / (Atom::one() + Atom::num(*k) * &e)),
    );
    let phase = vk_parse!("-𝑖/4096").unwrap();
    let input = &phase * &sum;
    assert!(
        RationalPart::new(&input, epsilon).is_none(),
        "one part admits it: the test does not exercise the split"
    );
    let prepared = RationalSeries::new(&input, epsilon).expect("split per spectator monomial");
    assert_eq!(prepared.valuation(), Some(0));
    let order = 3;
    let actual = prepared.series(order).unwrap().to_atom(epsilon);
    let expected = Atom::add_many(terms.iter().map(|(monomial, k)| {
        let mut series = Atom::Zero;
        for n in 0..=order {
            series += (Atom::num(-*k) * &e).pow(Atom::num(n));
        }
        &phase * monomial * &pi2 * series
    }));
    assert!((actual - expected).expand().is_zero());
}

// Depth-dependent declines of the exact plan (Q L4 entry 1621 at depth 7).
//
// The declined input of that entry is not readable here, so the factors
// below are production-SHAPED stand-ins: the rational factor (prefactor x
// reduction coefficient) of one master monomial of a four-loop R* remainder
// is a sum over colour x tensor structures of rational functions of
// ε = (4-d)/2 with IBP denominators, times a lower-loop factor that the
// caller truncates at its own depth. Both the number of spectator monomials
// and the numerator degree therefore grow with the depth, while the master
// monomial and the requested order of the term only set the work bound.

/// Structures of the stand-in for the Q L4 1621 factor: one whole exact part
/// admits it at depths 4 to 6 (894 to 1,639 summands) and declines it from
/// depth 7 on (2,086 summands, 1,120 spectator monomials), as reported for
/// that entry (passes at depth 4, declines at depth 7).
const STAND_IN_STRUCTURES: usize = 80;

/// The master monomial of Q L4 entry 1621 (two-point row, depth 7).
const Q_L4_1621_MONOMIAL: [(i64, i64); 5] = [(3, 1), (-2, 3), (1, -3), (-1, -2), (-4, -1)];

/// `Π PGamma(m)^p` (Gamma tags).
fn gamma_tags(tags: &[(i64, i64)]) -> Atom {
    Atom::mul_many(
        tags.iter()
            .map(|(m, p)| super::super::p_gamma_symbol().call(*m).pow(Atom::num(*p)))
            .collect::<Vec<_>>()
            .iter(),
    )
}

/// ε-linear factors `a + bε` of reductions to d = 4 - 2ε (d-3, d-1, 3d-10,
/// 2d-7, up to content), shared by the terms of one reduction.
const IBP_FACTORS: [(i64, i64); 4] = [(1, -2), (3, -2), (1, -3), (1, -4)];

/// A small deterministic generator for the stand-ins.
struct Lcg(u64);

impl Lcg {
    fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % n
    }
}

/// Products of π², ζ(3), ζ(5), ζ(7) and ζ(9) of total weight `weight`.
fn zeta_monomials(weight: i64) -> Vec<Atom> {
    fn walk(weight: i64, from: usize, out: &mut Vec<Atom>, current: Atom) {
        if weight == 0 {
            out.push(current);
            return;
        }
        for (i, w) in [2i64, 3, 5, 7, 9].into_iter().enumerate().skip(from) {
            if w <= weight {
                let factor = if w == 2 {
                    Atom::var(Symbol::PI).pow(Atom::num(2))
                } else {
                    symbolica::transcendental::zeta().call(w)
                };
                walk(weight - w, i, out, &current * factor);
            }
        }
    }
    let mut out = Vec::new();
    walk(weight, 0, &mut out, Atom::one());
    out
}

/// The stand-in rational factor at caller depth `depth` (comment above):
/// `structures` colour x tensor x coupling structures (ca^a cf^b nf^c with
/// a+b+c = 3, `transverse` or `longitudinal`, a power product of `gs2_i` and
/// `gs_minus`), each with one IBP term `(c0 + c1 ε) ε^-k / Π (a + bε)^p`,
/// times the lower-loop factor `Σ_(w=0)^(depth+1) ε^w Σ c_z z` over the zeta
/// monomials z of weight w, expanded, with the phase and normalisation of a
/// remainder.
fn r_star_like_factor(depth: i64, structures: usize, seed: u64) -> Atom {
    let epsilon = vk_symbol!("ε");
    let e = Atom::var(epsilon);
    let mut rng = Lcg(seed);
    let symbol = |name: &str| Atom::var(vk_symbol!(name));
    let (ca, cf, nf) = (symbol("ca"), symbol("cf"), symbol("nf"));
    let tensors = [symbol("transverse"), symbol("longitudinal")];
    let (g2, gm) = (symbol("gs2_i"), symbol("gs_minus"));
    let couplings = [
        Atom::one(),
        g2.clone(),
        g2.pow(Atom::num(2)),
        gm.pow(Atom::num(2)),
        &g2 * gm.pow(Atom::num(2)),
    ];
    let mut colours = Vec::new();
    for a in 0..=3i64 {
        for b in 0..=3 - a {
            colours
                .push(ca.pow(Atom::num(a)) * cf.pow(Atom::num(b)) * nf.pow(Atom::num(3 - a - b)));
        }
    }
    let mut base = Vec::new();
    for s in 0..structures {
        let structure = &colours[s % colours.len()]
            * &tensors[(s / colours.len()) % tensors.len()]
            * &couplings[(s / (colours.len() * tensors.len())) % couplings.len()];
        let pole = rng.below(4) as i64;
        let numerator =
            Atom::num(1 + rng.below(7) as i64) + Atom::num(rng.below(9) as i64 - 4) * &e;
        let mut denominator = Atom::one();
        for (a, b) in IBP_FACTORS {
            let power = rng.below(3) as i64;
            denominator *= (Atom::num(a) + Atom::num(b) * &e).pow(Atom::num(power));
        }
        base.push(structure * numerator * e.pow(Atom::num(-pole)) / denominator);
    }
    let lower = Atom::add_many(
        (0..=depth + 1)
            .flat_map(|w| {
                zeta_monomials(w)
                    .into_iter()
                    .map(|z| Atom::num(1 + rng.below(5) as i64) * z * e.pow(Atom::num(w)))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
            .iter(),
    );
    let normalisation = vk_parse!("-𝑖*𝜋^(-8)/4096").unwrap();
    (normalisation * Atom::add_many(base.iter()) * lower).expand()
}

/// The term of one master monomial as the production route computes it,
/// with the counts of exact and generic rational series.
fn master_term(
    coefficient: &Atom,
    monomial: &Atom,
    order: i64,
) -> (Option<Laurent>, super::SeriesCounts) {
    master_term_in(coefficient, monomial, vk_symbol!("ε"), order)
}

/// [`master_term`] in the caller's ε.
fn master_term_in(
    coefficient: &Atom,
    monomial: &Atom,
    epsilon: Symbol,
    order: i64,
) -> (Option<Laurent>, super::SeriesCounts) {
    let mut counts = super::SeriesCounts::default();
    let term = super::monomial_term(
        &Atom::one(),
        monomial,
        coefficient,
        &super::master_symbols(),
        epsilon,
        order,
        &mut counts,
        &super::masters::p_master_expansion,
    )
    .unwrap();
    (term, counts)
}

/// Order of an N-loop result at requester depth D: through ε^(D-N-1).
fn four_loop_order(depth: i64) -> i64 {
    depth - 5
}

/// Q L4 entry 1621 (two-point row): its master monomial with an R*-shaped
/// factor takes the exact series at depth 4, and at depth 7, where the whole
/// factor exceeds one exact part's support bound and origin/p_integral fell
/// back to the generic series (`relative_valuation`, exact plan available:
/// false).
#[test]
fn q_l4_1621_master_monomial_takes_the_exact_series_at_depth_4_and_7() {
    crate::Vakint::initialize_vakint_symbols();
    let monomial = gamma_tags(&Q_L4_1621_MONOMIAL);
    for depth in [4, 7] {
        let coefficient = r_star_like_factor(depth, STAND_IN_STRUCTURES, 1621);
        let (term, counts) = master_term(&coefficient, &monomial, four_loop_order(depth));
        assert!(term.is_some(), "depth {depth}: the term vanished");
        assert_eq!(
            (counts.exact, counts.generic),
            (1, 0),
            "depth {depth}: the rational factor took the generic series"
        );
    }
}

/// [`master_term`] through the generic series (the reference route).
fn generic_master_term(coefficient: &Atom, monomial: &Atom, order: i64) -> Option<Laurent> {
    generic_master_term_in(coefficient, monomial, vk_symbol!("ε"), order)
}

/// [`generic_master_term`] in the caller's ε.
fn generic_master_term_in(
    coefficient: &Atom,
    monomial: &Atom,
    epsilon: Symbol,
    order: i64,
) -> Option<Laurent> {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            super::FORCE_GENERIC_SERIES.with(|force| force.set(false));
        }
    }
    super::FORCE_GENERIC_SERIES.with(|force| force.set(true));
    let _reset = Reset;
    let (term, counts) = master_term_in(coefficient, monomial, epsilon, order);
    assert_eq!(
        (counts.exact, counts.generic),
        (0, 1),
        "the reference was not generic"
    );
    term
}

/// Master monomials of four-loop terms: the one of Q L4 entry 1621, other
/// Gamma-tag products, and products with `PG(x1, x2)` (whose Γ form adds
/// rational ε factors to the term's rational part).
fn four_loop_master_monomials() -> Vec<Atom> {
    let tags: [&[(i64, i64)]; 16] = [
        &Q_L4_1621_MONOMIAL,
        &[(-2, -4), (-1, 8), (1, 4)],
        &[(-2, -3), (-1, 6), (1, 3)],
        &[(3, 1), (-3, 1), (-1, -2)],
        &[(2, 1), (-2, 2), (1, -2), (-1, -1)],
        &[(4, 1), (-4, -1), (-2, 2), (1, -2)],
        &[(-3, 2), (1, -1), (-1, -1), (3, 1)],
        &[(5, -1), (-2, 3), (2, 1)],
        &[(1, -4)],
        &[(-1, 4), (-2, -2)],
        &[(3, 2), (-2, -1), (-1, -2), (1, -2)],
        &[(-4, 1), (2, -1), (-2, 1)],
        &[(2, 2), (-1, -3)],
        &[(-5, 1), (5, -1), (1, 1), (-1, -1)],
        &[(3, -1), (-2, 4), (1, -3)],
        &[(-3, -1), (-1, 3), (4, 1), (-2, -2)],
    ];
    let pg = |x1: i64, x2: i64| super::super::pg_symbol().call_args([Atom::num(x1), Atom::num(x2)]);
    let mut monomials = tags.iter().map(|t| gamma_tags(t)).collect::<Vec<_>>();
    monomials.push(pg(1, 1) * gamma_tags(&[(-2, -1), (-1, 2)]));
    monomials.push(pg(0, 1) * pg(1, 1));
    monomials.push(pg(1, 2) * gamma_tags(&Q_L4_1621_MONOMIAL));
    monomials.push(pg(2, 2) * gamma_tags(&[(1, -2)]));
    monomials
}

/// For R*-shaped factors that the whole exact part declines (the family of
/// Q L4 entry 1621), the split exact route equals the generic series through
/// the requested order on 20 master monomials (depths 7 to 9), and the
/// coefficients are identical atoms.
#[test]
fn split_exact_series_equals_the_generic_series_on_twenty_master_monomials() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let monomials = four_loop_master_monomials();
    assert_eq!(monomials.len(), 20);
    for (index, monomial) in monomials.iter().enumerate() {
        let depth = 7 + (index % 3) as i64;
        let coefficient = r_star_like_factor(depth, STAND_IN_STRUCTURES, 1621 + index as u64);
        assert!(
            RationalPart::new(&coefficient, epsilon).is_none(),
            "monomial {index}: the whole part admits the factor; the test does not exercise the split"
        );
        let order = four_loop_order(depth);
        let (exact, counts) = master_term(&coefficient, monomial, order);
        assert_eq!(
            (counts.exact, counts.generic),
            (1, 0),
            "monomial {index}: {monomial}"
        );
        let generic = generic_master_term(&coefficient, monomial, order);
        let (exact, generic) = (exact.unwrap(), generic.unwrap());
        assert_eq!(exact.valuation, generic.valuation, "monomial {index}");
        assert_eq!(
            exact.coefficients.len(),
            generic.coefficients.len(),
            "monomial {index}"
        );
        for (k, (a, b)) in exact
            .coefficients
            .iter()
            .zip(&generic.coefficients)
            .enumerate()
        {
            assert!(
                (a - b).expand().is_zero(),
                "monomial {index} ({monomial}), ε^{}: exact and generic differ",
                exact.valuation + k as i64
            );
            assert_eq!(
                a,
                b,
                "monomial {index}, ε^{}: equal values, different atoms",
                exact.valuation + k as i64
            );
        }
    }
}

/// `Σ_i m_i/(1 + k_i ε)` for `n` spectator monomials `m_i = a^(i%40) b^(i/40)`
/// over the denominators `1 + kε`, `k = 1..=denominators`, and its closed form
/// `Σ_i m_i Σ_(j<=order) (-k_i ε)^j`.
fn geometric_family(n: i64, denominators: i64, order: i64) -> (Atom, Atom) {
    let e = Atom::var(vk_symbol!("ε"));
    let (a, b) = (vk_parse!("a").unwrap(), vk_parse!("b").unwrap());
    let terms = (0..n)
        .map(|i| {
            (
                a.pow(Atom::num(i % 40)) * b.pow(Atom::num(i / 40)),
                i % denominators + 1,
            )
        })
        .collect::<Vec<_>>();
    let input = Atom::add_many(
        terms
            .iter()
            .map(|(m, k)| m / (Atom::one() + Atom::num(*k) * &e))
            .collect::<Vec<_>>()
            .iter(),
    );
    let closed = Atom::add_many(
        terms
            .iter()
            .map(|(m, k)| {
                m * Atom::add_many(
                    (0..=order)
                        .map(|j| (Atom::num(-*k) * &e).pow(Atom::num(j)))
                        .collect::<Vec<_>>()
                        .iter(),
                )
            })
            .collect::<Vec<_>>()
            .iter(),
    );
    (input, closed.expand())
}

/// One whole part admits the input but its series at a deep order exceeds
/// the work bound (`length² × support`, the depth-dependent decline with an
/// exact plan available): the series is then taken per spectator monomial,
/// equal to the closed form; a shallow order keeps the whole part's series.
#[test]
fn a_whole_part_beyond_its_work_bound_is_split_per_spectator_monomial() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let deep = 20;
    let (input, closed) = geometric_family(1000, 3, deep);
    let whole = RationalPart::new(&input, epsilon).expect("one part admits the input");
    assert!(
        whole.series(deep).is_none(),
        "the whole series must exceed its work bound"
    );
    let prepared = RationalSeries::new(&input, epsilon).unwrap();
    assert_eq!(prepared.valuation(), Some(0));
    let actual = prepared
        .series(deep)
        .expect("split series")
        .to_atom(epsilon);
    assert!((actual - &closed).expand().is_zero());
    let shallow = 2;
    let whole_series = whole.series(shallow).unwrap();
    let series = prepared.series(shallow).unwrap();
    assert_eq!(series.valuation, whole_series.valuation);
    assert_eq!(series.coefficients, whole_series.coefficients);
}

/// The split of a whole part keeps the whole part's (semantic) valuation:
/// a leading coefficient that only a proof makes zero (`sin² + cos² - 1`,
/// spread over three spectator monomials) is dropped, as the whole part does.
#[test]
fn a_split_whole_part_keeps_its_semantic_valuation() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let e = Atom::var(epsilon);
    let deep = 20;
    let (family, closed) = geometric_family(1000, 3, deep);
    let identity = vk_parse!("sin(x)^2+cos(x)^2-1").unwrap();
    let input = (&identity * e.pow(Atom::num(-3)) + family).expand();
    let whole = RationalPart::new(&input, epsilon).expect("one part admits the input");
    assert_eq!(whole.valuation(), Some(0));
    assert!(whole.series(deep).is_none());
    let prepared = RationalSeries::new(&input, epsilon).unwrap();
    assert_eq!(prepared.valuation(), Some(0));
    let series = prepared.series(deep).unwrap();
    assert_eq!(series.valuation, 0);
    assert_eq!(series.coefficients.len(), deep as usize + 1);
    assert!((series.to_atom(epsilon) - &closed).expand().is_zero());
}

/// Timing of both routes on the Q L4 1621 master monomial (ignored; run
/// alone, on one pinned core, with `--nocapture`): one line per depth with
/// the stand-in's size and the seconds of the whole exact part (if it admits
/// the factor), the split exact plan, the production term (exact) and the
/// generic term.
#[test]
#[ignore]
fn timing_q_l4_1621_master_monomial_both_routes() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let monomial = gamma_tags(&Q_L4_1621_MONOMIAL);
    let depths = std::env::var("VK_SERIES_DEPTHS").unwrap_or_else(|_| "4,7".to_owned());
    let structures = std::env::var("VK_SERIES_STRUCTURES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(STAND_IN_STRUCTURES);
    println!(
        "TIMING depth\torder\tsummands\tbytes\twhole_admits\twhole_s\tsplit_s\tterm_exact_s\tterm_generic_s\tequal"
    );
    for depth in depths.split(',').map(|d| d.trim().parse::<i64>().unwrap()) {
        let coefficient = r_star_like_factor(depth, structures, 1621);
        let summands = match coefficient.as_view() {
            symbolica::atom::AtomView::Add(sum) => sum.get_nargs(),
            _ => 1,
        };
        let order = four_loop_order(depth);
        let started = std::time::Instant::now();
        let whole = RationalPart::new(&coefficient, epsilon);
        let whole_admits = whole.is_some();
        if let Some(whole) = &whole {
            let _ = whole.series(order);
        }
        let whole_s = started.elapsed().as_secs_f64();
        let started = std::time::Instant::now();
        let split = RationalSeries::new(&coefficient, epsilon).map(|s| s.series(order));
        let split_s = started.elapsed().as_secs_f64();
        assert!(split.is_some_and(|s| s.is_some()));
        let started = std::time::Instant::now();
        let (exact, counts) = master_term(&coefficient, &monomial, order);
        let term_exact_s = started.elapsed().as_secs_f64();
        assert_eq!((counts.exact, counts.generic), (1, 0));
        let started = std::time::Instant::now();
        let generic = generic_master_term(&coefficient, &monomial, order);
        let term_generic_s = started.elapsed().as_secs_f64();
        let equal = (exact.unwrap().to_atom(epsilon) - generic.unwrap().to_atom(epsilon))
            .expand()
            .is_zero();
        println!(
            "TIMING {depth}\t{order}\t{summands}\t{}\t{whole_admits}\t{whole_s:.3}\t{split_s:.3}\t{term_exact_s:.3}\t{term_generic_s:.3}\t{equal}",
            coefficient.as_view().get_byte_size()
        );
    }
}

/// Sizes of the stand-in factors (ignored; `--nocapture`): summands, distinct
/// spectator monomials, whether one whole exact part admits the factor and
/// whether the exact plan (whole or split) does, for `VK_SERIES_STRUCTURES`
/// (comma list) and depths 4 to 9.
#[test]
#[ignore]
fn probe_stand_in_sizes() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let structures =
        std::env::var("VK_SERIES_STRUCTURES").unwrap_or_else(|_| "20,40,80".to_owned());
    println!("SIZE structures\tdepth\tsummands\tkeys\tbytes\twhole\tplan\tplan_s");
    for structures in structures
        .split(',')
        .map(|s| s.trim().parse::<usize>().unwrap())
    {
        for depth in 4..=9 {
            let factor = r_star_like_factor(depth, structures, 1621);
            let summands = match factor.as_view() {
                symbolica::atom::AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
                view => vec![view],
            };
            let keys = summands
                .iter()
                .map(|summand| {
                    let factors = match summand {
                        symbolica::atom::AtomView::Mul(product) => {
                            product.iter().collect::<Vec<_>>()
                        }
                        view => vec![*view],
                    };
                    Atom::mul_many(
                        factors
                            .into_iter()
                            .filter(|f| {
                                !matches!(f, symbolica::atom::AtomView::Num(_))
                                    && !f.contains_symbol(epsilon)
                            })
                            .collect::<Vec<_>>()
                            .iter(),
                    )
                })
                .collect::<std::collections::HashSet<_>>()
                .len();
            let whole = RationalPart::new(&factor, epsilon).is_some();
            let started = std::time::Instant::now();
            let plan = RationalSeries::new(&factor, epsilon)
                .and_then(|p| p.series(four_loop_order(depth)))
                .is_some();
            println!(
                "SIZE {structures}\t{depth}\t{}\t{keys}\t{}\t{whole}\t{plan}\t{:.3}",
                summands.len(),
                factor.as_view().get_byte_size(),
                started.elapsed().as_secs_f64()
            );
        }
    }
}

// The spectator split of a large input (the depth-7 record of Q L4 entry
// 1621, two-point row): the split checked the spectator coordinates and
// denominators of the WHOLE input before splitting it, a node walk capped at
// 2^20 visits. The recorded input (12.9 MB, 28,762 summands in 470 spectator
// groups) exceeds that walk although every group passes it, so it took the
// generic series. Every part runs the same check on its own expression.

/// The master monomial of that record: `PGamma(4) PGamma(-2)² PGamma(-3) /
/// (PGamma(1)² PGamma(2) PGamma(-1)² PGamma(-5))`.
const Q_L4_1621_RECORD_MONOMIAL: [(i64, i64); 7] = [
    (4, 1),
    (-2, 2),
    (-3, 1),
    (1, -2),
    (2, -1),
    (-1, -2),
    (-5, -1),
];

/// The order of the master term of that record (`order=5` at depth 7).
const Q_L4_1621_RECORD_ORDER: i64 = 5;

/// The ε-linear denominators `a + bε` of the summands of that record.
const RECORD_DENOMINATORS: [(i64, i64); 6] = [(-1, 4), (-1, 5), (-2, 5), (-3, 5), (1, 2), (1, 1)];

/// Visits of a node walk of `expression` that, like the spectator check of
/// the exact series, does not enter variables and functions.
fn spectator_walk_visits(expression: &Atom) -> usize {
    let mut visits = 0;
    expression.visitor(&mut |part| {
        visits += 1;
        !matches!(
            part,
            symbolica::atom::AtomView::Var(_) | symbolica::atom::AtomView::Fun(_)
        )
    });
    visits
}

/// A factor shaped like that record: `groups` distinct spectator monomials
/// (colour × ζ values × tensor × row power × flavour × the couplings
/// `gs⁶ gs_minus²`), each with `per_group` summands `c ε^p / Π D_j` of
/// distinct `(p, subset of RECORD_DENOMINATORS)`, `p` from -7 to 1, all times
/// the record's outer factors `𝑖 ε³ / (16 π² (1 - 2ε))`.
fn record_like_factor(groups: usize, per_group: usize, seed: u64) -> Atom {
    assert!(per_group <= 9 * 64, "(p, subset) pairs repeat");
    let e = Atom::var(vk_symbol!("ε"));
    let mut rng = Lcg(seed);
    let symbol = |name: &str| Atom::var(vk_symbol!(name));
    let (i2r, ca, cr) = (symbol("I2R"), symbol("cA"), symbol("cR"));
    let tensors = [symbol("transverse"), symbol("longitudinal")];
    let row = symbol("row_line_power");
    let flavour = symbol("skeleton_flavour_0");
    let couplings = symbol("gs").pow(Atom::num(6)) * symbol("gs_minus").pow(Atom::num(2));
    let mut keys = Vec::with_capacity(groups);
    'keys: for (power, flavours) in [(1i64, 1i64), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)] {
        for tensor in &tensors {
            for weight in 0..=9 {
                for zeta in zeta_monomials(weight) {
                    for a in 0..=3i64 {
                        for b in 0..=3 - a {
                            if keys.len() == groups {
                                break 'keys;
                            }
                            keys.push(
                                i2r.pow(Atom::num(a))
                                    * ca.pow(Atom::num(b))
                                    * cr.pow(Atom::num(3 - a - b))
                                    * &zeta
                                    * tensor
                                    * row.pow(Atom::num(power))
                                    * flavour.pow(Atom::num(flavours))
                                    * &couplings,
                            );
                        }
                    }
                }
            }
        }
    }
    assert_eq!(
        keys.len(),
        groups,
        "not enough distinct spectator monomials"
    );
    let denominators = RECORD_DENOMINATORS.map(|(a, b)| Atom::num(a) + Atom::num(b) * &e);
    let mut summands = Vec::with_capacity(groups * per_group);
    for (g, key) in keys.iter().enumerate() {
        for s in 0..per_group {
            let p = (s % 9) as i64 - 7;
            let subset = (7 * s + g) % 64;
            let sign = if rng.below(2) == 0 { 1 } else { -1 };
            let c =
                Atom::num(sign * (1 + rng.below(97) as i64)) / Atom::num(1 + rng.below(30) as i64);
            let mut summand = key * c * e.pow(Atom::num(p));
            for (j, denominator) in denominators.iter().enumerate() {
                if (subset >> j) & 1 == 1 {
                    summand /= denominator;
                }
            }
            summands.push(summand);
        }
    }
    vk_parse!("𝑖/16*ε^3*𝜋^(-2)/(1-2*ε)").unwrap() * Atom::add_many(summands.iter())
}

/// The coefficient of `ε^power` in `series` (zero outside its terms).
fn coefficient_at(series: &Laurent, power: i64) -> Atom {
    usize::try_from(power - series.valuation)
        .ok()
        .and_then(|i| series.coefficients.get(i).cloned())
        .unwrap_or(Atom::Zero)
}

/// The master term of `monomial` with `coefficient` takes the exact series,
/// and equals the generic series coefficient by coefficient through
/// `ε^order`.
fn assert_exact_term_equals_generic(
    coefficient: &Atom,
    monomial: &Atom,
    epsilon: Symbol,
    order: i64,
) {
    let (exact, counts) = master_term_in(coefficient, monomial, epsilon, order);
    assert_eq!(
        (counts.exact, counts.generic),
        (1, 0),
        "order {order}: the rational factor took the generic series"
    );
    let generic = generic_master_term_in(coefficient, monomial, epsilon, order);
    let (exact, generic) = (exact.expect("exact term"), generic.expect("generic term"));
    assert_eq!(
        exact.valuation, generic.valuation,
        "order {order}: valuation"
    );
    assert!(exact.order() >= order && generic.order() >= order);
    for power in exact.valuation..=order {
        assert!(
            (coefficient_at(&exact, power) - coefficient_at(&generic, power))
                .expand()
                .is_zero(),
            "order {order}, ε^{power}: exact and generic differ"
        );
    }
}

/// Many spectator groups, each inside the 2^20-visit walk of the spectator
/// check, the whole input beyond it (record-shaped: 470 groups × 64
/// summands): the split checks each group on its own, so the master term of
/// the record's monomial takes the exact series, equal to the generic series
/// coefficient by coefficient through ε^7 (the record's depth-7 order is
/// ε^5).
#[test]
fn spectator_groups_are_checked_per_group_not_as_one_whole_input() {
    crate::Vakint::initialize_vakint_symbols();
    let epsilon = vk_symbol!("ε");
    let factor = record_like_factor(470, 64, 1621);
    assert!(
        spectator_walk_visits(&factor) > 1 << 20,
        "the whole input must exceed the walk of the spectator check"
    );
    assert!(RationalSeries::validate_spectators(&factor, epsilon).is_none());
    assert!(
        RationalPart::new(&factor, epsilon).is_none(),
        "one part admits the factor: the test does not exercise the split"
    );
    let monomial = gamma_tags(&Q_L4_1621_RECORD_MONOMIAL);
    assert_exact_term_equals_generic(&factor, &monomial, epsilon, Q_L4_1621_RECORD_ORDER + 2);
}

/// One record of `VAKINT_SERIES_DECLINE_LOG`: the header `DECLINE phase=…
/// order=… exact_plan=…`, then `monomial=…` and `expression=…` in canonical
/// form.
struct DeclinedSeries {
    order: i64,
    monomial: Atom,
    expression: Atom,
    epsilon: Symbol,
}

/// The first record of the file named by `VAKINT_SERIES_FIXTURE`, in the ε
/// named by `VAKINT_SERIES_FIXTURE_EPSILON` (default `gammalooprs::ε`, the
/// caller of the Q L4 1621 record).
fn declined_series_fixture() -> DeclinedSeries {
    let path = std::env::var("VAKINT_SERIES_FIXTURE")
        .expect("VAKINT_SERIES_FIXTURE must name a record of VAKINT_SERIES_DECLINE_LOG");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{path}: {error}"));
    let mut lines = text.lines();
    let order = lines
        .next()
        .and_then(|header| {
            header
                .split_whitespace()
                .find_map(|field| field.strip_prefix("order="))
        })
        .and_then(|order| order.parse().ok())
        .expect("a header with order=");
    let parse = |text: &str| {
        Atom::parse(
            text,
            "vakint",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let mut field = |name: &str| {
        parse(
            lines
                .next()
                .and_then(|line| line.strip_prefix(name))
                .unwrap_or_else(|| panic!("a line {name}…")),
        )
    };
    let monomial = field("monomial=");
    let expression = field("expression=");
    let epsilon = std::env::var("VAKINT_SERIES_FIXTURE_EPSILON")
        .unwrap_or_else(|_| "gammalooprs::ε".to_owned());
    DeclinedSeries {
        order,
        monomial,
        expression,
        epsilon: parse(&epsilon).get_symbol().unwrap(),
    }
}

/// The recorded input of a declined exact series (ignored: the record is not
/// in the repository; `VAKINT_SERIES_FIXTURE` names it, e.g. the depth-7
/// record of Q L4 entry 1621, 12.9 MB, sha256 5f9eee5e540eb465…): its master
/// term takes the exact series at the recorded order, equal to the generic
/// series coefficient by coefficient.
#[test]
#[ignore]
fn a_recorded_declined_input_takes_the_exact_series() {
    crate::Vakint::initialize_vakint_symbols();
    let record = declined_series_fixture();
    println!(
        "FIXTURE order {} bytes {} spectator walk {} visits",
        record.order,
        record.expression.as_view().get_byte_size(),
        spectator_walk_visits(&record.expression)
    );
    assert_exact_term_equals_generic(
        &record.expression,
        &record.monomial,
        record.epsilon,
        record.order,
    );
}

/// Timing of one route on the record of `VAKINT_SERIES_FIXTURE` (ignored;
/// run it alone with `--nocapture` under `time -v` for the peak memory):
/// `VAKINT_SERIES_ROUTE` = `exact`, `generic`, or `parse` (reading the record
/// only, the baseline of the other two).
#[test]
#[ignore]
fn timing_a_recorded_declined_input() {
    crate::Vakint::initialize_vakint_symbols();
    let started = std::time::Instant::now();
    let record = declined_series_fixture();
    let parse_s = started.elapsed().as_secs_f64();
    let route = std::env::var("VAKINT_SERIES_ROUTE").unwrap_or_else(|_| "exact".to_owned());
    let started = std::time::Instant::now();
    let term = match route.as_str() {
        "parse" => None,
        "exact" => {
            let (term, counts) = master_term_in(
                &record.expression,
                &record.monomial,
                record.epsilon,
                record.order,
            );
            assert_eq!(
                (counts.exact, counts.generic),
                (1, 0),
                "the exact series was declined"
            );
            term
        }
        "generic" => generic_master_term_in(
            &record.expression,
            &record.monomial,
            record.epsilon,
            record.order,
        ),
        other => panic!("VAKINT_SERIES_ROUTE {other}: exact, generic or parse"),
    };
    println!(
        "TIMING route {route} order {} parse_s {parse_s:.3} route_s {:.3} valuation {:?} bytes {:?}",
        record.order,
        started.elapsed().as_secs_f64(),
        term.as_ref().map(|t| t.valuation),
        term.as_ref()
            .map(|t| t.to_atom(record.epsilon).as_view().get_byte_size())
    );
}
