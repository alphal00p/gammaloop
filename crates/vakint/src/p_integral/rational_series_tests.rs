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
            colours.push(ca.pow(Atom::num(a)) * cf.pow(Atom::num(b)) * nf.pow(Atom::num(3 - a - b)));
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
    let epsilon = vk_symbol!("ε");
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
