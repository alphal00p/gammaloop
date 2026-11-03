//! The p-finalizer, the p-master table and the `MZV` constant, checked
//! against closed forms. No FORM is run here.

use std::collections::HashMap;

use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::float::{Complex, RealLike};
use symbolica::function;
use symbolica::transcendental::{gamma, zeta};

use super::super::masters::{mzv, p_master_expansion, p_master_table};
use super::super::{ReducedPropagatorValue, p_master_symbol, pg_symbol};
use super::{finalize, pg_value};
use crate::rustred_evaluation::RustRedMasterMode;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{LoopNormalizationFactor, Vakint, VakintError, VakintSettings};

pub(crate) fn settings(normalization: LoopNormalizationFactor, terms: i64) -> VakintSettings {
    Vakint::initialize_vakint_symbols();
    VakintSettings {
        integral_normalization_factor: normalization,
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        run_time_decimal_precision: 40,
        ..VakintSettings::default()
    }
}

/// Parse after installing vakint's tensor-tagged momentum symbols.
fn parse(text: &str) -> Atom {
    Vakint::initialize_vakint_symbols();
    vk_parse!(text).unwrap()
}

fn epsilon(settings: &VakintSettings) -> Atom {
    Atom::var(vk_symbol!(&settings.epsilon_symbol))
}

/// Laurent coefficients `ε^first..=ε^last` of `expression` at
/// `Q² = μ² = 1` (`dot(p,p) = -1`, every logarithm zero).
pub(crate) fn coefficients_at_unit_scale(
    expression: &Atom,
    settings: &VakintSettings,
    head: &Atom,
    first: i64,
    last: i64,
) -> Vec<Atom> {
    let e = epsilon(settings);
    let shifted = (expression
        .replace(function!(symbolica::atom::Symbol::LOG, Atom::var(S.x_)).to_pattern())
        .with(Atom::Zero.to_pattern())
        .replace(function!(S.dot, head, head).to_pattern())
        .with(Atom::num(-1).to_pattern())
        * e.pow(Atom::num(-first)))
    .expand();
    let mut coefficients = vec![Atom::Zero; (last - first + 1) as usize];
    for (monomial, coefficient) in shifted.coefficient_list::<i8>(std::slice::from_ref(&e)) {
        let power = if monomial == Atom::num(1) {
            0
        } else if monomial == e {
            1
        } else {
            let symbolica::atom::AtomView::Pow(power) = monomial.as_view() else {
                panic!("unexpected monomial {monomial}");
            };
            crate::get_integer_from_atom(power.get_base_exp().1).unwrap()
        };
        assert!(power >= 0, "a pole below ε^{first}: {expression}");
        if let Some(slot) = coefficients.get_mut(power as usize) {
            *slot = coefficient.expand();
        }
    }
    coefficients
}

/// `G(a,b) = Γ(a+b-d/2)Γ(d/2-a)Γ(d/2-b)/(Γ(a)Γ(b)Γ(d-a-b))`, zero when an
/// index is a non-positive integer.
fn g(a: Atom, b: Atom, e: &Atom) -> Atom {
    let half_d = Atom::num(2) - e;
    let d = Atom::num(4) - Atom::num(2) * e;
    gamma().call(&a + &b - &half_d) * gamma().call(&half_d - &a) * gamma().call(&half_d - &b)
        / (gamma().call(a.clone()) * gamma().call(b.clone()) * gamma().call(d - a - b))
}

fn assert_series_equal(lhs: &Atom, rhs: &Atom, settings: &VakintSettings, depth: i64) {
    let e = vk_symbol!(&settings.epsilon_symbol);
    let difference = super::gamma_series(&(lhs - rhs), e, depth).unwrap();
    assert!(
        difference.is_zero(),
        "{lhs}\n!=\n{rhs}\ndifference: {difference}"
    );
}

#[test]
fn pg_matches_g_functions() {
    let settings = settings(LoopNormalizationFactor::pySecDec, 4);
    let e = epsilon(&settings);
    let g11 = g(Atom::num(1), Atom::num(1), &e);
    for (x1, x2) in [(0, 0), (1, 0), (2, 0), (3, 1), (2, 2), (4, 0)] {
        let expected = Atom::num(1 + x1 + x2)
            * g(
                Atom::num(1) + Atom::num(x1) * &e,
                Atom::num(1) + Atom::num(x2) * &e,
                &e,
            )
            / (&e * &g11);
        let value = pg_value(x1, x2, vk_symbol!(&settings.epsilon_symbol));
        assert_series_equal(&value, &expected, &settings, 6);
    }
    // GschemeConstants(0,0) = 1/ε exactly.
    assert_series_equal(
        &pg_value(0, 0, vk_symbol!(&settings.epsilon_symbol)),
        &(Atom::num(1) / &e),
        &settings,
        8,
    );
}

#[test]
fn p_master_table_holds_every_forcer_master() {
    let table = p_master_table().unwrap();
    let mut names = table.keys().cloned().collect::<Vec<_>>();
    names.sort();
    assert_eq!(
        names,
        [
            "bebe", "cross", "fastar2", "haha", "lala", "no", "no1", "no2", "no6", "nono",
            "nostar5", "nostar6", "t105", "t124", "t145", "t155"
        ]
    );
    let e = vk_symbol!("ε");
    // Master(no) = 20 ζ5 + ε(50 ζ6 - 80 ζ5 + 68 ζ3²) + ... through ε^7.
    let no = p_master_expansion("no", e).unwrap();
    assert_eq!(
        no.coefficient(Atom::var(e).pow(Atom::num(1)).as_view())
            .expand(),
        (Atom::num(50) * zeta().call(6) - Atom::num(80) * zeta().call(5)
            + Atom::num(68) * zeta().call(3).pow(Atom::num(2)))
        .expand()
    );
    let tail = Atom::var(e).pow(Atom::num(8))
        * vk_symbol!("Oep").call((8, function!(p_master_symbol(), Atom::var(vk_symbol!("no")))));
    assert!(
        (no.coefficient(Atom::var(e).pow(Atom::num(8)).as_view()) * Atom::var(e).pow(Atom::num(8))
            - tail)
            .expand()
            .is_zero()
    );
    // The weight-twelve depth-four MZV of the deepest coefficients.
    let t105 = p_master_expansion("t105", e).unwrap();
    assert!(
        t105.coefficient(Atom::var(e).pow(Atom::num(9)).as_view())
            .contains_symbol(vk_symbol!("MZV"))
    );
}

fn numerical(expression: &Atom, settings: &VakintSettings) -> f64 {
    let value =
        crate::precision::materialize(expression.as_view(), settings, &HashMap::new(), None)
            .unwrap();
    value
        .evaluate_with_prec(
            &ahash::HashMap::<Atom, Complex<symbolica::domains::float::Float>>::default(),
            settings.get_binary_precision(),
        )
        .unwrap_or_else(|error| panic!("{value}: {error}"))
        .re
        .to_f64()
}

/// Both sides evaluated separately, at the relative tolerance of `f64`.
fn assert_numerically_equal(lhs: &Atom, rhs: &Atom, settings: &VakintSettings) {
    let (a, b) = (numerical(lhs, settings), numerical(rhs, settings));
    assert!(
        (a - b).abs() <= 1e-15 * b.abs(),
        "{lhs} = {a} != {rhs} = {b}"
    );
}

#[test]
fn mzv_constants_follow_the_datamine_convention() {
    let settings = settings(LoopNormalizationFactor::pySecDec, 4);
    // ζ(5,3) = Σ_{m>n} 1/(m^5 n^3).
    assert_numerically_equal(
        &mzv(&[5, 3]),
        &vk_parse!("0.037707672984847544011304782293659`40").unwrap(),
        &settings,
    );
    // Stuffle: ζ(5)ζ(3) = ζ(5,3) + ζ(3,5) + ζ(8).
    assert_numerically_equal(
        &(zeta().call(5) * zeta().call(3)),
        &(mzv(&[5, 3]) + mzv(&[3, 5]) + zeta().call(8)),
        &settings,
    );
    // Depth one is the Riemann zeta value.
    assert_numerically_equal(&mzv(&[7]), &zeta().call(7), &settings);
    // The weight-twelve depth-four value of Masters.prc evaluates.
    let deep = numerical(&mzv(&[6, 4, 1, 1]), &settings);
    assert!(deep > 0.0 && deep < 1e-2, "{deep}");
}

/// The value at `Q² = 1`, w = 2L - Σa + r, of one reduced part.
fn value(scale_power: i64, reduced: Atom) -> ReducedPropagatorValue {
    ReducedPropagatorValue {
        scale_power,
        reduced,
    }
}

fn pg(x1: i64, x2: i64) -> Atom {
    pg_symbol().call((x1, x2))
}

#[test]
fn g_scheme_normalization_returns_the_forcer_value() {
    // The scalar kite: Forcer gives ((1-2ε)/ε) PG(0,0)² + ((3ε-1)/ε) PG(0,0)PG(1,0),
    // whose G-scheme value is 6ζ3 + O(ε). With five lines the reduced part is
    // -(this), at w = -1.
    let settings = settings(LoopNormalizationFactor::GScheme, 3);
    let e = epsilon(&settings);
    let kite = (Atom::num(1) - Atom::num(2) * &e) / &e * pg(0, 0).pow(Atom::num(2))
        + (Atom::num(3) * &e - Atom::num(1)) / &e * pg(0, 0) * pg(1, 0);
    let head = parse("p(3)");
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let result = finalize(
            &[value(-1, -kite.clone())],
            2,
            &head,
            &settings,
            mode,
            "test",
        )
        .unwrap();
        let coefficients = coefficients_at_unit_scale(&result, &settings, &head, -2, 0);
        assert!(
            coefficients[0].is_zero() && coefficients[1].is_zero(),
            "{result}"
        );
        // I = (Q²)^w (Q²/μ²)^(-2ε) (-C_G) at Q² = μ² = 1.
        assert_eq!(coefficients[2], (Atom::num(-6) * zeta().call(3)).expand());
    }
    // The finite G-scheme value of the three-loop tops LA and NO is 20ζ5.
    let settings = super::tests::settings(LoopNormalizationFactor::GScheme, 4);
    let no = function!(p_master_symbol(), Atom::var(vk_symbol!("no")));
    let result = finalize(
        &[value(-2, no)],
        3,
        &head,
        &settings,
        RustRedMasterMode::Symbolic,
        "test",
    )
    .unwrap();
    let coefficients = coefficients_at_unit_scale(&result, &settings, &head, -3, 0);
    assert_eq!(coefficients[3], (Atom::num(20) * zeta().call(5)).expand());
}

#[test]
fn normalizations_share_the_measure_and_keep_the_scale_logarithm() {
    // The one-loop bubble, reduced = PG(0,0) at w = 0: in pySecDec
    // normalization I = εG(1,1)·PG(0,0)·(Q²/μ²)^(-ε)(μ²)^(-ε) = G(1,1)(Q²)^(-ε).
    let head = parse("p(1)");
    let settings = settings(LoopNormalizationFactor::pySecDec, 4);
    let e = epsilon(&settings);
    let retained = finalize(
        &[value(0, pg(0, 0))],
        1,
        &head,
        &settings,
        RustRedMasterMode::Retain,
        "test",
    )
    .unwrap();
    assert!(retained.contains_symbol(pg_symbol()), "{retained}");
    let symbolic = finalize(
        &[value(0, pg(0, 0))],
        1,
        &head,
        &settings,
        RustRedMasterMode::Symbolic,
        "test",
    )
    .unwrap();
    let log_q = function!(
        symbolica::atom::Symbol::LOG,
        -function!(S.dot, &head, &head) / Atom::var(vk_symbol!("mursq"))
    );
    let log_mu = function!(symbolica::atom::Symbol::LOG, Atom::var(vk_symbol!("mursq")));
    // Opaque logarithms: Symbolica would merge log(Q²/μ²) + log(μ²).
    let (l_q, l_mu) = (Atom::var(vk_symbol!("Lq")), Atom::var(vk_symbol!("Lmu")));
    let opaque = symbolic
        .replace(log_q.to_pattern())
        .with(l_q.to_pattern())
        .replace(log_mu.to_pattern())
        .with(l_mu.to_pattern());
    let bubble = g(Atom::num(1), Atom::num(1), &e)
        * function!(symbolica::atom::Symbol::EXP, -&e * (&l_q + &l_mu));
    assert_series_equal(&opaque, &bubble, &settings, 2);
    // MSbar cancels log(mursq) exactly.
    let msbar = finalize(
        &[value(0, pg(0, 0))],
        1,
        &head,
        &super::tests::settings(LoopNormalizationFactor::MSbar, 4),
        RustRedMasterMode::Symbolic,
        "test",
    )
    .unwrap();
    let contains = |expression: &Atom, part: &Atom| {
        let mut found = false;
        expression.visitor(&mut |view| {
            found |= view == part.as_view();
            !found
        });
        found
    };
    assert!(!contains(&msbar, &log_mu), "{msbar}");
    assert!(contains(&msbar, &log_q), "{msbar}");
}

#[test]
fn p_master_tails_fail_closed() {
    // NO is known through ε^7. An ε^-1 coefficient needs ε^(target+1).
    let head = parse("p(1)");
    let no = function!(p_master_symbol(), Atom::var(vk_symbol!("no")));
    let evaluate = |terms: i64| {
        let settings = settings(LoopNormalizationFactor::GScheme, terms);
        let reduced = &no / epsilon(&settings);
        finalize(
            &[value(-2, reduced)],
            3,
            &head,
            &settings,
            RustRedMasterMode::Symbolic,
            "Forcer",
        )
    };
    // Ten terms from ε^-3 end at ε^6 and need NO through ε^7.
    assert!(!evaluate(10).unwrap().is_zero());
    match evaluate(11) {
        Err(VakintError::MissingMasterOrders(missing)) => {
            assert_eq!(missing.requested_order, Some(7));
            assert_eq!(missing.masters.len(), 1, "{missing}");
            assert!(missing.masters[0].master.contains("no"), "{missing}");
            assert_eq!(missing.masters[0].first_missing_order, 8, "{missing}");
        }
        other => panic!("expected missing master orders, got {other:?}"),
    }
}

#[test]
fn g_scheme_normalization_validates() {
    let settings = settings(LoopNormalizationFactor::GScheme, 4);
    let (normalization, _, _) = settings
        .integral_normalization_factor
        .validate(&settings)
        .unwrap();
    assert!(normalization.contains_symbol(symbolica::transcendental::gamma()));
}
