//! D6 (helper VH): p-integrals with open external Lorentz indices.
//!
//! A p-integral numerator with open loop indices is reduced by the same FeynKit
//! tensor reducer as a vacuum numerator, with the external momentum declared as
//! an external vector; the scalar coefficients go through the p lowering and
//! RustRed, and the result keeps the open indices on loop-free tensors (metrics,
//! `p(1, mu)`). The caller applies any projector afterwards.
//!
//! The acceptance identity used throughout ("closed case"): for every closing
//! tensor `T` built from metrics and `p(1)` over the open indices, the open
//! result contracted with `T` equals the existing scalar path on the input
//! contracted with `T`, order by order in `ε` through the last order the scalar
//! result carries. Every evaluation here is RustRed alone with its built-in
//! 1L-3L library; no FORM executable exists at the configured path.

#[path = "support/open_indices.rs"]
mod open_indices;

use open_indices::*;
use symbolica::atom::{Atom, AtomCore, Symbol};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, TensorReductionMethod, Vakint, VakintError, VakintSettings, vakint_parse,
    vakint_symbol,
};

const MISSING_FORM: &str = "/d6-open-index-p-integrals-must-not-call-form";

/// The bubble with distinct powers (no automorphism).
const BUBBLE_12: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,2))";
const BUBBLE_11: &str = "topo(P1L(p(1),1,1))";
const KITE: &str = "topo(P2L(p(1),1,1,1,1,1))";
const KITE_POWERS: &str = "topo(P2L(p(1),1,2,1,1,1))";
const LADDER_3L: &str = "topo(P3L_LA(p(1),1,1,1,1,1,1,1,1))";

fn engine() -> Vakint {
    Vakint::new().unwrap()
}

fn settings(terms: i64) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: MISSING_FORM.into(),
        integral_normalization_factor: LoopNormalizationFactor::pySecDec,
        number_of_terms_in_epsilon_expansion: terms,
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

fn parse(text: &str) -> Atom {
    vakint_parse!(text).unwrap_or_else(|error| panic!("{text}: {error}"))
}

fn evaluate(engine: &Vakint, settings: &VakintSettings, input: &str) -> Atom {
    engine
        .evaluate(settings, parse(input).as_view())
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

fn epsilon(settings: &VakintSettings) -> Symbol {
    vakint_symbol!(settings.epsilon_symbol.as_str())
}

/// Every closing tensor of rank `indices.len()` built from metrics and `p(1)`:
/// a set of indices carried by `p(1)` and a metric pairing of the others.
fn closing_tensors(indices: &[&str]) -> Vec<String> {
    fn pairings(points: &[&str]) -> Vec<Vec<(String, String)>> {
        if points.is_empty() {
            return vec![vec![]];
        }
        let mut out = Vec::new();
        for i in 1..points.len() {
            let rest = points[1..i]
                .iter()
                .chain(&points[i + 1..])
                .copied()
                .collect::<Vec<_>>();
            for mut pairing in pairings(&rest) {
                pairing.push((points[0].to_owned(), points[i].to_owned()));
                out.push(pairing);
            }
        }
        out
    }
    let n = indices.len();
    let mut out = Vec::new();
    for mask in 0..(1u32 << n) {
        let carried = (0..n)
            .filter(|i| mask & (1 << i) != 0)
            .map(|i| indices[i])
            .collect::<Vec<_>>();
        let rest = (0..n)
            .filter(|i| mask & (1 << i) == 0)
            .map(|i| indices[i])
            .collect::<Vec<_>>();
        if rest.len() % 2 == 1 {
            continue;
        }
        for pairing in pairings(&rest) {
            let factors = carried
                .iter()
                .map(|index| format!("p(1,{index})"))
                .chain(pairing.iter().map(|(a, b)| format!("g({a},{b})")))
                .collect::<Vec<_>>();
            out.push(if factors.is_empty() {
                "1".to_owned()
            } else {
                factors.join("*")
            });
        }
    }
    out
}

/// The closed-case identity of `numerator × topology` for every closing
/// tensor over `indices`.
fn assert_all_closed_cases(
    engine: &Vakint,
    settings: &VakintSettings,
    numerator: &str,
    topology: &str,
    indices: &[&str],
) {
    let open = evaluate(engine, settings, &format!("({numerator})*{topology}"));
    let rendered = open.to_canonical_string();
    assert!(!rendered.contains("topo("), "unreduced result: {rendered}");
    for closing in closing_tensors(indices) {
        let scalar = evaluate(
            engine,
            settings,
            &format!("({numerator})*({closing})*{topology}"),
        );
        assert_closed_case(
            &open,
            &parse(&closing),
            &scalar,
            epsilon(settings),
            &format!("{numerator} * {topology}"),
        );
    }
}

/// Rank one, one loop: `∫ k^μ / (k² (k-p)²) = p^μ/2 ∫ 1 / (k² (k-p)²)`
/// (exact by the symmetry `k -> p - k` of the bubble).
#[test]
fn the_rank_one_bubble_is_half_the_external_momentum() {
    let engine = engine();
    let settings = settings(3);
    let open = evaluate(&engine, &settings, &format!("k(1,1)*{BUBBLE_11}"));
    let scalar = evaluate(&engine, &settings, BUBBLE_11);
    assert_closed_case(
        &open,
        &parse("p(1,1)"),
        &(scalar * parse("dot(p(1),p(1))") / Atom::num(2)),
        epsilon(&settings),
        "rank-one bubble",
    );
}

/// Ranks 1 to 3 on the one-loop bubble (powers 1, 2: no symmetry), every
/// closing tensor.
#[test]
fn one_loop_ranks_one_to_three_equal_the_scalar_path() {
    let engine = engine();
    let settings = settings(3);
    for (numerator, indices) in [
        ("k(1,1)", &["1"][..]),
        ("k(1,1)*k(1,2)", &["1", "2"][..]),
        ("k(1,1)*k(1,2)*k(1,3)", &["1", "2", "3"][..]),
        ("k(1,1)*k(1,2)*dot(k(1),p(1))", &["1", "2"][..]),
    ] {
        assert_all_closed_cases(&engine, &settings, numerator, BUBBLE_12, indices);
    }
}

/// Two loops: mixed and repeated loop momenta, odd and even ranks.
#[test]
fn two_loop_tensors_equal_the_scalar_path() {
    let engine = engine();
    let settings = settings(3);
    for (numerator, topology, indices) in [
        ("k(1,1)", KITE, &["1"][..]),
        ("k(2,1)*k(1,2)", KITE_POWERS, &["1", "2"][..]),
        ("k(2,1)*k(2,2)*k(1,3)", KITE, &["1", "2", "3"][..]),
    ] {
        assert_all_closed_cases(&engine, &settings, numerator, topology, indices);
    }
}

/// Three loops, rank two.
#[test]
fn a_three_loop_rank_two_tensor_equals_the_scalar_path() {
    let engine = engine();
    let settings = settings(2);
    assert_all_closed_cases(&engine, &settings, "k(1,1)*k(3,2)", LADDER_3L, &["1", "2"]);
}

/// A loop momentum contracted with another vector `q = p(2)` (a soft leg of
/// an N-point function, kept): the reducer treats `q` like an open index.
#[test]
fn a_loop_momentum_contracted_with_another_vector_is_reduced() {
    let engine = engine();
    let settings = settings(3);
    let open = evaluate(&engine, &settings, &format!("dot(k(1),p(2))*k(1,1)*{BUBBLE_12}"));
    let rendered = open.to_canonical_string();
    assert!(!rendered.contains("topo("), "unreduced result: {rendered}");
    // Closing q with p(1) gives back the scalar path.
    let replaced = open.replace(parse("p(2,x_)").to_pattern()).with(parse("p(1,x_)").to_pattern());
    let replaced = replaced
        .replace(parse("dot(p(2),y_)").to_pattern())
        .with(parse("dot(p(1),y_)").to_pattern())
        .replace(parse("dot(y_,p(2))").to_pattern())
        .with(parse("dot(y_,p(1))").to_pattern());
    let scalar = evaluate(&engine, &settings, &format!("dot(k(1),p(1))*k(1,1)*p(1,1)*{BUBBLE_12}"));
    assert_closed_case(&replaced, &parse("p(1,1)"), &scalar, epsilon(&settings), "foreign vector");
}

/// The isotropic FORM projector cannot reduce a p-integral: refused by name.
#[test]
fn the_alphaloop_projector_is_refused_for_tensor_p_integrals() {
    let settings = VakintSettings {
        tensor_reduction_method: TensorReductionMethod::AlphaLoop,
        ..settings(3)
    };
    match engine().evaluate(&settings, parse(&format!("k(1,1)*{BUBBLE_12}")).as_view()) {
        Err(VakintError::InvalidNumerator(message)) => assert!(
            message.contains("needs the FeynKit tensor reduction"),
            "{message}"
        ),
        other => panic!("expected the AlphaLoop refusal, got {other:?}"),
    }
}

/// The DEFAULT configuration on a host without FORM: RustRed first in every
/// order, the FeynKit reduction, and open indices on a p-integral and on a
/// vacuum integral evaluated without FORM.
#[test]
fn the_default_configuration_without_form_runs_open_indices() {
    assert!(matches!(
        EvaluationOrder::all().0.first(),
        Some(EvaluationMethod::RustRed(_))
    ));
    let order = EvaluationOrder::available(MISSING_FORM, "/d6-missing-python");
    assert!(
        matches!(order.0.as_slice(), [EvaluationMethod::RustRed(_)]),
        "the default order without FORM and python is {order:?}"
    );
    let settings = VakintSettings {
        evaluation_order: order,
        form_exe_path: MISSING_FORM.into(),
        python_exe_path: "/d6-missing-python".into(),
        ..VakintSettings::default()
    };
    assert_eq!(settings.tensor_reduction_method, TensorReductionMethod::FeynKit);
    let engine = engine();
    for input in [
        format!("k(1,1)*k(1,2)*{BUBBLE_11}"),
        format!("k(1,1)*{KITE}"),
        "k(1,1)*k(1,2)*topo(prop(1,edge(1,1),k(1),muvsq,2))".to_owned(),
    ] {
        let value = evaluate(&engine, &settings, &input);
        let rendered = value.to_canonical_string();
        assert!(!rendered.contains("topo("), "{input}: unreduced result {rendered}");
        assert!(
            !value.is_zero() && (rendered.contains("g(") || rendered.contains("p(1,")),
            "{input}: no open-index tensor in {rendered}"
        );
    }
}

/// Scope (D6 note, section 2b): metrics contracting loop indices with open
/// ones, Spenso slots, and a numerator-only vector `p(2)` open and contracted.
#[test]
fn the_vacuum_input_scope_is_accepted() {
    let engine = engine();
    let settings = settings(3);
    // Metrics carry symbolic indices: vakint's dot-notation pass reads an
    // integer-labelled metric `g(1,3)` as a vector `g(1)` with index 3 (a
    // pre-existing ambiguity of the notation; it fails with an error).
    for (numerator, indices) in [
        ("g(mu1,mu2)*k(1,mu2)*k(1,mu3)", &["mu1", "mu3"][..]),
        ("g(mu2,mu3)*k(1,mu2)*k(1,mu3)*k(1,mu1)", &["mu1"][..]),
        ("k(1,1)*k(1,2)*dot(k(1),p(2))", &["1", "2"][..]),
        ("k(1,1)*p(2,2)", &["1", "2"][..]),
    ] {
        assert_all_closed_cases(&engine, &settings, numerator, BUBBLE_12, indices);
    }
    // Spenso slots: the same result as raw indices, closed.
    let slotted = evaluate(
        &engine,
        &settings,
        &format!(
            "k(1,spenso::mink(4-2*ε,1))*k(1,spenso::mink(4-2*ε,2))*{BUBBLE_12}"
        ),
    );
    for closing in closing_tensors(&["1", "2"]) {
        let scalar = evaluate(
            &engine,
            &settings,
            &format!("k(1,1)*k(1,2)*({closing})*{BUBBLE_12}"),
        );
        assert_closed_case(&slotted, &parse(&closing), &scalar, epsilon(&settings), "Spenso slots");
    }
    // An opaque gamma-like tensor carrying the open slot is bridged.
    let gamma = evaluate(
        &engine,
        &settings,
        &format!(
            "k(1,spenso::mink(4,mu))*spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))*{BUBBLE_12}"
        ),
    );
    let rendered = gamma.to_canonical_string();
    assert!(
        !rendered.contains("topo(") && rendered.contains("gamma"),
        "gamma slot: {rendered}"
    );
}

/// Scope edges that stay refused.
#[test]
fn the_scope_edges_are_refused() {
    let engine = engine();
    let settings = settings(3);
    let error = |input: String| {
        engine
            .evaluate(&settings, parse(&input).as_view())
            .map(|value| format!("accepted as {value}"))
            .unwrap_or_else(|error| error.to_string())
    };
    let bare = error(format!("k(1)*{BUBBLE_12}"));
    assert!(bare.contains("loop momentum k(1) outside a scalar product"), "{bare}");
    for input in [
        format!("dot(k(1),p(2))^-1*k(1,1)*{BUBBLE_12}"),
        format!("k(1,1)*k(1,1)*k(1,1)*p(1,2)*{BUBBLE_12}"),
    ] {
        let message = error(input.clone());
        assert!(!message.starts_with("accepted"), "{input}: {message}");
    }
}

/// Analytic control at two loops (exact by symmetry): in the sunset
/// `P2L(p, 1,1,0,0,1)` (lines `k1`, `k2`, `k1+k2+p`) the three lines are
/// permuted by its symmetries, so `∫ k1^μ = ∫ k2^μ = -∫ (k1+k2+p)^μ`, hence
/// `∫ k1^μ / (k1² k2² (k1+k2+p)²) = -(p^μ/3) × sunset`.
#[test]
fn the_rank_one_sunset_is_minus_a_third_of_the_external_momentum() {
    let engine = engine();
    let settings = settings(3);
    let sunset = "topo(P2L(p(1),1,1,0,0,1))";
    let open = evaluate(&engine, &settings, &format!("k(1,1)*{sunset}"));
    let scalar = evaluate(&engine, &settings, sunset);
    assert_closed_case(
        &open,
        &parse("p(1,1)"),
        &(scalar * parse("dot(p(1),p(1))") / Atom::num(-3)),
        epsilon(&settings),
        "rank-one sunset",
    );
}

/// Analytic controls at two loops (exact by the left-right reflection of the
/// symmetric kite, which exchanges lines 1 <-> 4 and 3 <-> 5 of
/// `P2L(p, 1,1,1,1,1)`): `k1.p = (D4 - D1 - p²)/2` and
/// `k2.p = (D5 - D3 - D4 + D1)/2`, the differences integrate to zero, so
/// `∫ k1^μ = -(p^μ/2) × kite` and `∫ k2^μ = 0`.
#[test]
fn the_rank_one_symmetric_kite_follows_from_its_reflection() {
    let engine = engine();
    let settings = settings(3);
    let kite = evaluate(&engine, &settings, KITE);
    let first = evaluate(&engine, &settings, &format!("k(1,1)*{KITE}"));
    assert_closed_case(
        &first,
        &parse("p(1,1)"),
        &(&kite * parse("dot(p(1),p(1))") / Atom::num(-2)),
        epsilon(&settings),
        "kite k1",
    );
    let second = evaluate(&engine, &settings, &format!("k(2,1)*{KITE}"));
    assert!(second.expand().is_zero(), "kite k2: {second}");
}
