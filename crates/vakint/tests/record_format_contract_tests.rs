//! The CONTRACT of the record format (`WH-RECORD-FORMAT.md` v1.4, "v1.4 §n",
//! and its corrections, "corr. n") through the PUBLIC entries on real
//! records (the vendored fixtures; the producer's fixtures under
//! `VAKINT_TEST_WH_FIXTURES` for the `#[ignore]`d ones): what comes back,
//! linearity and scaling identities, the same function written differently,
//! the masses, text quirks, and the validation rules of marked records on
//! every field and entry.

use std::path::PathBuf;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use vakint::p_integral::factorized::{
    ClosureValue, FactorizedKind, evaluate_factorized_closure, evaluate_factorized_ports,
};
use vakint::p_integral::{FactorizedError, FactorizedNumerator};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};

#[path = "support/wh_fixtures.rs"]
mod wh_fixtures;

const R: &str = "gammalooprs::{scalar}::rstar_reference_mass_squared";
const MAUX: &str = "gammalooprs::{scalar}::rstar_mass_squared";
const DIM: &str = "gammalooprs::{}::dim";
const EPS: &str = "gammalooprs::{}::ε";

fn init() {
    Vakint::initialize_vakint_symbols();
    let _ = symbolica::try_parse!(
        "spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(1,1))*spenso::{spenso::representation,spenso::self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(2,1))*spenso::{symmetric,linear,real,spenso::tensor}::g(x,y)*gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q(1,z)*gammalooprs::{scalar}::rstar_reference_mass_squared*gammalooprs::{scalar}::rstar_mass_squared"
    )
    .unwrap();
}

/// A fully qualified expression (the checks write every name with its
/// namespace: the test crate's default namespace is never used).
fn atom(text: &str) -> Atom {
    init();
    symbolica::try_parse!(text).unwrap_or_else(|e| panic!("{e}: {text}"))
}

fn eps() -> Symbol {
    init();
    symbolica::symbol!("gammalooprs::ε")
}

fn engine() -> Vakint {
    init();
    let mut vakint = Vakint::new().unwrap();
    if let Some(library) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") {
        vakint
            .load_p_library_directory(std::path::Path::new(&library), false)
            .unwrap();
    }
    vakint
}

/// The R* site settings (v1.4 §7 "settings": eps `gammalooprs::ε`, scale
/// `gammalooprs::μᵣ²`), RustRed only in `mode`, FORM absent.
fn site(mode: RustRedMasterMode) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        form_exe_path: "/nonexistent/helper-vh/form".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 3,
        use_dot_product_notation: true,
        n_cores: 1,
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        ..Default::default()
    }
}

fn symbolic() -> VakintSettings {
    site(RustRedMasterMode::Symbolic)
}

fn fixture_text(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src/p_integral/factorized/fixtures")
        .join(name);
    wh_fixtures::read_fixture(&path)
}

fn record(name: &str) -> FactorizedNumerator {
    init();
    FactorizedNumerator::from_text(&fixture_text(&format!("{name}.record.txt")))
        .unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn closure(
    v: &FactorizedNumerator,
    settings: &VakintSettings,
    argument: &Atom,
) -> Result<ClosureValue, FactorizedError> {
    evaluate_factorized_closure(
        &engine(),
        settings,
        std::slice::from_ref(v),
        std::slice::from_ref(argument),
    )
}

fn value(v: &FactorizedNumerator, settings: &VakintSettings) -> Atom {
    closure(v, settings, &v.prefactor.total())
        .unwrap_or_else(|e| panic!("{}: {e}", v.label))
        .value
}

fn is_zero(a: &Atom) -> bool {
    a.expand().is_zero()
}

/// The largest ε order present (the order the entry truncates at).
fn top_order(a: &Atom) -> i64 {
    *wh_fixtures::epsilon_coefficients(a, eps())
        .keys()
        .max()
        .expect("a nonzero value")
}

/// `a` truncated after ε^order.
fn truncate(a: &Atom, order: i64) -> Atom {
    let e = Atom::var(eps());
    let mut out = Atom::Zero;
    for (n, c) in wh_fixtures::epsilon_coefficients(a, eps()) {
        if n <= order {
            out += c * e.pow(Atom::num(n));
        }
    }
    out
}

fn assert_same(got: &Atom, expected: &Atom, what: &str) {
    let difference = (got - expected).expand();
    assert!(
        difference.is_zero(),
        "{what}:\n got      {got}\n expected {expected}\n difference {difference}"
    );
}

fn detail(e: &FactorizedError) -> String {
    e.to_string()
}

// ===========================================================================
// P1: the promises of v1.4 §7, sentence by sentence
// ===========================================================================

/// v1.4 §7 "value: Laurent series in eps through depth, exact rationals" and
/// "settings: the call takes the caller's VakintSettings (... evaluation
/// order)": in Numerical mode a marked K1 record evaluates, through the same
/// orders, with the Symbolic coefficients evaluated numerically (here to
/// 1e-12 relative).
#[test]
fn numerical_mode_k1_root_equals_symbolic_numerically() {
    let v = record("wh-k1-root");
    let exact = value(&v, &symbolic());
    let numerical = value(&v, &site(RustRedMasterMode::Numerical));
    let e = wh_fixtures::epsilon_coefficients(&exact, eps());
    let n = wh_fixtures::epsilon_coefficients(&numerical, eps());
    assert_eq!(
        e.keys().collect::<Vec<_>>(),
        n.keys().collect::<Vec<_>>(),
        "orders differ\nexact {exact}\nnumerical {numerical}"
    );
    let pi = atom("𝜋");
    let pi_value = atom("3.14159265358979323846264338327950288419716939937510582097494459");
    for (order, exact_c) in &e {
        let numerical_c = &n[order];
        let difference = (exact_c
            .replace(pi.clone())
            .with(pi_value.clone())
            .to_float(60)
            - numerical_c
                .replace(pi.clone())
                .with(pi_value.clone())
                .to_float(60))
        .expand();
        // every remaining term is a float times symbols: tiny relative to 1
        let terms = match difference.as_view() {
            AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
            _ if difference.is_zero() => Vec::new(),
            _ => vec![difference.clone()],
        };
        for term in terms {
            let mut magnitude = 1.0f64;
            if let AtomView::Mul(m) = term.as_view() {
                for f in m.iter() {
                    if let AtomView::Num(_) = f {
                        magnitude = parse_magnitude(&f.to_owned().to_string());
                    }
                }
            } else if let AtomView::Num(_) = term.as_view() {
                magnitude = parse_magnitude(&term.to_string());
            }
            assert!(
                magnitude.abs() < 1e-12,
                "eps^{order}: numerical differs from symbolic by {term}\nexact {exact_c}\nnumerical {numerical_c}"
            );
        }
    }
}

/// The largest modulus among the float literals of a printed number (a
/// real or complex float `a+b𝑖`); 1 when it has none (an exact number).
fn parse_magnitude(text: &str) -> f64 {
    let mut best: Option<f64> = None;
    let mut token = String::new();
    let flush = |token: &mut String, best: &mut Option<f64>| {
        if token.contains('.') || token.contains('e') {
            if let Ok(x) = token.parse::<f64>() {
                *best = Some(best.map_or(x.abs(), |b: f64| b.max(x.abs())));
            }
        }
        token.clear();
    };
    let chars = text.chars().collect::<Vec<_>>();
    for (i, c) in chars.iter().enumerate() {
        let exponent_sign = (*c == '-' || *c == '+') && i > 0 && chars[i - 1] == 'e';
        if c.is_ascii_digit() || *c == '.' || *c == 'e' || exponent_sign {
            token.push(*c);
        } else {
            flush(&mut token, &mut best);
        }
    }
    flush(&mut token, &mut best);
    best.unwrap_or(1.0)
}

/// v1.4 §7 "depth: ... Under master mode Retain there is no eps
/// truncation": a K1 record WITHOUT insertions evaluates in Retain; the
/// retained value, re-finalised by vakint's own re-finaliser in Symbolic,
/// equals the Symbolic value (the retained value is the exact object).
#[test]
fn retain_k1_root_without_insertions_refinalises_to_symbolic() {
    let v = record("wh-k1-root");
    let exact = value(&v, &symbolic());
    let top = top_order(&exact);
    let retained = closure(&v, &site(RustRedMasterMode::Retain), &v.prefactor.total())
        .unwrap_or_else(|e| panic!("Retain: {e}"))
        .value;
    eprintln!("CONTRACT RETAIN k1-root value: {retained}");
    // the retained value is exact in ε (Γ functions, powers M^(-2ε)): its
    // Laurent expansion through the Symbolic value's top order equals it
    let series = retained
        .series(
            eps(),
            Atom::Zero.as_view(),
            symbolica::domains::rational::Rational::from(top),
        )
        .map(|s| s.to_atom().expand());
    match series {
        Ok(series) => {
            eprintln!("CONTRACT RETAIN expanded: {series}");
            assert_same(
                &truncate(&series, top),
                &exact,
                "Retain expanded vs Symbolic",
            );
        }
        Err(e) => panic!("the Retain value cannot be expanded by symbolica: {e}"),
    }
}

/// v1.4 §7 "UNKNOWN masters (C8): UNKNOWN(master_id, k) eps^k atoms inside
/// the coefficients, linear, stable master_id; ... allow-unknowns is OPT-IN"
/// and corr. 4 ("UNKNOWN(master, k)*ε^-2" allowed): the UNKNOWN atom vakint
/// writes (`vakint::UNKNOWN`, truncation_guard.rs) in an insertion
/// coefficient of a marked record is carried LINEARLY: the value equals the
/// value of the record with coefficient c times (1 + UNKNOWN(..)) = value(c)
/// + UNKNOWN(..) value(c) (allow-unknowns ON, as the producer's runs that
/// need it).
#[test]
fn unknown_from_vakint_in_an_insertion_is_carried_linearly() {
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1") };
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    let unknown = atom("vakint::UNKNOWN(vakint::ggppss_int_53,0)");
    let mut v = base.clone();
    v.insertions[0].coefficient = &v.insertions[0].coefficient * (Atom::num(1) + &unknown);
    let got = closure(&v, &symbolic(), &v.prefactor.total());
    match got {
        Ok(got) => assert_same(
            &got.value,
            &(&plain * (Atom::num(1) + &unknown)),
            "UNKNOWN linear",
        ),
        Err(e) => panic!(
            "the insertion with vakint's UNKNOWN atom is refused: {}",
            detail(&e)
        ),
    }
}

/// The same with the caller's own UNKNOWN symbol (`gammalooprs::UNKNOWN`,
/// outside vakint's namespace): carried linearly.
#[test]
fn unknown_caller_symbol_in_an_insertion_is_carried_linearly() {
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    let unknown = atom("gammalooprs::UNKNOWN(gammalooprs::m53,0)");
    let mut v = base.clone();
    v.insertions[0].coefficient = &v.insertions[0].coefficient * (Atom::num(1) + &unknown);
    let got = closure(&v, &symbolic(), &v.prefactor.total()).unwrap();
    assert_same(
        &got.value,
        &(&plain * (Atom::num(1) + &unknown)),
        "caller UNKNOWN linear",
    );
}

/// UNKNOWN in the PREFACTOR ARGUMENT, multiplying an order the depth CUTS:
/// `x * ε^2` (x a caller scalar) times the K1 value (orders ε^-2, ε^-1)
/// lands entirely beyond the truncation of the returned product, whose
/// orders are fixed by the record's depth (v1.4 §7 "depth": truncated at
/// eps^(depth-L-1)). Expected: the value is the product truncated there,
/// i.e. ZERO for an argument ε^2 (the argument rides inside the boundary
/// coefficients, before the truncation). The check records which: exact
/// product (no truncation of the argument's orders) or truncated.
#[test]
fn an_argument_power_beyond_depth_is_truncated() {
    let v = record("wh-k1-root");
    let plain = value(&v, &symbolic());
    let top = top_order(&plain);
    let x = atom(&format!("gammalooprs::x*{EPS}^2"));
    let got = closure(&v, &symbolic(), &(v.prefactor.total() * &x))
        .unwrap()
        .value;
    eprintln!("CONTRACT ARGUMENT eps^2: plain top order {top}; got {got}");
    assert_same(
        &truncate(&(&plain * &x), top),
        &got,
        "argument x eps^2 truncated",
    );
}

/// v1.4 §7 "zero: exact 0 for a scaleless remainder (as their Built::Zero)":
/// a closure record whose massless lines form a sub-loop that does not
/// depend on the closing line (after the Taylor step a massless bubble at
/// zero external momentum: scaleless) returns EXACTLY 0.
#[test]
fn a_scaleless_remainder_is_exactly_zero() {
    let text = scaleless_record_text();
    let v = FactorizedNumerator::from_text(&text).unwrap_or_else(|e| panic!("{e}"));
    let got = closure(&v, &symbolic(), &v.prefactor.total())
        .unwrap_or_else(|e| panic!("refused: {}", detail(&e)));
    assert!(got.value.is_zero(), "not an exact zero: {}", got.value);
    assert!(
        got.vakint_value.is_zero(),
        "not an exact zero: {}",
        got.vakint_value
    );
}

/// Two loops: a massless bubble of lines 1 (k1+q) and 2 (k1) between
/// vertices 1 and 2, and the closing line 3 (k2) a self-loop at vertex 1:
/// after the Taylor step the bubble has no external momentum (scaleless).
fn scaleless_record_text() -> String {
    let q = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
    let mink = "spenso::{spenso::representation,spenso::self_dual}::mink";
    let hedge = "gammalooprs::{spenso::index}::hedge";
    let g = "spenso::{symmetric,linear,real,spenso::tensor}::g";
    let lines = vec![
        "schema\tpint6.factorized-numerator.v3".to_string(),
        "label\ta scaleless remainder".into(),
        "kind\tclosure\texternals=conserved".into(),
        format!("dimension\t{DIM}"),
        "momenta\tloops=2\texternals=1".into(),
        "loop\t1\tedge=2".into(),
        "loop\t2\tedge=3".into(),
        "external\t1\tedge=0".into(),
        "edges\t5".into(),
        "edge\t0\tsource=ext\ttarget=1\thedges=-,0\tparticle=B\tline=scalar\tpower=1\tmomentum=0,0,1\t0".into(),
        "edge_factor\t0\t1".into(),
        "edge\t1\tsource=1\ttarget=2\thedges=2,3\tparticle=s\tline=scalar\tpower=1\tmomentum=1,0,1\t0".into(),
        "edge_factor\t1\t1".into(),
        "edge\t2\tsource=2\ttarget=1\thedges=4,5\tparticle=s\tline=scalar\tpower=1\tmomentum=1,0,0\t0".into(),
        "edge_factor\t2\t1".into(),
        format!("edge\t3\tsource=1\ttarget=1\thedges=6,7\tparticle=s\tline=scalar\tpower=1\tmomentum=0,1,0\t{R}"),
        "edge_factor\t3\t1".into(),
        "edge\t4\tsource=2\ttarget=ext\thedges=1,-\tparticle=B\tline=scalar\tpower=1\tmomentum=0,0,1\t0".into(),
        "edge_factor\t4\t1".into(),
        "vertices\t2".into(),
        format!("vertex\t1\trule=V\t{q}(1,{mink}({DIM},{hedge}(0,1)))"),
        format!("vertex\t2\trule=V\t{q}(1,{mink}({DIM},{hedge}(1,1)))"),
        "ports\t2".into(),
        format!("port\t{mink}({DIM},{hedge}(0,1))"),
        format!("port\t{mink}({DIM},{hedge}(1,1))"),
        format!("projector\tp2_power=-1\t{g}({mink}({DIM},{hedge}(0,1)),{mink}({DIM},{hedge}(1,1)))"),
        "open\t0".into(),
        "prefactor\tloop_sign=1\tnf_power=0".into(),
        "prefactor_overall\t1".into(),
        "prefactor_scalars\tUFO::{}::gs^2".into(),
        "prefactor_colour\t1".into(),
        "colour_branch\tnone".into(),
        format!("closure\tedge=3\tpower=1\t{R}"),
        format!("taylor\tdegree=2\tcompensate_mass=1\tscaled_edges=0\torders=all\t{MAUX}"),
        "softs\t0".into(),
        "original_cycles\tfermion=0\tghost=0".into(),
        "depth\t2".into(),
        "insertions\t0".into(),
    ];
    let n = lines.len() + 1;
    let mut text = lines.join("\n");
    text.push_str(&format!("\nlines\t{n}\nend\n"));
    text
}

/// v1.4 §7 "weights: ... never n_f, loop sign, 1/AutG": for a record WITH
/// `original_cycles` (the producer's) the record cannot carry n_f
/// (`nf_power=0` is required) and the returned value contains no flavour
/// symbol and no vakint symbol.
#[test]
fn never_nf_with_original_cycles() {
    let v = record("wh-k1-root");
    let mut nf = v.clone();
    nf.prefactor.nf_power = 1;
    assert!(
        nf.validate().is_err(),
        "nf_power=1 accepted with original_cycles"
    );
    let got = value(&v, &symbolic());
    let mut bad = None;
    got.as_view().visitor(&mut |part| {
        let name = match part {
            AtomView::Var(x) => x.get_symbol().get_name().to_string(),
            AtomView::Fun(f) => f.get_symbol().get_name().to_string(),
            _ => return true,
        };
        if name.ends_with("::nf") || name.starts_with("vakint::") {
            bad = Some(name);
        }
        bad.is_none()
    });
    assert_eq!(bad, None, "the value carries {bad:?}: {got}");
}

/// The same promise for a closure record WITHOUT `original_cycles` (allowed
/// by the code: the v1 standalone weights): with a closed fermion loop the
/// record must say nf_power=1, and `Prefactor::total()` (the documented
/// argument, closure.rs) writes `vakint::nf`: the returned value then
/// carries n_f, in vakint's reserved namespace. Expected by v1.4 §7: never
/// n_f in the value (and no vakint:: symbol, corr. 1). The check records it.
#[test]
fn never_nf_without_original_cycles() {
    let base = record("wh-k1-root");
    let mut v = base.clone();
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.original_cycles = None;
    }
    v.prefactor.loop_sign = -1; // the carried graph's ghost loop (v1 standalone rule)
    match v.validate() {
        Ok(()) => {}
        Err(e) => {
            eprintln!("CONTRACT NEVER-NF: without original_cycles the record is refused: {e}");
            return;
        }
    }
    let got = value(&v, &symbolic());
    assert!(
        !got.to_string().contains("nf"),
        "the value carries n_f: {got}"
    );
}

// ===========================================================================
// P2: linearity and scaling identities on real records
// ===========================================================================

/// A fresh caller scalar.
fn z() -> Atom {
    atom("caller::z")
}

/// One closure record of each kind (vendored): K1 root, K1 root rc, K2 child,
/// K2 child with an insertion, K3 hs, K3 Dirac, K4 quotient, K5 surrogate.
const CLOSURE_KINDS: [&str; 8] = [
    "wh-k1-root",
    "wh-k1-root-rc",
    "wh-k2-child",
    "wh-k2-child-ins",
    "wh-k3-hs",
    "wh-k3-dirac",
    "wh-k4-quotient",
    "wh-k5-surrogate",
];

/// (a) z in the prefactor ARGUMENT scales the value by z exactly (Symbolic),
/// every kind.
#[test]
fn z_in_the_argument_scales_every_kind() {
    for name in CLOSURE_KINDS {
        let v = record(name);
        let plain = value(&v, &symbolic());
        let got = closure(&v, &symbolic(), &(v.prefactor.total() * z()))
            .unwrap()
            .value;
        assert_same(&got, &(&plain * z()), &format!("{name}: argument z"));
    }
}

/// (a) z in prefactor_scalars, in an edge or vertex FACTOR (the first one
/// that is not 1) and in an INSERTION coefficient scales the value by z.
#[test]
fn z_in_a_field_scales_every_kind() {
    for name in CLOSURE_KINDS {
        let v = record(name);
        let plain = value(&v, &symbolic());
        let expected = &plain * z();
        let mut s = v.clone();
        s.prefactor.scalars = &s.prefactor.scalars * z();
        assert_same(
            &value(&s, &symbolic()),
            &expected,
            &format!("{name}: prefactor_scalars z"),
        );
        let mut f = v.clone();
        let vertex = f
            .vertices
            .iter_mut()
            .find(|x| x.factor != Atom::num(1))
            .expect("a vertex factor");
        vertex.factor = &vertex.factor * z();
        assert_same(
            &value(&f, &symbolic()),
            &expected,
            &format!("{name}: vertex factor z"),
        );
        if !v.insertions.is_empty() {
            let mut i = v.clone();
            i.insertions[0].coefficient = &i.insertions[0].coefficient * z();
            assert_same(
                &value(&i, &symbolic()),
                &expected,
                &format!("{name}: insertion z"),
            );
        }
    }
}

/// (a) 1/(1+zε) in the argument / in prefactor_scalars / in an insertion
/// coefficient (pole order 0, no depth raise): the value is the series
/// product truncated at the record's top order (v1.4 §7 "Laurent series in
/// eps through depth"), every kind (Symbolic).
#[test]
fn one_over_one_plus_z_eps_is_the_truncated_product() {
    let factor = atom(&format!("(1+caller::z*{EPS})^-1"));
    for name in CLOSURE_KINDS {
        let v = record(name);
        let plain = value(&v, &symbolic());
        let top = top_order(&plain);
        let series = |a: &Atom| -> Atom {
            // a * (1 - z eps + z^2 eps^2 - ...) through `top`
            let mut s = Atom::Zero;
            let e = Atom::var(eps());
            for k in 0..=(top
                - wh_fixtures::epsilon_coefficients(a, eps())
                    .keys()
                    .min()
                    .copied()
                    .unwrap()
                + 1)
            {
                s += (Atom::num(-1) * z() * &e).pow(Atom::num(k));
            }
            truncate(&(a * s).expand(), top)
        };
        let expected = series(&plain);
        let got = closure(&v, &symbolic(), &(v.prefactor.total() * &factor))
            .unwrap()
            .value;
        assert_same(
            &truncate(&got, top),
            &expected,
            &format!("{name}: argument 1/(1+z eps)"),
        );
        assert_eq!(
            top_order(&got),
            top,
            "{name}: argument: the returned orders changed"
        );
        let mut s = v.clone();
        s.prefactor.scalars = &s.prefactor.scalars * &factor;
        let got = value(&s, &symbolic());
        assert_same(
            &truncate(&got, top),
            &expected,
            &format!("{name}: prefactor_scalars 1/(1+z eps)"),
        );
        if !v.insertions.is_empty() {
            let mut i = v.clone();
            i.insertions[0].coefficient = &i.insertions[0].coefficient * &factor;
            let got = value(&i, &symbolic());
            assert_same(
                &truncate(&got, top),
                &expected,
                &format!("{name}: insertion 1/(1+z eps)"),
            );
        }
    }
}

/// (a) in RETAIN, a fresh z in the argument and in a factor scales the
/// retained value by z exactly (records without insertions: no truncation
/// is involved, v1.4 §7).
#[test]
fn z_scales_in_retain() {
    for name in ["wh-k1-root", "wh-k2-child", "wh-k3-hs"] {
        let v = record(name);
        let retain = site(RustRedMasterMode::Retain);
        let plain = value(&v, &retain);
        let got = closure(&v, &retain, &(v.prefactor.total() * z()))
            .unwrap()
            .value;
        assert_same(&got, &(&plain * z()), &format!("{name}: Retain argument z"));
        let mut f = v.clone();
        let vertex = f
            .vertices
            .iter_mut()
            .find(|x| x.factor != Atom::num(1))
            .unwrap();
        vertex.factor = &vertex.factor * z();
        assert_same(
            &value(&f, &retain),
            &(&plain * z()),
            &format!("{name}: Retain factor z"),
        );
    }
}

/// (a) in RETAIN, 1/(1+zε) in the ARGUMENT of a record without insertions:
/// the retained value is exact in ε, so the value with the factor equals
/// the value times the factor exactly (checked at rational points of ε and
/// z, where both sides are exact).
#[test]
fn one_over_one_plus_z_eps_is_exact_in_retain_without_insertions() {
    let retain = site(RustRedMasterMode::Retain);
    let factor = atom(&format!("(1+caller::z*{EPS})^-1"));
    for name in ["wh-k1-root", "wh-k2-child"] {
        let v = record(name);
        let plain = value(&v, &retain);
        let got = closure(&v, &retain, &(v.prefactor.total() * &factor))
            .unwrap()
            .value;
        let at = |a: &Atom| {
            a.replace(Atom::var(eps()))
                .with(atom("3/7"))
                .replace(z())
                .with(atom("5/11"))
                .expand()
        };
        assert_same(
            &at(&got),
            &at(&(&plain * &factor)),
            &format!("{name}: Retain 1/(1+z eps) at rational points"),
        );
    }
}

/// (b) one record evaluated with half its argument, twice in one call, sums
/// to the record's value (the entry sums its records, v1.4 §7 via the
/// closure entry's ClosureValue).
#[test]
fn two_halves_sum_to_the_whole() {
    let v = record("wh-k2-child-ins");
    let whole = value(&v, &symbolic());
    let half = v.prefactor.total() * atom("1/2");
    let got = evaluate_factorized_closure(
        &engine(),
        &symbolic(),
        &[v.clone(), v.clone()],
        &[half.clone(), half],
    )
    .unwrap();
    assert_same(&got.value, &whole, "two halves");
}

/// (d) an insertion with coefficient 1 (pole order 0) equals the record
/// with the insertion's factor merged into its vertex (v1.4 §5: the
/// insertion is a local factor at its vertex; the contracted route of the
/// assembler's tests).
#[test]
fn a_unit_insertion_equals_its_factor_on_the_vertex() {
    let base = record("wh-k2-child-ins");
    let mut unit = base.clone();
    unit.insertions[0].coefficient = Atom::num(1);
    let mut merged = base.clone();
    let insertion = merged.insertions.remove(0);
    let vertex = merged
        .vertices
        .iter_mut()
        .find(|x| x.id == insertion.vertex)
        .unwrap();
    vertex.factor = &vertex.factor * &insertion.factor;
    vertex.rule = "V_merged".into();
    let a = value(&unit, &symbolic());
    let b = value(&merged, &symbolic());
    let top = top_order(&a).min(top_order(&b));
    assert_same(
        &truncate(&a, top),
        &truncate(&b, top),
        "unit insertion vs merged factor",
    );
}

/// (e) records in another order in one call: the per-record values and the
/// sum are unchanged.
#[test]
fn record_order_is_irrelevant() {
    let a = record("wh-k1-root");
    let b = record("wh-k2-child");
    let ab = evaluate_factorized_closure(
        &engine(),
        &symbolic(),
        &[a.clone(), b.clone()],
        &[a.prefactor.total(), b.prefactor.total()],
    )
    .unwrap();
    let ba = evaluate_factorized_closure(
        &engine(),
        &symbolic(),
        &[b.clone(), a.clone()],
        &[b.prefactor.total(), a.prefactor.total()],
    )
    .unwrap();
    assert_same(&ab.value, &ba.value, "sum");
    assert_same(&ab.records[0].value, &ba.records[1].value, "record a");
}

/// (e) two insertions in another order (the two poled children of the
/// assembler's L2 test, built from k1-root's vertices is not available
/// here: the producer's k5-stress records carry four; vendored: none with
/// two) -- covered by the external check below.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn insertion_order_is_irrelevant_k5_stress() {
    let dir = wh_fixtures::fixture_dir("k5-stress");
    let records = wh_fixtures::record_dirs(&dir);
    let input = wh_fixtures::read_fixture(&records[0].join("input.txt"));
    let v = FactorizedNumerator::from_text(&input).unwrap();
    assert!(v.insertions.len() >= 2);
    let mut r = v.clone();
    r.insertions.reverse();
    // the format orders insertions by vertex: another order is refused
    // (never a different value)
    match r.validate() {
        Err(FactorizedError::Insertion { detail, .. }) => {
            eprintln!("CONTRACT insertion order: refused: {detail}")
        }
        other => panic!("reversed insertions: {other:?}"),
    }
}

/// (e) relabelling the record's DUMMY indices (a contracted pair in a vertex
/// factor, `edge(e, k)` labels) and renumbering an INTERNAL edge: the value
/// is invariant.
#[test]
fn relabelling_dummies_and_internal_edges_is_irrelevant() {
    let v = record("wh-k1-root");
    let plain = value(&v, &symbolic());
    // the edge ids: the last internal edge 5 (loop edge of loop 2) -> 50
    let text = v.to_text().unwrap();
    let renumbered = text
        .replace("edge\t5\t", "edge\t50\t")
        .replace("edge_factor\t5\t", "edge_factor\t50\t")
        .replace("loop\t2\tedge=5", "loop\t2\tedge=50")
        .replace(
            "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q(5,",
            "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q(50,",
        );
    assert_ne!(renumbered, text);
    let r = FactorizedNumerator::from_text(&renumbered).unwrap_or_else(|e| panic!("{e}"));
    assert_same(&value(&r, &symbolic()), &plain, "internal edge 5 -> 50");
    // dummy labels: every edge(e,k) label k -> k + 100
    let relabelled = text.replace("index-label:e}::edge(", "index-label:e}::edge(9,");
    if relabelled != text {
        let r = FactorizedNumerator::from_text(&relabelled);
        match r {
            Ok(r) => assert_same(&value(&r, &symbolic()), &plain, "dummy labels"),
            Err(e) => eprintln!("CONTRACT relabel: refused: {e}"),
        }
    }
}

// ===========================================================================
// P3: the same function written differently
// ===========================================================================

/// An insertion coefficient c ε^-1 written as -2 c/(D-4) and as c (D-4+2ε+1)/ε
/// (all equal after D -> 4-2ε, v1.4 §5 / corr. 4): the same value.
#[test]
fn equal_coefficients_give_equal_values() {
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    let c = base.insertions[0].coefficient.clone();
    let c_over_eps = (&c * Atom::var(eps())).expand(); // c without its eps^-1
    for form in [
        format!("-2*({DIM}-4)^-1"),
        format!("({DIM}-4+2*{EPS}+1)*{EPS}^-1"),
        format!("(1+{EPS})^2*{EPS}^-1*(1+2*{EPS}+{EPS}^2)^-1"),
        format!("({EPS}*(1+{EPS})^-2)^-1*(1+{EPS})^-2"),
    ] {
        let mut v = base.clone();
        v.insertions[0].coefficient = &c_over_eps * atom(&form);
        match closure(&v, &symbolic(), &v.prefactor.total()) {
            Ok(got) => assert_same(&got.value, &plain, &format!("coefficient form {form}")),
            Err(e) => panic!("coefficient form {form}: refused: {}", detail(&e)),
        }
    }
}

/// The same forms of a pole-free factor 1 in the prefactor ARGUMENT
/// ((D-4+2ε)/... is undefined; here forms of 1): equal values.
#[test]
fn equal_arguments_give_equal_values() {
    let v = record("wh-k1-root");
    let plain = value(&v, &symbolic());
    for form in [
        format!("(-2*{EPS})*({DIM}-4)^-1"),
        format!("(1+{EPS})*(1+{EPS})^-1"),
        format!("({DIM}-2)*(2-2*{EPS})^-1"),
    ] {
        let got = closure(&v, &symbolic(), &(v.prefactor.total() * atom(&form)))
            .unwrap_or_else(|e| panic!("{form}: {e}"))
            .value;
        assert_same(&got, &plain, &format!("argument form {form}"));
    }
}

/// The text form and the typed form of every vendored record: to_text ->
/// from_text gives the same value.
#[test]
fn text_round_trip_gives_the_same_value() {
    for name in ["wh-k1-root", "wh-k2-child-ins", "wh-k3-hs"] {
        let v = record(name);
        let r = FactorizedNumerator::from_text(&v.to_text().unwrap()).unwrap();
        assert_same(&value(&r, &symbolic()), &value(&v, &symbolic()), name);
    }
}

/// vakint's epsilon `vakint::ε` in a record evaluated with the settings'
/// epsilon `gammalooprs::ε`: refused (corr. 3, "accepted only when it is the
/// settings' epsilon"), never a different value.
#[test]
fn vakint_epsilon_under_another_settings_epsilon_is_refused() {
    let base = record("wh-k2-child-ins");
    let mut v = base.clone();
    v.insertions[0].coefficient = v.insertions[0]
        .coefficient
        .replace(Atom::var(eps()))
        .with(atom("vakint::ε"));
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("CONTRACT vakint::ε: refused: {}", detail(&e)),
        Ok(got) => panic!("vakint::ε accepted under gammalooprs::ε: {}", got.value),
    }
}

// ===========================================================================
// P4: the auxiliary and reference masses
// ===========================================================================

/// The reference mass R inside a FUNCTION of a vertex factor (`log(R)`):
/// v1.4 §4.1 (R -> (1-λ²) M_aux² "everywhere in the integrand") and the
/// namespace rule (no vakint symbol in the value, corr. 1): either refused
/// or a value free of vakint's internal symbols (the Taylor parameter).
#[test]
fn the_reference_mass_inside_a_function() {
    let base = record("wh-k1-root");
    let mut v = base.clone();
    let vertex = v
        .vertices
        .iter_mut()
        .find(|x| x.factor != Atom::num(1))
        .unwrap();
    vertex.factor = &vertex.factor * atom(&format!("log({R})"));
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("CONTRACT log(R): refused: {}", detail(&e)),
        Ok(got) => {
            let text = got.value.to_string();
            assert!(
                !text.contains("vakint::") && !text.contains("closure_lambda"),
                "log(R) in a factor: the value carries vakint's internal symbol: {text}"
            );
        }
    }
}

/// An ε-dependent auxiliary mass M_aux² (1+ε): v1.4 §7 ("M_aux^2 powers and
/// log(M_aux^2) as split_auxiliary_mass expects") presumes a plain symbol;
/// expected: refused, or a value equal to the plain-symbol value with the
/// substitution M_aux² -> M_aux² (1+ε) re-expanded (never a different one).
#[test]
fn an_epsilon_dependent_auxiliary_mass() {
    let base = record("wh-k1-root");
    let plain = value(&base, &symbolic());
    let mut v = base.clone();
    let shifted = atom(&format!("{MAUX}*(1+{EPS})"));
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = shifted.clone();
    }
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("CONTRACT M_aux(1+eps): refused: {}", detail(&e)),
        Ok(got) => {
            eprintln!("CONTRACT M_aux(1+eps): accepted: {}", got.value);
            let top = top_order(&plain);
            // the plain value with M -> M(1+eps): log(M(1+eps)) = log M + log(1+eps)
            let m = atom(MAUX);
            let substituted = plain
                .replace(atom(&format!("log({MAUX})")))
                .with(atom(&format!("log({MAUX})+{EPS}-{EPS}^2/2+{EPS}^3/3")))
                .replace(m.clone())
                .with(shifted.clone());
            assert_same(
                &truncate(&got.value, top),
                &truncate(&substituted.expand(), top),
                "M_aux(1+eps)",
            );
        }
    }
}

/// An auxiliary mass that is ZERO only after D -> 4-2ε (M_aux² (D-4+2ε)):
/// expected: refused (a zero auxiliary mass is refused, `the auxiliary mass
/// is zero`), never evaluated.
#[test]
fn an_auxiliary_mass_zero_after_the_substitution() {
    let mut v = record("wh-k1-root");
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = atom(&format!("{MAUX}*({DIM}-4+2*{EPS})"));
    }
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!(
            "CONTRACT M_aux zero after D->4-2eps: refused: {}",
            detail(&e)
        ),
        Ok(got) => panic!(
            "an auxiliary mass zero after D -> 4-2eps is accepted: {}",
            got.value
        ),
    }
}

/// A D-dependent auxiliary mass M_aux² (D-2)/2 (not zero): expected refused
/// or the substituted value; recorded.
#[test]
fn a_dimension_dependent_auxiliary_mass() {
    let mut v = record("wh-k1-root");
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = atom(&format!("{MAUX}*({DIM}-2)/2"));
    }
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("CONTRACT M_aux (D-2)/2: refused: {}", detail(&e)),
        Ok(got) => panic!("a D-dependent auxiliary mass is accepted: {}", got.value),
    }
}

/// A reference mass written as a product (2 R, not a symbol) is refused
/// (`the reference mass ... is not a symbol`); a reference mass that is a
/// different symbol than the closing edge's is refused.
#[test]
fn reference_mass_not_a_symbol_or_two_masses() {
    let base = record("wh-k1-root");
    let text = base.to_text().unwrap();
    let doubled = text.replace(R, &format!("(2*{R})"));
    match FactorizedNumerator::from_text(&doubled) {
        Err(e) => eprintln!("CONTRACT 2R: refused: {e}"),
        Ok(v) => panic!("2R as the reference mass accepted: {:?}", v.kind),
    }
    let mut v = base.clone();
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.mass_squared = atom("gammalooprs::{scalar}::other_mass");
    }
    assert!(
        v.validate().is_err(),
        "a reference mass other than the closing edge's accepted"
    );
}

// ===========================================================================
// P5: text quirks of from_text / parse_text for marked records
// ===========================================================================

/// The outcome of a text variant: the same value as the clean form, a
/// refusal, or a DIFFERENT value (a defect).
fn text_variant(clean: &str, variant: &str, what: &str) -> Result<(), String> {
    let base = FactorizedNumerator::from_text(clean).unwrap();
    let plain = value(&base, &symbolic());
    match FactorizedNumerator::from_text(variant) {
        Err(e) => {
            eprintln!("CONTRACT {what}: refused: {e}");
            Ok(())
        }
        Ok(v) => match closure(&v, &symbolic(), &v.prefactor.total()) {
            Err(e) => {
                eprintln!("CONTRACT {what}: read, refused by the entry: {e}");
                Ok(())
            }
            Ok(got) => {
                if (&got.value - &plain).expand().is_zero() {
                    eprintln!("CONTRACT {what}: accepted, same value");
                    Ok(())
                } else {
                    Err(format!(
                        "{what}: ACCEPTED-BUT-DIFFERENT\n got {}\n clean {plain}",
                        got.value
                    ))
                }
            }
        },
    }
}

#[test]
fn line_ends_bom_trailing_whitespace_final_newline() {
    let clean = fixture_text("wh-k2-child-ins.record.txt");
    let mut mismatches = Vec::new();
    for (what, variant) in [
        ("CRLF", clean.replace('\n', "\r\n")),
        ("BOM", format!("\u{feff}{clean}")),
        ("trailing spaces", clean.replace('\n', "  \n")),
        ("trailing tabs", clean.replace('\n', "\t\n")),
        ("no final newline", clean.trim_end_matches('\n').to_string()),
        ("an extra empty final line", format!("{clean}\n")),
    ] {
        if let Err(f) = text_variant(&clean, &variant, what) {
            mismatches.push(f);
        }
    }
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n\n"));
}

/// Other Unicode spellings of ε in the insertion coefficient (U+03F5 ϵ,
/// U+1D700 𝜀, U+1D716 𝜖): not the settings' epsilon `gammalooprs::ε`
/// (U+03B5): refused (v1.4 §8 "an eps symbol other than the settings'
/// epsilon_symbol: Insertion / Evaluation"), never a constant.
#[test]
fn other_spellings_of_epsilon() {
    let clean = fixture_text("wh-k2-child-ins.record.txt");
    let mut mismatches = Vec::new();
    for spelling in ['\u{03f5}', '\u{1d700}', '\u{1d716}'] {
        let line = clean
            .lines()
            .find(|l| l.starts_with("insertion\t"))
            .unwrap();
        let changed = line.replace(
            "gammalooprs::{}::ε",
            &format!("gammalooprs::{{}}::{spelling}"),
        );
        assert_ne!(changed, line);
        let variant = clean.replace(line, &changed);
        if let Err(f) = text_variant(
            &clean,
            &variant,
            &format!("epsilon spelled U+{:04X}", spelling as u32),
        ) {
            mismatches.push(f);
        }
    }
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n\n"));
}

/// Other spellings of 𝜋 / 𝑖: `π` (U+03C0) and `i` are not Symbolica's
/// built-ins (corr. 3): read into vakint:: and refused, never a different
/// constant.
#[test]
fn other_spellings_of_pi_and_i() {
    let clean = fixture_text("wh-k2-child-ins.record.txt");
    let line = clean
        .lines()
        .find(|l| l.starts_with("insertion\t"))
        .unwrap();
    let mut mismatches = Vec::new();
    for (what, from, to) in [
        ("π U+03C0", "𝜋", "π"),
        ("bare i", "𝑖", "i"),
        ("pi alias", "𝜋", "pi"),
    ] {
        let variant = clean.replace(line, &line.replace(from, to));
        if let Err(f) = text_variant(&clean, &variant, what) {
            mismatches.push(f);
        }
    }
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n\n"));
}

/// Duplicated / reordered sections and numeric-field spellings: refused, or
/// accepted with the same value.
#[test]
fn sections_and_numbers() {
    let clean = fixture_text("wh-k1-root.record.txt");
    let mut mismatches = Vec::new();
    let depth_line = "depth\t2";
    assert!(clean.contains(depth_line));
    let variants = vec![
        ("depth +2", clean.replace(depth_line, "depth\t+2")),
        ("depth 02", clean.replace(depth_line, "depth\t02")),
        ("depth 2e0", clean.replace(depth_line, "depth\t2e0")),
        (
            "power=+1",
            clean.replacen("power=1\tmomentum=1,0,0", "power=+1\tmomentum=1,0,0", 1),
        ),
        (
            "power=01",
            clean.replacen("power=1\tmomentum=1,0,0", "power=01\tmomentum=1,0,0", 1),
        ),
        (
            "momentum +1",
            clean.replacen("momentum=1,0,0\t", "momentum=+1,0,0\t", 1),
        ),
        (
            "duplicated depth line",
            clean.replace(depth_line, &format!("{depth_line}\n{depth_line}")),
        ),
        ("prefactor_colour before prefactor_scalars", {
            let s = clean
                .lines()
                .find(|l| l.starts_with("prefactor_scalars\t"))
                .unwrap();
            let c = clean
                .lines()
                .find(|l| l.starts_with("prefactor_colour\t"))
                .unwrap();
            clean.replace(&format!("{s}\n{c}"), &format!("{c}\n{s}"))
        }),
        ("a huge integer literal in prefactor_scalars", {
            let s = clean
                .lines()
                .find(|l| l.starts_with("prefactor_scalars\t"))
                .unwrap();
            clean.replace(s, &format!("{s}*100000000000000000000000000000000000000/100000000000000000000000000000000000000"))
        }),
        ("a leading + in prefactor_scalars", {
            let s = clean
                .lines()
                .find(|l| l.starts_with("prefactor_scalars\t"))
                .unwrap();
            clean.replace(s, &s.replace("prefactor_scalars\t", "prefactor_scalars\t+"))
        }),
    ];
    for (what, variant) in variants {
        if let Err(f) = text_variant(&clean, &variant, what) {
            mismatches.push(f);
        }
    }
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n\n"));
}

/// a non-ASCII name in a namespace-qualified caller symbol is refused at
/// the entry (outside ASCII only the settings' epsilon and scale and
/// Symbolica's built-ins), never read as a constant; the ASCII twin is
/// accepted.
#[test]
fn non_ascii_names() {
    let clean = fixture_text("wh-k1-root.record.txt");
    let variant = clean.replace("UFO::{}::gs_minus", "UFO::{}::gs_mïnüs");
    assert_ne!(variant, clean);
    let v = FactorizedNumerator::from_text(&variant).unwrap();
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(FactorizedError::Graph { detail }) => assert!(
            detail.starts_with("a non-ASCII name UFO::gs_mïnüs in vertex "),
            "{detail}"
        ),
        other => panic!("non-ASCII name: {other:?}"),
    }
    let ascii = FactorizedNumerator::from_text(
        &clean.replace("UFO::{}::gs_minus", "UFO::{}::gs_minus_ascii"),
    )
    .unwrap();
    closure(&ascii, &symbolic(), &ascii.prefactor.total()).unwrap();
}

/// an epsilon of the right spelling in ANOTHER namespace (`caller::ε`) in
/// a vertex factor is not the settings' epsilon: refused at the entry; the
/// settings' own `gammalooprs::ε` there is accepted (and the value differs
/// from the plain one only by that factor).
#[test]
fn an_epsilon_of_another_namespace_is_refused() {
    let base = record("wh-k1-root");
    let mut v = base.clone();
    let vertex = v
        .vertices
        .iter_mut()
        .find(|x| x.factor != Atom::num(1))
        .unwrap();
    vertex.factor = &vertex.factor * atom("1+caller::ε");
    let id = vertex.id;
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(FactorizedError::Graph { detail }) => assert!(
            detail.starts_with(&format!(
                "a non-ASCII name caller::ε in vertex {id} factor: "
            )),
            "{detail}"
        ),
        other => panic!("caller::ε: {other:?}"),
    }
    let mut ok = base.clone();
    let vertex = ok
        .vertices
        .iter_mut()
        .find(|x| x.factor != Atom::num(1))
        .unwrap();
    vertex.factor = &vertex.factor * atom(&format!("1+{EPS}"));
    closure(&ok, &symbolic(), &ok.prefactor.total()).unwrap();
}

/// The direct / ports entry for a marked K5 direct record: z in the argument
/// scales the value (P2 for the direct kind).
#[test]
fn direct_record_z_in_the_argument() {
    let v = record("wh-k5-direct");
    let vakint = engine();
    // the record's depth is the number of terms
    let settings = VakintSettings {
        number_of_terms_in_epsilon_expansion: v.record.as_ref().unwrap().depth.unwrap() as i64,
        ..symbolic()
    };
    let plain = evaluate_factorized_ports(&vakint, &settings, &v, &v.prefactor.total())
        .unwrap()
        .value;
    let got = evaluate_factorized_ports(&vakint, &settings, &v, &(v.prefactor.total() * z()))
        .unwrap()
        .value;
    assert_same(&got, &(&plain * z()), "k5-direct argument z");
    let mut i = v.clone();
    i.insertions[0].coefficient = &i.insertions[0].coefficient * z();
    let got = evaluate_factorized_ports(&vakint, &settings, &i, &i.prefactor.total())
        .unwrap()
        .value;
    assert_same(&got, &(&plain * z()), "k5-direct insertion z");
}

/// P3: equal forms of 1 multiplied into a VERTEX FACTOR (FactorRules admits
/// tensor-free scalars rational in D and ε): the same value.
#[test]
fn equal_factor_forms_give_equal_values() {
    let v = record("wh-k1-root");
    let plain = value(&v, &symbolic());
    for form in [
        format!("(1+{EPS})*(1+{EPS})^-1"),
        format!("({DIM}-2)*(2-2*{EPS})^-1"),
        format!("(-2*{EPS})*({DIM}-4)^-1"),
    ] {
        let mut f = v.clone();
        let vertex = f
            .vertices
            .iter_mut()
            .find(|x| x.factor != Atom::num(1))
            .unwrap();
        vertex.factor = &vertex.factor * atom(&form);
        match closure(&f, &symbolic(), &f.prefactor.total()) {
            Ok(got) => assert_same(&got.value, &plain, &format!("factor form {form}")),
            Err(e) => panic!("factor form {form}: refused: {}", detail(&e)),
        }
    }
}

/// v1.4 §7 "value: ... in the caller's notation" + corr. 1: no vakint
/// symbol in the returned value of any kind.
#[test]
fn every_kind_returns_the_callers_notation() {
    for name in CLOSURE_KINDS {
        let got = value(&record(name), &symbolic());
        let mut bad = None;
        got.as_view().visitor(&mut |part| {
            let name = match part {
                AtomView::Var(x) => x.get_symbol().get_name().to_string(),
                AtomView::Fun(f) => f.get_symbol().get_name().to_string(),
                _ => return true,
            };
            if name.starts_with("vakint::") && name != "vakint::ε" {
                bad = Some(name);
            }
            bad.is_none()
        });
        assert_eq!(bad, None, "{name}: the value carries a vakint symbol");
    }
}

/// v1.4 §7 "linear spectators (C8): ... beta_projection::transverse ...
/// returned linear, never expanded into": the K1 value is homogeneous of
/// degree 1 in the marker.
#[test]
fn linear_spectator_marker() {
    let got = value(&record("wh-k1-root"), &symbolic());
    let marker = atom("beta_projection::{}::transverse");
    let doubled = got.replace(marker.clone()).with(Atom::num(2) * &marker);
    assert_same(&doubled, &(&got * Atom::num(2)), "marker linear");
    let without = got.replace(marker).with(Atom::Zero);
    assert!(is_zero(&without), "a term without the marker: {without}");
}

/// P2 for K6 (the producer's k6-norm records, direct, through the ports
/// entry): z in the argument scales the value.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn k6_z_in_the_argument() {
    let dir = wh_fixtures::fixture_dir("k6-norm");
    let mut compared = 0;
    for path in std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().path()) {
        let input = if path.is_dir() {
            path.join("input.txt")
        } else {
            path.clone()
        };
        if !input.is_file() {
            continue;
        }
        let text = wh_fixtures::read_fixture(&input);
        if !text.starts_with("schema\tpint6.factorized-numerator.v3") {
            continue;
        }
        let v = FactorizedNumerator::from_text(&text)
            .unwrap_or_else(|e| panic!("{}: {e}", input.display()));
        let vakint = engine();
        let plain = evaluate_factorized_ports(&vakint, &symbolic(), &v, &v.prefactor.total())
            .unwrap_or_else(|e| panic!("{}: {e}", input.display()))
            .value;
        let got = evaluate_factorized_ports(&vakint, &symbolic(), &v, &(v.prefactor.total() * z()))
            .unwrap()
            .value;
        assert_same(&got, &(&plain * z()), &input.display().to_string());
        compared += 1;
    }
    assert!(compared > 0, "no k6 record");
    eprintln!("CONTRACT K6 compared {compared}");
}

/// R under a NEGATIVE or a non-integer power in a factor
/// (`1/R`, `sqrt(R)`: FactorRules admits tensor-free scalars): refused, or a
/// value free of vakint's internal symbols.
#[test]
fn the_reference_mass_under_a_power() {
    let base = record("wh-k1-root");
    let mut leaks = Vec::new();
    for form in [
        format!("{R}^-1"),
        format!("{R}^(1/2)"),
        format!("({R}+{MAUX})^-1"),
    ] {
        let mut v = base.clone();
        let vertex = v
            .vertices
            .iter_mut()
            .find(|x| x.factor != Atom::num(1))
            .unwrap();
        vertex.factor = &vertex.factor * atom(&form);
        match closure(&v, &symbolic(), &v.prefactor.total()) {
            Err(e) => eprintln!("CONTRACT {form}: refused: {}", detail(&e)),
            Ok(got) => {
                let text = got.value.to_string();
                eprintln!("CONTRACT {form}: accepted: {text}");
                if text.contains("closure_lambda") {
                    leaks.push(form);
                }
            }
        }
    }
    assert!(
        leaks.is_empty(),
        "vakint::closure_lambda in the value for {leaks:?}"
    );
}

/// P3: the same record evaluated with the settings' epsilon under ANOTHER
/// NAME (vakint's default `ε` = `vakint::ε`, the record's ε respelled to
/// it): the same value up to the renaming of ε.
#[test]
fn the_settings_epsilon_under_another_name() {
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    let ours = atom("vakint::ε");
    let mut v = base.clone();
    v.insertions[0].coefficient = v.insertions[0]
        .coefficient
        .replace(Atom::var(eps()))
        .with(ours.clone());
    let settings = VakintSettings {
        epsilon_symbol: "ε".into(),
        ..symbolic()
    };
    let got = closure(&v, &settings, &v.prefactor.total())
        .unwrap_or_else(|e| panic!("{}", detail(&e)))
        .value
        .replace(ours)
        .with(Atom::var(eps()));
    assert_same(&got, &plain, "epsilon renamed");
}

/// P2 (c): two insertions with coefficients a and b, and the same record
/// with a*b on the first and 1 on the second: the same value (the
/// coefficients are scalars of the integrand; the depth raise is the sum of
/// the pole orders either way). A producer k5-stress record (four
/// insertions), the smallest one.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn moving_a_coefficient_between_insertions_k5_stress() {
    let dir = wh_fixtures::fixture_dir("k5-stress");
    let mut records = wh_fixtures::record_dirs(&dir)
        .into_iter()
        .map(|d| wh_fixtures::read_fixture(&d.join("input.txt")))
        .map(|t| FactorizedNumerator::from_text(&t).unwrap())
        .filter(|v| v.loops() <= 4)
        .collect::<Vec<_>>();
    records.sort_by_key(|v| v.to_text().unwrap().len());
    let v = records.remove(0);
    assert!(v.insertions.len() >= 2);
    let plain = value(&v, &symbolic());
    let mut moved = v.clone();
    let b = moved.insertions[1].coefficient.clone();
    moved.insertions[0].coefficient = &moved.insertions[0].coefficient * &b;
    moved.insertions[1].coefficient = Atom::num(1);
    assert_same(&value(&moved, &symbolic()), &plain, "a*b on one insertion");
}

/// P1 Numerical mode with an INSERTION (k2-child-ins): the same orders and
/// the Symbolic coefficients numerically.
#[test]
fn numerical_mode_with_an_insertion() {
    let v = record("wh-k2-child-ins");
    let exact = value(&v, &symbolic());
    let numerical = value(&v, &site(RustRedMasterMode::Numerical));
    let e = wh_fixtures::epsilon_coefficients(&exact, eps());
    let n = wh_fixtures::epsilon_coefficients(&numerical, eps());
    assert_eq!(
        e.keys().collect::<Vec<_>>(),
        n.keys().collect::<Vec<_>>(),
        "orders"
    );
    let pi = atom("𝜋");
    let pi_value = atom("3.14159265358979323846264338327950288419716939937510582097494459");
    for (order, exact_c) in &e {
        let difference = (exact_c
            .replace(pi.clone())
            .with(pi_value.clone())
            .to_float(60)
            - n[order]
                .replace(pi.clone())
                .with(pi_value.clone())
                .to_float(60))
        .expand();
        let terms = match difference.as_view() {
            AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
            _ if difference.is_zero() => Vec::new(),
            _ => vec![difference.clone()],
        };
        for term in terms {
            let magnitude = parse_magnitude(&term.to_string());
            assert!(magnitude < 1e-12, "eps^{order}: {term}");
        }
    }
}

// ===========================================================================
// INVENTORY for the result backstop: the symbols a real result contains that
// are NOT among its inputs (record atoms, prefactor argument, settings).
// ===========================================================================

fn symbols_of(a: &Atom, out: &mut std::collections::BTreeSet<String>) {
    a.as_view().visitor(&mut |part| {
        match part {
            AtomView::Var(x) => {
                out.insert(x.get_symbol().get_name().to_string());
            }
            AtomView::Fun(f) => {
                out.insert(f.get_symbol().get_name().to_string());
            }
            _ => {}
        }
        true
    });
}

fn input_symbols(v: &FactorizedNumerator, argument: &Atom) -> std::collections::BTreeSet<String> {
    let mut s = std::collections::BTreeSet::new();
    for e in &v.edges {
        symbols_of(&e.factor, &mut s);
        symbols_of(&e.mass_squared, &mut s);
    }
    for x in &v.vertices {
        symbols_of(&x.factor, &mut s);
    }
    for p in &v.ports {
        symbols_of(p, &mut s);
    }
    if let Some(p) = &v.projector {
        symbols_of(&p.tensor, &mut s);
    }
    symbols_of(&v.prefactor.overall, &mut s);
    symbols_of(&v.prefactor.scalars, &mut s);
    symbols_of(&v.prefactor.colour, &mut s);
    for i in &v.insertions {
        symbols_of(&i.coefficient, &mut s);
        symbols_of(&i.factor, &mut s);
    }
    if let FactorizedKind::Closure(spec) = &v.kind {
        symbols_of(&spec.mass_squared, &mut s);
        symbols_of(&spec.auxiliary_mass_squared, &mut s);
    }
    symbols_of(argument, &mut s);
    s.insert("gammalooprs::ε".into());
    s.insert("gammalooprs::μᵣ²".into());
    s
}

fn inventory(
    records: &[FactorizedNumerator],
    closure_kind: bool,
) -> std::collections::BTreeMap<String, std::collections::BTreeSet<String>> {
    let vakint = engine();
    let mut out = std::collections::BTreeMap::<String, std::collections::BTreeSet<String>>::new();
    for (mode_name, mode) in [
        ("symbolic", RustRedMasterMode::Symbolic),
        ("numerical", RustRedMasterMode::Numerical),
        ("retain", RustRedMasterMode::Retain),
    ] {
        let settings = site(mode);
        for v in records {
            let argument = v.prefactor.total();
            let result = if closure_kind {
                evaluate_factorized_closure(
                    &vakint,
                    &settings,
                    std::slice::from_ref(v),
                    std::slice::from_ref(&argument),
                )
                .map(|r| r.value)
            } else {
                evaluate_factorized_ports(&vakint, &settings, v, &argument).map(|r| r.value)
            };
            let result = match result {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("INVENTORY {mode_name} {}: ERROR {e}", v.label);
                    continue;
                }
            };
            let inputs = input_symbols(v, &argument);
            let mut found = std::collections::BTreeSet::new();
            symbols_of(&result, &mut found);
            for name in found.difference(&inputs) {
                out.entry(format!("{mode_name} {name}"))
                    .or_default()
                    .insert(v.label.split_whitespace().next().unwrap_or("").to_string());
            }
            if inputs.contains("gammalooprs::dim") && found.contains("gammalooprs::dim") {
                out.entry(format!("{mode_name} DIMENSION gammalooprs::dim"))
                    .or_default()
                    .insert(v.label.clone());
            }
        }
    }
    for (k, v) in &out {
        eprintln!(
            "INVENTORY {k} <- {} records, e.g. {:?}",
            v.len(),
            v.iter().next()
        );
    }
    out
}

#[test]
fn inventory_vendored_records() {
    let closures = CLOSURE_KINDS.iter().map(|n| record(n)).collect::<Vec<_>>();
    inventory(&closures, true);
    inventory(&[record("wh-k5-direct")], false);
}

#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn inventory_producer_records() {
    let root = wh_fixtures::fixture_root();
    let mut closures = Vec::new();
    let mut directs = Vec::new();
    let only = std::env::var("VAKINT_TEST_INVENTORY_SETS")
        .unwrap_or_else(|_| "k1-quick,k5-stress,k6-norm".into());
    for set in only.split(',') {
        let dir = root.join(set);
        let mut inputs = Vec::new();
        for d in wh_fixtures::record_dirs(&dir) {
            inputs.push(d.join("input.txt"));
        }
        if set == "k6-norm" {
            for e in std::fs::read_dir(&dir).unwrap() {
                let p = e.unwrap().path();
                if p.is_file() && p.extension().is_some_and(|x| x == "txt") {
                    inputs.push(p);
                }
            }
        }
        for p in inputs {
            let text = wh_fixtures::read_fixture(&p);
            if !text.starts_with("schema\tpint6.factorized-numerator.v3") {
                continue;
            }
            let v = FactorizedNumerator::from_text(&text).unwrap();
            if matches!(v.kind, FactorizedKind::Closure(_)) {
                closures.push(v);
            } else {
                directs.push(v);
            }
        }
    }
    eprintln!(
        "INVENTORY producer: {} closure, {} direct",
        closures.len(),
        directs.len()
    );
    inventory(&closures, true);
    inventory(&directs, false);
}

/// THE RESULT CHECK ALONE catches  and  through the public closure
/// entry (with no validation rule for the masses' positions): an internal
/// `Evaluation` error naming the symbol, prefix `result check: `.
#[test]
fn the_result_check_catches_f1_and_f3_without_a_validation_rule() {
    let expect = |v: &FactorizedNumerator, what: &str, needle: &str| match closure(
        v,
        &symbolic(),
        &v.prefactor.total(),
    ) {
        Err(FactorizedError::Evaluation { detail }) => {
            assert!(
                detail.starts_with("result check: ") && detail.contains(needle),
                "{what}: {detail}"
            )
        }
        // a later validation rule may refuse first: then not Evaluation
        Err(e) => eprintln!("{what}: refused before the result check: {}", detail(&e)),
        Ok(got) => panic!("{what}: accepted: {}", got.value),
    };
    let base = record("wh-k1-root");
    for form in [format!("log({R})"), format!("{R}^-1"), format!("{R}^(1/2)")] {
        let mut v = base.clone();
        let vertex = v
            .vertices
            .iter_mut()
            .find(|x| x.factor != Atom::num(1))
            .unwrap();
        vertex.factor = &vertex.factor * atom(&form);
        expect(&v, &form, "vakint::closure_lambda");
    }
    let mut v = base.clone();
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = atom(&format!("{MAUX}*({DIM}-2)/2"));
    }
    expect(
        &v,
        "M_aux (D-2)/2",
        "the dimension gammalooprs::dim outside a representation slot",
    );
}

// ===========================================================================
// the masses of a marked closure record by position and form
// ===========================================================================

fn closure_refusal(v: &FactorizedNumerator, prefix: &str, what: &str) {
    match v.validate() {
        Err(FactorizedError::Closure { detail }) => {
            assert!(
                detail.starts_with(prefix),
                "{what}: {detail}\n(expected prefix {prefix:?})"
            )
        }
        other => panic!("{what}: expected Closure {prefix:?}, got {other:?}"),
    }
}

/// The reference mass in a factor or the projector (polynomially or in a
/// function), the auxiliary mass in a factor, the projector or a prefactor
/// field: refused at validation with the field named; the record without
/// it is accepted.
#[test]
fn the_masses_occur_only_in_their_positions() {
    let base = record("wh-k1-root");
    base.validate().unwrap();
    let vertex_id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    let in_vertex = |extra: &str| {
        let mut v = base.clone();
        let vertex = v.vertices.iter_mut().find(|x| x.id == vertex_id).unwrap();
        vertex.factor = &vertex.factor * atom(extra);
        v
    };
    let r_prefix = format!(
        "the reference mass {} in vertex {vertex_id} factor: it may occur only as the closure mass and as the closing edge's mass",
        atom(R)
    );
    for form in [
        R.to_string(),
        format!("log({R})"),
        format!("{R}^-1"),
        format!("{R}^(1/2)"),
        format!("({R}+{MAUX})^-1"),
    ] {
        closure_refusal(&in_vertex(&form), &r_prefix, &form);
    }
    let m_prefix = format!(
        "the auxiliary mass {} in vertex {vertex_id} factor: it may occur only as the taylor line's mass",
        atom(MAUX)
    );
    closure_refusal(&in_vertex(MAUX), &m_prefix, "M_aux in a factor");
    let mut p = base.clone();
    let projector = p.projector.as_mut().unwrap();
    projector.tensor = &projector.tensor * atom(R);
    closure_refusal(
        &p,
        &format!("the reference mass {} in projector: ", atom(R)),
        "R in the projector",
    );
    let mut s = base.clone();
    s.prefactor.scalars = &s.prefactor.scalars * atom(MAUX);
    closure_refusal(
        &s,
        &format!("the auxiliary mass {} in prefactor_scalars: ", atom(MAUX)),
        "M_aux in prefactor_scalars",
    );
    // the prefactor ARGUMENT
    match closure(&base, &symbolic(), &(base.prefactor.total() * atom(MAUX))) {
        Err(FactorizedError::NotScalar { location, detail }) => {
            assert_eq!(location, "prefactor argument");
            assert!(detail.starts_with("the auxiliary mass "), "{detail}");
        }
        other => panic!("M_aux in the argument: {other:?}"),
    }
}

/// The auxiliary mass and an edge mass must be plain symbols (or 0): an
/// expression in ε or the dimension, a polynomial that is zero
/// (^2-m^2+6m-9), a product: refused at validation.
#[test]
fn masses_are_plain_symbols() {
    let base = record("wh-k1-root");
    let zero_polynomial = "(gammalooprs::m-3)^2-gammalooprs::m^2+6*gammalooprs::m-9";
    for form in [
        format!("{MAUX}*({DIM}-2)/2"),
        format!("{MAUX}*(1+{EPS})"),
        format!("{MAUX}*({DIM}-4+2*{EPS})"),
        zero_polynomial.to_string(),
        format!("2*{MAUX}"),
    ] {
        let mut v = base.clone();
        if let FactorizedKind::Closure(spec) = &mut v.kind {
            spec.auxiliary_mass_squared = atom(&form);
        }
        if atom(&form).expand().is_zero() && atom(&form).is_zero() {
            continue; // a literal zero is the zero-mass rule's
        }
        closure_refusal(&v, "the auxiliary mass ", &form);
        let detail = v.validate().unwrap_err().to_string();
        assert!(detail.contains("is not a plain symbol"), "{form}: {detail}");
    }
    // an edge mass written as the zero polynomial (it was classified massive)
    let mut v = base.clone();
    let edge = v
        .edges
        .iter_mut()
        .find(|e| !e.is_external() && e.mass_squared.is_zero())
        .unwrap();
    edge.mass_squared = atom(zero_polynomial);
    let id = edge.id;
    closure_refusal(
        &v,
        &format!("edge {id}: mass "),
        "zero polynomial edge mass",
    );
    assert!(
        v.validate()
            .unwrap_err()
            .to_string()
            .contains("is neither 0 nor a plain symbol")
    );
}

/// Every producer v3 record (all 75 under the wh-fixtures root, runs/
/// excluded) and every vendored one is ACCEPTED by validation (from_text
/// validates): the grammar rules refuse no real record.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn every_real_record_validates() {
    let mut files = Vec::new();
    let mut stack = vec![wh_fixtures::fixture_root()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            if p.is_dir() {
                if name != "runs" {
                    stack.push(p);
                }
            } else if name == "input.txt" || (name.contains("record") && name.ends_with(".txt")) {
                files.push(p);
            }
        }
    }
    let vendored =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/p_integral/factorized/fixtures");
    for e in std::fs::read_dir(&vendored).unwrap() {
        let p = e.unwrap().path();
        if p.file_name()
            .unwrap()
            .to_string_lossy()
            .ends_with(".record.txt")
        {
            files.push(p);
        }
    }
    init();
    let mut accepted = 0;
    let mut refused = Vec::new();
    for p in files {
        let text = wh_fixtures::read_fixture(&p);
        if !text.starts_with("schema\tpint6.factorized-numerator.v3") {
            continue;
        }
        match FactorizedNumerator::from_text(&text) {
            Ok(_) => accepted += 1,
            Err(e) => refused.push(format!("{}: {e}", p.display())),
        }
    }
    assert!(
        refused.is_empty(),
        "refused real records:\n{}",
        refused.join("\n")
    );
    wh_fixtures::assert_compared("REAL_RECORDS_VALIDATE", accepted, Some(85));
}

// ===========================================================================
// opaque caller scalars depend on nothing the evaluation transforms
// ===========================================================================

/// A function other than the structural heads, or a non-integer power, whose
/// argument / base / exponent depends on ε, the dimension, a mass or a
/// momentum, in a factor or a prefactor field: refused at validation
/// (`Graph`, prefix `a non-polynomial dependence in <field>: `); the same
/// function of a caller scalar is accepted.
#[test]
fn opaque_functions_and_powers_are_free_of_eps_dimension_masses_momenta() {
    let base = record("wh-k1-root");
    let vertex_id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    let with = |extra: &str| {
        let mut v = base.clone();
        let vertex = v.vertices.iter_mut().find(|x| x.id == vertex_id).unwrap();
        vertex.factor = &vertex.factor * atom(extra);
        v
    };
    let dependence = format!("a non-polynomial dependence in vertex {vertex_id} factor: ");
    let exponent = format!("a non-integer exponent in vertex {vertex_id} factor: ");
    for (bad, prefix) in [
        (format!("log(1+{EPS})"), &dependence),
        (format!("caller::f({DIM})"), &dependence),
        (format!("sin({DIM}*caller::x)"), &dependence),
        (format!("(1+{EPS})^(1/2)"), &exponent),
        (format!("caller::x^{EPS}"), &exponent),
    ] {
        match with(&bad).validate() {
            Err(FactorizedError::Graph { detail }) => {
                assert!(detail.starts_with(prefix.as_str()), "{bad}: {detail}")
            }
            other => panic!("{bad}: {other:?}"),
        }
    }
    for good in [
        format!("(1+{EPS})"),
        "caller::f(caller::x)".to_string(),
        "log(caller::x)".to_string(),
        "symbolica::zeta(3)".to_string(),
        format!("({DIM}-1)^-2"),
    ] {
        with(&good)
            .validate()
            .unwrap_or_else(|e| panic!("{good}: {e}"));
    }
    // a prefactor field and the prefactor ARGUMENT
    let mut s = base.clone();
    s.prefactor.scalars = &s.prefactor.scalars * atom(&format!("exp({EPS})"));
    match s.validate() {
        Err(FactorizedError::Graph { detail }) => assert!(
            detail.starts_with("a non-integer exponent in prefactor_scalars: the exponent "),
            "{detail}"
        ),
        other => panic!("exp(eps) in prefactor_scalars: {other:?}"),
    }
    match closure(
        &base,
        &symbolic(),
        &(base.prefactor.total() * atom(&format!("log({DIM})"))),
    ) {
        Err(FactorizedError::Graph { detail }) => assert!(
            detail.starts_with("a non-polynomial dependence in prefactor argument: "),
            "{detail}"
        ),
        other => panic!("log(D) in the argument: {other:?}"),
    }
    closure(
        &base,
        &symbolic(),
        &(base.prefactor.total() * atom("log(caller::x)")),
    )
    .unwrap();
}

// ===========================================================================
// The one vakint atom a marked record may carry: UNKNOWN(master, k)
// ===========================================================================

/// With allow-unknowns OFF (the default), a caller-supplied vakint UNKNOWN
/// atom reaches the result and vakint's unknown-master policy refuses the
/// result (`Evaluation`, the unknown-master text), as for an unknown master
/// of the integral itself: never a value.
#[test]
fn unknown_from_vakint_with_allow_unknowns_off_is_refused_by_the_policy() {
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::remove_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS") };
    let mut v = record("wh-k2-child-ins");
    v.insertions[0].coefficient = &v.insertions[0].coefficient
        * (Atom::num(1) + atom("vakint::UNKNOWN(vakint::ggppss_int_53,0)"));
    v.validate().unwrap();
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(FactorizedError::Evaluation { detail }) => assert!(
            detail.contains("expected master coefficient") && detail.contains("UNKNOWN"),
            "{detail}"
        ),
        other => panic!("{other:?}"),
    }
}

/// A caller UNKNOWN multiplying an order the record's depth CUTS is dropped
/// with that order (the returned product is truncated at ε^(depth-L-1)), as
/// any term there.
#[test]
fn unknown_at_a_cut_order_is_dropped_with_the_order() {
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1") };
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    let mut v = base.clone();
    v.insertions[0].coefficient = &v.insertions[0].coefficient
        * (Atom::num(1) + atom(&format!("vakint::UNKNOWN(vakint::ggppss_int_53,0)*{EPS}^3")));
    assert_same(
        &value(&v, &symbolic()),
        &plain,
        "UNKNOWN x eps^3 beyond the depth",
    );
}

/// Only the exact shape UNKNOWN(<master symbol>, <integer>) of vakint's own
/// head is the exception: lookalike heads, another arity, a non-integer
/// order stay reserved symbols; the caller's own `user::UNKNOWN` is a caller
/// scalar.
#[test]
fn only_vakint_unknown_in_its_exact_shape_is_allowed() {
    let base = record("wh-k2-child-ins");
    for bad in [
        "vakint::UNKNOWN2(vakint::m53,0)",
        "vakint::unknown(vakint::m53,0)",
        "vakint::UNKNOWN(vakint::m53,1/2)",
        "vakint::UNKNOWN(vakint::m53,0,1)",
        "vakint::UNKNOWN(vakint::m53+1,0)",
    ] {
        let mut v = base.clone();
        v.insertions[0].coefficient = &v.insertions[0].coefficient * atom(bad);
        match v.validate() {
            Err(FactorizedError::Graph { detail }) => {
                assert!(
                    detail.starts_with("reserved symbol vakint::"),
                    "{bad}: {detail}"
                )
            }
            other => panic!("{bad}: {other:?}"),
        }
    }
    for good in [
        "vakint::UNKNOWN(vakint::ggppss_int_53,0)",
        "vakint::UNKNOWN(vakint::ggppss_int_53,-2)",
        "user::UNKNOWN(user::m53,0)",
    ] {
        let mut v = base.clone();
        v.insertions[0].coefficient = &v.insertions[0].coefficient * atom(good);
        v.validate().unwrap_or_else(|e| panic!("{good}: {e}"));
    }
}

// ===========================================================================
// every exponent a caller writes is an integer literal
// ===========================================================================

/// P(t) = (t-3)^2 - t^2 + 6t - 9, identically zero but not a literal zero.
fn zero_polynomial(t: &str) -> String {
    format!("(({t}-3)^2-{t}^2+6*{t}-9)")
}

/// inputs: on the real K2 record, ε^-1 (P(x)+1)^(-sqrt(2)),
/// ε^-1 (P(x)+1)^(-π) and ε^-1 (P(x)+1)^(P(y)-1) in the insertion
/// coefficient (the last exponent is identically -1 but not WRITTEN as an
/// integer; the base is identically 1, so the zero test of a base passes
/// it): refused at validation, `Insertion`, prefix `a non-integer exponent in
/// insertion <v> coefficient: the exponent `; nothing is expanded to decide.
#[test]
fn exponents_identically_integer_but_not_written_so_are_refused() {
    let base = record("wh-k2-child-ins");
    let vertex = base.insertions[0].vertex;
    let p = format!("({}+1)", zero_polynomial("caller::x"));
    for exponent in [
        "-2^(1/2)".to_string(),
        "-𝜋".to_string(),
        format!("{}-1", zero_polynomial("caller::y")),
        "(caller::y-3)^2-caller::y^2+6*caller::y-10".to_string(),
    ] {
        let mut v = base.clone();
        v.insertions[0].coefficient = atom(&format!("{EPS}^-1*{p}^({exponent})"));
        match v.validate() {
            Err(FactorizedError::Insertion { vertex: at, detail }) => {
                assert_eq!(at, vertex);
                assert!(
                    detail.starts_with(&format!(
                        "a non-integer exponent in insertion {vertex} coefficient: the exponent "
                    )),
                    "{exponent}: {detail}"
                );
            }
            other => panic!("{exponent}: {other:?}"),
        }
    }
}

/// In EVERY caller field (edge and vertex factor, projector, prefactor
/// fields, insertion coefficient and factor, the entry's prefactor
/// argument): x^(1/2) and x^y of a caller scalar are refused; x^-3, and
/// x^(2-3) which the parser normalises to x^-1, are accepted.
#[test]
fn exponents_are_integer_literals_in_every_field() {
    let base = record("wh-k2-child-ins");
    base.validate().unwrap();
    let edge_id = base
        .edges
        .iter()
        .find(|e| !e.is_external() && e.factor != Atom::num(1))
        .unwrap()
        .id;
    let vertex_id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    let fields: Vec<(&str, Box<dyn Fn(&mut FactorizedNumerator, &Atom)>)> = vec![
        (
            "edge factor",
            Box::new(move |v, x| {
                let e = v.edges.iter_mut().find(|e| e.id == edge_id).unwrap();
                e.factor = &e.factor * x;
            }),
        ),
        (
            "vertex factor",
            Box::new(move |v, x| {
                let e = v.vertices.iter_mut().find(|e| e.id == vertex_id).unwrap();
                e.factor = &e.factor * x;
            }),
        ),
        (
            "prefactor_scalars",
            Box::new(|v, x| v.prefactor.scalars = &v.prefactor.scalars * x),
        ),
        (
            "prefactor_colour",
            Box::new(|v, x| v.prefactor.colour = &v.prefactor.colour * x),
        ),
        (
            "insertion coefficient",
            Box::new(|v, x| v.insertions[0].coefficient = &v.insertions[0].coefficient * x),
        ),
        (
            "insertion factor",
            Box::new(|v, x| v.insertions[0].factor = &v.insertions[0].factor * x),
        ),
    ];
    for (field, put) in &fields {
        for bad in ["caller::x^(1/2)", "caller::x^caller::y"] {
            let mut v = base.clone();
            put(&mut v, &atom(bad));
            let detail = v
                .validate()
                .expect_err(&format!("{field}: {bad} accepted"))
                .to_string();
            assert!(
                detail.contains("a non-integer exponent in "),
                "{field} {bad}: {detail}"
            );
        }
        for good in ["caller::x^-3", "caller::x^(2-3)"] {
            let mut v = base.clone();
            put(&mut v, &atom(good));
            v.validate()
                .unwrap_or_else(|e| panic!("{field} {good}: {e}"));
        }
    }
    // the projector (on a record with one) and the prefactor argument
    let root = record("wh-k1-root");
    let mut v = root.clone();
    let projector = v.projector.as_mut().unwrap();
    projector.tensor = &projector.tensor * atom("caller::x^(1/2)");
    assert!(
        v.validate()
            .unwrap_err()
            .to_string()
            .contains("a non-integer exponent in projector: ")
    );
    match closure(
        &root,
        &symbolic(),
        &(root.prefactor.total() * atom("caller::x^caller::y")),
    ) {
        Err(FactorizedError::Graph { detail }) => assert!(
            detail.starts_with("a non-integer exponent in prefactor argument: "),
            "{detail}"
        ),
        other => panic!("argument x^y: {other:?}"),
    }
    closure(
        &root,
        &symbolic(),
        &(root.prefactor.total() * atom("caller::x^-3")),
    )
    .unwrap();
}

/// a v3 schema line with trailing spaces or a tab is refused with
/// advice that names v3 (`TextSyntax`), not as an unknown schema of a reader
/// that "knows v1"; a v1 text keeps its earlier error.
#[test]
fn a_v3_schema_with_trailing_whitespace_names_v3() {
    let clean = fixture_text("wh-k1-root.record.txt");
    for pad in ["  ", "\t"] {
        let variant = clean.replacen(
            "pint6.factorized-numerator.v3",
            &format!("pint6.factorized-numerator.v3{pad}"),
            1,
        );
        match FactorizedNumerator::from_text(&variant) {
            Err(e @ FactorizedError::TextSchema { line: 1, .. }) => {
                let text = e.to_string();
                assert!(text.contains("pint6.factorized-numerator.v3, written exactly, without trailing whitespace"), "{text}")
            }
            other => panic!("{pad:?}: {other:?}"),
        }
    }
    FactorizedNumerator::from_text(&clean).unwrap();
}

// ===========================================================================
// ROUND 2: the allow-list rules attacked, every rule on every field
// ===========================================================================

/// The outcome of a marked record through its public entry (the closure
/// entry for a closure record, the ports entry for a direct one): Ok, or the
/// refusal text.
fn through_entry(
    v: &FactorizedNumerator,
    argument: &Atom,
    settings: &VakintSettings,
) -> Result<Atom, String> {
    if let Err(e) = v.validate() {
        return Err(format!("validate: {e}"));
    }
    if matches!(v.kind, FactorizedKind::Closure(_)) {
        evaluate_factorized_closure(
            &engine(),
            settings,
            std::slice::from_ref(v),
            std::slice::from_ref(argument),
        )
        .map(|r| r.value)
        .map_err(|e| format!("entry: {e}"))
    } else {
        evaluate_factorized_ports(&engine(), settings, v, argument)
            .map(|r| r.value)
            .map_err(|e| format!("entry: {e}"))
    }
}

/// One table row: rule x field -> refused?
type Setter = Box<dyn Fn(&mut FactorizedNumerator, &Atom)>;

fn field_setters(v: &FactorizedNumerator) -> Vec<(&'static str, Setter)> {
    let mut out: Vec<(&'static str, Setter)> = Vec::new();
    if let Some(e) = v
        .edges
        .iter()
        .find(|e| !e.is_external() && e.factor != Atom::num(1))
    {
        let id = e.id;
        out.push((
            "edge factor",
            Box::new(move |v, x| {
                let e = v.edges.iter_mut().find(|e| e.id == id).unwrap();
                e.factor = &e.factor * x;
            }),
        ));
    }
    if let Some(x) = v.vertices.iter().find(|x| x.factor != Atom::num(1)) {
        let id = x.id;
        out.push((
            "vertex factor",
            Box::new(move |v, x| {
                let e = v.vertices.iter_mut().find(|e| e.id == id).unwrap();
                e.factor = &e.factor * x;
            }),
        ));
    }
    if v.projector.is_some() {
        out.push((
            "projector",
            Box::new(|v, x| {
                let p = v.projector.as_mut().unwrap();
                p.tensor = &p.tensor * x;
            }),
        ));
    }
    out.push((
        "prefactor_overall",
        Box::new(|v, x| v.prefactor.overall = &v.prefactor.overall * x),
    ));
    out.push((
        "prefactor_scalars",
        Box::new(|v, x| v.prefactor.scalars = &v.prefactor.scalars * x),
    ));
    out.push((
        "prefactor_colour",
        Box::new(|v, x| v.prefactor.colour = &v.prefactor.colour * x),
    ));
    out.push((
        "weights",
        Box::new(|v, x| {
            let w = ("extra".to_string(), false, x.clone());
            match &mut v.kind {
                FactorizedKind::Closure(spec) => spec.weights.push(w),
                _ => v.record.as_mut().unwrap().weights.push(w),
            }
        }),
    ));
    if !v.insertions.is_empty() {
        out.push((
            "insertion coefficient",
            Box::new(|v, x| v.insertions[0].coefficient = &v.insertions[0].coefficient * x),
        ));
        out.push((
            "insertion factor",
            Box::new(|v, x| v.insertions[0].factor = &v.insertions[0].factor * x),
        ));
    }
    out
}

/// EVERY RULE ON EVERY FIELD: each rule's violating atom in each caller
/// field of a closure record with a projector (k1-root), one with an
/// insertion (k2-child-ins) and a v3 direct record (k5-direct), and in the
/// prefactor argument: refused (any variant), never a value. The table is
/// printed; a "no" is a defect.
#[test]
fn every_rule_on_every_field() {
    let settings = symbolic();
    let rules: Vec<(&str, String)> = vec![
        ("reserved namespace", "vakint::closure_lambda".into()),
        ("wildcard name", "caller::x_".into()),
        ("non-ASCII name", "caller::𝜀".into()),
        ("reference mass", R.into()),
        ("auxiliary mass", MAUX.into()),
        ("opaque function of D", format!("caller::f({DIM})")),
        ("opaque function of eps", format!("caller::f({EPS})")),
        ("non-integer exponent", "caller::x^(1/2)".into()),
        ("symbolic exponent", "caller::x^caller::y".into()),
        (
            "UNKNOWN bad shape",
            "vakint::UNKNOWN(vakint::m53,1/2)".into(),
        ),
        ("UNKNOWN with eps", format!("vakint::UNKNOWN({EPS},0)")),
    ];
    let mut misses = Vec::new();
    let mut table = Vec::new();
    for name in ["wh-k1-root", "wh-k2-child-ins", "wh-k5-direct"] {
        let base = record(name);
        let closure_kind = matches!(base.kind, FactorizedKind::Closure(_));
        // CONTROL: a harmless caller scalar in each field must be accepted,
        // else the column says nothing about the rules
        let mut controls = Vec::new();
        for (field, put) in field_setters(&base) {
            let mut v = base.clone();
            put(&mut v, &atom("caller::harmless"));
            let ok = through_entry(&v, &v.prefactor.total(), &settings).is_ok();
            controls.push(format!(
                "{field}: {}",
                if ok {
                    "control ok"
                } else {
                    "CONTROL REFUSED (column uninformative)"
                }
            ));
        }
        eprintln!("CONTRACT {name} | CONTROL | {}", controls.join(" | "));
        for (rule, text) in &rules {
            if !closure_kind && (*rule == "reference mass" || *rule == "auxiliary mass") {
                continue; // a direct record has no reference / auxiliary mass
            }
            let x = atom(text);
            let mut row = vec![format!("{name} | {rule}")];
            for (field, put) in field_setters(&base) {
                let mut v = base.clone();
                put(&mut v, &x);
                let refused = through_entry(&v, &v.prefactor.total(), &settings).is_err();
                row.push(format!("{field}: {}", if refused { "yes" } else { "NO" }));
                if !refused {
                    misses.push(format!("{name} {rule} in {field}"));
                }
            }
            let refused = through_entry(&base, &(base.prefactor.total() * &x), &settings).is_err();
            row.push(format!("argument: {}", if refused { "yes" } else { "NO" }));
            if !refused {
                misses.push(format!("{name} {rule} in the prefactor argument"));
            }
            table.push(row.join(" | "));
        }
    }
    for row in &table {
        eprintln!("CONTRACT {row}");
    }
    assert!(misses.is_empty(), "accepted:\n{}", misses.join("\n"));
}

/// BYPASSES of the function / exponent rules: sqrt and exp of ε (powers
/// after the parser), caller heads named like Symbolica's or like the
/// structural heads in another namespace with ε or D inside: refused.
#[test]
fn function_and_exponent_bypasses() {
    let base = record("wh-k1-root");
    let id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    let mut accepted = Vec::new();
    for bad in [
        format!("sqrt(1+{EPS})"),
        format!("exp({EPS})"),
        format!("caller::sqrt(1+{EPS})"),
        format!("caller::root({DIM},2)"),
        format!("caller::mink({DIM},{EPS})"),
        format!("caller::Q(1,{EPS})"),
        format!("caller::dot({EPS},caller::x)"),
        format!("caller::g({DIM},{DIM})"),
        format!("caller::gamma({EPS})"),
    ] {
        let mut v = base.clone();
        let vertex = v.vertices.iter_mut().find(|x| x.id == id).unwrap();
        vertex.factor = &vertex.factor * atom(&bad);
        if through_entry(&v, &v.prefactor.total(), &symbolic()).is_ok() {
            accepted.push(bad);
        }
    }
    assert!(accepted.is_empty(), "accepted: {accepted:?}");
}

/// A caller symbol NAMED like the masses in another namespace is a
/// different symbol: accepted as an opaque scalar and carried linearly
/// (no substitution reaches it).
#[test]
fn mass_names_in_another_namespace_are_plain_scalars() {
    let base = record("wh-k1-root");
    let plain = value(&base, &symbolic());
    for name in [
        "caller::rstar_reference_mass_squared",
        "caller::rstar_mass_squared",
    ] {
        let mut v = base.clone();
        v.prefactor.scalars = &v.prefactor.scalars * atom(name);
        assert_same(&value(&v, &symbolic()), &(&plain * atom(name)), name);
    }
}

/// Non-ASCII spellings: ε with a combining mark, ε followed by a zero-width
/// space, an ASCII name carrying a zero-width joiner or a combining mark,
/// in an insertion coefficient (typed) and through the text form: refused.
#[test]
fn combining_marks_and_zero_width_characters() {
    let base = record("wh-k2-child-ins");
    let mut accepted = Vec::new();
    for name in [
        "gammalooprs::ε\u{0301}",
        "gammalooprs::ε\u{200b}",
        "UFO::gs\u{200d}",
        "UFO::gs\u{0308}",
        "gammalooprs::\u{0065}\u{0301}",
    ] {
        init();
        let symbol = symbolica::symbol!(name);
        let mut v = base.clone();
        v.insertions[0].coefficient =
            &v.insertions[0].coefficient * (Atom::num(1) + Atom::var(symbol));
        if through_entry(&v, &v.prefactor.total(), &symbolic()).is_ok() {
            accepted.push(format!("typed {name:?}"));
        }
        // the text form
        let text = fixture_text("wh-k2-child-ins.record.txt");
        let line = text.lines().find(|l| l.starts_with("insertion\t")).unwrap();
        let variant = text.replace(line, &format!("{line}*(1+{name})"));
        if let Ok(r) = FactorizedNumerator::from_text(&variant)
            && through_entry(&r, &r.prefactor.total(), &symbolic()).is_ok()
        {
            accepted.push(format!("text {name:?}"));
        }
    }
    assert!(accepted.is_empty(), "accepted: {accepted:?}");
}

/// The settings' epsilon SET to a lookalike (`gammalooprs::𝜀`) and the record
/// written in it consistently: an epsilon of another name, the same value up
/// to the renaming; the record in the real `gammalooprs::ε` under that
/// setting is refused.
#[test]
fn a_lookalike_settings_epsilon() {
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    init();
    let look = Atom::var(symbolica::symbol!("gammalooprs::𝜀"));
    let settings = VakintSettings {
        epsilon_symbol: "gammalooprs::𝜀".into(),
        ..symbolic()
    };
    let mut v = base.clone();
    v.insertions[0].coefficient = v.insertions[0]
        .coefficient
        .replace(Atom::var(eps()))
        .with(look.clone());
    let got = through_entry(&v, &v.prefactor.total(), &settings).unwrap_or_else(|e| panic!("{e}"));
    assert_same(
        &got.replace(look).with(Atom::var(eps())),
        &plain,
        "lookalike settings epsilon",
    );
    assert!(through_entry(&base, &base.prefactor.total(), &settings).is_err());
}

/// The TYPED path with atoms the text form does not write: a float exponent
/// (refused); an integer exponent 10^30 on an ε- or D-dependent base in a
/// prefactor field (must be refused or evaluated within seconds: a value
/// computed by expanding (3-2ε)^(10^30) would not finish).
#[test]
fn typed_exponents() {
    let base = record("wh-k1-root");
    let mut v = base.clone();
    let vertex = v
        .vertices
        .iter_mut()
        .find(|x| x.factor != Atom::num(1))
        .unwrap();
    vertex.factor = &vertex.factor * atom("caller::x^0.5");
    assert!(
        through_entry(&v, &v.prefactor.total(), &symbolic()).is_err(),
        "float exponent accepted"
    );
    // 10^30 on (-1+D): run in a thread with a time limit
    let mut v = base.clone();
    v.prefactor.scalars =
        &v.prefactor.scalars * atom(&format!("(-1+{DIM})^1000000000000000000000000000000"));
    let (send, receive) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let r = through_entry(&v, &v.prefactor.total(), &symbolic()).map(|_| ());
        let _ = send.send(r);
    });
    match receive.recv_timeout(std::time::Duration::from_secs(60)) {
        Ok(Err(e)) => eprintln!("CONTRACT (D-1)^(10^30): refused: {e}"),
        Ok(Ok(())) => panic!("(D-1)^(10^30) in prefactor_scalars: accepted and evaluated"),
        Err(_) => panic!(
            "(D-1)^(10^30) in prefactor_scalars: no answer within 60 s (accepted, the evaluation does not finish)"
        ),
    }
}

/// Forbidden content INSIDE the arguments of STRUCTURAL heads (a slot
/// index, a momentum, dot, gamma): refused (by the exponent / function
/// rules or by the index rules), never a value.
#[test]
fn inside_structural_heads() {
    let base = record("wh-k1-root");
    let id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    let mink = "spenso::{spenso::representation,spenso::self_dual}::mink";
    let q = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
    let mut accepted = Vec::new();
    for bad in [
        format!("{mink}({DIM},caller::i^(1/2))*{mink}({DIM},caller::i^(1/2))"),
        format!("{mink}({DIM},log({EPS}))*{mink}({DIM},log({EPS}))"),
        format!("spenso::dot({q}(2,{mink}({DIM})),{q}(2,{mink}({DIM})))^(1/2)"),
        format!("spenso::dot({q}(2,{mink}({DIM})),log({EPS}))"),
        format!("{q}(2,{mink}({DIM},caller::j))*{q}(2,{mink}({DIM},caller::j))*{EPS}^(1/2)"),
    ] {
        let mut v = base.clone();
        let vertex = v.vertices.iter_mut().find(|x| x.id == id).unwrap();
        vertex.factor = &vertex.factor * atom(&bad);
        if through_entry(&v, &v.prefactor.total(), &symbolic()).is_ok() {
            accepted.push(bad);
        }
    }
    assert!(accepted.is_empty(), "accepted: {accepted:?}");
}

/// vakint::UNKNOWN misused: squared, in a denominator, inside a slot, with
/// ε or a mass in its arguments. Expected by C8 ("linear"): refused except a
/// linear occurrence; the check records which are accepted.
#[test]
fn unknown_misused() {
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1") };
    let base = record("wh-k2-child-ins");
    let u = "vakint::UNKNOWN(vakint::ggppss_int_53,0)";
    let mut accepted = Vec::new();
    for bad in [
        format!("{u}^2"),
        format!("{u}^-1"),
        format!("vakint::UNKNOWN({EPS},0)"),
        format!("vakint::UNKNOWN({R},0)"),
        format!(
            "vakint::UNKNOWN(vakint::ggppss_int_53,0)*vakint::UNKNOWN(vakint::ggppss_int_54,1)"
        ),
    ] {
        let mut v = base.clone();
        v.insertions[0].coefficient = &v.insertions[0].coefficient * (Atom::num(1) + atom(&bad));
        if through_entry(&v, &v.prefactor.total(), &symbolic()).is_ok() {
            accepted.push(bad);
        }
    }
    // inside a slot of a vertex factor
    let mut v = base.clone();
    let vertex = v
        .vertices
        .iter_mut()
        .find(|x| x.factor != Atom::num(1))
        .unwrap();
    vertex.factor = &vertex.factor
        * atom(&format!(
            "spenso::{{spenso::representation,spenso::self_dual}}::mink({DIM},{u})"
        ));
    if through_entry(&v, &v.prefactor.total(), &symbolic()).is_ok() {
        accepted.push("inside a slot".into());
    }
    eprintln!("CONTRACT UNKNOWN accepted: {accepted:?}");
    assert!(accepted.is_empty(), "accepted: {accepted:?}");
}

// ===========================================================================
// The three checks owed from round 1
// ===========================================================================

/// v1.4 §7 "open spinor chains only between the two declared bispinor
/// ports (OQ-6)": the k3-dirac value has no bispinor slot other than its
/// declared bispinor ports (here: none, the projector closes the quark
/// line).
#[test]
fn open_spinor_chains_only_between_declared_ports_k3_dirac() {
    let v = record("wh-k3-dirac");
    let got = value(&v, &symbolic());
    let declared = v
        .ports
        .iter()
        .filter(|p| p.to_string().contains("bis("))
        .map(|p| p.to_canonical_string())
        .collect::<std::collections::BTreeSet<_>>();
    let mut found = std::collections::BTreeSet::new();
    got.as_view().visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && f.get_symbol().get_name() == "spenso::bis"
        {
            found.insert(part.to_owned().to_canonical_string());
            return false;
        }
        true
    });
    let stray = found.difference(&declared).collect::<Vec<_>>();
    eprintln!(
        "CONTRACT OPEN SPINOR k3-dirac: declared bispinor ports {declared:?}, bispinor slots in the value {found:?}"
    );
    assert!(
        stray.is_empty(),
        "bispinor slots outside the declared ports: {stray:?}"
    );
}

/// Numerical mode on the VACUUM entry for a MARKED value (a two-loop
/// massive sunrise, v3 direct record): the same orders as Symbolic, the
/// Symbolic coefficients numerically.
#[test]
fn numerical_mode_on_the_vacuum_entry() {
    use vakint::p_integral::factorized::{
        DirectRecord, Endpoint, FactorizedEdge, FactorizedVertex, LineKind, Prefactor,
        evaluate_factorized_vacuum,
    };
    init();
    let line =
        |id: usize, s: usize, t: usize, h: (usize, usize), momentum: [i64; 2]| FactorizedEdge {
            id,
            source: Endpoint::Vertex(s),
            target: Endpoint::Vertex(t),
            source_hedge: Some(h.0),
            target_hedge: Some(h.1),
            particle: "x".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: atom("caller::m2"),
            momentum: momentum.to_vec(),
            factor: Atom::num(1),
        };
    let v = FactorizedNumerator {
        label: "vacuum sunrise".into(),
        kind: FactorizedKind::Direct,
        dimension: symbolica::symbol!("gammalooprs::dim"),
        loop_edges: vec![Some(1), Some(2)],
        external_edges: Vec::new(),
        edges: vec![
            line(1, 1, 2, (1, 2), [1, 0]),
            line(2, 1, 2, (3, 4), [0, 1]),
            line(3, 2, 1, (5, 6), [1, 1]),
        ],
        vertices: vec![
            FactorizedVertex {
                id: 1,
                rule: "V".into(),
                factor: Atom::num(1),
            },
            FactorizedVertex {
                id: 2,
                rule: "V".into(),
                factor: Atom::num(1),
            },
        ],
        ports: Vec::new(),
        projector: None,
        prefactor: Prefactor {
            scalars: atom("UFO::gs^2"),
            ..Prefactor::default()
        },
        colour_branch: None,
        insertions: Vec::new(),
        record: Some(DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
    };
    let vakint = engine();
    let run = |mode| {
        evaluate_factorized_vacuum(&vakint, &site(mode), &v, &v.prefactor.total())
            .unwrap_or_else(|e| panic!("{mode:?}: {e}"))
            .value
    };
    let exact = run(RustRedMasterMode::Symbolic);
    let numerical = run(RustRedMasterMode::Numerical);
    eprintln!("CONTRACT VACUUM symbolic {exact}\nP1 VACUUM numerical {numerical}");
    let e = wh_fixtures::epsilon_coefficients(&exact, eps());
    let n = wh_fixtures::epsilon_coefficients(&numerical, eps());
    assert_eq!(
        e.keys().collect::<Vec<_>>(),
        n.keys().collect::<Vec<_>>(),
        "orders"
    );
    let pi = atom("𝜋");
    let pi_value = atom("3.14159265358979323846264338327950288419716939937510582097494459");
    let mut skipped = Vec::new();
    for (order, c) in &e {
        // vakint's transcendental constants (HI, HR, ...) stay symbolic in
        // the Symbolic value: such a coefficient is not compared here
        let text = c.to_string();
        if text.contains("HI(") || text.contains("HR(") {
            skipped.push(*order);
            continue;
        }
        let difference = (c.replace(pi.clone()).with(pi_value.clone()).to_float(60)
            - n[order]
                .replace(pi.clone())
                .with(pi_value.clone())
                .to_float(60))
        .expand();
        let terms = match difference.as_view() {
            AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
            _ if difference.is_zero() => Vec::new(),
            _ => vec![difference.clone()],
        };
        for term in terms {
            assert!(
                parse_magnitude(&term.to_string()) < 1e-12,
                "eps^{order}: {term}"
            );
        }
    }
    eprintln!(
        "CONTRACT VACUUM numerical = symbolic at every order but {skipped:?} (vakint constants HI / HR there)"
    );
    assert!(skipped.len() < e.len(), "no coefficient compared");
}

/// UNKNOWN masters produced by the BACKEND through the closure entry: the
/// five-loop K1 record L5-60 rr1 evaluated three Laurent terms deeper than
/// its depth needs master coefficients the catalogue does not know: with
/// allow-unknowns OFF the entry refuses (`Evaluation`, the unknown-master
/// text); ON, the value carries `vakint::UNKNOWN(..)` linearly and the
/// result check lets it through.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES + VAKINT_TEST_P_LIBRARY_1_4"]
fn unknown_masters_from_the_backend_through_the_closure_entry() {
    wh_fixtures::required_env("VAKINT_TEST_P_LIBRARY_1_4", "five loops");
    let dir = wh_fixtures::fixture_dir("k1-l5").join("L5-60").join("rr1");
    let mut v =
        FactorizedNumerator::from_text(&wh_fixtures::read_fixture(&dir.join("input.txt"))).unwrap();
    let FactorizedKind::Closure(spec) = &mut v.kind else {
        unreachable!()
    };
    let depth = spec.depth.unwrap();
    let extra: usize = std::env::var("VAKINT_TEST_EXTRA_DEPTH")
        .ok()
        .and_then(|x| x.parse().ok())
        .unwrap_or(6);
    spec.depth = Some(depth + extra);
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::remove_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS") };
    let off = closure(&v, &symbolic(), &v.prefactor.total());
    eprintln!(
        "CONTRACT BACKEND UNKNOWN off: {:?}",
        off.as_ref()
            .map(|r| r.value.to_string().len())
            .map_err(|e| e.to_string())
    );
    let mut settings = symbolic();
    settings.allow_unknown_master_coefficients = true;
    let on = closure(&v, &settings, &v.prefactor.total());
    match (&off, &on) {
        (Err(FactorizedError::Evaluation { detail }), Ok(on)) => {
            assert!(detail.contains("expected master coefficient"), "{detail}");
            assert!(
                on.value.to_string().contains("UNKNOWN"),
                "no UNKNOWN in {}",
                on.value
            );
        }
        (Ok(_), Ok(on)) => {
            // depth + 3 still within the known coefficients: say so
            assert!(!on.value.to_string().contains("UNKNOWN"));
            panic!(
                "depth + {extra} needs no unknown master coefficient on L5-60 rr1: choose a deeper depth"
            );
        }
        other => panic!("{other:?}"),
    }
}

/// The UNKNOWN exception may not carry vakint's own names: the real k1-root
/// (L2-22) with prefactor_scalars x vakint::UNKNOWN(vakint::closure_tensor_0, 0)
/// (the closure entry's slot placeholder inside the atom), and the same with
/// the Taylor parameter, a core symbol or a wildcard name as the master:
/// refused at validation with allow-unknowns on and off; a catalogue-like
/// master name is accepted.
#[test]
fn unknown_with_an_internal_name_as_master_is_refused() {
    let base = record("wh-k1-root");
    for allow in [true, false] {
        // SAFETY: nextest runs each test in its own process
        unsafe {
            if allow {
                std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1")
            } else {
                std::env::remove_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS")
            }
        }
        for master in [
            "vakint::closure_tensor_0",
            "vakint::closure_lambda",
            "vakint::d",
            "vakint::ep",
            "vakint::slash",
            "UFO::x_",
            "UFO::x",
        ] {
            let mut v = base.clone();
            v.prefactor.scalars =
                &v.prefactor.scalars * atom(&format!("vakint::UNKNOWN({master},0)"));
            match through_entry(&v, &v.prefactor.total(), &symbolic()) {
                Err(e) => assert!(
                    e.contains("reserved symbol vakint::UNKNOWN"),
                    "{master}: {e}"
                ),
                Ok(got) => panic!("UNKNOWN({master}, 0) accepted (allow={allow}): {got}"),
            }
        }
    }
    // a catalogue master name OUTSIDE an insertion: refused (data atoms
    // stand only in insertion coefficients and factors)
    unsafe { std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1") };
    let mut v = base.clone();
    let u = atom("vakint::UNKNOWN(vakint::ggppss_int_53,-2)");
    v.prefactor.scalars = &v.prefactor.scalars * &u;
    assert!(through_entry(&v, &v.prefactor.total(), &symbolic()).is_err());
    // an order beyond the bound
    let mut v = base.clone();
    v.prefactor.scalars = &v.prefactor.scalars * atom("vakint::UNKNOWN(vakint::ggppss_int_53,65)");
    assert!(v.validate().is_err(), "an UNKNOWN order 65 accepted");
}

// ===========================================================================
// The vakint data atoms (UNKNOWN, MZV): shape, position, linearity, data
// ===========================================================================

/// One insertion coefficient of the real k2-child-ins times `extra`.
fn k2_coefficient_times(extra: &str) -> FactorizedNumerator {
    let mut v = record("wh-k2-child-ins");
    v.insertions[0].coefficient = &v.insertions[0].coefficient * atom(extra);
    v
}

/// The two data atoms in their exact shape in an insertion coefficient:
/// accepted and carried linearly (the value of c·A is A times the value of c);
/// UNKNOWN · ε^k raises the pole order by k.
#[test]
fn data_atoms_in_an_insertion_are_carried_linearly() {
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1") };
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    for data in [
        "vakint::MZV(5,3)",
        "vakint::MZV(3,5,3)",
        "vakint::UNKNOWN(vakint::ggppss_int_53,0)",
    ] {
        let v = k2_coefficient_times(data);
        assert_same(&value(&v, &symbolic()), &(&plain * atom(data)), data);
    }
    let order = base.insertion_pole_order(eps()).unwrap();
    for k in [1i64, 2, 3] {
        let v = k2_coefficient_times(&format!(
            "vakint::UNKNOWN(vakint::ggppss_int_53,{k})*{EPS}^(-{k})"
        ));
        assert_eq!(
            v.insertion_pole_order(eps()).unwrap(),
            order + k as u64,
            "UNKNOWN eps^-{k}"
        );
    }
}

/// The shapes: arguments are literals of the stated type; nothing else
/// inside (the placeholder hole must not reopen through MZV).
#[test]
fn data_atoms_have_their_exact_shape() {
    for bad in [
        "vakint::MZV(vakint::closure_tensor_0)",
        "vakint::MZV(caller::x)",
        "vakint::MZV(5/2,3)",
        "vakint::MZV(0,3)",
        "vakint::MZV(65)",
        "vakint::MZV()",
        "vakint::UNKNOWN(vakint::m53,0)",
        "vakint::UNKNOWN(vakint::ggppss_int_x,0)",
        "vakint::UNKNOWN(vakint::closure_tensor_0,0)",
        "vakint::UNKNOWN(vakint::ggppss_int_53,caller::k)",
    ] {
        let detail = k2_coefficient_times(bad)
            .validate()
            .expect_err(bad)
            .to_string();
        assert!(
            detail.contains("reserved symbol vakint::"),
            "{bad}: {detail}"
        );
    }
    for good in [
        "vakint::MZV(5,3)",
        "vakint::UNKNOWN(vakint::ggppss_int_53,0)",
        "vakint::UNKNOWN(vakint::x_t277_1_1,2)",
    ] {
        k2_coefficient_times(good)
            .validate()
            .unwrap_or_else(|e| panic!("{good}: {e}"));
    }
}

/// The positions: only in an insertion's coefficient and factor, as a
/// factor of a coefficient: refused in a vertex factor, the prefactor fields,
/// the prefactor argument, a mass, inside a function argument, as an
/// exponent; accepted as a factor of the insertion coefficient.
#[test]
fn data_atoms_only_in_insertions_and_not_inside_functions_or_exponents() {
    let base = record("wh-k2-child-ins");
    for data in [
        "vakint::MZV(5,3)",
        "vakint::UNKNOWN(vakint::ggppss_int_53,0)",
    ] {
        let x = atom(data);
        let mut f = base.clone();
        let vertex = f
            .vertices
            .iter_mut()
            .find(|v| v.factor != Atom::num(1))
            .unwrap();
        vertex.factor = &vertex.factor * &x;
        assert!(f.validate().is_err(), "{data} in a vertex factor");
        let mut s = base.clone();
        s.prefactor.scalars = &s.prefactor.scalars * &x;
        assert!(s.validate().is_err(), "{data} in prefactor_scalars");
        assert!(
            through_entry(&base, &(base.prefactor.total() * &x), &symbolic()).is_err(),
            "{data} in the prefactor argument"
        );
        let mut m = base.clone();
        if let FactorizedKind::Closure(spec) = &mut m.kind {
            spec.auxiliary_mass_squared = x.clone();
        }
        assert!(m.validate().is_err(), "{data} as the auxiliary mass");
        for misplaced in [
            format!("caller::f({data})"),
            format!("caller::x^{data}"),
            format!("log({data})"),
        ] {
            let detail = k2_coefficient_times(&misplaced)
                .validate()
                .expect_err(&misplaced)
                .to_string();
            assert!(
                detail.contains("a misplaced vakint data atom")
                    || detail.contains("a non-integer exponent"),
                "{misplaced}: {detail}"
            );
        }
        k2_coefficient_times(data)
            .validate()
            .unwrap_or_else(|e| panic!("{data}: {e}"));
    }
}

/// The data atoms are derived from the code and the data: every vakint
/// symbol in the master coefficients of a catalogue (its `exact` values) is
/// one of VAKINT_DATA_ATOMS (MZV, UNKNOWN with its master), every other
/// constant is Symbolica's (zeta, pi). The built-in catalogue (default suite).
#[test]
fn the_builtin_catalogue_carries_only_listed_data_atoms() {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/p_master/p-master-catalogue-5l.json");
    catalogue_constants(&path);
}

/// The same for the SHIPPED catalogue of the data package
/// (`$VAKINT_TEST_P_LIBRARY_1_5/p-master-catalogue.json`, e.g.
/// ggppss-5L-v4-rc4a, which carries UNKNOWN entries).
#[test]
#[ignore = "external data: VAKINT_TEST_P_LIBRARY_1_5"]
fn the_shipped_catalogue_carries_only_listed_data_atoms() {
    let dir = wh_fixtures::required_env(
        "VAKINT_TEST_P_LIBRARY_1_5",
        "the shipped p master catalogue",
    );
    catalogue_constants(&dir.join("p-master-catalogue.json"));
}

/// Every function head and bare name in a catalogue's master coefficients:
/// only `zeta`, `MZV`, `UNKNOWN` (with a master name and an integer), `pi`.
fn catalogue_constants(path: &std::path::Path) {
    let text = wh_fixtures::read_fixture(path);
    let json: serde_json::Value = serde_json::from_str(&text).unwrap();
    let mut seen = std::collections::BTreeMap::<String, usize>::new();
    let mut coefficients = 0;
    for master in json["masters"].as_array().unwrap() {
        let Some(map) = master["coefficients"].as_object() else {
            continue;
        };
        for entry in map.values() {
            let Some(exact) = entry["exact"].as_str() else {
                continue;
            };
            coefficients += 1;
            let mut chars = exact.char_indices().peekable();
            while let Some((i, c)) = chars.next() {
                if c.is_ascii_alphabetic() || c == '_' {
                    let mut end = i + c.len_utf8();
                    while let Some((j, d)) = chars.peek() {
                        if d.is_ascii_alphanumeric() || *d == '_' {
                            end = j + d.len_utf8();
                            chars.next();
                        } else {
                            break;
                        }
                    }
                    let name = &exact[i..end];
                    let head = exact[end..].starts_with('(');
                    let kind = match (name, head) {
                        ("zeta", true) => "zeta (symbolica::zeta after parse_constant)".to_string(),
                        ("MZV", true) => "MZV (vakint data atom)".to_string(),
                        ("UNKNOWN", true) => "UNKNOWN (vakint data atom)".to_string(),
                        ("pi", false) => "pi (symbolica::𝜋 after parse_constant)".to_string(),
                        (n, false) if n.starts_with("ggppss_int_") || n.starts_with("x_t") => {
                            "an UNKNOWN master name".to_string()
                        }
                        (n, h) => format!("UNLISTED {n} head={h}"),
                    };
                    *seen.entry(kind).or_default() += 1;
                }
            }
        }
    }
    eprintln!(
        "CATALOGUE {}: {coefficients} exact coefficients, constants {seen:?}",
        path.display()
    );
    assert!(coefficients > 0);
    let unlisted = seen
        .keys()
        .filter(|k| k.starts_with("UNLISTED"))
        .collect::<Vec<_>>();
    assert!(
        unlisted.is_empty(),
        "constants not in the list: {unlisted:?}"
    );
}

/// A MARKED two-loop massive vacuum sunrise (v3 direct record) for the
/// vacuum entry; `mass` is the lines' mass, `factor` the first vertex's.
fn marked_vacuum_sunrise(mass: &str, factor: &str) -> FactorizedNumerator {
    use vakint::p_integral::factorized::{
        DirectRecord, Endpoint, FactorizedEdge, FactorizedVertex, LineKind, Prefactor,
    };
    init();
    let line =
        |id: usize, s: usize, t: usize, h: (usize, usize), momentum: [i64; 2]| FactorizedEdge {
            id,
            source: Endpoint::Vertex(s),
            target: Endpoint::Vertex(t),
            source_hedge: Some(h.0),
            target_hedge: Some(h.1),
            particle: "x".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: atom(mass),
            momentum: momentum.to_vec(),
            factor: Atom::num(1),
        };
    FactorizedNumerator {
        label: "marked vacuum sunrise".into(),
        kind: FactorizedKind::Direct,
        dimension: symbolica::symbol!("gammalooprs::dim"),
        loop_edges: vec![Some(1), Some(2)],
        external_edges: Vec::new(),
        edges: vec![
            line(1, 1, 2, (1, 2), [1, 0]),
            line(2, 1, 2, (3, 4), [0, 1]),
            line(3, 2, 1, (5, 6), [1, 1]),
        ],
        vertices: vec![
            FactorizedVertex {
                id: 1,
                rule: "V".into(),
                factor: atom(factor),
            },
            FactorizedVertex {
                id: 2,
                rule: "V".into(),
                factor: Atom::num(1),
            },
        ],
        ports: Vec::new(),
        projector: None,
        prefactor: Prefactor {
            scalars: atom("UFO::gs^2"),
            ..Prefactor::default()
        },
        colour_branch: None,
        insertions: Vec::new(),
        record: Some(DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
    }
}

/// The vacuum entry applies a MARKED value's atom rules (it never calls
/// `validate`): a vertex factor x log(ε), a fractional or symbolic exponent,
/// an edge mass that is not a plain symbol: refused; the plain record
/// evaluates.
#[test]
fn the_vacuum_entry_applies_the_marked_atom_rules() {
    use vakint::p_integral::factorized::evaluate_factorized_vacuum;
    let vakint = engine();
    let run = |v: &FactorizedNumerator| {
        evaluate_factorized_vacuum(&vakint, &symbolic(), v, &v.prefactor.total()).map(|r| r.value)
    };
    run(&marked_vacuum_sunrise("caller::m2", "1")).unwrap();
    for (mass, factor, prefix) in [
        (
            "caller::m2",
            format!("log({EPS})"),
            "a non-polynomial dependence in vertex 1 factor: ",
        ),
        (
            "caller::m2",
            "caller::x^(1/2)".to_string(),
            "a non-integer exponent in vertex 1 factor: ",
        ),
        (
            "caller::m2",
            "caller::x^caller::y".to_string(),
            "a non-integer exponent in vertex 1 factor: ",
        ),
        ("caller::m2^(1/2)", "1".to_string(), "edge 1: mass "),
        ("log(caller::m2)", "1".to_string(), "edge 1: mass "),
    ] {
        match run(&marked_vacuum_sunrise(mass, &factor)) {
            Err(FactorizedError::Graph { detail }) => {
                assert!(detail.starts_with(prefix), "{mass} {factor}: {detail}")
            }
            other => panic!("{mass} {factor}: {other:?}"),
        }
    }
}

/// Outside ASCII only the settings' epsilon and scale and Symbolica's 𝜋 and
/// 𝑒: a caller symbol in Symbolica's namespace with a non-ASCII name
/// (`symbolica::ε`, a lookalike of the epsilon) is refused in a factor; a
/// non-ASCII DIMENSION name is refused; 𝜋 and 𝑒 are accepted.
#[test]
fn non_ascii_names_only_the_settings_and_symbolica_constants() {
    let base = record("wh-k1-root");
    let id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    init();
    let lookalike = Atom::var(symbolica::symbol!("symbolica::ε"));
    let mut v = base.clone();
    let vertex = v.vertices.iter_mut().find(|x| x.id == id).unwrap();
    vertex.factor = &vertex.factor * lookalike.pow(Atom::num(-1));
    match through_entry(&v, &v.prefactor.total(), &symbolic()) {
        Err(e) => assert!(
            e.contains(&format!(
                "a non-ASCII name symbolica::ε in vertex {id} factor"
            )),
            "{e}"
        ),
        Ok(got) => panic!("symbolica::ε accepted: {got}"),
    }
    for ok in ["𝜋^-2", "𝑒^2"] {
        let mut v = base.clone();
        let vertex = v.vertices.iter_mut().find(|x| x.id == id).unwrap();
        vertex.factor = &vertex.factor * atom(ok);
        through_entry(&v, &v.prefactor.total(), &symbolic())
            .unwrap_or_else(|e| panic!("{ok}: {e}"));
    }
    // a non-ASCII dimension name, used consistently
    let text =
        fixture_text("wh-k1-root.record.txt").replace("gammalooprs::{}::dim", "gammalooprs::{}::𝒟");
    let v = FactorizedNumerator::from_text(&text).unwrap_or_else(|e| panic!("{e}"));
    match through_entry(&v, &v.prefactor.total(), &symbolic()) {
        Err(e) => assert!(e.contains("a non-ASCII name gammalooprs::𝒟 in "), "{e}"),
        Ok(got) => panic!("a non-ASCII dimension accepted: {got}"),
    }
}

/// The masses by ROLE: an auxiliary mass with another name than the
/// producer's (`caller::maux`) inside an insertion coefficient, and the
/// reference mass in prefactor_overall of a record without original_cycles:
/// refused; the literal 0 stays the massless lines' mass (accepted), a mass
/// written as a zero polynomial is refused.
#[test]
fn the_masses_are_identified_by_their_role() {
    let base = record("wh-k2-child-ins");
    let mut v = base.clone();
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = atom("caller::maux");
    }
    v.validate().unwrap();
    let mut bad = v.clone();
    bad.insertions[0].coefficient = &bad.insertions[0].coefficient * atom("caller::maux");
    match bad.validate() {
        Err(FactorizedError::Insertion { detail, .. }) => assert!(
            detail.starts_with(&format!(
                "the auxiliary mass {} in insertion ",
                atom("caller::maux")
            )),
            "{detail}"
        ),
        other => panic!("{other:?}"),
    }
    // R in prefactor_overall (no original_cycles: overall is free)
    let mut r = record("wh-k1-root");
    if let FactorizedKind::Closure(spec) = &mut r.kind {
        spec.original_cycles = None;
    }
    r.prefactor.loop_sign = -1;
    r.validate().unwrap();
    let mut bad = r.clone();
    bad.prefactor.overall = &bad.prefactor.overall * atom(R);
    closure_refusal(
        &bad,
        &format!("the reference mass {} in prefactor_overall: ", atom(R)),
        "R in overall",
    );
    // the literal 0 on massless lines: every real record has it (accepted)
    assert!(base.edges.iter().any(|e| e.mass_squared.is_zero()));
    base.validate().unwrap();
}

/// Integer exponents are bounded (|n| <= 64): (1+x)^65, x^-65 and D^1000000
/// in a vertex factor or the prefactor argument are refused; x^64 is
/// accepted. A bare ε, dimension or mass as a DIRECT argument of a
/// structural head (hedge(ε, 1), mink(D, ε), g(.., R)) is refused.
#[test]
fn exponents_are_bounded_and_structural_arguments_are_checked() {
    let base = record("wh-k1-root");
    let id = base
        .vertices
        .iter()
        .find(|x| x.factor != Atom::num(1))
        .unwrap()
        .id;
    let with = |extra: &str| {
        let mut v = base.clone();
        let vertex = v.vertices.iter_mut().find(|x| x.id == id).unwrap();
        vertex.factor = &vertex.factor * atom(extra);
        v
    };
    for bad in [
        "(1+caller::x)^65",
        "caller::x^-65",
        &format!("{DIM}^1000000"),
    ] {
        let detail = with(bad).validate().expect_err(bad).to_string();
        assert!(detail.contains("exceeds 64 in modulus"), "{bad}: {detail}");
    }
    with("caller::x^64").validate().unwrap();
    match through_entry(
        &base,
        &(base.prefactor.total() * atom("(1+caller::x)^1000000000")),
        &symbolic(),
    ) {
        Err(e) => assert!(e.contains("exceeds 64 in modulus"), "{e}"),
        Ok(_) => panic!("(1+x)^(10^9) in the argument accepted"),
    }
    let mink = "spenso::{spenso::representation,spenso::self_dual}::mink";
    let hedge = "gammalooprs::{spenso::index}::hedge";
    for bad in [
        format!("{mink}({DIM},{hedge}({EPS},1))*{mink}({DIM},{hedge}({EPS},1))"),
        format!("{mink}({DIM},{EPS})*{mink}({DIM},{EPS})"),
        format!("{mink}({DIM},{R})*{mink}({DIM},{R})"),
    ] {
        assert!(with(&bad).validate().is_err(), "{bad} accepted");
    }
}

/// MARKED three-loop massive vacuum values through the vacuum entry in
/// Symbolic and Numerical mode at 3 and 6 terms: the basketball (four lines
/// of one mass between two vertices) and the all-massive mercedes (the
/// tetrahedron, six lines), whose Symbolic value at 6 terms carries vakint's
/// three-loop vacuum constant D6bar1: no result-check refusal (vakint's
/// vacuum constants are output, not internal names).
#[test]
fn a_marked_three_loop_massive_tadpole_through_the_vacuum_entry() {
    use vakint::p_integral::factorized::{
        DirectRecord, Endpoint, FactorizedEdge, FactorizedVertex, LineKind, Prefactor,
        evaluate_factorized_vacuum,
    };
    init();
    let line =
        |id: usize, s: usize, t: usize, h: (usize, usize), momentum: [i64; 3]| FactorizedEdge {
            id,
            source: Endpoint::Vertex(s),
            target: Endpoint::Vertex(t),
            source_hedge: Some(h.0),
            target_hedge: Some(h.1),
            particle: "x".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: atom("caller::m2"),
            momentum: momentum.to_vec(),
            factor: Atom::num(1),
        };
    let v = FactorizedNumerator {
        label: "marked three-loop basketball".into(),
        kind: FactorizedKind::Direct,
        dimension: symbolica::symbol!("gammalooprs::dim"),
        loop_edges: vec![Some(1), Some(2), Some(3)],
        external_edges: Vec::new(),
        edges: vec![
            line(1, 1, 2, (1, 2), [1, 0, 0]),
            line(2, 1, 2, (3, 4), [0, 1, 0]),
            line(3, 1, 2, (5, 6), [0, 0, 1]),
            line(4, 2, 1, (7, 8), [1, 1, 1]),
        ],
        vertices: vec![
            FactorizedVertex {
                id: 1,
                rule: "V".into(),
                factor: Atom::num(1),
            },
            FactorizedVertex {
                id: 2,
                rule: "V".into(),
                factor: Atom::num(1),
            },
        ],
        ports: Vec::new(),
        projector: None,
        prefactor: Prefactor {
            scalars: atom("UFO::gs^2"),
            ..Prefactor::default()
        },
        colour_branch: None,
        insertions: Vec::new(),
        record: Some(DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
    };
    // the three-loop all-massive tetrahedron (six lines)
    let mut mercedes = v.clone();
    mercedes.label = "marked three-loop mercedes".into();
    mercedes.edges = vec![
        line(1, 1, 2, (1, 2), [1, 0, 0]),
        line(2, 2, 3, (3, 4), [0, 1, 0]),
        line(3, 3, 1, (5, 6), [0, 0, 1]),
        line(4, 4, 1, (7, 8), [1, 0, -1]),
        line(5, 2, 4, (9, 10), [1, -1, 0]),
        line(6, 4, 3, (11, 12), [0, -1, 1]),
    ];
    mercedes.vertices = (1..=4)
        .map(|id| FactorizedVertex {
            id,
            rule: "V".into(),
            factor: Atom::num(1),
        })
        .collect();
    let vakint = engine();
    let mut failures = Vec::new();
    let mut d6bar1 = false;
    for v in [&v, &mercedes] {
        for (mode, terms) in [
            (RustRedMasterMode::Symbolic, 3),
            (RustRedMasterMode::Numerical, 3),
            (RustRedMasterMode::Symbolic, 6),
            (RustRedMasterMode::Numerical, 6),
        ] {
            let settings = VakintSettings {
                number_of_terms_in_epsilon_expansion: terms,
                ..site(mode)
            };
            match evaluate_factorized_vacuum(&vakint, &settings, v, &v.prefactor.total()) {
                Ok(r) => {
                    eprintln!(
                        "TADPOLE3 {} {mode:?} terms={terms}: ok {}",
                        v.label, r.value
                    );
                    d6bar1 |= r.value.to_canonical_string().contains("D6bar1");
                }
                Err(e) => {
                    eprintln!("TADPOLE3 {} {mode:?} terms={terms}: ERROR {e}", v.label);
                    failures.push(format!("{} {mode:?} terms={terms}: {e}", v.label));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:?}");
    // the witness is live: a value did carry the constant
    assert!(d6bar1, "no value carried D6bar1");
}

/// Every unique v3 record under the roots of `VAKINT_TEST_RECORD_ROOTS`
/// (colon-separated; files `*.record.txt` and `input.txt`) through
/// validation only (`from_text`). The roots this split is stated for:
/// `/common/build-cache/ben/consol:/common/ben/beta/results/consolidation/wh-fixtures`
/// (6,240 files, 2,633 unique v3 records). EXACTLY 2,611 are accepted and 22
/// refused, every one a Closure refusal "every line is massive: not a
/// closure" (one-loop remainders whose only line left is the closing massive
/// line: a massive tadpole, the producer's designed decline
/// vakint:all_massive, dumped before the vakint call); any other refusal, or
/// another count, fails.
#[test]
#[ignore = "external data: VAKINT_TEST_RECORD_ROOTS"]
fn every_production_record_validates() {
    const ALL_MASSIVE: &str = "every line is massive: not a closure";
    let roots =
        wh_fixtures::required_env("VAKINT_TEST_RECORD_ROOTS", "the producer's record roots");
    let mut files = Vec::new();
    for root in roots.to_string_lossy().split(':') {
        let mut stack = vec![PathBuf::from(root)];
        while let Some(d) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&d) else {
                continue;
            };
            for e in entries {
                let p = e.unwrap().path();
                let name = p.file_name().unwrap().to_string_lossy().to_string();
                if p.is_dir() {
                    stack.push(p);
                } else if name == "input.txt" || name.ends_with(".record.txt") {
                    files.push(p);
                }
            }
        }
    }
    init();
    let started = std::time::Instant::now();
    let mut seen = std::collections::BTreeSet::new();
    let (mut accepted, mut refused) = (
        0usize,
        std::collections::BTreeMap::<String, Vec<String>>::new(),
    );
    for p in &files {
        let Ok(text) = std::fs::read_to_string(p) else {
            continue;
        };
        if !text.starts_with("schema\tpint6.factorized-numerator.v3") || !seen.insert(text.clone())
        {
            continue;
        }
        match FactorizedNumerator::from_text(&text) {
            Ok(_) => accepted += 1,
            Err(e) => {
                let detail = e.to_string();
                let expected = matches!(&e, FactorizedError::Closure { detail }
                    if detail.starts_with(ALL_MASSIVE));
                let prefix = if expected {
                    ALL_MASSIVE.to_string()
                } else {
                    format!("UNEXPECTED {detail}")
                };
                refused
                    .entry(prefix)
                    .or_default()
                    .push(p.display().to_string());
            }
        }
    }
    eprintln!(
        "PRODUCTION files={} unique_v3={} accepted={accepted} refused={} seconds={:.1}",
        files.len(),
        seen.len(),
        refused.values().map(Vec::len).sum::<usize>(),
        started.elapsed().as_secs_f64()
    );
    for (prefix, paths) in &refused {
        eprintln!(
            "PRODUCTION REFUSED {} x {prefix}: e.g. {}",
            paths.len(),
            paths[0]
        );
    }
    let refused_total = refused.values().map(Vec::len).sum::<usize>();
    assert_eq!(
        refused.keys().collect::<Vec<_>>(),
        vec![ALL_MASSIVE],
        "a refusal other than the all-massive tadpole: {refused:?}"
    );
    assert_eq!(
        (seen.len(), accepted, refused_total),
        (2633, 2611, 22),
        "(unique v3, accepted, refused)"
    );
}

/// The producer's plain caller scalars stay opaque coefficients: UFO
/// couplings, beta_projection::transverse, uv_isotropic_channel_n,
/// insertion_closure_metric / _longitudinal in an insertion coefficient and
/// in a vertex factor validate and evaluate (input symbols for the result
/// check), the value linear in each.
#[test]
fn the_producers_plain_scalars_are_accepted() {
    let base = record("wh-k2-child-ins");
    let plain = value(&base, &symbolic());
    for name in [
        "UFO::gs_minus",
        "UFO::gs2_i",
        "beta_projection::transverse",
        "gammalooprs::uv_isotropic_channel_1",
        "gammalooprs::insertion_closure_metric",
        "gammalooprs::insertion_closure_longitudinal",
    ] {
        let x = atom(name);
        let mut c = base.clone();
        c.insertions[0].coefficient = &c.insertions[0].coefficient * &x;
        assert_same(
            &value(&c, &symbolic()),
            &(&plain * &x),
            &format!("{name} in the coefficient"),
        );
        let mut f = base.clone();
        let vertex = f
            .vertices
            .iter_mut()
            .find(|v| v.factor != Atom::num(1))
            .unwrap();
        vertex.factor = &vertex.factor * &x;
        assert_same(
            &value(&f, &symbolic()),
            &(&plain * &x),
            &format!("{name} in a vertex factor"),
        );
    }
}

/// EVERY RULE ON EVERY FIELD AND ENTRY, typed and through the text form: the
/// closure entry (k1-root, k2-child-ins), the ports and the direct entry
/// (k5-direct), the vacuum entry (a marked vacuum sunrise); each rule's
/// violating atom in each caller field and in the prefactor argument is
/// refused (with a control per field that a harmless scalar is accepted).
#[test]
fn every_rule_on_every_field_and_entry() {
    use vakint::p_integral::factorized::{evaluate_factorized, evaluate_factorized_vacuum};
    let settings = symbolic();
    let vakint = engine();
    #[derive(Clone, Copy, Debug)]
    enum Entry {
        Closure,
        Ports,
        Direct,
        Vacuum,
    }
    let run = |entry: Entry, v: &FactorizedNumerator, arg: &Atom| -> Result<(), String> {
        match entry {
            Entry::Closure => evaluate_factorized_closure(
                &vakint,
                &settings,
                std::slice::from_ref(v),
                std::slice::from_ref(arg),
            )
            .map(|_| ()),
            Entry::Ports => evaluate_factorized_ports(&vakint, &settings, v, arg).map(|_| ()),
            Entry::Direct => evaluate_factorized(&vakint, &settings, v, arg).map(|_| ()),
            Entry::Vacuum => evaluate_factorized_vacuum(&vakint, &settings, v, arg).map(|_| ()),
        }
        .map_err(|e| e.to_string())
    };
    let rules: Vec<(&str, String)> = vec![
        ("reserved namespace", "vakint::closure_lambda".into()),
        ("wildcard name", "caller::x_".into()),
        ("non-ASCII name", "caller::𝜀".into()),
        ("opaque function of D", format!("caller::f({DIM})")),
        ("opaque function of eps", format!("caller::f({EPS})")),
        ("non-integer exponent", "caller::x^(1/2)".into()),
        ("unbounded exponent", "caller::x^65".into()),
        ("data atom outside an insertion", "vakint::MZV(5,3)".into()),
    ];
    let cases: Vec<(&str, FactorizedNumerator, Vec<Entry>)> = vec![
        ("k1-root", record("wh-k1-root"), vec![Entry::Closure]),
        (
            "k2-child-ins",
            record("wh-k2-child-ins"),
            vec![Entry::Closure],
        ),
        (
            "k5-direct",
            record("wh-k5-direct"),
            vec![Entry::Ports, Entry::Direct],
        ),
        (
            "vacuum sunrise",
            marked_vacuum_sunrise("caller::m2", "caller::g"),
            vec![Entry::Vacuum],
        ),
    ];
    let mut misses = Vec::new();
    for (name, base, entries) in &cases {
        for entry in entries {
            for text_path in [false, true] {
                let through = |v: &FactorizedNumerator, arg: &Atom| -> Result<(), String> {
                    let v = if text_path {
                        let text = v.to_text().map_err(|e| e.to_string())?;
                        FactorizedNumerator::parse_text(&text).map_err(|e| e.to_string())?
                    } else {
                        v.clone()
                    };
                    run(*entry, &v, arg)
                };
                for (field, put) in field_setters(base) {
                    let mut c = base.clone();
                    put(&mut c, &atom("caller::harmless"));
                    let informative = through(&c, &c.prefactor.total()).is_ok();
                    for (rule, text) in &rules {
                        if *rule == "data atom outside an insertion"
                            && field.starts_with("insertion")
                        {
                            continue;
                        }
                        let mut v = base.clone();
                        put(&mut v, &atom(text));
                        if through(&v, &v.prefactor.total()).is_ok() {
                            misses.push(format!(
                                "{name} {entry:?} text={text_path} {rule} in {field}"
                            ));
                        }
                    }
                    eprintln!(
                        "CONTRACT {name} {entry:?} text={text_path} {field}: control {}",
                        if informative {
                            "ok"
                        } else {
                            "refused (uninformative)"
                        }
                    );
                }
                for (rule, text) in &rules {
                    if through(base, &(base.prefactor.total() * atom(text))).is_ok() {
                        misses.push(format!(
                            "{name} {entry:?} text={text_path} {rule} in the prefactor argument"
                        ));
                    }
                }
            }
        }
    }
    assert!(misses.is_empty(), "accepted:\n{}", misses.join("\n"));
}

/// The audit's inputs per rule (six each for G1, G2, G3, the exponents and
/// UNKNOWN) not already a cell of `every_rule_on_every_field_and_entry`: each
/// is REFUSED (at parse, validation or evaluation; never a panic), or, where
/// it is a legitimate spelling, gives the SAME value as a control that differs
/// only by a renaming, or (marked OBSERVE) is reported without an expectation.
#[test]
fn the_audit_inputs_per_rule() {
    use vakint::p_integral::factorized::evaluate_factorized_vacuum;
    // SAFETY: nextest runs each test in its own process
    unsafe { std::env::set_var("VAKINT_ALLOW_UNKNOWN_MASTER_COEFFICIENTS", "1") };
    let vakint = engine();
    let settings = symbolic();
    let k2_text = fixture_text("wh-k2-child-ins.record.txt");
    let k1_text = fixture_text("wh-k1-root.record.txt");
    type Input = Result<(FactorizedNumerator, Atom, bool), String>;
    let try_atom = |t: &str| -> Result<Atom, String> {
        init();
        std::panic::catch_unwind(|| symbolica::try_parse!(t))
            .map_err(|_| format!("parse panicked: {t}"))?
            .map_err(|e| format!("parse: {e}"))
    };
    let parse = |t: &str| FactorizedNumerator::from_text(t).map_err(|e| e.to_string());
    let plain = |v: FactorizedNumerator| -> Input {
        let a = v.prefactor.total();
        Ok((v, a, false))
    };
    let in_field = |text: &str, field: &str, x: &str| -> Input {
        let mut v = parse(text)?;
        let x = try_atom(x)?;
        let (_, put) = field_setters(&v)
            .into_iter()
            .find(|(f, _)| *f == field)
            .unwrap_or_else(|| panic!("no field {field}"));
        put(&mut v, &x);
        plain(v)
    };
    let in_argument = |text: &str, x: &str| -> Input {
        let v = parse(text)?;
        let a = v.prefactor.total() * try_atom(x)?;
        Ok((v, a, false))
    };
    let vacuum = |mass: &str, factor: &str| -> Input {
        let v = marked_vacuum_sunrise(mass, factor);
        let a = v.prefactor.total();
        Ok((v, a, true))
    };
    let run = |input: Input| -> Result<Atom, String> {
        let (v, a, is_vacuum) = input?;
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if is_vacuum {
                evaluate_factorized_vacuum(&vakint, &settings, &v, &a).map(|r| r.value)
            } else {
                evaluate_factorized_closure(
                    &vakint,
                    &settings,
                    std::slice::from_ref(&v),
                    std::slice::from_ref(&a),
                )
                .map(|r| r.value)
            }
        }));
        match outcome {
            Ok(r) => r.map_err(|e| e.to_string()),
            Err(_) => Err("PANIC".into()),
        }
    };
    enum Expect {
        Refused,
        // the same value as the control input after the renaming (pattern, image)
        Same(Input, &'static str, &'static str),
        Observe,
    }
    let edge5_r = |to: &str| {
        k2_text.replacen(
            &format!("momentum=1,0,0\t{R}"),
            &format!("momentum=1,0,0\t{to}"),
            1,
        )
    };
    // the text form of a closure record requires `original_cycles`: the
    // record without it is typed (the spec field set to None)
    let no_cycles = |x: &str| -> Input {
        let mut v = parse(&k2_text)?;
        if let FactorizedKind::Closure(spec) = &mut v.kind {
            spec.original_cycles = None;
        }
        v.prefactor.overall = &v.prefactor.overall * try_atom(x)?;
        plain(v)
    };
    let ok_master = "vakint::ggppss_int_1";
    let rows: Vec<(&str, Input, Expect)> = vec![
        // G1
        (
            "G1-1 vertex x symbolica::ε^-1",
            in_field(&k2_text, "vertex factor", "symbolica::ε^-1"),
            Expect::Refused,
        ),
        (
            "G1-2 UFO::g<ZWSP>s in a vertex factor",
            in_field(&k2_text, "vertex factor", "UFO::g\u{200B}s"),
            Expect::Refused,
        ),
        (
            "G1-3 ε<U+0301> in an insertion coefficient",
            in_field(
                &k2_text,
                "insertion coefficient",
                "gammalooprs::ε\u{0301}^-1",
            ),
            Expect::Refused,
        ),
        (
            "G1-4 the dimension named gammalooprs::𝒟",
            parse(&k2_text.replace(DIM, "gammalooprs::{}::𝒟")).and_then(plain),
            Expect::Refused,
        ),
        (
            "G1-5 the settings' μ as a function head μᵣ²(UFO::x)",
            in_field(&k2_text, "vertex factor", "gammalooprs::μᵣ²(UFO::x)"),
            Expect::Same(
                in_field(&k2_text, "vertex factor", "caller::h(UFO::x)"),
                "caller::h(UFO::x)",
                "gammalooprs::μᵣ²(UFO::x)",
            ),
        ),
        (
            "G1-6 UNKNOWN(UFO::μ, 0) in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                "vakint::UNKNOWN(UFO::μ,0)",
            ),
            Expect::Refused,
        ),
        // G2
        (
            "G2-1 the auxiliary mass named UFO::maux, x eps^-1 in an insertion",
            in_field(
                &k2_text.replace(MAUX, "UFO::maux"),
                "insertion coefficient",
                "UFO::maux*gammalooprs::ε^-1",
            ),
            Expect::Refused,
        ),
        (
            "G2-1c control: the auxiliary mass named UFO::maux",
            parse(&k2_text.replace(MAUX, "UFO::maux")).and_then(plain),
            Expect::Observe,
        ),
        (
            "G2-2 R in prefactor_overall without original_cycles (typed)",
            no_cycles(R),
            Expect::Refused,
        ),
        (
            "G2-2c control: no original_cycles (typed), overall x caller::c",
            no_cycles("caller::c"),
            Expect::Observe,
        ),
        (
            "G2-3 R in a weight",
            in_field(&k2_text, "weights", R),
            Expect::Refused,
        ),
        (
            "G2-4a vacuum edge mass UFO::m^(1/2)",
            vacuum("UFO::m^(1/2)", "1"),
            Expect::Refused,
        ),
        (
            "G2-4b vacuum edge mass log(UFO::m)",
            vacuum("log(UFO::m)", "1"),
            Expect::Refused,
        ),
        (
            "G2-5 a massless line with mass 0.0",
            parse(&k2_text.replacen("momentum=-1,-1,0\t0", "momentum=-1,-1,0\t0.0", 1))
                .and_then(plain),
            Expect::Observe,
        ),
        (
            "G2-6a the closing edge's mass written R^1",
            parse(&edge5_r(&format!("{R}^1"))).and_then(plain),
            Expect::Same(
                parse(&k2_text).and_then(plain),
                "caller::none",
                "caller::none",
            ),
        ),
        (
            "G2-6b the closing edge's mass written 1*R",
            parse(&edge5_r(&format!("1*{R}"))).and_then(plain),
            Expect::Same(
                parse(&k2_text).and_then(plain),
                "caller::none",
                "caller::none",
            ),
        ),
        // G3
        (
            "G3-1 hedge(ε, 1) in a vertex factor",
            in_field(
                &k2_text,
                "vertex factor",
                &format!("gammalooprs::{{spenso::index}}::hedge({EPS},1)"),
            ),
            Expect::Refused,
        ),
        (
            "G3-2 g(mink(D, a), R) in a vertex factor",
            in_field(
                &k2_text,
                "vertex factor",
                &format!(
                    "spenso::{{symmetric,linear,real,spenso::tensor}}::g(spenso::{{spenso::representation,spenso::self_dual}}::mink({DIM},gammalooprs::{{spenso::index}}::hedge(90,1)),{R})"
                ),
            ),
            Expect::Refused,
        ),
        (
            "G3-3 UFO::f(UFO::g(ε)) in the prefactor argument",
            in_argument(&k2_text, &format!("UFO::f(UFO::g({EPS}))")),
            Expect::Refused,
        ),
        (
            "G3-4 user::Q(1, ε) in a vertex factor",
            in_field(&k2_text, "vertex factor", &format!("user::Q(1,{EPS})")),
            Expect::Refused,
        ),
        (
            "G3-5 vacuum vertex factor x log(ε)",
            vacuum("caller::m2", &format!("log({EPS})")),
            Expect::Refused,
        ),
        (
            "G3-6 UFO::f(μᵣ²) in a vertex factor",
            in_field(&k2_text, "vertex factor", "UFO::f(gammalooprs::μᵣ²)"),
            Expect::Same(
                in_field(&k2_text, "vertex factor", "UFO::f(caller::mu)"),
                "caller::mu",
                "gammalooprs::μᵣ²",
            ),
        ),
        // exponents
        (
            "E-1 typed x^2.0 (a float exponent)",
            in_field(&k2_text, "vertex factor", "caller::x^2.0"),
            Expect::Refused,
        ),
        (
            "E-2 (1+UFO::x)^1000000000 in the prefactor argument",
            in_argument(&k2_text, "(1+UFO::x)^1000000000"),
            Expect::Refused,
        ),
        (
            "E-3 dim^1000000 in a vertex factor",
            in_field(&k2_text, "vertex factor", &format!("{DIM}^1000000")),
            Expect::Refused,
        ),
        (
            "E-4 typed (x-x)^0 in a vertex factor",
            parse(&k2_text).and_then(|mut v| {
                let x = Atom::var(symbolica::symbol!("caller::x"));
                let zero = &x - &x;
                let p = std::panic::catch_unwind(|| zero.pow(Atom::num(0)))
                    .map_err(|_| "pow panicked".to_string())?;
                let (_, put) = field_setters(&v)
                    .into_iter()
                    .find(|(f, _)| *f == "vertex factor")
                    .unwrap();
                put(&mut v, &p);
                plain(v)
            }),
            Expect::Observe,
        ),
        (
            "E-5 ε^(-10^30) in a vertex factor",
            in_field(
                &k2_text,
                "vertex factor",
                &format!("{EPS}^(-1000000000000000000000000000000)"),
            ),
            Expect::Refused,
        ),
        (
            "E-6 2^(1/2) in a projector",
            in_field(&k1_text, "projector", "2^(1/2)"),
            Expect::Refused,
        ),
        // UNKNOWN (allow-unknowns ON for the whole test)
        (
            "U-2 UNKNOWN(ε, 1) in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                &format!("vakint::UNKNOWN({EPS},1)"),
            ),
            Expect::Refused,
        ),
        (
            "U-3a UNKNOWN(UFO::m, 0)^2 in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                "vakint::UNKNOWN(UFO::m,0)^2",
            ),
            Expect::Refused,
        ),
        (
            "U-3b UNKNOWN(<master>, 0)^2 in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                &format!("vakint::UNKNOWN({ok_master},0)^2"),
            ),
            Expect::Refused,
        ),
        (
            "U-4a 1/UNKNOWN(UFO::m, 0) in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                "1/vakint::UNKNOWN(UFO::m,0)",
            ),
            Expect::Refused,
        ),
        (
            "U-4b 1/UNKNOWN(<master>, 0) in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                &format!("1/vakint::UNKNOWN({ok_master},0)"),
            ),
            Expect::Refused,
        ),
        (
            "U-5a mink(D, UNKNOWN(UFO::m, 0)) as an index in a vertex factor",
            in_field(
                &k2_text,
                "vertex factor",
                &format!(
                    "spenso::{{spenso::representation,spenso::self_dual}}::mink({DIM},vakint::UNKNOWN(UFO::m,0))"
                ),
            ),
            Expect::Refused,
        ),
        (
            "U-5b mink(D, UNKNOWN(<master>, 0)) as an index in an insertion factor",
            in_field(
                &k2_text,
                "insertion factor",
                &format!(
                    "spenso::{{spenso::representation,spenso::self_dual}}::mink({DIM},vakint::UNKNOWN({ok_master},0))"
                ),
            ),
            Expect::Refused,
        ),
        (
            "U-6 UNKNOWN(UFO::x_, 0) in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                "vakint::UNKNOWN(UFO::x_,0)",
            ),
            Expect::Refused,
        ),
        (
            "U-c control: UNKNOWN(<master>, 0) in an insertion",
            in_field(
                &k2_text,
                "insertion coefficient",
                &format!("vakint::UNKNOWN({ok_master},0)"),
            ),
            Expect::Observe,
        ),
    ];
    let mut wrong = Vec::new();
    for (name, input, expect) in rows {
        let started = std::time::Instant::now();
        let got = run(input);
        let seconds = started.elapsed().as_secs_f64();
        let shown = match &got {
            Ok(v) => format!("ACCEPTED ({} chars)", v.to_string().len()),
            Err(e) => format!("REFUSED {}", e.chars().take(160).collect::<String>()),
        };
        eprintln!("CONTRACT AUDIT {name}: {shown} [{seconds:.2}s]");
        if matches!(&got, Err(e) if e == "PANIC") {
            wrong.push(format!("{name}: a panic"));
            continue;
        }
        match expect {
            Expect::Refused if got.is_ok() => wrong.push(format!("{name}: accepted")),
            Expect::Same(control, pattern, image) => match (&got, run(control)) {
                (Ok(g), Ok(c)) => {
                    let c = c.replace(atom(pattern)).with(atom(image));
                    if !(g - &c).expand().is_zero() {
                        wrong.push(format!("{name}: a different value from the control"));
                    }
                }
                (g, c) => wrong.push(format!(
                    "{name}: {:?} vs control {:?}",
                    g.as_ref().err(),
                    c.err()
                )),
            },
            _ => {}
        }
        if seconds > 30.0 {
            wrong.push(format!("{name}: {seconds:.1} s"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
