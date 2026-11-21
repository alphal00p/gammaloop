//! PROBES of the agreed record format against the code (tests only). Each
//! probe states the EXPECTED behaviour from the agreed text
//! (`WH-RECORD-FORMAT.md` v1.4, "v1.4 §n", and `WH-RECORD-FORMAT-CORRECTIONS.md`,
//! "corr. n") and asserts it through the PUBLIC entries on real records (the
//! vendored fixtures; the producer's fixtures under `VAKINT_TEST_WH_FIXTURES`
//! for the `#[ignore]`d ones). A probe that documents a FINDING (the code does
//! something else) is `#[ignore = "finding: ..."]`, so that the suite stays
//! green; the finding is described in the lane's probe report.

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

/// A fully qualified expression (the probes write every name with its
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
fn p1_numerical_mode_k1_root_equals_symbolic_numerically() {
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
fn p1_retain_k1_root_without_insertions_refinalises_to_symbolic() {
    let v = record("wh-k1-root");
    let exact = value(&v, &symbolic());
    let top = top_order(&exact);
    let retained = closure(&v, &site(RustRedMasterMode::Retain), &v.prefactor.total())
        .unwrap_or_else(|e| panic!("Retain: {e}"))
        .value;
    eprintln!("P1 RETAIN k1-root value: {retained}");
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
            eprintln!("P1 RETAIN expanded: {series}");
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
#[ignore = "finding F-4: vakint's own UNKNOWN atom in an insertion is refused as a reserved symbol"]
fn p1_unknown_from_vakint_in_an_insertion_is_carried_linearly() {
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
fn p1_unknown_caller_symbol_in_an_insertion_is_carried_linearly() {
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
/// coefficients, before the truncation). The probe records which: exact
/// product (no truncation of the argument's orders) or truncated.
#[test]
fn p1_an_argument_power_beyond_depth_is_truncated() {
    let v = record("wh-k1-root");
    let plain = value(&v, &symbolic());
    let top = top_order(&plain);
    let x = atom(&format!("gammalooprs::x*{EPS}^2"));
    let got = closure(&v, &symbolic(), &(v.prefactor.total() * &x))
        .unwrap()
        .value;
    eprintln!("P1 ARGUMENT eps^2: plain top order {top}; got {got}");
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
fn p1_a_scaleless_remainder_is_exactly_zero() {
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
        "label\tprobe scaleless remainder".into(),
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
fn p1_never_nf_with_original_cycles() {
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
/// n_f in the value (and no vakint:: symbol, corr. 1). The probe records it.
#[test]
fn p1_never_nf_without_original_cycles() {
    let base = record("wh-k1-root");
    let mut v = base.clone();
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.original_cycles = None;
    }
    v.prefactor.loop_sign = -1; // the carried graph's ghost loop (v1 standalone rule)
    match v.validate() {
        Ok(()) => {}
        Err(e) => {
            eprintln!("P1 NEVER-NF: without original_cycles the record is refused: {e}");
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
    atom("probe::z")
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
fn p2_z_in_the_argument_scales_every_kind() {
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
fn p2_z_in_a_field_scales_every_kind() {
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
fn p2_one_over_one_plus_z_eps_is_the_truncated_product() {
    let factor = atom(&format!("(1+probe::z*{EPS})^-1"));
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
fn p2_z_scales_in_retain() {
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
fn p2_one_over_one_plus_z_eps_is_exact_in_retain_without_insertions() {
    let retain = site(RustRedMasterMode::Retain);
    let factor = atom(&format!("(1+probe::z*{EPS})^-1"));
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
fn p2_two_halves_sum_to_the_whole() {
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
fn p2_a_unit_insertion_equals_its_factor_on_the_vertex() {
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
fn p2_record_order_is_irrelevant() {
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
/// two) -- covered by the external probe below.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn p2_insertion_order_is_irrelevant_k5_stress() {
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
            eprintln!("P2 insertion order: refused: {detail}")
        }
        other => panic!("reversed insertions: {other:?}"),
    }
}

/// (e) relabelling the record's DUMMY indices (a contracted pair in a vertex
/// factor, `edge(e, k)` labels) and renumbering an INTERNAL edge: the value
/// is invariant.
#[test]
fn p2_relabelling_dummies_and_internal_edges_is_irrelevant() {
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
            Err(e) => eprintln!("P2 relabel: refused: {e}"),
        }
    }
}

// ===========================================================================
// P3: the same function written differently
// ===========================================================================

/// An insertion coefficient c ε^-1 written as -2 c/(D-4) and as c (D-4+2ε+1)/ε
/// (all equal after D -> 4-2ε, v1.4 §5 / corr. 4): the same value.
#[test]
fn p3_equal_coefficients_give_equal_values() {
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
fn p3_equal_arguments_give_equal_values() {
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
fn p3_text_round_trip_gives_the_same_value() {
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
fn p3_vakint_epsilon_under_another_settings_epsilon_is_refused() {
    let base = record("wh-k2-child-ins");
    let mut v = base.clone();
    v.insertions[0].coefficient = v.insertions[0]
        .coefficient
        .replace(Atom::var(eps()))
        .with(atom("vakint::ε"));
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("P3 vakint::ε: refused: {}", detail(&e)),
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
fn p4_the_reference_mass_inside_a_function() {
    let base = record("wh-k1-root");
    let mut v = base.clone();
    let vertex = v
        .vertices
        .iter_mut()
        .find(|x| x.factor != Atom::num(1))
        .unwrap();
    vertex.factor = &vertex.factor * atom(&format!("log({R})"));
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("P4 log(R): refused: {}", detail(&e)),
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
fn p4_an_epsilon_dependent_auxiliary_mass() {
    let base = record("wh-k1-root");
    let plain = value(&base, &symbolic());
    let mut v = base.clone();
    let shifted = atom(&format!("{MAUX}*(1+{EPS})"));
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = shifted.clone();
    }
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("P4 M_aux(1+eps): refused: {}", detail(&e)),
        Ok(got) => {
            eprintln!("P4 M_aux(1+eps): accepted: {}", got.value);
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
fn p4_an_auxiliary_mass_zero_after_the_substitution() {
    let mut v = record("wh-k1-root");
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = atom(&format!("{MAUX}*({DIM}-4+2*{EPS})"));
    }
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("P4 M_aux zero after D->4-2eps: refused: {}", detail(&e)),
        Ok(got) => panic!(
            "an auxiliary mass zero after D -> 4-2eps is accepted: {}",
            got.value
        ),
    }
}

/// A D-dependent auxiliary mass M_aux² (D-2)/2 (not zero): expected refused
/// or the substituted value; recorded.
#[test]
fn p4_a_dimension_dependent_auxiliary_mass() {
    let mut v = record("wh-k1-root");
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        spec.auxiliary_mass_squared = atom(&format!("{MAUX}*({DIM}-2)/2"));
    }
    match closure(&v, &symbolic(), &v.prefactor.total()) {
        Err(e) => eprintln!("P4 M_aux (D-2)/2: refused: {}", detail(&e)),
        Ok(got) => panic!("a D-dependent auxiliary mass is accepted: {}", got.value),
    }
}

/// A reference mass written as a product (2 R, not a symbol) is refused
/// (`the reference mass ... is not a symbol`); a reference mass that is a
/// different symbol than the closing edge's is refused.
#[test]
fn p4_reference_mass_not_a_symbol_or_two_masses() {
    let base = record("wh-k1-root");
    let text = base.to_text().unwrap();
    let doubled = text.replace(R, &format!("(2*{R})"));
    match FactorizedNumerator::from_text(&doubled) {
        Err(e) => eprintln!("P4 2R: refused: {e}"),
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
/// refusal, or a DIFFERENT value (a finding).
fn text_variant(clean: &str, variant: &str, what: &str) -> Result<(), String> {
    let base = FactorizedNumerator::from_text(clean).unwrap();
    let plain = value(&base, &symbolic());
    match FactorizedNumerator::from_text(variant) {
        Err(e) => {
            eprintln!("P5 {what}: refused: {e}");
            Ok(())
        }
        Ok(v) => match closure(&v, &symbolic(), &v.prefactor.total()) {
            Err(e) => {
                eprintln!("P5 {what}: read, refused by the entry: {e}");
                Ok(())
            }
            Ok(got) => {
                if (&got.value - &plain).expand().is_zero() {
                    eprintln!("P5 {what}: accepted, same value");
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
fn p5_line_ends_bom_trailing_whitespace_final_newline() {
    let clean = fixture_text("wh-k2-child-ins.record.txt");
    let mut findings = Vec::new();
    for (what, variant) in [
        ("CRLF", clean.replace('\n', "\r\n")),
        ("BOM", format!("\u{feff}{clean}")),
        ("trailing spaces", clean.replace('\n', "  \n")),
        ("trailing tabs", clean.replace('\n', "\t\n")),
        ("no final newline", clean.trim_end_matches('\n').to_string()),
        ("an extra empty final line", format!("{clean}\n")),
    ] {
        if let Err(f) = text_variant(&clean, &variant, what) {
            findings.push(f);
        }
    }
    assert!(findings.is_empty(), "{}", findings.join("\n\n"));
}

/// Other Unicode spellings of ε in the insertion coefficient (U+03F5 ϵ,
/// U+1D700 𝜀, U+1D716 𝜖): not the settings' epsilon `gammalooprs::ε`
/// (U+03B5): refused (v1.4 §8 "an eps symbol other than the settings'
/// epsilon_symbol: Insertion / Evaluation"), never a constant.
#[test]
#[ignore = "finding F-2: epsilon spelled U+1D700 / U+1D716 is accepted as a constant"]
fn p5_other_spellings_of_epsilon() {
    let clean = fixture_text("wh-k2-child-ins.record.txt");
    let mut findings = Vec::new();
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
            findings.push(f);
        }
    }
    assert!(findings.is_empty(), "{}", findings.join("\n\n"));
}

/// Other spellings of 𝜋 / 𝑖: `π` (U+03C0) and `i` are not Symbolica's
/// built-ins (corr. 3): read into vakint:: and refused, never a different
/// constant.
#[test]
fn p5_other_spellings_of_pi_and_i() {
    let clean = fixture_text("wh-k2-child-ins.record.txt");
    let line = clean
        .lines()
        .find(|l| l.starts_with("insertion\t"))
        .unwrap();
    let mut findings = Vec::new();
    for (what, from, to) in [
        ("π U+03C0", "𝜋", "π"),
        ("bare i", "𝑖", "i"),
        ("pi alias", "𝜋", "pi"),
    ] {
        let variant = clean.replace(line, &line.replace(from, to));
        if let Err(f) = text_variant(&clean, &variant, what) {
            findings.push(f);
        }
    }
    assert!(findings.is_empty(), "{}", findings.join("\n\n"));
}

/// Duplicated / reordered sections and numeric-field spellings: refused, or
/// accepted with the same value.
#[test]
fn p5_sections_and_numbers() {
    let clean = fixture_text("wh-k1-root.record.txt");
    let mut findings = Vec::new();
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
            findings.push(f);
        }
    }
    assert!(findings.is_empty(), "{}", findings.join("\n\n"));
}

/// Non-ASCII names in a namespace-qualified caller symbol: accepted with
/// the same value up to the renaming.
#[test]
fn p5_non_ascii_names() {
    let clean = fixture_text("wh-k1-root.record.txt");
    let variant = clean.replace("UFO::{}::gs_minus", "UFO::{}::gs_mïnüs");
    assert_ne!(variant, clean);
    let base = FactorizedNumerator::from_text(&clean).unwrap();
    let v = FactorizedNumerator::from_text(&variant).unwrap();
    let renamed = value(&base, &symbolic())
        .replace(atom("UFO::gs_minus"))
        .with(atom("UFO::gs_mïnüs"));
    assert_same(&value(&v, &symbolic()), &renamed, "non-ASCII name");
}

/// The direct / ports entry for a marked K5 direct record: z in the argument
/// scales the value (P2 for the direct kind).
#[test]
fn p2_direct_record_z_in_the_argument() {
    let v = record("wh-k5-direct");
    let vakint = engine();
    let plain = evaluate_factorized_ports(&vakint, &symbolic(), &v, &v.prefactor.total())
        .unwrap()
        .value;
    let got = evaluate_factorized_ports(&vakint, &symbolic(), &v, &(v.prefactor.total() * z()))
        .unwrap()
        .value;
    assert_same(&got, &(&plain * z()), "k5-direct argument z");
    let mut i = v.clone();
    i.insertions[0].coefficient = &i.insertions[0].coefficient * z();
    let got = evaluate_factorized_ports(&vakint, &symbolic(), &i, &i.prefactor.total())
        .unwrap()
        .value;
    assert_same(&got, &(&plain * z()), "k5-direct insertion z");
}

/// P3: equal forms of 1 multiplied into a VERTEX FACTOR (FactorRules admits
/// tensor-free scalars rational in D and ε): the same value.
#[test]
fn p3_equal_factor_forms_give_equal_values() {
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
fn p1_every_kind_returns_the_callers_notation() {
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
fn p1_linear_spectator_marker() {
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
fn p2_k6_z_in_the_argument() {
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
    eprintln!("P2 K6 compared {compared}");
}

/// F-1 extended: R under a NEGATIVE or a non-integer power in a factor
/// (`1/R`, `sqrt(R)`: FactorRules admits tensor-free scalars): refused, or a
/// value free of vakint's internal symbols.
#[test]
fn p4_the_reference_mass_under_a_power() {
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
            Err(e) => eprintln!("P4 {form}: refused: {}", detail(&e)),
            Ok(got) => {
                let text = got.value.to_string();
                eprintln!("P4 {form}: accepted: {text}");
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
fn p3_the_settings_epsilon_under_another_name() {
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
fn p2_moving_a_coefficient_between_insertions_k5_stress() {
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
fn p1_numerical_mode_with_an_insertion() {
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
    let only = std::env::var("PROBE_INVENTORY_SETS")
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

/// THE RESULT CHECK ALONE catches F-1 and F-3 through the public closure
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
