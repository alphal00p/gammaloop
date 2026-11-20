//! Counterterm insertions ([`Insertion`]) against the CONTRACTED route.
//!
//! Fixtures are built from the two in-tree L1 records (background-field gluon
//! and quark bubbles): one internal line is replaced by an insertion
//! surrogate of `n` scalar bubbles in series (n = 1, 2, 3: L2, L3, L4), its
//! coefficient vertex carrying `coefficient × factor` with the factor's slots
//! on the two inner boundary half-edges (one at each end of the chain), as the
//! R* insertion surrogates have them; and an R* child shape: a contracted
//! vertex (`ct:`) whose factor is a vertex polynomial with ε poles.
//!
//! The reference is the contracted route of the caller: the same coefficient ×
//! factor multiplied by hand into the coefficient vertex's own factor (what
//! the surrogate's `with_coefficient` does), contracted by idenso, mapped and
//! evaluated by `Vakint::evaluate` with no contract check. Values must be
//! IDENTICAL (exact difference zero).
use super::*;
use crate::utils::vakint_macros::vk_parse;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};

const GLUON: &str = include_str!("fixtures/l1-gluon-bubble.txt");
const QUARK: &str = include_str!("fixtures/l1-quark-bubble.txt");
const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const BIS: &str = "spenso::{spenso::representation,spenso::self_dual}::bis";
const METRIC: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
const GAMMA: &str = "spenso::{linear,spenso::tensor}::gamma";
const MOMENTUM: &str = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const EDGE: &str = "gammalooprs::{spenso::index,spenso::index-label:e}::edge";
const DIM: &str = "gammalooprs::{}::dim";
const EPS: &str = "vakint::{}::ε";

fn atom(text: &str) -> Atom {
    crate::Vakint::initialize_vakint_symbols();
    vk_parse!(text).unwrap_or_else(|e| panic!("{e}: {text}"))
}

fn read(text: &str) -> FactorizedNumerator {
    crate::Vakint::initialize_vakint_symbols();
    FactorizedNumerator::from_text(text).unwrap_or_else(|e| panic!("{e}\n{text}"))
}

fn mink(h: usize) -> String {
    format!("{MINK}({DIM},{HEDGE}({h},1))")
}

fn bis(h: usize) -> String {
    format!("{BIS}(4,{HEDGE}({h},1))")
}

fn dummy(e: usize, k: usize) -> String {
    format!("{MINK}({DIM},{EDGE}({e},{k}))")
}

fn q(e: usize, slot: &str) -> String {
    format!("{MOMENTUM}({e},{slot})")
}

fn metric(a: &str, b: &str) -> String {
    format!("{METRIC}({a},{b})")
}

fn gamma(row: &str, column: &str, mu: &str) -> String {
    format!("{GAMMA}({row},{column},{mu})")
}

/// The vertex of the coefficient (the chain's first vertex).
const C0: usize = 10;
/// The two inner boundary half-edges of the surrogate.
const H_IN: usize = 40;
const H_OUT: usize = 41;
/// The edge leaving the chain.
const OUT: usize = 30;

/// `value` with internal edge 3 replaced by `bubbles` scalar bubbles in
/// series: 3 = (source of 3) -> C0, the chain C0 -> .. -> C0+n, OUT = C0+n ->
/// (target of 3); the coefficient vertex C0 has rule `ins:se` and factor 1.
/// `boundary(e_in, e_out)` gives the two boundary edges' new factors.
fn surrogate(
    value: &FactorizedNumerator,
    bubbles: usize,
    boundary: impl Fn(&FactorizedEdge) -> (Atom, Atom),
) -> FactorizedNumerator {
    let mut out = value.clone();
    let old_loops = value.loops();
    let loops = old_loops + bubbles;
    let widen = |m: &[i64]| {
        let mut w = m[..old_loops].to_vec();
        w.extend(std::iter::repeat_n(0, bubbles));
        w.extend_from_slice(&m[old_loops..]);
        w
    };
    for edge in &mut out.edges {
        edge.momentum = widen(&edge.momentum);
    }
    let three = out.edges.iter().position(|e| e.id == 3).unwrap();
    let original = out.edges.remove(three);
    let (into, from) = boundary(&original);
    let m3 = original.momentum.clone();
    let mut first = original.clone();
    first.target = Endpoint::Vertex(C0);
    first.target_hedge = Some(H_IN);
    first.factor = into;
    let mut last = original.clone();
    last.id = OUT;
    last.source = Endpoint::Vertex(C0 + bubbles);
    last.source_hedge = Some(H_OUT);
    last.factor = from;
    out.edges.push(first);
    out.edges.push(last);
    for i in 0..bubbles {
        let mut b = vec![0; loops + 1];
        b[old_loops + i] = 1;
        let rest = m3.iter().zip(&b).map(|(m, b)| m - b).collect::<Vec<_>>();
        for (k, momentum) in [(0, b), (1, rest)] {
            out.edges.push(FactorizedEdge {
                id: 50 + 2 * i + k,
                source: Endpoint::Vertex(C0 + i),
                target: Endpoint::Vertex(C0 + i + 1),
                source_hedge: Some(60 + 4 * i + 2 * k),
                target_hedge: Some(61 + 4 * i + 2 * k),
                particle: "s".into(),
                line: LineKind::Scalar,
                power: 1,
                mass_squared: Atom::Zero,
                momentum,
                factor: Atom::num(1),
            });
        }
        out.loop_edges.push(Some(50 + 2 * i));
    }
    for i in 0..=bubbles {
        out.vertices.push(FactorizedVertex {
            id: C0 + i,
            rule: if i == 0 {
                "ins:se".into()
            } else {
                format!("scalar_insertion_{i}")
            },
            factor: Atom::num(1),
        });
    }
    out.edges.sort_by_key(|e| e.id);
    out.vertices.sort_by_key(|v| v.id);
    out.label = format!("{} + {bubbles}-bubble surrogate on edge 3", value.label);
    out
}

/// The gluon bubble with a gluon self-energy surrogate on edge 3:
/// `X = a g(in, out) Q(3)^2 + b Q(3, in) Q(OUT, out)`, `a = (D-2)/(D-1)`,
/// `b = -1`, coefficient `ε^-1 + 1/3`.
fn gluon_surrogate(bubbles: usize) -> FactorizedNumerator {
    let base = read(GLUON);
    let mut value = surrogate(&base, bubbles, |e| {
        let (s, t) = (e.source_hedge.unwrap(), e.target_hedge.unwrap());
        (
            atom(&metric(&mink(s), &mink(H_IN))),
            atom(&metric(&mink(H_OUT), &mink(t))),
        )
    });
    let mu = dummy(3, 7);
    let factor = format!(
        "({DIM}-2)/({DIM}-1)*{}*{}^2-{}*{}",
        metric(&mink(H_IN), &mink(H_OUT)),
        q(3, &mu),
        q(3, &mink(H_IN)),
        q(OUT, &mink(H_OUT))
    );
    value.insertions.push(Insertion {
        vertex: C0,
        hedges: vec![H_IN, H_OUT],
        cycles: (0, 0),
        coefficient: atom(&format!("{EPS}^-1+1/3")),
        factor: atom(&factor),
    });
    value
}

/// The quark bubble with a quark self-energy surrogate on edge 3:
/// `X = gamma(out, in, mu) Q(3, mu)` (the opened slash), coefficient
/// `(D-4)/(D-3) + ε^-1`; the Dirac line passes through the chain.
fn quark_surrogate(bubbles: usize) -> FactorizedNumerator {
    let base = read(QUARK);
    let mut value = surrogate(&base, bubbles, |e| {
        let s = e.source_hedge.unwrap();
        let t = e.target_hedge.unwrap();
        let slash = |edge: usize, row: usize, column: usize| {
            let mu = dummy(edge, 1);
            atom(&format!(
                "{}*{}",
                q(edge, &mu),
                gamma(&bis(row), &bis(column), &mu)
            ))
        };
        (slash(3, H_IN, s), slash(OUT, t, H_OUT))
    });
    let mu = dummy(3, 7);
    value.insertions.push(Insertion {
        vertex: C0,
        hedges: vec![H_IN, H_OUT],
        cycles: (0, 0),
        coefficient: atom(&format!("({DIM}-4)/({DIM}-3)+{EPS}^-1")),
        factor: atom(&format!(
            "{}*{}",
            gamma(&bis(H_OUT), &bis(H_IN), &mu),
            q(3, &mu)
        )),
    });
    value
}

/// [`gluon_surrogate`] with vertex 1 also a contracted child (`ct:`): its
/// vertex polynomial (momentum of edge 3 renamed to OUT, the edge now at
/// vertex 1) becomes an insertion with coefficient `ε^-2 - 2 ε^-1`.
fn gluon_surrogate_and_child(bubbles: usize) -> FactorizedNumerator {
    let mut value = gluon_surrogate(bubbles);
    let vertex = value.vertices.iter_mut().find(|v| v.id == 1).unwrap();
    let factor = vertex
        .factor
        .to_canonical_string()
        .replace(&format!("{MOMENTUM}(3,"), &format!("{MOMENTUM}({OUT},"));
    vertex.factor = Atom::num(1);
    vertex.rule = "ct:child".into();
    let hedges = value
        .local_hedges(Location::Vertex(1))
        .into_iter()
        .collect::<Vec<_>>();
    value.insertions.insert(
        0,
        Insertion {
            vertex: 1,
            hedges,
            cycles: (0, 0),
            coefficient: atom(&format!("{EPS}^-2-2*{EPS}^-1")),
            factor: atom(&factor),
        },
    );
    value.label = format!("{} + ct: child at vertex 1", value.label);
    value
}

/// The caller's contracted route: every insertion multiplied by hand into its
/// vertex's own factor, the value contracted by idenso and evaluated by
/// `Vakint::evaluate`, no contract check (the surrogate vertex is not local).
fn contracted_route(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
) -> (Atom, f64) {
    let mut hand = value.clone();
    for insertion in std::mem::take(&mut hand.insertions) {
        let vertex = hand
            .vertices
            .iter_mut()
            .find(|v| v.id == insertion.vertex)
            .unwrap();
        vertex.factor = &vertex.factor * &insertion.coefficient * &insertion.factor;
    }
    let started = std::time::Instant::now();
    let (input, _) =
        super::evaluate::vakint_input_unchecked(settings, &hand, &Atom::num(1)).unwrap();
    let value = vakint.evaluate(settings, input.as_view()).unwrap();
    (value, started.elapsed().as_secs_f64())
}

fn settings(terms: i64) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: "/nonexistent/helper-vh/form".into(),
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: terms,
        run_time_decimal_precision: 40,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

/// The built-in p library (one to three loops), or the 1-4 loop package of
/// `VAKINT_TEST_P_LIBRARY_1_4` (needed at four loops).
fn engine() -> (Vakint, bool) {
    let mut engine = Vakint::new().unwrap();
    let Some(directory) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") else {
        return (engine, false);
    };
    let library = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(directory.into()),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap();
    engine.inject_p_library(library).unwrap();
    (engine, true)
}

/// The insertion route equals the contracted route exactly; also the text
/// round trip of the record. Prints one timing line (median of 3 paired
/// repetitions, single core).
fn identical(vakint: &Vakint, settings: &VakintSettings, value: &FactorizedNumerator) {
    value
        .validate()
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    let text = value.to_text().unwrap();
    let again = FactorizedNumerator::from_text(&text).unwrap();
    assert_eq!(&again, value, "{}: text round trip", value.label);
    assert_eq!(again.to_text().unwrap(), text);
    let (mut ours, mut theirs) = (Vec::new(), Vec::new());
    let mut got = None;
    let mut reference = None;
    for _ in 0..3 {
        let started = std::time::Instant::now();
        let value_ours = evaluate_factorized(vakint, settings, value, &Atom::num(1))
            .unwrap_or_else(|e| panic!("{}: {e}", value.label));
        ours.push(started.elapsed().as_secs_f64());
        let (value_theirs, seconds) = contracted_route(vakint, settings, value);
        theirs.push(seconds);
        got = Some(value_ours.value);
        reference = Some(value_theirs);
    }
    let (got, reference) = (got.unwrap(), reference.unwrap());
    assert!(!reference.is_zero(), "{}: zero reference", value.label);
    assert!(
        (&got - &reference).expand().is_zero(),
        "{}:\ninsertion  {got}\ncontracted {reference}",
        value.label
    );
    let median = |v: &mut Vec<f64>| {
        v.sort_by(f64::total_cmp);
        v[1]
    };
    eprintln!(
        "INSERTION_IDENTITY label=\"{}\" loops={} identical=true uncontracted_s={:.4} contracted_s={:.4} value={}",
        value.label,
        value.loops(),
        median(&mut ours),
        median(&mut theirs),
        got
    );
}

#[test]
fn gluon_surrogates_at_two_three_and_four_loops_equal_the_contracted_route() {
    let (vakint, four) = engine();
    let settings = settings(3);
    for bubbles in [1, 2, 3] {
        if bubbles == 3 && !four {
            eprintln!("INSERTION_IDENTITY SKIPPED L4: VAKINT_TEST_P_LIBRARY_1_4 not set");
            continue;
        }
        identical(&vakint, &settings, &gluon_surrogate(bubbles));
    }
}

#[test]
fn quark_surrogates_with_a_dirac_line_through_the_chain_equal_the_contracted_route() {
    let (vakint, four) = engine();
    let settings = settings(3);
    for bubbles in [1, 2, 3] {
        if bubbles == 3 && !four {
            eprintln!("INSERTION_IDENTITY SKIPPED L4: VAKINT_TEST_P_LIBRARY_1_4 not set");
            continue;
        }
        let value = quark_surrogate(bubbles);
        // the fermion loop runs through the chain (one node for the flow)
        let loops = value.closed_loops().unwrap();
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].edges, vec![2, 3, OUT]);
        identical(&vakint, &settings, &value);
    }
}

#[test]
fn a_contracted_child_with_poles_and_a_surrogate_in_one_record_equal_the_contracted_route() {
    let (vakint, _) = engine();
    let settings = settings(3);
    for bubbles in [1, 2] {
        identical(&vakint, &settings, &gluon_surrogate_and_child(bubbles));
    }
}

/// With no insertion, the text form and every value are those of the
/// unchanged records (the fixtures keep `insertions 0`).
#[test]
fn records_without_insertions_are_unchanged() {
    for text in [GLUON, QUARK] {
        let value = read(text);
        assert!(value.insertions.is_empty());
        assert_eq!(value.to_text().unwrap(), text);
        let product = value.factor_product();
        let mut by_hand = Atom::num(1);
        for edge in &value.edges {
            by_hand *= &edge.factor;
        }
        for vertex in &value.vertices {
            by_hand *= &vertex.factor;
        }
        assert_eq!(product, by_hand);
    }
}

fn refused(value: &FactorizedNumerator) -> FactorizedError {
    value
        .validate()
        .expect_err(&format!("{} must be refused", value.label))
}

fn ins(value: &mut FactorizedNumerator) -> &mut Insertion {
    &mut value.insertions[0]
}

fn is_insertion(error: &FactorizedError) -> bool {
    matches!(error, FactorizedError::Insertion { .. })
}

#[test]
fn every_insertion_rule_has_its_refusal() {
    let good = gluon_surrogate(1);
    good.validate().unwrap();
    let with = |change: &dyn Fn(&mut FactorizedNumerator)| {
        let mut value = good.clone();
        change(&mut value);
        value
    };
    // structure
    for (what, value) in [
        ("missing vertex", with(&|v| ins(v).vertex = 999)),
        (
            "vertex factor not 1",
            with(&|v| v.vertices.iter_mut().find(|x| x.id == C0).unwrap().factor = Atom::num(2)),
        ),
        ("empty hedges", with(&|v| ins(v).hedges.clear())),
        (
            "unsorted hedges",
            with(&|v| ins(v).hedges = vec![H_OUT, H_IN]),
        ),
        (
            "unknown hedge",
            with(&|v| ins(v).hedges = vec![H_IN, H_OUT, 999]),
        ),
        (
            "no hedge at the vertex",
            with(&|v| {
                ins(v).hedges = vec![H_OUT];
            }),
        ),
        (
            "rule neither ct: nor ins:",
            with(&|v| v.vertices.iter_mut().find(|x| x.id == C0).unwrap().rule = "V_x".into()),
        ),
        (
            "ct: with hedges other than the incident ones",
            with(&|v| v.vertices.iter_mut().find(|x| x.id == C0).unwrap().rule = "ct:x".into()),
        ),
        (
            "two insertions on one vertex",
            with(&|v| {
                let copy = v.insertions[0].clone();
                v.insertions.push(copy);
            }),
        ),
        (
            "colour invariant in the coefficient",
            with(&|v| ins(v).coefficient = atom("spenso::cas(2,spenso::coad(spenso::dA))")),
        ),
        (
            "flavour symbol in the coefficient",
            with(&|v| ins(v).coefficient = atom("vakint::nf")),
        ),
        (
            "mUVexp in the factor",
            with(&|v| ins(v).factor = &ins(v).factor * atom("gammalooprs::mUVexp")),
        ),
        (
            "degree four",
            with(&|v| ins(v).factor = &ins(v).factor * atom(&format!("{}^2", q(3, &dummy(3, 8))))),
        ),
    ] {
        let error = refused(&value);
        assert!(is_insertion(&error), "{what}: {error}");
        eprintln!("INSERTION_REFUSAL {what}: {error}");
    }
    // five open slots: the contracted child at vertex 1 (three half-edges)
    // with two more spinor slots
    let mut five = gluon_surrogate_and_child(1);
    five.validate().unwrap();
    ins(&mut five).factor = &ins(&mut five).factor * atom(&metric(&bis(0), &bis(3)));
    let error = refused(&five);
    assert!(is_insertion(&error), "five slots: {error}");
    eprintln!("INSERTION_REFUSAL five open slots: {error}");
    // a non-scalar coefficient, a non-boundary momentum, a non-local slot
    let error = refused(&with(&|v| {
        ins(v).coefficient = atom(&format!("{}^2", q(3, &dummy(3, 8))))
    }));
    assert!(
        matches!(error, FactorizedError::NotScalar { .. }),
        "{error}"
    );
    let error = refused(&with(&|v| {
        ins(v).factor = atom(&format!("{}*{}", q(2, &mink(H_IN)), q(OUT, &mink(H_OUT))))
    }));
    assert!(
        matches!(error, FactorizedError::UnsupportedTensor { .. }),
        "{error}"
    );
    let error = refused(&with(&|v| {
        ins(v).factor = atom(&format!("{}*{}", q(3, &mink(H_IN)), q(OUT, &mink(4))))
    }));
    assert!(
        matches!(
            error,
            FactorizedError::IndexNotLocal { .. } | FactorizedError::IndexOverContracted { .. }
        ),
        "{error}"
    );
    // weights: a closed cycle inside the insertion changes the loop sign
    let error = refused(&with(&|v| ins(v).cycles = (1, 0)));
    assert!(matches!(error, FactorizedError::LoopSign { .. }), "{error}");
    let mut balanced = with(&|v| ins(v).cycles = (1, 0));
    balanced.prefactor.loop_sign = -1;
    let error = refused(&balanced);
    assert!(matches!(error, FactorizedError::NfPower { .. }), "{error}");
    balanced.prefactor.nf_power = 1;
    balanced.validate().unwrap();
    // a second epsilon at evaluation
    let wrong = with(&|v| ins(v).coefficient = atom("gammalooprs::ε^-1"));
    wrong.validate().unwrap();
    let error = super::evaluate::vakint_input(&settings(2), &wrong, &Atom::num(1)).unwrap_err();
    assert!(is_insertion(&error), "{error}");
}

#[test]
fn damaged_insertion_lines_are_refused_by_the_reader() {
    let text = gluon_surrogate(1).to_text().unwrap();
    assert!(text.contains("insertions\t1\ninsertion\t10\thedges=40,41\tcycles=0,0\t"));
    let fewer = text.replace("insertions\t1", "insertions\t2");
    assert!(matches!(
        FactorizedNumerator::parse_text(&fewer),
        Err(FactorizedError::TextSectionCount { .. })
    ));
    let mismatch = text.replace("insertion_factor\t10\t", "insertion_factor\t11\t");
    assert!(matches!(
        FactorizedNumerator::parse_text(&mismatch),
        Err(FactorizedError::TextSyntax { .. })
    ));
    let cycles = text.replace("cycles=0,0", "cycles=0");
    assert!(matches!(
        FactorizedNumerator::parse_text(&cycles),
        Err(FactorizedError::TextSyntax { .. })
    ));
    let missing = text.replace("\tcycles=0,0", "");
    assert!(FactorizedNumerator::parse_text(&missing).is_err());
}

/// The vacuum entry: a unit-mass two-loop sunrise whose vertex 1 is a
/// contracted child (`ct:`) carrying `(ε^-1 + 2) × (k1 + (D-2) p)^mu` with an
/// open port, against the same value with the insertion multiplied into the
/// vertex factor by hand.
#[test]
fn a_contracted_child_in_the_vacuum_entry_equals_the_contracted_route() {
    let mink_port = |h: usize| mink(h);
    let edge = |id: usize, s: usize, t: usize, h: (usize, usize), m: [i64; 3], factor: Atom| {
        FactorizedEdge {
            id,
            source: Endpoint::Vertex(s),
            target: Endpoint::Vertex(t),
            source_hedge: Some(h.0),
            target_hedge: Some(h.1),
            particle: "x".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: atom("muvsq"),
            momentum: m.to_vec(),
            factor,
        }
    };
    let leg = FactorizedEdge {
        id: 20,
        source: Endpoint::External,
        target: Endpoint::Vertex(1),
        source_hedge: None,
        target_hedge: Some(120),
        particle: "x".into(),
        line: LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: vec![0, 0, 1],
        factor: Atom::num(1),
    };
    // degree <= 2 in the boundary momenta (the default limit)
    let polynomial = atom(&format!(
        "{}+({DIM}-2)*{}",
        q(1, &mink_port(10)),
        q(20, &mink_port(10)),
    ));
    let coefficient = atom(&format!("{EPS}^-1+2"));
    let mut value = FactorizedNumerator {
        label: "2L sunrise, ct: child at vertex 1".into(),
        kind: FactorizedKind::Direct,
        dimension: FactorizedNumerator::default_dimension(),
        loop_edges: vec![Some(1), Some(2)],
        external_edges: vec![Some(20)],
        edges: vec![
            edge(1, 1, 2, (1, 2), [1, 0, 0], Atom::num(1)),
            edge(2, 1, 2, (3, 4), [0, 1, 0], Atom::num(1)),
            edge(3, 2, 1, (5, 6), [1, 1, 0], atom(&q(3, &mink_port(11)))),
            leg,
        ],
        vertices: vec![
            FactorizedVertex {
                id: 1,
                rule: "ct:child".into(),
                factor: Atom::num(1),
            },
            FactorizedVertex {
                id: 2,
                rule: "test".into(),
                factor: Atom::num(1),
            },
        ],
        ports: vec![atom(&mink_port(10)), atom(&mink_port(11))],
        projector: None,
        prefactor: Prefactor::default(),
        colour_branch: None,
        insertions: Vec::new(),
        record: None,
    };
    value.insertions.push(Insertion {
        vertex: 1,
        hedges: vec![1, 3, 6, 120],
        cycles: (0, 0),
        coefficient: coefficient.clone(),
        factor: polynomial.clone(),
    });
    let mut hand = value.clone();
    hand.insertions.clear();
    hand.vertices[0].factor = &coefficient * &polynomial;
    let vakint = Vakint::new().unwrap();
    let mut settings = settings(3);
    settings.integral_normalization_factor = LoopNormalizationFactor::MSbar;
    let got = evaluate_factorized_vacuum(&vakint, &settings, &value, &Atom::num(1)).unwrap();
    let reference = evaluate_factorized_vacuum(&vakint, &settings, &hand, &Atom::num(1)).unwrap();
    assert!(!reference.value.is_zero());
    assert!(
        (&got.value - &reference.value).expand().is_zero(),
        "insertion {}\ncontracted {}",
        got.value,
        reference.value
    );
    eprintln!(
        "INSERTION_IDENTITY label=\"{}\" loops=2 identical=true value={}",
        value.label, got.value
    );
    // a non-local slot is refused by the vacuum entry too
    let mut bad = value.clone();
    bad.insertions[0].factor = atom(&q(1, &mink(2)));
    let error = evaluate_factorized_vacuum(&vakint, &settings, &bad, &Atom::num(1)).unwrap_err();
    assert!(
        matches!(
            error,
            FactorizedError::IndexNotLocal { .. } | FactorizedError::Insertion { .. }
        ),
        "{error}"
    );
}

/// The beta's real record with two contracted children (L3 graph 185,
/// remainder 3, sector 0; written by their d5-caller dump; its closure section
/// removed and nf_power set to the original diagram's, so that it is a
/// direct-kind value): their insertion lines (vertices listed 2 then 0, a
/// squared indexed momentum, ε poles in `gammalooprs::ε`) pass every insertion
/// rule, and the quark cycle inside insertion 2 and the ghost cycle inside
/// insertion 0 are what the weights check counts.
#[test]
fn the_beta_threads_real_insertion_lines_pass_the_rules() {
    const REAL: &str = include_str!("fixtures/l3-185-insertions.txt");
    let value = read(REAL);
    assert_eq!(
        value
            .insertions
            .iter()
            .map(|i| i.vertex)
            .collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert_eq!(value.insertions[0].cycles, (0, 1));
    assert_eq!(value.insertions[1].cycles, (1, 0));
    assert!(value.closed_loops().unwrap().is_empty());
    let text = value.to_text().unwrap();
    assert_eq!(FactorizedNumerator::from_text(&text).unwrap(), value);
    let mut wrong = value.clone();
    wrong.prefactor.nf_power = 0;
    assert!(matches!(
        wrong.validate(),
        Err(FactorizedError::NfPower { .. })
    ));
    let mut wrong = value.clone();
    wrong.insertions[1].hedges = vec![9];
    assert!(is_insertion(&wrong.validate().unwrap_err()));
}

/// Ruling OQ-23: the imaginary unit and the UFO couplings may sit INSIDE edge,
/// vertex and insertion factors (prefactor scalars left at 1): accepted, and
/// the value is that of the contracted route.
#[test]
fn couplings_and_the_imaginary_unit_inside_factors_are_accepted() {
    let (vakint, _) = engine();
    let settings = settings(2);
    let mut value = gluon_surrogate(1);
    value.prefactor.scalars = Atom::num(1);
    let i = Atom::i();
    for edge in value.edges.iter_mut().filter(|e| e.id == 2) {
        edge.factor = &edge.factor * &i * atom("UFO::gs_minus");
    }
    for vertex in value.vertices.iter_mut().filter(|v| v.id == 0) {
        vertex.factor = &vertex.factor * atom("-1/2*UFO::gs_i");
    }
    ins(&mut value).coefficient = &ins(&mut value).coefficient * &i * atom("UFO::gs_minus^2");
    ins(&mut value).factor = &ins(&mut value).factor * atom("UFO::gs_i") * &i;
    value.label = format!("{} (i and couplings inside the factors)", value.label);
    identical(&vakint, &settings, &value);
}

/// `spenso::dot(Q(e, mink(D)), Q(f, mink(D)))` (how the beta's surrogate
/// factors write p²) is a scalar product of two boundary momenta: accepted,
/// degree 2, the same value as the indexed square.
#[test]
fn a_dot_of_boundary_momenta_inside_an_insertion_is_accepted() {
    let (vakint, _) = engine();
    let settings = settings(2);
    let mut value = gluon_surrogate(1);
    let compact = format!("{MINK}({DIM})");
    let dot = format!("spenso::dot({},{})", q(3, &compact), q(OUT, &compact));
    ins(&mut value).factor = atom(&format!(
        "({DIM}-2)/({DIM}-1)*{}*{dot}-{}*{}",
        metric(&mink(H_IN), &mink(H_OUT)),
        q(3, &mink(H_IN)),
        q(OUT, &mink(H_OUT))
    ));
    value.label = format!("{} (p² as a dot)", value.label);
    identical(&vakint, &settings, &value);
    // degree: dot counts two momenta
    let mut high = value.clone();
    ins(&mut high).factor = &ins(&mut high).factor * atom(&dot);
    assert!(is_insertion(&refused(&high)));
    // a dot with a non-boundary momentum
    let mut bad = value.clone();
    ins(&mut bad).factor = atom(&format!(
        "{}*spenso::dot({},{})",
        metric(&mink(H_IN), &mink(H_OUT)),
        q(2, &compact),
        q(3, &compact)
    ));
    assert!(matches!(
        refused(&bad),
        FactorizedError::UnsupportedTensor { .. }
    ));
}

/// OQ-28: the reference / closing mass inside an insertion is refused (a
/// massive-reference-mass child stays on the caller's path), whether named by
/// R*'s symbols or by an edge's mass.
#[test]
fn a_reference_mass_inside_an_insertion_is_refused() {
    for name in [
        "gammalooprs::rstar_reference_mass_squared",
        "gammalooprs::rstar_mass_squared",
    ] {
        let mut value = gluon_surrogate(1);
        ins(&mut value).coefficient = &ins(&mut value).coefficient * atom(name);
        let error = refused(&value);
        assert!(is_insertion(&error), "{error}");
        assert!(error.to_string().contains("mass symbol"), "{error}");
        let mut value = gluon_surrogate(1);
        ins(&mut value).factor = &ins(&mut value).factor * atom(&format!("1+{name}"));
        assert!(is_insertion(&refused(&value)));
    }
    // an edge's own mass symbol
    let mut value = gluon_surrogate(1);
    value
        .edges
        .iter_mut()
        .find(|e| e.id == 2)
        .unwrap()
        .mass_squared = atom("helper_test_m2");
    ins(&mut value).coefficient = atom("helper_test_m2");
    assert!(is_insertion(&refused(&value)));
}

/// The pole order the closure entry raises its depth by.
#[test]
fn the_pole_order_of_the_insertions() {
    let epsilon = atom(EPS);
    let AtomView::Var(epsilon) = epsilon.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    assert_eq!(read(GLUON).insertion_pole_order(epsilon).unwrap(), 0);
    assert_eq!(gluon_surrogate(1).insertion_pole_order(epsilon).unwrap(), 1);
    assert_eq!(quark_surrogate(1).insertion_pole_order(epsilon).unwrap(), 1);
    // SUM over the insertions (the product multiplies the record): the
    // child's eps^-2 and the surrogate's eps^-1 give 3, not the maximum 2
    assert_eq!(
        gluon_surrogate_and_child(1)
            .insertion_pole_order(epsilon)
            .unwrap(),
        3
    );
    // two eps^-1 insertions: raise 2
    let mut two = gluon_surrogate_and_child(1);
    two.insertions[0].coefficient = atom(&format!("{EPS}^-1"));
    assert_eq!(two.insertion_pole_order(epsilon).unwrap(), 2);
    // a pole only inside the factor (coefficient 1), and a regular insertion
    let mut value = gluon_surrogate(1);
    let factor = &ins(&mut value).factor * atom(&format!("{EPS}^-3+{DIM}"));
    ins(&mut value).factor = factor;
    ins(&mut value).coefficient = Atom::num(1);
    assert_eq!(value.insertion_pole_order(epsilon).unwrap(), 3);
    ins(&mut value).factor = atom(&format!(
        "{}*(1/(1-{EPS}))",
        metric(&mink(H_IN), &mink(H_OUT))
    ));
    assert_eq!(value.insertion_pole_order(epsilon).unwrap(), 0);
}

/// A compact slash `gamma(a, b, Q(x, mink(D)))` in an insertion factor (the
/// producers' K3 Dirac reassembly inside a quotient's insertion) is opened on
/// read, exactly: the record equals the one written with the opened slash.
#[test]
fn a_compact_slash_in_an_insertion_is_opened_on_read() {
    let (vakint, _) = engine();
    let settings = settings(2);
    let opened = quark_surrogate(1);
    let text = opened.to_text().unwrap();
    let factor = one_line(&opened.insertions[0].factor);
    let compact = atom(&gamma(
        &bis(H_OUT),
        &bis(H_IN),
        &q(3, &format!("{MINK}({DIM})")),
    ));
    let compact_text = text.replace(&factor, &compact.to_canonical_string());
    assert_ne!(compact_text, text);
    let read = FactorizedNumerator::from_text(&compact_text).unwrap();
    assert_ne!(read.insertions[0].factor, opened.insertions[0].factor);
    let a = evaluate_factorized(&vakint, &settings, &read, &Atom::num(1)).unwrap();
    let b = evaluate_factorized(&vakint, &settings, &opened, &Atom::num(1)).unwrap();
    assert!(!b.value.is_zero());
    assert!(
        (&a.value - &b.value).expand().is_zero(),
        "{}\n{}",
        a.value,
        b.value
    );
    // two compact slashes in one factor get distinct dummies
    let other = atom(&gamma(
        &bis(H_IN),
        &bis(H_OUT),
        &q(OUT, &format!("{MINK}({DIM})")),
    ));
    let twice = super::expand_compact_slash(&(&compact * &other), C0 as i64);
    assert!(
        !twice
            .to_canonical_string()
            .contains(&format!("{MINK}({DIM})")),
        "{twice}"
    );
}

fn one_line(atom: &Atom) -> String {
    atom.to_canonical_string()
}

/// Format v1.4 section 5: zeta(n) and pi are scalars in an insertion's
/// coefficient and factor; symbol dummies (`gammalooprs::dummy(n)`,
/// `insertion_slash_<n>`) are slots; bispinor hedge slots inside an
/// insertion factor (a Dirac decorated vertex) are local slots.
#[test]
fn zeta_pi_symbol_dummies_and_spinor_slots_inside_insertions() {
    let (vakint, _) = engine();
    let settings = settings(2);
    let mut value = gluon_surrogate(1);
    ins(&mut value).coefficient = &ins(&mut value).coefficient * atom("zeta(3)*𝜋^-2");
    let dummy = format!("{MINK}({DIM},gammalooprs::dummy(5))");
    let slash = format!("{MINK}({DIM},gammalooprs::insertion_slash_3)");
    ins(&mut value).factor = atom(&format!(
        "(zeta(5)+{DIM})*{}*{}*{}-{}*{}*{}*{}",
        metric(&mink(H_IN), &mink(H_OUT)),
        q(3, &dummy),
        q(3, &dummy),
        q(3, &mink(H_IN)),
        q(OUT, &mink(H_OUT)),
        metric(&slash, &slash),
        atom("1/gammalooprs::dim").to_canonical_string()
    ));
    value.label = format!("{} (zeta, pi, symbol dummies)", value.label);
    identical(&vakint, &settings, &value);
    // a Dirac insertion with spinor hedge slots: the quark surrogate's factor
    // written with a symbol dummy
    let mut quark = quark_surrogate(1);
    ins(&mut quark).factor = atom(&format!(
        "{}*{}",
        gamma(&bis(H_OUT), &bis(H_IN), &slash),
        q(3, &slash)
    ));
    quark.label = format!("{} (spinor slots, symbol dummy)", quark.label);
    identical(&vakint, &settings, &quark);
}
