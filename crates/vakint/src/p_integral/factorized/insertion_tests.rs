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

#[path = "../../../tests/support/wh_fixtures.rs"]
mod wh_fixtures;

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

/// The INDEPENDENT reference of the insertion identities: it shares no code
/// with the factorised entry beyond idenso's contraction and vakint's
/// evaluation. Written here: the product of every edge and vertex factor, of
/// each insertion's coefficient x factor and of the projector (plain
/// multiplication); idenso contracts it (`contract_with_idenso`, the base's
/// oracle); every scalar product of two edge momenta `g(Q(a, mink(D)), Q(b,
/// mink(D)))` / `dot(..)` becomes vakint's `dot` of the edges' momenta over
/// `k(i)`, `p(1)`; the topology; `D -> 4 - 2 eps`; `Vakint::evaluate`.
/// Valid for closed two-point values (a projector closing every port).
fn independent_reference(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
) -> Atom {
    use crate::symbols::S;
    let mut product = Atom::num(1);
    for edge in &value.edges {
        product *= &edge.factor;
    }
    for vertex in &value.vertices {
        product *= &vertex.factor;
    }
    for insertion in &value.insertions {
        product = &product * &insertion.coefficient * &insertion.factor;
    }
    let projector = value.projector.as_ref().expect("a closed value");
    product *= &projector.tensor;
    let contracted = contract_with_idenso(&product, false).unwrap();
    let loops = value.loops();
    let momentum = |id: i64| -> Atom {
        let edge = value
            .edges
            .iter()
            .find(|e| e.id as i64 == id)
            .unwrap_or_else(|| panic!("edge {id}"));
        let mut out = Atom::Zero;
        for (i, c) in edge.momentum.iter().enumerate() {
            if *c != 0 {
                let vector = if i < loops {
                    S.k.call(i as i64 + 1)
                } else {
                    S.p.call(i as i64 - loops as i64 + 1)
                };
                out += Atom::num(*c) * vector;
            }
        }
        out
    };
    let compact = |view: AtomView<'_>| -> Option<i64> {
        let AtomView::Fun(q) = view else { return None };
        if q.get_symbol().get_name() != "gammalooprs::Q" || q.get_nargs() != 2 {
            return None;
        }
        let mut args = q.iter();
        let e = i64::try_from(args.next()?).ok()?;
        matches!(args.next()?, AtomView::Fun(r) if r.get_nargs() == 1).then_some(e)
    };
    let mapped = contracted.replace_map(|part, _, out| {
        let AtomView::Fun(f) = part else { return };
        let symbol = f.get_symbol();
        let name = symbol.get_name();
        if (name == "spenso::g" || name == "spenso::dot") && f.get_nargs() == 2 {
            let args = f.iter().collect::<Vec<_>>();
            if let (Some(a), Some(b)) = (compact(args[0]), compact(args[1])) {
                **out = S.dot.call((&momentum(a), &momentum(b)));
            }
        }
    });
    let mut leftover = false;
    mapped.visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && matches!(
                f.get_symbol().get_name(),
                "gammalooprs::Q" | "spenso::g" | "spenso::gamma"
            )
        {
            leftover = true;
        }
        !leftover
    });
    assert!(!leftover, "{}: tensor notation left: {mapped}", value.label);
    let mut topology = Atom::num(1);
    for edge in value.edges.iter().filter(|e| !e.is_external()) {
        let (Endpoint::Vertex(s), Endpoint::Vertex(t)) = (edge.source, edge.target) else {
            unreachable!()
        };
        topology *= S.prop.call((
            edge.id as i64,
            S.edge.call((s as i64, t as i64)),
            &momentum(edge.id as i64),
            Atom::Zero,
            Atom::num(edge.power),
        ));
    }
    let head = S.p.call(1);
    let epsilon = Atom::var(crate::utils::vakint_macros::vk_symbol!(
        &settings.epsilon_symbol
    ));
    let input = (mapped
        * S.dot
            .call((&head, &head))
            .pow(Atom::num(projector.p2_power))
        * S.topo.call(topology))
    .replace(Atom::var(value.dimension))
    .with(Atom::num(4) - Atom::num(2) * epsilon);
    vakint.evaluate(settings, input.as_view()).unwrap()
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

/// The entry equals the INDEPENDENT reference ([`independent_reference`])
/// exactly; also the text round trip of the record. Prints one timing line
/// (median of 3 paired repetitions, single core) against the hand-contracted
/// route (the insertions multiplied into their vertex factors, evaluated
/// without the contract checks: a timing baseline, not a reference).
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
    let (got, _contracted) = (got.unwrap(), reference.unwrap());
    let reference = independent_reference(vakint, settings, value);
    assert!(!reference.is_zero(), "{}: zero reference", value.label);
    assert!(
        (&got - &reference).expand().is_zero(),
        "{}:\ninsertion   {got}\nindependent {reference}",
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
fn gluon_surrogates_at_two_and_three_loops_equal_the_independent_reference() {
    let (vakint, _) = engine();
    let settings = settings(3);
    for bubbles in [1, 2] {
        identical(&vakint, &settings, &gluon_surrogate(bubbles));
    }
}

/// The 1-4 loop package (`VAKINT_TEST_P_LIBRARY_1_4`, outside the
/// repository) is REQUIRED: without it the test fails.
fn four_loop_engine() -> Vakint {
    let (vakint, four) = engine();
    assert!(four, "four loops need VAKINT_TEST_P_LIBRARY_1_4");
    vakint
}

#[test]
#[ignore = "needs the 1-4 loop package of VAKINT_TEST_P_LIBRARY_1_4"]
fn four_loop_surrogates_equal_the_independent_reference() {
    let vakint = four_loop_engine();
    let settings = settings(3);
    identical(&vakint, &settings, &gluon_surrogate(3));
    identical(&vakint, &settings, &quark_surrogate(3));
}

#[test]
fn quark_surrogates_with_a_dirac_line_through_the_chain_equal_the_independent_reference() {
    let (vakint, _) = engine();
    let settings = settings(3);
    for bubbles in [1, 2] {
        let value = quark_surrogate(bubbles);
        // the fermion loop runs through the chain (one node for the flow)
        let loops = value.closed_loops().unwrap();
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].edges, vec![2, 3, OUT]);
        identical(&vakint, &settings, &value);
    }
}

#[test]
fn a_contracted_child_with_poles_and_a_surrogate_in_one_record_equal_the_independent_reference() {
    let (vakint, _) = engine();
    let settings = settings(3);
    for bubbles in [1, 2] {
        identical(&vakint, &settings, &gluon_surrogate_and_child(bubbles));
    }
}

/// Two LOCAL contracted children (`ct:`) with poles on the L1 gluon bubble:
/// against the BASE's entry on the UNMARKED record (no insertion), each
/// coefficient x factor multiplied into its vertex's own factor by plain
/// multiplication here, through `evaluate_factorized` of a value without any
/// new-feature marker.
#[test]
fn local_contracted_children_equal_the_base_entry_on_the_unmarked_record() {
    let (vakint, _) = engine();
    let settings = settings(3);
    let unmarked = read(GLUON);
    let mut marked = unmarked.clone();
    let mut reference_record = unmarked.clone();
    for (vertex, coefficient) in [
        (0usize, atom(&format!("{EPS}^-1"))),
        (1, atom(&format!("3*{EPS}^-1+1/2"))),
    ] {
        let hedges = marked
            .local_hedges(Location::Vertex(vertex))
            .into_iter()
            .collect::<Vec<_>>();
        let v = marked.vertices.iter_mut().find(|v| v.id == vertex).unwrap();
        let factor = std::mem::replace(&mut v.factor, Atom::num(1));
        v.rule = format!("ct:child{vertex}");
        let r = reference_record
            .vertices
            .iter_mut()
            .find(|v| v.id == vertex)
            .unwrap();
        r.factor = &coefficient * &factor;
        marked.insertions.push(Insertion {
            vertex,
            hedges,
            cycles: (0, 0),
            coefficient,
            factor,
        });
    }
    marked.validate().unwrap();
    assert!(reference_record.insertions.is_empty() && reference_record.record.is_none());
    let got = evaluate_factorized(&vakint, &settings, &marked, &Atom::num(1)).unwrap();
    let reference =
        evaluate_factorized(&vakint, &settings, &reference_record, &Atom::num(1)).unwrap();
    assert!(!reference.value.is_zero());
    assert!(
        (&got.value - &reference.value).expand().is_zero(),
        "marked {}\nunmarked {}",
        got.value,
        reference.value
    );
    // and the independent reference agrees
    let independent = independent_reference(&vakint, &settings, &marked);
    assert!((&got.value - &independent).expand().is_zero());
}

fn refused(value: &FactorizedNumerator) -> FactorizedError {
    value
        .validate()
        .expect_err(&format!("{} must be refused", value.label))
}

fn ins(value: &mut FactorizedNumerator) -> &mut Insertion {
    &mut value.insertions[0]
}

fn c0(v: &mut FactorizedNumerator) -> &mut FactorizedVertex {
    v.vertices.iter_mut().find(|x| x.id == C0).unwrap()
}

/// `error` is an `Insertion` refusal at `vertex` whose detail starts with
/// `prefix` (the frozen REFUSALS list).
fn insertion_refusal(error: &FactorizedError, vertex: usize, prefix: &str, what: &str) {
    match error {
        FactorizedError::Insertion { vertex: v, detail } => {
            assert_eq!(*v, vertex, "{what}: {error}");
            assert!(
                detail.starts_with(prefix),
                "{what}: {detail:?} does not start with {prefix:?}"
            );
        }
        other => panic!("{what}: not an Insertion refusal: {other:?}"),
    }
}

/// Every structural insertion rule, each on an input that the PRISTINE
/// surrogate record (accepted, asserted first) turns into by ONE change: the
/// refusal's variant and detail prefix are asserted (the frozen REFUSALS
/// list), so the test fails if the rule is removed (the next rule that might
/// catch the input has another prefix or none).
#[test]
fn every_insertion_rule_has_its_refusal() {
    let good = gluon_surrogate(1);
    good.validate().expect("the pristine record is accepted");
    let with = |change: &dyn Fn(&mut FactorizedNumerator)| {
        let mut value = good.clone();
        change(&mut value);
        value
    };
    // the chain's second vertex holds H_OUT and the target half-edge 61 of the
    // first bubble line; 60 is that line's source half-edge, at C0
    for (what, value, vertex, prefix) in [
        (
            "missing vertex",
            with(&|v| ins(v).vertex = 999),
            999,
            "the vertex does not exist",
        ),
        (
            "vertex factor not 1",
            with(&|v| c0(v).factor = Atom::num(2)),
            C0,
            "the vertex factor must be 1",
        ),
        (
            "empty hedges",
            with(&|v| ins(v).hedges.clear()),
            C0,
            "hedges must be non-empty, sorted and distinct",
        ),
        (
            "unsorted hedges",
            with(&|v| ins(v).hedges = vec![H_OUT, H_IN]),
            C0,
            "hedges must be non-empty, sorted and distinct",
        ),
        (
            "unknown hedge",
            with(&|v| ins(v).hedges = vec![H_IN, 999]),
            C0,
            "half-edge 999 is not in the graph",
        ),
        (
            "no hedge at the vertex",
            with(&|v| ins(v).hedges = vec![H_OUT, 61]),
            C0,
            "no boundary half-edge sits at the vertex",
        ),
        (
            "rule neither ct: nor ins:",
            with(&|v| c0(v).rule = "V_x".into()),
            C0,
            "the vertex rule V_x is neither",
        ),
        (
            "ct: with hedges other than the incident ones",
            with(&|v| c0(v).rule = "ct:x".into()),
            C0,
            "a contracted vertex (rule ct:x) has hedges",
        ),
        (
            "two insertions on one vertex",
            with(&|v| {
                let copy = v.insertions[0].clone();
                v.insertions.push(copy);
            }),
            C0,
            "insertions are not sorted by vertex, or two share a vertex",
        ),
        (
            "degree four",
            with(&|v| ins(v).factor = &ins(v).factor * atom(&format!("{}^2", q(3, &dummy(3, 8))))),
            C0,
            "the factor has degree 4 in the boundary momenta",
        ),
        (
            "a mass symbol",
            with(&|v| {
                ins(v).coefficient =
                    &ins(v).coefficient * atom("gammalooprs::rstar_reference_mass_squared")
            }),
            C0,
            "contains a mass symbol ",
        ),
    ] {
        insertion_refusal(&refused(&value), vertex, prefix, what);
        eprintln!("INSERTION_REFUSAL {what}: {}", refused(&value));
    }
    // five open slots: the contracted child at vertex 1 (three half-edges)
    // with two more spinor slots; its pristine record is accepted
    let mut five = gluon_surrogate_and_child(1);
    five.validate().expect("pristine child record");
    ins(&mut five).factor = &ins(&mut five).factor * atom(&metric(&bis(0), &bis(3)));
    insertion_refusal(
        &refused(&five),
        1,
        "the factor has 5 open slots, at most 4 are supported",
        "five slots",
    );
    // colour, flavour and mUVexp inside an insertion: the format's variants
    match refused(&with(&|v| {
        ins(v).coefficient = atom("spenso::cas(2,spenso::coad(spenso::dA))")
    })) {
        FactorizedError::ColourInFactor { location, tensor } => {
            assert_eq!(location, format!("insertion at vertex {C0} coefficient"));
            assert!(tensor.starts_with("the colour object "), "{tensor}");
        }
        other => panic!("colour: {other:?}"),
    }
    match refused(&with(&|v| ins(v).coefficient = atom("beta::nf"))) {
        FactorizedError::Graph { detail } => assert!(
            detail.starts_with(&format!(
                "insertion at vertex {C0} coefficient contains the flavour symbol "
            )),
            "{detail}"
        ),
        other => panic!("flavour: {other:?}"),
    }
    match refused(&with(&|v| {
        ins(v).factor = &ins(v).factor * atom("gammalooprs::mUVexp")
    })) {
        FactorizedError::Graph { detail } => assert!(
            detail.starts_with(&format!(
                "insertion at vertex {C0} factor contains the auxiliary mass symbol "
            )),
            "{detail}"
        ),
        other => panic!("mUVexp: {other:?}"),
    }
    // a momentum in the coefficient, a non-boundary momentum, a non-local slot
    match refused(&with(&|v| {
        ins(v).coefficient = atom(&format!("{}^2", q(3, &dummy(3, 8))))
    })) {
        FactorizedError::NotScalar { location, detail } => {
            assert_eq!(location, format!("insertion at vertex {C0} coefficient"));
            assert!(
                detail.starts_with("must be free of indices and momenta"),
                "{detail}"
            );
        }
        other => panic!("momentum in the coefficient: {other:?}"),
    }
    match refused(&with(&|v| {
        ins(v).factor = atom(&format!("{}*{}", q(2, &mink(H_IN)), q(OUT, &mink(H_OUT))))
    })) {
        FactorizedError::UnsupportedTensor { location, tensor } => {
            assert_eq!(location, format!("insertion at vertex {C0}"));
            assert!(
                tensor.ends_with("(not a boundary momentum of the insertion)"),
                "{tensor}"
            );
        }
        other => panic!("non-boundary momentum: {other:?}"),
    }
    match refused(&with(&|v| {
        ins(v).factor = atom(&format!("{}*{}", q(3, &mink(H_IN)), q(OUT, &mink(4))))
    })) {
        FactorizedError::IndexNotLocal { location, .. } => {
            assert_eq!(location, format!("insertion at vertex {C0}"))
        }
        other => panic!("non-local slot: {other:?}"),
    }
    // weights (standalone direct record): a closed cycle inside the insertion
    let error = refused(&with(&|v| ins(v).cycles = (1, 0)));
    assert!(
        matches!(error, FactorizedError::LoopSign { found: 1, loops: 1 }),
        "{error}"
    );
    let mut balanced = with(&|v| ins(v).cycles = (1, 0));
    balanced.prefactor.loop_sign = -1;
    let error = refused(&balanced);
    assert!(
        matches!(error, FactorizedError::NfPower { found: 0, loops: 1 }),
        "{error}"
    );
    balanced.prefactor.nf_power = 1;
    balanced
        .validate()
        .expect("the balanced weights are accepted");
    // a second epsilon at evaluation (the record itself validates)
    let wrong = with(&|v| ins(v).coefficient = atom("gammalooprs::ε^-1"));
    wrong.validate().unwrap();
    let error = super::evaluate::vakint_input(&settings(2), &wrong, &Atom::num(1)).unwrap_err();
    insertion_refusal(
        &error,
        C0,
        "the symbol gammalooprs::ε is not the evaluation's epsilon",
        "second eps",
    );
    let _ = super::evaluate::vakint_input(&settings(2), &good, &Atom::num(1))
        .expect("the pristine record's input");
    // a fractional power of eps: refused at validation (not a Laurent
    // polynomial); the pole order's own rule still answers for a value that
    // was not validated (the public insertion_pole_order)
    let half = with(&|v| ins(v).coefficient = atom(&format!("{EPS}^(1/2)")));
    insertion_refusal(
        &half.validate().unwrap_err(),
        C0,
        "the dependence on ε is not a Laurent polynomial",
        "eps^(1/2) at validation",
    );
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    insertion_refusal(
        &half.insertion_pole_order(epsilon).unwrap_err(),
        C0,
        "a fractional power of ε",
        "eps^(1/2)",
    );
    good.insertion_pole_order(epsilon)
        .expect("the pristine pole order");
}

/// The pole order's own Laurent expansion can fail: that is an INTERNAL
/// failure, `Evaluation` with the vertex and symbolica's reason. Both inputs
/// that make it fail are refused by validation first: an essential
/// singularity (exp(1/ε): not a Laurent polynomial) and an integer power
/// beyond i64 (ε^(-10^30), which symbolica does not expand: an ε power out of
/// range). Only the public `insertion_pole_order` on a value that was not
/// validated reaches it. The same record with a pole has order 1.
#[test]
fn a_failing_laurent_series_of_an_insertion_is_an_evaluation_failure() {
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    let mut essential = gluon_surrogate(1);
    ins(&mut essential).coefficient = atom(&format!("exp(1/{EPS})"));
    insertion_refusal(
        &essential.validate().unwrap_err(),
        C0,
        "the dependence on ε is not a Laurent polynomial",
        "exp(1/eps) at validation",
    );
    let error = essential.insertion_pole_order(epsilon).unwrap_err();
    let prefix = format!("insertion at vertex {C0}: the Laurent series in ε failed: ");
    match &error {
        FactorizedError::Evaluation { detail } => assert!(
            detail.starts_with(&prefix) && detail.len() > prefix.len(),
            "{detail}"
        ),
        other => panic!("expected Evaluation, got {other:?}"),
    }
    let mut huge = gluon_surrogate(1);
    ins(&mut huge).coefficient = atom(&format!("{EPS}^(-10^30)"));
    insertion_refusal(
        &huge.validate().unwrap_err(),
        C0,
        "an ε power out of range",
        "eps^(-10^30) at validation",
    );
    match &huge.insertion_pole_order(epsilon).unwrap_err() {
        FactorizedError::Evaluation { detail } => {
            assert!(detail.starts_with(&prefix), "{detail}")
        }
        other => panic!("expected Evaluation, got {other:?}"),
    }
    let mut pole = gluon_surrogate(1);
    ins(&mut pole).coefficient = atom(&format!("{EPS}^-1"));
    assert_eq!(pole.insertion_pole_order(epsilon).unwrap(), 1);
}

/// The degree count is not the pole order: a negative power of a sum whose
/// terms cancel at ε = 0, ((1-ε)^-1 - 1 - ε)^-n (its base is ε²/(1-ε)),
/// counts n but has pole order 2n. With n = 5 (count 5, pole order 10) it is
/// accepted with its exact order; times ε^-55 (count 60, within the
/// per-insertion rule with the factor's 2) its exact pole order is 65, which
/// the depth raise refuses (`an ε power out of range`, the total above 64).
/// Checked at the depth raise itself (`insertion_depth_raise`, what the
/// closure entry calls before any reduction). Kept small: symbolica's series
/// of such an unsimplified power takes ~1 s at n = 10 and ~50 s at n = 20.
#[test]
fn a_cancelling_sum_under_a_negative_power_keeps_its_exact_pole_order() {
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    let good = gluon_surrogate(1);
    let with = |coefficient: String| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(&coefficient);
        value
    };
    let cancelling = format!("((1-{EPS})^-1-1-{EPS})^-5");
    let twin = with(cancelling.clone());
    twin.validate().expect("count 5: accepted");
    assert_eq!(twin.insertion_depth_raise(epsilon).unwrap(), 10);
    let witness = with(format!("{EPS}^-55*{cancelling}"));
    witness
        .validate()
        .expect("count 60: accepted by the per-insertion rule");
    assert_eq!(witness.insertion_pole_order(epsilon).unwrap(), 65);
    insertion_refusal(
        &witness.insertion_depth_raise(epsilon).unwrap_err(),
        C0,
        "an ε power out of range: the insertions' pole orders add up to 65, above 64",
        "the total pole order",
    );
}

/// An insertion that is UNDEFINED after D -> 4-2ε (a negative power of a base
/// that vanishes identically then) is refused (Insertion, "a denominator
/// vanishes identically after D -> 4-2ε"), wherever the zero is: (D-4+2ε)^-1;
/// (x-x+D-4+2ε)^-1 (the parser cancels x-x: the same); a base that is zero
/// only after expansion, ((1+ε)^2-1-2ε-ε^2)^-1; nested,
/// ((D-4+2ε)^-1+1)^-1; in the factor. A zero NUMERATOR is defined (the
/// insertion contributes zero): D-4+2ε and (D-4+2ε)·ε^-1 are accepted with
/// pole order 0; an ordinary pole (a base vanishing only at ε = 0), 1/(D-4),
/// is accepted with pole order 1. A 0/0 cannot be written: the parser (and
/// Atom arithmetic alike) cancels (D-4+2ε)/(D-4+2ε) to 1.
#[test]
fn an_insertion_undefined_after_the_substitution_is_refused() {
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    let good = gluon_surrogate(1);
    let with_coefficient = |c: &str| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(c);
        value
    };
    let zero = format!("({DIM}-4+2*{EPS})");
    let prefix = "a denominator vanishes identically";
    for (what, value) in [
        ("1/(D-4+2eps)", with_coefficient(&format!("{zero}^-1"))),
        (
            "(x-x+D-4+2eps)^-1",
            with_coefficient(&format!("(witness::x-witness::x+{DIM}-4+2*{EPS})^-1")),
        ),
        (
            "zero only after expansion",
            with_coefficient(&format!("((1+{EPS})^2-1-2*{EPS}-{EPS}^2)^-1")),
        ),
        ("nested", with_coefficient(&format!("({zero}^-1+1)^-1"))),
        ("in the factor", {
            let mut value = good.clone();
            let factor = ins(&mut value).factor.clone();
            ins(&mut value).factor = factor * atom(&format!("{zero}^-1"));
            value
        }),
    ] {
        insertion_refusal(&value.validate().unwrap_err(), C0, prefix, what);
    }
    for (what, coefficient, order) in [
        ("a zero numerator", zero.clone(), 0),
        ("zero times a pole", format!("{zero}*{EPS}^-1"), 0),
        ("an ordinary pole", format!("1/({DIM}-4)"), 1),
    ] {
        let value = with_coefficient(&coefficient);
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what} is defined: {e}"));
        assert_eq!(
            value.insertion_pole_order(epsilon).unwrap(),
            order,
            "{what}"
        );
    }
    assert_eq!(
        atom(&format!("{zero}/{zero}")),
        Atom::num(1),
        "the parser cancels 0/0"
    );
    assert_eq!(
        atom(&zero) * atom(&format!("{zero}^-1")),
        Atom::num(1),
        "and so does Atom arithmetic"
    );
}

/// The zero test of a denominator does not depend on what the base contains:
/// identically zero bases are refused (P(t) = (t-3)^2 - t^2 + 6t - 9 for
/// t = D, a caller scalar x, a number 2 - 2 spelled as x - x + 2 - 2 inside a
/// sum; in the factor too), their twins with base 1 ((t-3)^2 - t^2 + 6t - 8)
/// are accepted, as are 1/(D-4) (a pole), 1/x and 1/(x-3) (zero only at one
/// value of x). A large base under a negative power, ((x+y+z+w)^40)^-1, is
/// decided by the cheap exact pre-filter (no expansion).
#[test]
fn a_vanishing_denominator_is_refused_whatever_it_depends_on() {
    let good = gluon_surrogate(1);
    let with_coefficient = |c: &str| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(c);
        value
    };
    let p = |t: &str, k: i64| format!("(({t}-3)^2-{t}^2+6*{t}-{k})");
    let prefix = "a denominator vanishes identically";
    for (what, value) in [
        ("1/P(D)", with_coefficient(&format!("1/{}", p(DIM, 9)))),
        (
            "1/P(x)",
            with_coefficient(&format!("1/{}", p("witness::x", 9))),
        ),
        (
            "1/(x - x + P(y))",
            with_coefficient(&format!("1/(witness::x-witness::x+{})", p("witness::y", 9))),
        ),
        ("1/P(D) in the factor", {
            let mut value = good.clone();
            let factor = ins(&mut value).factor.clone();
            ins(&mut value).factor = factor * atom(&format!("1/{}", p(DIM, 9)));
            value
        }),
    ] {
        insertion_refusal(&value.validate().unwrap_err(), C0, prefix, what);
    }
    for (what, coefficient) in [
        ("base 1 in D", format!("1/{}", p(DIM, 8))),
        ("base 1 in x", format!("1/{}", p("witness::x", 8))),
        ("a pole", format!("1/({DIM}-4)")),
        ("1/x", "1/witness::x".to_string()),
        ("1/(x-3)", "1/(witness::x-3)".to_string()),
    ] {
        with_coefficient(&coefficient)
            .validate()
            .unwrap_or_else(|e| panic!("{what} is defined: {e}"));
    }
    let started = std::time::Instant::now();
    with_coefficient("((witness::x+witness::y+witness::z+witness::w)^40)^-1")
        .validate()
        .expect("a large nonzero base");
    eprintln!("LARGE BASE validation: {:?}", started.elapsed());
}

/// Exponents and function arguments, at validation: a vanishing base P
/// (identically zero) under the integer 2 is accepted (the coefficient is 0:
/// the insertion contributes zero, pole order 0); the zero test refuses P^y
/// (a symbol), x^(1/Q) (Q identically zero inside the exponent: the exponent
/// is walked) and f(1/P) (a caller function's argument is walked); the forms
/// it passes, x^(-sqrt(2)) and P^(1/2), are refused by the integer-exponent
/// rule. The zero test alone
/// ([`vanishing_denominator`], the rule's second line): a vanishing base
/// under -sqrt(2), -π, a symbol, or inside an exponent (x^(1/Q)) is found;
/// x^(-sqrt(2)), P^(1/2), P^2 and (P+1)^(-sqrt(2)) are not.
#[test]
fn exponents_and_function_arguments_are_walked_for_vanishing_bases() {
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    let good = gluon_surrogate(1);
    let with_coefficient = |c: &str| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(c);
        value
    };
    let p = |t: &str| format!("(({t}-3)^2-{t}^2+6*{t}-9)");
    let (x, q) = ("witness::x", "witness::q");
    let value = with_coefficient(&format!("{}^2", p(x)));
    value.validate().unwrap_or_else(|e| panic!("P^2: {e}"));
    assert_eq!(value.insertion_pole_order(epsilon).unwrap(), 0, "P^2");
    for (what, coefficient, prefix) in [
        (
            "x^(-sqrt(2))",
            format!("{x}^(-2^(1/2))"),
            "a non-integer exponent",
        ),
        (
            "P^(1/2)",
            format!("{}^(1/2)", p(x)),
            "a non-integer exponent",
        ),
        (
            "P^y",
            format!("{}^witness::y", p(x)),
            "a denominator vanishes identically",
        ),
        (
            "x^(1/Q)",
            format!("{x}^(1/{})", p(q)),
            "a denominator vanishes identically",
        ),
        (
            "f(1/P)",
            format!("witness::f(1/{})", p(x)),
            "a denominator vanishes identically",
        ),
    ] {
        insertion_refusal(
            &with_coefficient(&coefficient).validate().unwrap_err(),
            C0,
            prefix,
            what,
        );
    }
    let dimension = good.dimension;
    for (what, expression, found) in [
        ("P^(-sqrt(2))", format!("{}^(-2^(1/2))", p(x)), true),
        ("P^(-π)", format!("{}^(-𝜋)", p(x)), true),
        ("P^y", format!("{}^witness::y", p(x)), true),
        ("x^(1/Q)", format!("{x}^(1/{})", p(q)), true),
        ("x^(-sqrt(2))", format!("{x}^(-2^(1/2))"), false),
        ("P^(1/2)", format!("{}^(1/2)", p(x)), false),
        ("P^2", format!("{}^2", p(x)), false),
        (
            "(P+1)^(-sqrt(2))",
            format!("({}+1)^(-2^(1/2))", p(x)),
            false,
        ),
    ] {
        let expression = atom(&expression);
        assert_eq!(
            vanishing_denominator(expression.as_view(), dimension, epsilon).is_some(),
            found,
            "{what}"
        );
    }
}

/// Validation recognises ε by its name; an epsilon setting of another name
/// (`user::ep`) gets the pole order's rules where the settings are known, the
/// entries' input check (here `vakint_input`, the direct and ports path): the
/// same variant and prefixes; the twin with a pole is accepted.
#[test]
fn an_epsilon_setting_of_another_name_gets_the_pole_order_rules_at_the_entry() {
    let mut user = settings(2);
    user.epsilon_symbol = "user::ep".into();
    let good = gluon_surrogate(1);
    let with = |coefficient: &str| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(coefficient);
        value
    };
    // log(user::ep): a function head that is not allowed, refused by name at
    // validation already (whatever the settings' epsilon)
    insertion_refusal(
        &with("log(user::ep)").validate().unwrap_err(),
        C0,
        &format!(
            "a function head that is not allowed in insertion {C0} coefficient: symbolica::log"
        ),
        "log(user::ep)",
    );
    for (what, coefficient, prefix) in [(
        "user::ep^-40 (1+user::ep)^-30 (degree 70, every exponent within 64)",
        "user::ep^-40*(1+user::ep)^-30",
        "an ε power out of range",
    )] {
        let value = with(coefficient);
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what}: validation does not know the settings: {e}"));
        let error = super::evaluate::vakint_input(&user, &value, &Atom::num(1)).unwrap_err();
        insertion_refusal(&error, C0, prefix, what);
    }
    let _accepted = super::evaluate::vakint_input(&user, &with("user::ep^-1"), &Atom::num(1))
        .expect("user::ep^-1 is a Laurent coefficient");
}

/// The pole order is bounded: a bound on |the powers of ε and of the
/// dimension| of coefficient x factor (degrees: sums the largest, products
/// added, integer powers times |n|) may reach at most
/// INSERTION_MAX_EPSILON_POWER (64), the coefficient's and the factor's bounds
/// added; beyond, and for an exponent beyond i64, the record is refused at
/// validation, `Insertion` "an ε power out of range". At the cap the record is
/// accepted with its pole order (at most the bound).
#[test]
fn an_epsilon_power_beyond_the_cap_is_refused_at_validation() {
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    assert_eq!(super::INSERTION_MAX_EPSILON_POWER, 64);
    let good = gluon_surrogate(1);
    let with = |coefficient: &str, extra: &str| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(coefficient);
        let factor = ins(&mut value).factor.clone();
        ins(&mut value).factor = factor * atom(extra);
        value
    };
    // the surrogate's own factor, (D-2)/(D-1) x tensors, has the bound 2
    // (degree 1 in D twice), so 62 is the coefficient's power at the cap
    let dimension = good.dimension;
    let depends = |s: Symbol| super::is_epsilon_symbol(s) || s == dimension;
    assert_eq!(
        super::epsilon_power_bound(ins(&mut good.clone()).factor.as_view(), &depends),
        Some(2)
    );
    let prefix = "an ε power out of range";
    for (what, value) in [
        (
            "eps^-63 (one above the cap)",
            with(&format!("{EPS}^-63"), "1"),
        ),
        (
            "eps^63 (one above the cap)",
            with(&format!("{EPS}^63"), "1"),
        ),
        (
            "eps^-31 x (D-4)^-32 (the bounds add, one above)",
            with(&format!("{EPS}^-31"), &format!("({DIM}-4)^-32")),
        ),
        (
            "eps^(-10^30) (beyond i64)",
            with(&format!("{EPS}^(-10^30)"), "1"),
        ),
        (
            "(1+eps)^(2^62) (a power of a sum)",
            with(&format!("(1+{EPS})^4611686018427387904"), "1"),
        ),
    ] {
        insertion_refusal(&value.validate().unwrap_err(), C0, prefix, what);
    }
    for (what, value, order) in [
        ("eps^-62 (at the cap)", with(&format!("{EPS}^-62"), "1"), 62),
        ("eps^62 (at the cap)", with(&format!("{EPS}^62"), "1"), 0),
        (
            "eps^-31 x (D-4)^-31 (at the cap)",
            with(&format!("{EPS}^-31"), &format!("({DIM}-4)^-31")),
            62,
        ),
    ] {
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what} is at the cap: {e}"));
        assert_eq!(
            value.insertion_pole_order(epsilon).unwrap(),
            order,
            "{what}"
        );
    }
}

/// An insertion's coefficient and factor are Laurent expressions in ε: ε and
/// the dimension (which the pole order sets to 4 - 2ε) may occur only through
/// sums, products and integer powers. Anything
/// else has no pole order and is refused at validation, `Insertion` with the
/// prefix `the dependence on ε is not a Laurent polynomial`: log ε, exp(1/ε),
/// √ε, ε^x, ε inside the factor's function, the dimension inside a log (in
/// the coefficient and in the factor) or under a square root (factor). The
/// dimension as a representation slot's dimension argument is no dependence.
/// The defect-free twins are accepted with their pole order: negative integer
/// powers of sums (1/(1+ε), 1/((D-1)(D-2)), 1/(D-4)^2), a caller's scalar
/// function of no ε, a fractional power of another scalar.
#[test]
fn a_non_laurent_dependence_on_epsilon_is_refused_at_validation() {
    let eps_atom = atom(EPS);
    let AtomView::Var(epsilon) = eps_atom.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    let good = gluon_surrogate(1);
    let with_coefficient = |c: &str| {
        let mut value = good.clone();
        ins(&mut value).coefficient = atom(c);
        value
    };
    // the factor times `extra`, coefficient 1 (so the pole order is the
    // factor's own)
    let with_factor = |extra: &str| {
        let mut value = good.clone();
        let factor = ins(&mut value).factor.clone();
        ins(&mut value).factor = factor * atom(extra);
        ins(&mut value).coefficient = Atom::num(1);
        value
    };
    let prefix = "the dependence on ε is not a Laurent polynomial";
    for (what, value) in [
        ("log(eps)", with_coefficient(&format!("log({EPS})"))),
        ("exp(1/eps)", with_coefficient(&format!("exp(1/{EPS})"))),
        ("sqrt(eps)", with_coefficient(&format!("sqrt({EPS})"))),
        ("eps^x", with_coefficient(&format!("{EPS}^witness::x"))),
        ("eps in an exponent", with_coefficient(&format!("2^{EPS}"))),
        (
            "log of the dimension in the coefficient",
            with_coefficient(&format!("log({DIM}-4)")),
        ),
        (
            "eps in a function of the factor",
            with_factor(&format!("exp({EPS})")),
        ),
        (
            "log of the dimension in the factor",
            with_factor(&format!("log({DIM}-4)")),
        ),
        (
            "a square root of the dimension in the factor",
            with_factor(&format!("({DIM}-2)^(1/2)")),
        ),
    ] {
        let error = value.validate().unwrap_err();
        insertion_refusal(&error, C0, prefix, what);
    }
    // a fractional power of a scalar that is neither eps nor the dimension is
    // not a dependence on eps; the R* reference mass is refused in an
    // insertion by the mass rule (whatever its power), not by this one
    insertion_refusal(
        &with_coefficient(&format!(
            "gammalooprs::rstar_reference_mass_squared^(1/2)*{EPS}^-1"
        ))
        .validate()
        .unwrap_err(),
        C0,
        "contains a mass symbol ",
        "the reference mass^(1/2)",
    );
    for (what, coefficient, order) in [
        ("eps^-1", format!("{EPS}^-1"), 1),
        ("3 eps^-1 + 1/2", format!("3*{EPS}^-1+1/2"), 1),
        ("eps^-2 pi^2", format!("{EPS}^-2*𝜋^2"), 2),
        ("zeta(3) eps", format!("zeta(3)*{EPS}"), 0),
        ("a ratio in eps", format!("(1+{EPS})^-1*{EPS}^-1"), 1),
        (
            "a ratio in the dimension",
            format!("({DIM}-4)/({DIM}-3)"),
            0,
        ),
        (
            "a negative power of a sum in eps",
            format!("1/(1+{EPS})"),
            0,
        ),
        (
            "a caller's plain scalar of no eps (a caller FUNCTION is refused by the head allow-list, below)",
            format!("gammalooprs::master7*{EPS}^-2"),
            2,
        ),
        (
            "an integer power of another scalar (a fractional one is refused: exponents are integer literals)",
            format!("witness::m2^3*{EPS}^-1"),
            1,
        ),
    ] {
        let value = with_coefficient(&coefficient);
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what} is a Laurent coefficient: {e}"));
        assert_eq!(
            value.insertion_pole_order(epsilon).unwrap(),
            order,
            "{what}"
        );
    }
    insertion_refusal(
        &with_coefficient(&format!(
            "gammalooprs::UNKNOWN(gammalooprs::master7,3)*{EPS}^-2"
        ))
        .validate()
        .unwrap_err(),
        C0,
        &format!(
            "a function head that is not allowed in insertion {C0} coefficient: gammalooprs::UNKNOWN"
        ),
        "a caller's function of no eps",
    );
    for (what, extra, order) in [
        ("the dimension only in slots", "1".to_string(), 0),
        ("a polynomial in the dimension", format!("{DIM}-2"), 0),
        (
            "a projector normalisation",
            format!("1/(({DIM}-1)*({DIM}-2))"),
            0,
        ),
        ("a pole in the dimension", format!("1/({DIM}-4)"), 1),
        (
            "a double pole in the dimension",
            format!("1/({DIM}-4)^2"),
            2,
        ),
    ] {
        let value = with_factor(&extra);
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what} is a Laurent factor: {e}"));
        assert_eq!(
            value.insertion_pole_order(epsilon).unwrap(),
            order,
            "{what}"
        );
    }
}

/// Integers a marked direct value carries are checked before use: an edge or
/// vertex id above i64::MAX (written as Q(id, ..) and edge(source, target):
/// it wrapped to a negative id) is refused (Graph, "id .."); a momentum
/// coefficient above 2^31 in absolute value (multiplied by ±1 and summed per
/// vertex: it could wrap) is refused (MomentumBasis); a tensor power
/// above 2 of a slot-carrying factor (Q(3, ..)^(10^12): the rules copied its
/// slots n times before checking them) is refused at once as the
/// over-contraction it is; the settings' number of terms outside 1..=128 is
/// refused by the entries' input check (Graph, "the depth .. is out of
/// range"), 2 accepted.
#[test]
fn a_marked_value_refuses_integers_out_of_range() {
    let good = gluon_surrogate(1);
    good.validate().unwrap();
    let mut edge = good.clone();
    edge.edges[0].id = usize::MAX;
    let mut vertex = good.clone();
    vertex.vertices[0].id = usize::MAX;
    for (what, value) in [("edge id", edge), ("vertex id", vertex)] {
        match value.validate().unwrap_err() {
            FactorizedError::Graph { detail } => assert!(
                detail.starts_with(&format!("id {}:", usize::MAX)),
                "{what}: {detail}"
            ),
            other => panic!("{what}: expected the Graph id refusal, got {other:?}"),
        }
    }
    let mut coefficient = good.clone();
    coefficient.edges[0].momentum[0] = (1 << 31) + 1;
    match coefficient.validate().unwrap_err() {
        FactorizedError::MomentumBasis { detail } => assert!(
            detail.contains("momentum coefficient 2147483649 out of range"),
            "{detail}"
        ),
        other => panic!("expected the momentum-coefficient refusal, got {other:?}"),
    }
    let mut power = good.clone();
    let tensor = atom(&format!("{}^1000000000000", q(3, &mink(40))));
    let factor = &power.vertices[0].factor * &tensor;
    power.vertices[0].factor = factor;
    assert!(
        matches!(
            power.validate().unwrap_err(),
            FactorizedError::IndexOverContracted { .. }
        ),
        "a huge tensor power"
    );
    for terms in [-1i64, 0, 129] {
        let mut settings = settings(2);
        settings.number_of_terms_in_epsilon_expansion = terms;
        match super::evaluate::vakint_input(&settings, &good, &Atom::num(1)).unwrap_err() {
            FactorizedError::Graph { detail } => assert!(
                detail.starts_with(&format!("the depth {terms} is out of range")),
                "{terms}: {detail}"
            ),
            other => panic!("{terms}: expected the Graph depth refusal, got {other:?}"),
        }
    }
    let _accepted =
        super::evaluate::vakint_input(&settings(2), &good, &Atom::num(1)).expect("2 terms");
}

/// The surrogate coefficient vertex (`ins:`) carries exactly the two inner
/// boundary half-edges (format section 8): a third half-edge of the
/// vertex (the first bubble line's source, 60) is refused by THAT rule; the
/// pristine record (two hedges) is accepted.
#[test]
fn a_surrogate_vertex_needs_exactly_two_hedges() {
    let good = gluon_surrogate(1);
    good.validate().expect("the pristine record is accepted");
    let mut three = good.clone();
    ins(&mut three).hedges = vec![H_IN, H_OUT, 60];
    insertion_refusal(
        &refused(&three),
        C0,
        "a surrogate coefficient vertex (rule ins:se) has 3 hedges",
        "ins: three hedges",
    );
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
            mass_squared: atom("helper_test::muvsq"),
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
    insertion_refusal(
        &wrong.validate().unwrap_err(),
        2,
        "a contracted vertex (rule ct:edges-6-7) has hedges [9]",
        "real record, a hedge dropped",
    );
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
    insertion_refusal(
        &refused(&high),
        C0,
        "the factor has degree 4 in the boundary momenta",
        "dot x dot",
    );
    // a dot with a non-boundary momentum
    let mut bad = value.clone();
    ins(&mut bad).factor = atom(&format!(
        "{}*spenso::dot({},{})",
        metric(&mink(H_IN), &mink(H_OUT)),
        q(2, &compact),
        q(3, &compact)
    ));
    match refused(&bad) {
        FactorizedError::UnsupportedTensor { location, tensor } => {
            assert_eq!(location, format!("insertion at vertex {C0}"));
            assert!(
                tensor.ends_with("(not a boundary momentum of the insertion)"),
                "{tensor}"
            );
        }
        other => panic!("a dot with a non-boundary momentum: {other:?}"),
    }
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
        insertion_refusal(
            &refused(&value),
            C0,
            &format!("contains a mass symbol {name}"),
            name,
        );
        let mut value = gluon_surrogate(1);
        ins(&mut value).factor = &ins(&mut value).factor * atom(&format!("1+{name}"));
        insertion_refusal(
            &refused(&value),
            C0,
            &format!("contains a mass symbol {name}"),
            name,
        );
    }
    // an edge's own mass symbol
    let mut value = gluon_surrogate(1);
    value
        .edges
        .iter_mut()
        .find(|e| e.id == 2)
        .unwrap()
        .mass_squared = atom("helper_test::m2");
    ins(&mut value).coefficient = atom("helper_test::m2");
    insertion_refusal(
        &refused(&value),
        C0,
        "contains a mass symbol ",
        "an edge's mass",
    );
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

/// Weights convention (B): a direct record WITH `original_cycles` leaves the
/// closed-loop sign and n_f to the caller (loop_sign 1, nf_power 0), also when
/// an insertion holds a quark loop; its value carries neither the -1 nor n_f
/// and equals the contracted route. Without `original_cycles` the v1
/// standalone rule stays (the record's own sign and n_f, checked).
#[test]
fn direct_records_with_original_cycles_leave_the_weights_to_the_caller() {
    let (vakint, _) = engine();
    let settings = settings(2);
    let mut value = gluon_surrogate(1);
    ins(&mut value).cycles = (1, 0);
    // (B): every weight is the caller's (the symmetry factor too)
    value.prefactor.overall = Atom::num(1);
    value.record = Some(DirectRecord {
        original_cycles: Some((1, 0)),
        depth: None,
        weights: Vec::new(),
    });
    value.validate().unwrap();
    let prefactor = value.prefactor.total();
    let got = evaluate_factorized(&vakint, &settings, &value, &prefactor).unwrap();
    assert!(!got.value.contains_symbol(symbolica::symbol!(names::NF)));
    let reference = independent_reference(&vakint, &settings, &value);
    assert!(!reference.is_zero());
    assert!(
        (&got.value - &reference * &prefactor).expand().is_zero(),
        "{}\n{}",
        got.value,
        reference
    );
    // the record's own sign under (B) is refused, more cycles than the original too
    let mut signed = value.clone();
    signed.prefactor.loop_sign = -1;
    match signed.validate() {
        Err(FactorizedError::NotScalar { location, detail }) => {
            assert_eq!(location, "prefactor");
            assert!(
                detail.starts_with("with original_cycles the weights are the caller's"),
                "{detail}"
            );
        }
        other => panic!("a sign under (B): {other:?}"),
    }
    let mut more = value.clone();
    ins(&mut more).cycles = (2, 0);
    match more.validate() {
        Err(FactorizedError::Graph { detail }) => assert!(
            detail.starts_with("the record and its insertions carry "),
            "{detail}"
        ),
        other => panic!("more cycles than original: {other:?}"),
    }
    // v1 standalone (no original_cycles): the inner quark loop must be in the
    // record's own weights
    let mut standalone = value.clone();
    standalone.record = None;
    assert!(matches!(
        standalone.validate(),
        Err(FactorizedError::LoopSign { .. })
    ));
    standalone.prefactor.loop_sign = -1;
    standalone.prefactor.nf_power = 1;
    standalone.validate().unwrap();
}

/// THE BETA'S REAL K5 DIRECT RECORD (the surrogate of L2_17 at <= 4 loops:
/// kind direct v3, two open Lorentz ports, one `ins:` coefficient vertex, a
/// ghost loop through it, weights left to the caller): through the ports
/// entry, equal to the value their `integrate_physical_numerator` returns,
/// exactly through its valid_through (eps^-1) and no further; the ghost
/// loop's -1 is not in the value.
#[test]
fn the_beta_k5_direct_record_equals_its_contracted_value() {
    let record = include_str!("fixtures/wh-k5-direct.record.txt");
    let value_text = include_str!("fixtures/wh-k5-direct.value.txt");
    let value = FactorizedNumerator::from_text(record).unwrap_or_else(|e| panic!("{e}"));
    let expected = atom(wh_fixtures::value_line(value_text));
    let settings = VakintSettings {
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 2,
        ..settings(2)
    };
    let (vakint, _) = engine();
    let got = evaluate_factorized_ports(&vakint, &settings, &value, &value.prefactor.total())
        .unwrap_or_else(|e| panic!("{e}"));
    let epsilon = atom("gammalooprs::ε");
    let AtomView::Var(epsilon) = epsilon.as_view() else {
        unreachable!()
    };
    let epsilon = epsilon.get_symbol();
    // the shared comparator: exact through the file's valid_through (eps^-1)
    // and no further, signed, with a nonzero reference
    let through = wh_fixtures::valid_through(value_text);
    assert_eq!(through, -1);
    assert!(
        wh_fixtures::assert_equal_through(&got.value, &expected, epsilon, through, "k5-direct") > 0
    );
    // a wrong oracle (the lowest coefficient doubled) fails
    let coefficients = wh_fixtures::epsilon_coefficients(&expected, epsilon);
    let (low, c) = coefficients.iter().next().expect("a nonzero reference");
    let wrong = &expected + c * Atom::var(epsilon).pow(Atom::num(*low));
    assert!(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            wh_fixtures::assert_equal_through(&got.value, &wrong, epsilon, through, "k5-direct")
        }))
        .is_err(),
        "a wrong k5-direct oracle passed"
    );
    eprintln!("K5_DIRECT identical through eps^-1: {}", got.value);
}

/// Record format section 7: the value of a NEW kind does not depend on the
/// caller's tensor-reduction method or notation: a direct record with an
/// insertion and a ports record (the real k5-direct record through the ports
/// entry), each with every combination (closure records: the closure tests).
#[test]
fn new_kinds_do_not_depend_on_the_reduction_method_or_the_notation() {
    use crate::TensorReductionMethod::{AlphaLoop, FeynKit};
    let (vakint, _) = engine();
    let combinations = [
        (FeynKit, true),
        (FeynKit, false),
        (AlphaLoop, true),
        (AlphaLoop, false),
    ];
    let flipped = |base: &VakintSettings, method, dots| {
        let mut s = base.clone();
        s.tensor_reduction_method = method;
        s.use_dot_product_notation = dots;
        s
    };
    let same = |values: Vec<Atom>, what: &str| {
        assert!(!values[0].is_zero(), "{what}");
        for v in &values[1..] {
            assert!(
                (v - &values[0]).expand().is_zero(),
                "{what}: {v}\nvs {}",
                values[0]
            );
        }
    };
    let direct = gluon_surrogate(1);
    let base = settings(2);
    same(
        combinations
            .iter()
            .map(|(m, d)| {
                evaluate_factorized(&vakint, &flipped(&base, *m, *d), &direct, &Atom::num(1))
                    .unwrap()
                    .value
            })
            .collect(),
        "direct with an insertion",
    );
    let ports =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k5-direct.record.txt")).unwrap();
    let site = VakintSettings {
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        ..settings(2)
    };
    same(
        combinations
            .iter()
            .map(|(m, d)| {
                evaluate_factorized_ports(&vakint, &flipped(&site, *m, *d), &ports, &Atom::num(1))
                    .unwrap()
                    .value
            })
            .collect(),
        "ports record",
    );
}

/// The reserved `vakint::` namespace inside an insertion (marked value): a
/// `vakint::` symbol in the coefficient and a `vakint::` function head in the
/// factor are refused by the record walk (Graph, `reserved symbol <name> in
/// insertion <v> coefficient|factor: `); the same atoms in another namespace
/// are accepted.
#[test]
fn a_reserved_vakint_symbol_inside_an_insertion_is_refused() {
    let good = gluon_surrogate(1);
    good.validate().expect("the pristine record is accepted");
    let reserved = |value: &FactorizedNumerator, prefix: &str, what: &str| match refused(value) {
        FactorizedError::Graph { detail } => {
            assert!(
                detail.starts_with(prefix),
                "{what}: {detail:?} does not start with {prefix:?}"
            )
        }
        other => panic!("{what}: {other:?}"),
    };
    for (name, accepted) in [
        ("vakint::closure_tensor_0", "witness::closure_tensor_0"),
        ("vakint::helper_scalar", "witness::helper_scalar"),
    ] {
        let mut bad = good.clone();
        ins(&mut bad).coefficient = &ins(&mut bad).coefficient * atom(name);
        reserved(
            &bad,
            &format!(
                "reserved symbol {name} in insertion {C0} coefficient: the vakint:: namespace is reserved"
            ),
            name,
        );
        let mut twin = good.clone();
        ins(&mut twin).coefficient = &ins(&mut twin).coefficient * atom(accepted);
        twin.validate()
            .unwrap_or_else(|e| panic!("{accepted}: {e}"));
    }
    // a function head in the factor
    let mut bad = good.clone();
    ins(&mut bad).factor = &ins(&mut bad).factor * atom("vakint::helper_head(2)");
    reserved(
        &bad,
        &format!(
            "reserved symbol vakint::helper_head in insertion {C0} factor: the vakint:: namespace is reserved"
        ),
        "function head",
    );
    // the caller-namespaced twin: as a FUNCTION head refused by the head
    // allow-list, as a plain scalar accepted
    let mut twin = good.clone();
    ins(&mut twin).factor = &ins(&mut twin).factor * atom("witness::helper_head(2)");
    match refused(&twin) {
        FactorizedError::Insertion { detail, .. } => assert!(
            detail.starts_with(&format!(
                "a function head that is not allowed in insertion {C0} factor: witness::helper_head"
            )),
            "{detail}"
        ),
        other => panic!("witness::helper_head(2): {other:?}"),
    }
    let mut twin = good.clone();
    ins(&mut twin).factor = &ins(&mut twin).factor * atom("witness::helper_head");
    twin.validate()
        .unwrap_or_else(|e| panic!("witness::helper_head: {e}"));
}
