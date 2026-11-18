//! The vacuum input mode against vakint's existing vacuum evaluation of the
//! HAND-CONTRACTED input, exactly, at one, two and three loops (unit-mass
//! single-scale families), with open indices, an external vector and a Dirac
//! trace; and its CLOSURE back end (one massive line among massless ones)
//! against vakint's massive closure of the HAND-REDUCED, HAND-CUT input, at one,
//! two and three loops, and against the family back end where both apply.
//! Every test runs with RustRed as the only backend and FORM pointed at a path
//! that does not exist (no FORM anywhere).

use super::super::{FactorizedEdge, FactorizedVertex, LineKind, Prefactor, Projector};
use super::*;
use crate::utils::vakint_macros::vk_parse;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode,
};

const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const BIS: &str = "spenso::{spenso::representation,spenso::self_dual}::bis";
const METRIC: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
const GAMMA: &str = "spenso::{linear,spenso::tensor}::gamma";
const MOMENTUM: &str = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const EDGE: &str = "gammalooprs::{spenso::index,spenso::index-label:e}::edge";
const DIM: &str = "gammalooprs::{}::dim";

fn atom(text: &str) -> Atom {
    crate::Vakint::initialize_vakint_symbols();
    vk_parse!(text).unwrap_or_else(|e| panic!("{e}: {text}"))
}

fn hedge(h: usize) -> String {
    format!("{MINK}({DIM},{HEDGE}({h},1))")
}

fn index(e: usize, k: usize) -> String {
    format!("{MINK}({DIM},{EDGE}({e},{k}))")
}

fn q(e: usize, index: &str) -> String {
    format!("{MOMENTUM}({e},{index})")
}

fn gamma(out: usize, into: usize, mu: &str) -> String {
    format!("{GAMMA}({BIS}(4,{HEDGE}({out},1)),{BIS}(4,{HEDGE}({into},1)),{mu})")
}

#[allow(clippy::too_many_arguments)]
fn edge(
    id: usize,
    s: usize,
    t: usize,
    hedges: (usize, usize),
    momentum: &[i64],
    line: LineKind,
    power: i64,
    factor: Atom,
) -> FactorizedEdge {
    FactorizedEdge {
        id,
        source: Endpoint::Vertex(s),
        target: Endpoint::Vertex(t),
        source_hedge: Some(hedges.0),
        target_hedge: Some(hedges.1),
        particle: "x".into(),
        line,
        power,
        mass_squared: atom("muvsq"),
        momentum: momentum.to_vec(),
        factor,
    }
}

/// An external leg `id` into vertex `v` carrying `p(j)` (one-based) of `n`
/// external momenta, after `loops` loop entries.
fn leg(id: usize, v: usize, loops: usize, j: usize, n: usize) -> FactorizedEdge {
    let mut momentum = vec![0; loops + n];
    momentum[loops + j - 1] = 1;
    FactorizedEdge {
        id,
        source: Endpoint::External,
        target: Endpoint::Vertex(v),
        source_hedge: None,
        target_hedge: Some(100 + id),
        particle: "x".into(),
        line: LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum,
        factor: Atom::num(1),
    }
}

fn vertex(id: usize, factor: Atom) -> FactorizedVertex {
    FactorizedVertex {
        id,
        rule: "test".into(),
        factor,
    }
}

fn value(
    loop_edges: Vec<Option<usize>>,
    external_edges: Vec<Option<usize>>,
    edges: Vec<FactorizedEdge>,
    vertices: Vec<FactorizedVertex>,
    ports: Vec<Atom>,
) -> FactorizedNumerator {
    FactorizedNumerator {
        label: "vacuum test".into(),
        kind: FactorizedKind::Direct,
        dimension: FactorizedNumerator::default_dimension(),
        loop_edges,
        external_edges,
        edges,
        vertices,
        ports,
        projector: None,
        prefactor: Prefactor::default(),
        colour_branch: None,
        insertions: Vec::new(),
    }
}

/// RustRed only, FORM at a path that does not exist.
fn no_form_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: "/nonexistent/helper-vh/form".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 3,
        use_dot_product_notation: true,
        n_cores: 1,
        ..Default::default()
    }
}

/// The vacuum mode equals vakint's evaluation of `hand` (vakint notation:
/// ports as raw indices PORT + i, external vectors p(j)) times the topology.
fn check(value: &FactorizedNumerator, hand: &str) -> VacuumValue {
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let got = evaluate_factorized_vacuum(&vakint, &settings, value, &Atom::num(1))
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    assert_eq!(got.backend, VacuumBackend::Family, "{}", value.label);
    let topology = got.input.as_view().to_owned();
    let mut topologies = Vec::new();
    topology.as_view().visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && f.get_symbol() == S.topo
        {
            topologies.push(part.to_owned());
            return false;
        }
        true
    });
    let reference_input = atom(hand) * &topologies[0];
    let reference = vakint
        .evaluate(&settings, reference_input.as_view())
        .unwrap_or_else(|e| panic!("{}: reference: {e}", value.label));
    assert!(
        (&got.vakint_value - &reference).expand().is_zero(),
        "{}:\nvacuum mode {}\nhand        {}",
        value.label,
        got.vakint_value,
        reference
    );
    assert!(!reference.is_zero(), "{}: zero reference", value.label);
    eprintln!(
        "VACUUM_MODE label=\"{}\" contract_s={:.4} evaluate_s={:.4} value={}",
        value.label, got.contract_seconds, got.evaluate_seconds, got.value
    );
    got
}

/// One loop: `k^mu k^nu (k.p)^2` on one massive line of power 3, two open
/// ports and one external vector.
#[test]
fn one_loop_open_indices_and_an_external_vector() {
    let factor = format!(
        "{}*{}*{}*{}*{}*{}",
        q(1, &hedge(10)),
        q(1, &hedge(11)),
        q(1, &index(1, 8)),
        q(20, &index(1, 8)),
        q(1, &index(1, 9)),
        q(20, &index(1, 9))
    );
    let mut v = value(
        vec![Some(1)],
        vec![Some(20)],
        vec![
            edge(1, 0, 0, (1, 2), &[1, 0], LineKind::Scalar, 3, Atom::num(1)),
            leg(20, 0, 1, 1, 1),
        ],
        vec![vertex(0, atom(&factor))],
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "1L k^mu k^nu (k.p)^2".into();
    let got = check(&v, "k(1,9001)*k(1,9002)*dot(k(1),p(1))^2");
    // back in the caller's notation: the ports and the leg's momentum
    assert!(got.value.to_canonical_string().contains("hedge(10,1)"));
    assert!(!got.value.contains_symbol(S.p));
}

/// One loop with a Dirac trace: a massive fermion loop of two lines with the
/// two open ports on its vertices: Tr[gamma^mu k gamma^nu k] = 4 (2 k^mu k^nu
/// - k^2 g^{mu nu}) (D dimensions, Tr 1 = 4).
#[test]
fn one_loop_dirac_trace_with_open_ports() {
    let slash = |e: usize, hs: usize, ht: usize| {
        let mu = index(e, 1);
        atom(&format!("{}*{}", q(e, &mu), gamma(ht, hs, &mu)))
    };
    let mut v = value(
        vec![Some(1)],
        Vec::new(),
        vec![
            edge(1, 0, 1, (1, 2), &[1], LineKind::Fermion, 1, slash(1, 1, 2)),
            edge(2, 1, 0, (3, 4), &[1], LineKind::Fermion, 1, slash(2, 3, 4)),
        ],
        vec![
            vertex(0, atom(&gamma(1, 4, &hedge(10)))),
            vertex(1, atom(&gamma(3, 2, &hedge(11)))),
        ],
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "1L fermion loop trace".into();
    check(&v, "4*(2*k(1,9001)*k(1,9002)-g(9001,9002)*dot(k(1),k(1)))");
    // closed by the caller's projector g(mu, nu): a scalar
    v.projector = Some(Projector {
        tensor: atom(&format!("{METRIC}({},{})", hedge(10), hedge(11))),
        p2_power: 0,
    });
    v.label = "1L fermion loop trace, metric projector".into();
    check(&v, "4*(2-(4-2*ε))*dot(k(1),k(1))");
}

/// Two loops, the unit-mass sunrise: `k1^mu (k2.p)^2 (k1+k2)^nu` (even rank).
#[test]
fn two_loop_sunrise_open_indices() {
    let v0 = format!(
        "{}*{}*{}*{}*{}",
        q(1, &hedge(10)),
        q(2, &index(2, 8)),
        q(20, &index(2, 8)),
        q(2, &index(2, 9)),
        q(20, &index(2, 9))
    );
    let mut v = value(
        vec![Some(1), Some(2)],
        vec![Some(20)],
        vec![
            edge(
                1,
                1,
                2,
                (1, 2),
                &[1, 0, 0],
                LineKind::Scalar,
                1,
                Atom::num(1),
            ),
            edge(
                2,
                1,
                2,
                (3, 4),
                &[0, 1, 0],
                LineKind::Scalar,
                1,
                Atom::num(1),
            ),
            edge(
                3,
                2,
                1,
                (5, 6),
                &[1, 1, 0],
                LineKind::Scalar,
                1,
                atom(&q(3, &hedge(11))),
            ),
            leg(20, 1, 2, 1, 1),
        ],
        vec![vertex(1, atom(&v0))],
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "2L sunrise k1^mu (k2.p)^2 (k1+k2)^nu".into();
    check(&v, "k(1,9001)*dot(k(2),p(1))^2*(k(1,9002)+k(2,9002))");
}

/// Three loops, the unit-mass six-line parent: `k1^mu k2^nu k1.k1`.
#[test]
fn three_loop_parent_open_indices() {
    let routes: [(usize, usize, usize, [i64; 3]); 6] = [
        (1, 1, 2, [1, 0, 0]),
        (2, 2, 3, [0, 1, 0]),
        (3, 3, 1, [0, 0, 1]),
        (4, 1, 4, [-1, 0, 1]),
        (5, 2, 4, [1, -1, 0]),
        (6, 3, 4, [0, 1, -1]),
    ];
    let edges = routes
        .iter()
        .map(|(id, s, t, m)| {
            let factor = match id {
                1 => atom(&format!(
                    "{}*{}*{}*{METRIC}({},{})",
                    q(1, &hedge(10)),
                    q(1, &index(1, 8)),
                    q(1, &index(1, 9)),
                    index(1, 8),
                    index(1, 9)
                )),
                2 => atom(&q(2, &hedge(11))),
                _ => Atom::num(1),
            };
            edge(
                *id,
                *s,
                *t,
                (2 * id, 2 * id + 1),
                m,
                LineKind::Scalar,
                1,
                factor,
            )
        })
        .collect();
    let mut v = value(
        vec![Some(1), Some(2), Some(3)],
        Vec::new(),
        edges,
        Vec::new(),
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "3L parent k1^mu k2^nu k1.k1".into();
    check(&v, "k(1,9001)*k(2,9002)*dot(k(1),k(1))");
}

/// A massless line (the closure tests).
fn massless(mut edge: FactorizedEdge) -> FactorizedEdge {
    edge.mass_squared = Atom::Zero;
    edge
}

/// vakint's massive closure `{p(1), muvsq, power}` of a hand-cut, hand-reduced
/// input (vakint notation: the cut graph in `k(1..L-1)` and `p(1)`, the
/// numerator in epsilon).
fn closure_of(input: &str, power: i64) -> Atom {
    let vakint = Vakint::new().unwrap();
    MassiveClosure {
        momentum: S.p.call(1),
        mass_squared: atom("muvsq"),
        power,
    }
    .evaluate(&vakint, &no_form_settings(), atom(input).as_view())
    .unwrap_or_else(|e| panic!("hand closure {input}: {e}"))
}

/// The closure back end equals `tensor` times the hand closure of `hand`.
fn check_closure(value: &FactorizedNumerator, tensor: &str, hand: &str, power: i64) -> VacuumValue {
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let got = evaluate_factorized_vacuum(&vakint, &settings, value, &Atom::num(1))
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    assert!(
        matches!(got.backend, VacuumBackend::Closure(_)),
        "{}: {:?}",
        value.label,
        got.backend
    );
    let reference = atom(tensor) * closure_of(hand, power);
    assert!(
        (&got.vakint_value - &reference).expand().is_zero(),
        "{}:\nclosure back end {}\nhand             {}",
        value.label,
        got.vakint_value,
        reference
    );
    assert!(!reference.is_zero(), "{}: zero reference", value.label);
    eprintln!(
        "VACUUM_CLOSURE label=\"{}\" contract_s={:.4} evaluate_s={:.4} value={}",
        value.label, got.contract_seconds, got.evaluate_seconds, got.value
    );
    got
}

/// One loop, a massive and a massless line (the closure back end): `k^mu k^nu`
/// with two open ports is, by hand, `g^{mu nu}/D` times the closure of the
/// tree cut `1/p^2`. With `k.k` more in the numerator the massless line
/// cancels and the integral is vakint's FAMILY tadpole of `k^mu k^nu`: both back
/// ends give the same value (one normalisation, one truncation).
#[test]
fn one_loop_closure_massive_and_massless_lines_with_open_ports() {
    let mut v = value(
        vec![Some(1)],
        Vec::new(),
        vec![
            edge(1, 0, 1, (1, 2), &[1], LineKind::Scalar, 1, Atom::num(1)),
            massless(edge(
                2,
                1,
                0,
                (3, 4),
                &[1],
                LineKind::Scalar,
                1,
                Atom::num(1),
            )),
        ],
        vec![vertex(
            0,
            atom(&format!("{}*{}", q(1, &hedge(10)), q(1, &hedge(11)))),
        )],
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "1L closure k^mu k^nu".into();
    let got = check_closure(
        &v,
        "g(9001,9002)",
        "dot(p(1),p(1))/(4-2*ε)*dot(p(1),p(1))^(-1)",
        1,
    );
    let VacuumBackend::Closure(cut) = &got.backend else {
        unreachable!()
    };
    assert_eq!(cut.routing, vec![(atom("k(1)"), atom("p(1)"))]);
    assert_eq!(cut.topology, atom("dot(p(1),p(1))^(-1)"));
    // the massless line cancelled: the family tadpole of k^mu k^nu
    v.vertices[0].factor = atom(&format!(
        "{}*{}*{}*{}*{METRIC}({},{})",
        q(1, &hedge(10)),
        q(1, &hedge(11)),
        q(1, &index(1, 8)),
        q(1, &index(1, 9)),
        index(1, 8),
        index(1, 9)
    ));
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let closure = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1)).unwrap();
    // the family finaliser writes log(1/m^2) where the closure writes -log(m^2)
    let family = vakint
        .evaluate(
            &settings,
            atom("k(1,9001)*k(1,9002)*topo(prop(1,edge(1,1),k(1),muvsq,1))").as_view(),
        )
        .unwrap()
        .replace(atom("log(1/muvsq)"))
        .with(atom("-log(muvsq)"));
    assert!(
        (&closure.vakint_value - &family).expand().is_zero(),
        "closure back end {}\nfamily back end  {}",
        closure.vakint_value,
        family
    );
}

/// Two loops, the sunrise with ONE massive line (k1 + k2) and two massless
/// ones: cut by hand, k(2) -> p(1) - k(1); `k1^mu (k2.p)` (one port, one
/// external vector) reduces to `p^mu k1.k2 / D`, `k1^mu k2^nu` to
/// `g^{mu nu} k1.k2 / D`; and the value does not depend on the closing line's
/// orientation (another routing of the same integral).
#[test]
fn two_loop_closure_one_massive_line() {
    let lines = |closing: (usize, usize, [i64; 3]), factor: String, ports: Vec<Atom>| {
        value(
            vec![Some(1), Some(2)],
            vec![Some(20)],
            vec![
                massless(edge(
                    1,
                    1,
                    2,
                    (1, 2),
                    &[1, 0, 0],
                    LineKind::Scalar,
                    1,
                    Atom::num(1),
                )),
                massless(edge(
                    2,
                    1,
                    2,
                    (3, 4),
                    &[0, 1, 0],
                    LineKind::Scalar,
                    1,
                    Atom::num(1),
                )),
                edge(
                    3,
                    closing.0,
                    closing.1,
                    (5, 6),
                    &closing.2,
                    LineKind::Scalar,
                    1,
                    Atom::num(1),
                ),
                leg(20, 1, 2, 1, 1),
            ],
            vec![vertex(1, atom(&factor))],
            ports,
        )
    };
    let cut = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(1,2),-k(1)+p(1),0,1))";
    let reduced = format!("(dot(k(1),p(1))-dot(k(1),k(1)))/(4-2*ε)*{cut}");
    let one_port = format!(
        "{}*{}*{}",
        q(1, &hedge(10)),
        q(2, &index(2, 8)),
        q(20, &index(2, 8))
    );
    let mut v = lines((2, 1, [1, 1, 0]), one_port.clone(), vec![atom(&hedge(10))]);
    v.label = "2L closure k1^mu (k2.p)".into();
    let got = check_closure(&v, "p(1,9001)", &reduced, 1);
    let VacuumBackend::Closure(derived) = &got.backend else {
        unreachable!()
    };
    assert_eq!(derived.topology, atom(cut));
    // the closing line reversed: another routing, the same value
    let mut reversed = lines((1, 2, [-1, -1, 0]), one_port, vec![atom(&hedge(10))]);
    reversed.label = "2L closure k1^mu (k2.p), closing line reversed".into();
    let again = check_closure(&reversed, "p(1,9001)", &reduced, 1);
    assert!((&got.vakint_value - &again.vakint_value).expand().is_zero());
    let two_ports = format!("{}*{}", q(1, &hedge(10)), q(2, &hedge(11)));
    let mut v = lines(
        (2, 1, [1, 1, 0]),
        two_ports,
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "2L closure k1^mu k2^nu".into();
    check_closure(&v, "g(9001,9002)", &reduced, 1);
}

/// Three loops, the six-line parent with ONE massive line (6: k2 - k3) among
/// five massless ones: cut by hand, k(3) -> k(2) - p(1), a two-loop massless
/// p-integral; `k1^mu k2^nu k1.k1` reduces to `g^{mu nu} k1.k2 k1.k1 / D`.
#[test]
fn three_loop_closure_one_massive_line() {
    let routes: [(usize, usize, usize, [i64; 3]); 6] = [
        (1, 1, 2, [1, 0, 0]),
        (2, 2, 3, [0, 1, 0]),
        (3, 3, 1, [0, 0, 1]),
        (4, 1, 4, [-1, 0, 1]),
        (5, 2, 4, [1, -1, 0]),
        (6, 3, 4, [0, 1, -1]),
    ];
    let edges = routes
        .iter()
        .map(|(id, s, t, m)| {
            let factor = match id {
                1 => atom(&format!(
                    "{}*{}*{}*{METRIC}({},{})",
                    q(1, &hedge(10)),
                    q(1, &index(1, 8)),
                    q(1, &index(1, 9)),
                    index(1, 8),
                    index(1, 9)
                )),
                2 => atom(&q(2, &hedge(11))),
                _ => Atom::num(1),
            };
            let line = edge(
                *id,
                *s,
                *t,
                (2 * id, 2 * id + 1),
                m,
                LineKind::Scalar,
                1,
                factor,
            );
            if *id == 6 { line } else { massless(line) }
        })
        .collect();
    let mut v = value(
        vec![Some(1), Some(2), Some(3)],
        Vec::new(),
        edges,
        Vec::new(),
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = "3L closure k1^mu k2^nu k1.k1".into();
    let got = check_closure(
        &v,
        "g(9001,9002)",
        "dot(k(1),k(2))*dot(k(1),k(1))/(4-2*ε)*topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,3),k(2),0,1)*prop(3,edge(3,1),k(2)-p(1),0,1)*prop(4,edge(1,4),-k(1)+k(2)-p(1),0,1)*prop(5,edge(2,4),k(1)-k(2),0,1))",
        1,
    );
    let VacuumBackend::Closure(cut) = &got.backend else {
        unreachable!()
    };
    assert_eq!(
        cut.routing,
        vec![
            (atom("k(1)"), atom("k(1)")),
            (atom("k(2)"), atom("k(2)")),
            (atom("k(3)"), atom("k(2)-p(1)")),
        ]
    );
}

/// Three loops, two massless bubbles IN SERIES with the massive line (the
/// bare remainder of the beta's L3 diagram 185 at lambda^0, scalar numerator):
/// the cut keeps the lines that carry p alone (5, and 8 merged from two lines
/// in series) inside the two-point graph; by hand the cut as gammalooprs's
/// VacuumCut writes it (the two p-lines as one of power 3, other labels).
#[test]
fn three_loop_closure_two_bubbles_in_series() {
    let lines: [(usize, usize, usize, [i64; 3], i64); 8] = [
        (2, 1, 0, [1, 0, 0], 1),
        (3, 0, 1, [1, -1, 0], 1),
        (4, 0, 5, [0, 1, 0], 1),
        (5, 1, 3, [0, -1, 0], 1),
        (6, 3, 2, [0, 0, 1], 1),
        (7, 2, 3, [0, 1, 1], 1),
        (8, 2, 4, [0, -1, 0], 1),
        (9, 4, 5, [0, -1, 0], 1),
    ];
    let edges = lines
        .iter()
        .map(|(id, s, t, m, power)| {
            let line = edge(
                *id,
                *s,
                *t,
                (2 * id, 2 * id + 1),
                m,
                LineKind::Scalar,
                *power,
                Atom::num(1),
            );
            if *id == 4 { line } else { massless(line) }
        })
        .collect();
    let mut v = value(
        vec![Some(2), Some(4), Some(6)],
        Vec::new(),
        edges,
        Vec::new(),
        Vec::new(),
    );
    v.label = "3L closure, two bubbles in series".into();
    let got = check_closure(
        &v,
        "1",
        "topo(prop(1,edge(0,1),k(1),0,1)*prop(2,edge(2,0),p(1),0,3)*prop(3,edge(2,3),k(2),0,1)*prop(4,edge(1,0),-p(1)+k(1),0,1)*prop(6,edge(3,2),k(2)+p(1),0,1))",
        1,
    );
    let VacuumBackend::Closure(cut) = &got.backend else {
        unreachable!()
    };
    assert_eq!(
        cut.topology,
        atom(
            "topo(prop(2,edge(1,0),k(1),0,1)*prop(3,edge(0,1),k(1)-p(1),0,1)*prop(5,edge(1,3),-p(1),0,1)*prop(6,edge(3,2),k(2),0,1)*prop(7,edge(2,3),k(2)+p(1),0,1)*prop(8,edge(2,5),-p(1),0,2))"
        )
    );
}

/// A massive tadpole and a massless tadpole at one vertex: scaleless, zero;
/// two massive lines among massless ones: refused (neither a family nor one
/// closure).
#[test]
fn the_closure_back_end_is_zero_when_scaleless_and_refuses_two_massive_lines() {
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let v = value(
        vec![Some(1), Some(2)],
        Vec::new(),
        vec![
            edge(1, 0, 0, (1, 2), &[1, 0], LineKind::Scalar, 1, Atom::num(1)),
            massless(edge(
                2,
                0,
                0,
                (3, 4),
                &[0, 1],
                LineKind::Scalar,
                1,
                Atom::num(1),
            )),
        ],
        Vec::new(),
        Vec::new(),
    );
    let got = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1)).unwrap();
    assert_eq!(got.backend, VacuumBackend::Scaleless);
    assert!(got.vakint_value.is_zero());
    let v = value(
        vec![Some(1), Some(2)],
        Vec::new(),
        vec![
            edge(1, 1, 2, (1, 2), &[1, 0], LineKind::Scalar, 1, Atom::num(1)),
            edge(2, 1, 2, (3, 4), &[0, 1], LineKind::Scalar, 1, Atom::num(1)),
            massless(edge(
                3,
                2,
                1,
                (5, 6),
                &[1, 1],
                LineKind::Scalar,
                1,
                Atom::num(1),
            )),
        ],
        Vec::new(),
        Vec::new(),
    );
    let error = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1))
        .unwrap_err()
        .to_string();
    assert!(error.contains("2 massive lines among 3 lines"), "{error}");
}

#[test]
fn the_vacuum_mode_refuses_what_is_outside_its_scope() {
    // a propagator graph (an internal line carries the external momentum)
    let mut v = value(
        vec![Some(1)],
        vec![Some(20)],
        vec![
            edge(1, 0, 0, (1, 2), &[1, 1], LineKind::Scalar, 1, Atom::num(1)),
            leg(20, 0, 1, 1, 1),
        ],
        Vec::new(),
        Vec::new(),
    );
    assert!(contract_to_vakint_vacuum(&v).is_err());
    // two mass scales
    v.edges[0].momentum = vec![1, 0];
    let mut other = edge(2, 0, 0, (3, 4), &[1, 0], LineKind::Scalar, 1, Atom::num(1));
    other.mass_squared = atom("mother");
    v.edges.push(other);
    v.edges.sort_by_key(|e| e.id);
    assert!(contract_to_vakint_vacuum(&v).is_err());
}
