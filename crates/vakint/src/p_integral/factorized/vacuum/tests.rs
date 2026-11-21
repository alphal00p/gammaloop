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
        record: None,
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

/// The FLAG-FREE form of the reader rules (the beta's hard part as a vacuum
/// graph of its own lines): the boundary appears only in the factors, as slots
/// of half-edges the graph does not have (rule 1: open ports) and as `Q(x)` of
/// an edge that is not in the graph (rule 2: a constant external vector). The
/// leg form of `one_loop_open_indices_and_an_external_vector` and this form
/// give the same value, in vakint's and in the caller's notation; the
/// flag-free value validates with its ports and does not without them.
#[test]
fn the_flag_free_form_equals_the_leg_form() {
    let factor = format!(
        "{}*{}*{}*{}*{}*{}",
        q(1, &hedge(10)),
        q(1, &hedge(11)),
        q(1, &index(1, 8)),
        q(20, &index(1, 8)),
        q(1, &index(1, 9)),
        q(20, &index(1, 9))
    );
    let ports = vec![atom(&hedge(10)), atom(&hedge(11))];
    let leg_form = value(
        vec![Some(1)],
        vec![Some(20)],
        vec![
            edge(1, 0, 0, (1, 2), &[1, 0], LineKind::Scalar, 3, Atom::num(1)),
            leg(20, 0, 1, 1, 1),
        ],
        vec![vertex(0, atom(&factor))],
        ports.clone(),
    );
    let mut flag_free = value(
        vec![Some(1)],
        Vec::new(),
        vec![edge(
            1,
            0,
            0,
            (1, 2),
            &[1],
            LineKind::Scalar,
            3,
            Atom::num(1),
        )],
        vec![vertex(0, atom(&factor))],
        ports,
    );
    flag_free.label = "1L flag-free k^mu k^nu (k.Q20)^2".into();
    assert_eq!(flag_free.constant_vectors(), vec![20]);
    assert!(leg_form.constant_vectors().is_empty());
    flag_free
        .validate()
        .unwrap_or_else(|e| panic!("the flag-free form validates: {e}"));
    let mut undeclared = flag_free.clone();
    undeclared.ports.clear();
    assert!(undeclared.validate().is_err());
    // the FeynKit ingress derives exactly these ports (rule 1)
    #[cfg(feature = "feynkit-ingress")]
    assert_eq!(undeclared.open_lorentz_ports(), flag_free.ports);
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let leg = evaluate_factorized_vacuum(&vakint, &settings, &leg_form, &Atom::num(1)).unwrap();
    let free = evaluate_factorized_vacuum(&vakint, &settings, &flag_free, &Atom::num(1)).unwrap();
    assert!(!leg.vakint_value.is_zero());
    assert!((&leg.vakint_value - &free.vakint_value).expand().is_zero());
    assert!(
        (&leg.value - &free.value).expand().is_zero(),
        "caller notation: leg form {}\nflag-free {}",
        leg.value,
        free.value
    );
    assert!(!free.value.contains_symbol(S.p));
}

/// The size predictor (H2, D5-HARD-KIND 8.2): the exact number of expanded
/// Lorentz monomials when none merge, an upper bound when they do; a sum free
/// of Minkowski slots is a coefficient (counted once). The resource guard
/// refuses, naming D, when the contract-first bound exceeds the caller's
/// process memory budget, and records D otherwise.
#[test]
fn the_size_predictor_and_the_resource_guard() {
    let terms = |a: &Atom| match a.as_view() {
        AtomView::Add(sum) => sum.iter().count(),
        _ => 1,
    };
    let product = atom(&format!(
        "({}+{}+{})*({}+{})*{METRIC}({},{})",
        q(1, &hedge(10)),
        q(2, &hedge(10)),
        q(3, &hedge(10)),
        q(1, &hedge(11)),
        q(4, &hedge(11)),
        hedge(12),
        hedge(13)
    ));
    assert_eq!(distributed_lorentz_monomials(product.as_view()), 6.0);
    assert_eq!(terms(&product.expand()), 6);
    let scalar = atom(&format!("(x+y+z)*{}", q(1, &hedge(10))));
    assert_eq!(distributed_lorentz_monomials(scalar.as_view()), 1.0);
    let square = atom(&format!("({}+{})^2", q(1, &hedge(10)), q(2, &hedge(10))));
    assert_eq!(distributed_lorentz_monomials(square.as_view()), 4.0);
    assert!(terms(&square.expand()) <= 4);
    let v = value(
        vec![Some(1)],
        vec![Some(20)],
        vec![
            edge(1, 0, 0, (1, 2), &[1, 0], LineKind::Scalar, 3, Atom::num(1)),
            leg(20, 0, 1, 1, 1),
        ],
        vec![vertex(
            0,
            atom(&format!(
                "({}+{})*{}",
                q(1, &hedge(10)),
                q(20, &hedge(10)),
                q(1, &hedge(11))
            )),
        )],
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    let vakint = Vakint::new().unwrap();
    let mut settings = no_form_settings();
    let got = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1)).unwrap();
    assert_eq!(got.predicted_monomials, 2.0);
    settings.rustred_resources = crate::ResourceLimits {
        max_process_rss_bytes: Some(1_000_000),
        ..Default::default()
    };
    let refused = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1))
        .unwrap_err()
        .to_string();
    assert!(
        refused.contains("D = 2.000e0 Lorentz monomials"),
        "{refused}"
    );
}

/// THE STREAMING RULE (D5-HARD-KIND 9.4): one piece up to 1e5 predicted
/// Lorentz monomials, else the smallest power of two not below D / 1e5, at
/// most 32; nf0 #32's hard terms at D = 3.4e4 / 5.3e5 / 3.2e7 / 1.4e8 get
/// 1 / 8 / 32 / 32 pieces.
#[test]
fn the_streaming_rule() {
    assert_eq!(stream_pieces_for(0.0), 1);
    assert_eq!(stream_pieces_for(3.4e4), 1);
    assert_eq!(stream_pieces_for(1.0e5), 1);
    assert_eq!(stream_pieces_for(1.5e5), 2);
    assert_eq!(stream_pieces_for(5.3e5), 8);
    assert_eq!(stream_pieces_for(3.2e7), 32);
    assert_eq!(stream_pieces_for(1.4e8), 32);
    assert_eq!(stream_pieces_for(f64::MAX), 32);
}

/// A sunrise with sums in its numerator (two open ports, one contracted pair,
/// an external vector): `closing` (s, t, momentum) is the massive line in
/// the closure case, every line massive in the family case.
fn streamed_sunrise(closure: bool) -> FactorizedNumerator {
    let factor = format!(
        "({}+{}+{})*({}-{})*({}+{}+2*{})*({}-{})",
        q(1, &hedge(10)),
        q(2, &hedge(10)),
        q(20, &hedge(10)),
        q(1, &hedge(11)),
        q(2, &hedge(11)),
        q(1, &index(1, 8)),
        q(2, &index(1, 8)),
        q(20, &index(1, 8)),
        q(1, &index(1, 8)),
        q(20, &index(1, 8))
    );
    let line = |id: usize, s: usize, t: usize, hedges: (usize, usize), m: &[i64], massive: bool| {
        let e = edge(id, s, t, hedges, m, LineKind::Scalar, 1, Atom::num(1));
        if massive { e } else { massless(e) }
    };
    let mut v = value(
        vec![Some(1), Some(2)],
        vec![Some(20)],
        vec![
            line(1, 1, 2, (1, 2), &[1, 0, 0], !closure),
            line(2, 1, 2, (3, 4), &[0, 1, 0], !closure),
            line(3, 2, 1, (5, 6), &[1, 1, 0], true),
            leg(20, 1, 2, 1, 1),
        ],
        vec![vertex(1, atom(&factor))],
        vec![atom(&hedge(10)), atom(&hedge(11))],
    );
    v.label = format!("2L sunrise with sums, closure={closure}");
    v
}

/// The STREAMED contraction is exact and independent of the piece count: on
/// both back ends (family, closure) the values with 1, 8 and 32 pieces are
/// byte-identical, and equal to the unstreamed path; the pieces are really
/// contracted apart; Auto keeps the unstreamed path below the rule's threshold.
#[test]
fn streaming_is_exact_and_independent_of_the_piece_count() {
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    for closure in [false, true] {
        let v = streamed_sunrise(closure);
        let run = |streaming: Streaming| {
            evaluate_factorized_vacuum_streamed(&vakint, &settings, &v, &Atom::num(1), streaming)
                .unwrap_or_else(|e| panic!("{}: {streaming:?}: {e}", v.label))
        };
        let off = run(Streaming::Off);
        assert_eq!(off.pieces, 0);
        assert!(!off.vakint_value.is_zero(), "{}", v.label);
        let one = run(Streaming::Pieces(1));
        let eight = run(Streaming::Pieces(8));
        let many = run(Streaming::Pieces(32));
        assert_eq!((one.pieces, eight.pieces), (1, 8), "{}", v.label);
        assert!(many.pieces > 8, "{}: {} pieces", v.label, many.pieces);
        let text = one.vakint_value.to_canonical_string();
        assert_eq!(
            eight.vakint_value.to_canonical_string(),
            text,
            "{}",
            v.label
        );
        assert_eq!(many.vakint_value.to_canonical_string(), text, "{}", v.label);
        assert!(
            (&off.vakint_value - &one.vakint_value).expand().is_zero(),
            "{}: unstreamed {}\nstreamed {}",
            v.label,
            off.vakint_value,
            one.vakint_value
        );
        let auto = run(Streaming::Auto);
        assert_eq!(auto.pieces, 0, "{}: small input, unstreamed", v.label);
        assert_eq!(
            auto.vakint_value.to_canonical_string(),
            off.vakint_value.to_canonical_string()
        );
        eprintln!(
            "STREAMING label=\"{}\" pieces=1/8/{} contract_s={:.4}/{:.4}/{:.4} off_s={:.4}",
            v.label,
            many.pieces,
            one.contract_seconds,
            eight.contract_seconds,
            many.contract_seconds,
            off.contract_seconds
        );
    }
}

/// The caller's budget stops a streamed contraction cleanly (an error, not a
/// kill): here an elapsed budget of zero.
#[test]
fn the_streamed_contraction_fails_cleanly_under_the_budget() {
    let vakint = Vakint::new().unwrap();
    let mut settings = no_form_settings();
    settings.rustred_resources = crate::ResourceLimits {
        max_elapsed: Some(std::time::Duration::ZERO),
        ..Default::default()
    };
    let v = streamed_sunrise(true);
    let error = evaluate_factorized_vacuum_streamed(
        &vakint,
        &settings,
        &v,
        &Atom::num(1),
        Streaming::Pieces(8),
    )
    .unwrap_err()
    .to_string();
    assert!(error.contains("streamed contraction"), "{error}");
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
    // the closure back end reduces the integrand as ONE term read off its
    // product (`vacuum_term`): the same term as vakint's `split_integrals` and
    // the same reduced value as `Vakint::tensor_reduce` of the integrand
    let term = vacuum_term(&got.input)
        .unwrap()
        .unwrap_or_else(|| panic!("{}: the integrand is not one term", value.label));
    let split = VakintExpression::try_from(got.input.clone()).unwrap();
    assert_eq!(split.0.len(), 1, "{}", value.label);
    assert_eq!(split.0[0].integral, term.integral, "{}", value.label);
    assert!(
        (&split.0[0].numerator - &term.numerator).expand().is_zero(),
        "{}",
        value.label
    );
    let mut reducing = settings.clone();
    reducing.use_dot_product_notation = true;
    let by_vakint = vakint
        .tensor_reduce(&reducing, got.input.as_view())
        .unwrap();
    let direct = reduce_vacuum(&vakint, &settings, &got.input).unwrap();
    assert!(
        (&by_vakint - &direct).expand().is_zero(),
        "{}:\nVakint::tensor_reduce {by_vakint}\none term             {direct}",
        value.label
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

/// The two-loop sunrise of `two_loop_sunrise_open_indices` (unmarked).
fn sunrise_value() -> FactorizedNumerator {
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
    v.label = "2L sunrise".into();
    // a namespace-qualified mass (a marked record may not carry the bare
    // `muvsq`, read into vakint::)
    for e in &mut v.edges {
        if !e.is_external() {
            e.mass_squared = atom("helper_test::muvsq");
        }
    }
    v
}

/// the vacuum mode walks a MARKED value's own atoms (it never calls
/// `validate`): a reserved head inside an edge factor and a wildcard
/// dimension are refused; the marked value without the defect evaluates to
/// the unmarked one's value; an UNMARKED value is not walked (base).
#[test]
fn the_vacuum_mode_walks_a_marked_record() {
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let run =
        |v: &FactorizedNumerator| evaluate_factorized_vacuum(&vakint, &settings, v, &Atom::num(1));
    let unmarked = sunrise_value();
    let mut marked = unmarked.clone();
    marked.record = Some(super::super::DirectRecord {
        original_cycles: None,
        depth: None,
        weights: Vec::new(),
    });
    assert!(super::super::is_marked(&marked));
    assert!(
        (&run(&marked).unwrap().value - &run(&unmarked).unwrap().value)
            .expand()
            .is_zero()
    );
    let graph = |result: Result<VacuumValue>| match result {
        Err(FactorizedError::Graph { detail }) => detail,
        other => panic!("{other:?}"),
    };
    for (name, head) in [
        ("vakint::dot(vakint::k(1),vakint::k(1))", "vakint::dot"),
        ("vakint::topo(1)", "vakint::topo"),
    ] {
        let mut v = marked.clone();
        v.edges[0].factor = &v.edges[0].factor * atom(name);
        assert_eq!(
            graph(run(&v)),
            format!(
                "reserved symbol {head} in edge 1 factor: the vakint:: namespace is reserved for vakint's internal symbols; a name written without a namespace is read into vakint::"
            )
        );
        // unmarked: not walked (the base behaviour)
        let mut u = unmarked.clone();
        u.edges[0].factor = &u.edges[0].factor * atom(name);
        assert!(
            !matches!(run(&u), Err(FactorizedError::Graph { ref detail }) if detail.starts_with("reserved symbol")),
            "an unmarked value was walked"
        );
    }
    let mut v = marked.clone();
    v.dimension = symbolica::symbol!("UFO::D_");
    assert_eq!(
        graph(run(&v)),
        "a symbol ending in _ cannot be used as dimension: UFO::D_ is a pattern wildcard"
    );
    // the streamed variant and the pub input function walk it too
    let mut v = marked.clone();
    v.edges[0].factor = &v.edges[0].factor * atom("vakint::topo(1)");
    assert!(
        graph(evaluate_factorized_vacuum_streamed(
            &vakint,
            &settings,
            &v,
            &Atom::num(1),
            Streaming::Pieces(2)
        ))
        .starts_with("reserved symbol vakint::topo")
    );
    assert!(matches!(
        vakint_vacuum_input(&settings, &v, &Atom::num(1)),
        Err(FactorizedError::Graph { .. })
    ));
}

/// THE RESULT CHECK on MARKED massive vacuum values: the unit-mass two-loop
/// sunrise, the three-loop banana and the three-loop six-line parent, every
/// line massive (mass witness::msq), marked as v3 direct records, through the
/// vacuum entry in Symbolic and in Numerical mode: no result-check error (the
/// massive-vacuum constants are vakint's output vocabulary); the vakint
/// symbols of the values are printed (`VACUUM_SYMBOLS`).
#[test]
#[ignore = "massive two- and three-loop vacuum values in two modes"]
fn marked_massive_vacuum_values_pass_the_result_check() {
    /// A graph: its name, its lines (id, source, target, loop momentum
    /// coefficients) and its loop count.
    type Graph = (&'static str, Vec<(usize, usize, usize, Vec<i64>)>, usize);
    let graphs: [Graph; 3] = [
        (
            "2L sunrise",
            vec![
                (1, 1, 2, vec![1, 0]),
                (2, 1, 2, vec![0, 1]),
                (3, 2, 1, vec![1, 1]),
            ],
            2,
        ),
        (
            "3L banana",
            vec![
                (1, 1, 2, vec![1, 0, 0]),
                (2, 1, 2, vec![0, 1, 0]),
                (3, 1, 2, vec![0, 0, 1]),
                (4, 2, 1, vec![1, 1, 1]),
            ],
            2,
        ),
        (
            "3L six-line parent",
            vec![
                (1, 1, 2, vec![1, 0, 0]),
                (2, 2, 3, vec![0, 1, 0]),
                (3, 3, 1, vec![0, 0, 1]),
                (4, 1, 4, vec![-1, 0, 1]),
                (5, 2, 4, vec![1, -1, 0]),
                (6, 3, 4, vec![0, 1, -1]),
            ],
            4,
        ),
    ];
    let vakint = Vakint::new().unwrap();
    let mut seen = std::collections::BTreeSet::new();
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions { master_mode: mode },
            )]),
            ..no_form_settings()
        };
        for (name, routes, vertices) in &graphs {
            let loops = routes[0].3.len();
            let edges = routes
                .iter()
                .map(|(id, s, t, m)| {
                    let mut e = edge(
                        *id,
                        *s,
                        *t,
                        (2 * id, 2 * id + 1),
                        m,
                        LineKind::Scalar,
                        1,
                        Atom::num(1),
                    );
                    e.mass_squared = atom("witness::msq");
                    e
                })
                .collect();
            let mut v = value(
                (1..=loops).map(Some).collect(),
                Vec::new(),
                edges,
                (1..=*vertices).map(|i| vertex(i, Atom::num(1))).collect(),
                Vec::new(),
            );
            v.label = format!("{name}, every line massive");
            v.record = Some(super::super::DirectRecord {
                original_cycles: None,
                depth: None,
                weights: Vec::new(),
            });
            let got = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1))
                .unwrap_or_else(|e| panic!("{mode:?} {name}: {e}"));
            assert!(!got.value.is_zero(), "{mode:?} {name}");
            got.value.as_view().visitor(&mut |part| {
                let symbol = match part {
                    AtomView::Var(x) => Some(x.get_symbol()),
                    AtomView::Fun(f) => Some(f.get_symbol()),
                    _ => None,
                };
                if let Some(symbol) = symbol
                    && symbol.get_name().starts_with("vakint::")
                {
                    seen.insert(format!("{mode:?} {}", symbol.get_name()));
                }
                true
            });
        }
    }
    eprintln!("VACUUM_SYMBOLS {seen:?}");
}

/// The BOUNDS reach the vacuum entry, typed and through the text form (the
/// vacuum mode never calls `validate`): a marked two-loop massive sunrise
/// with its record depth usize::MAX, an edge id or a vertex id of
/// i64::MAX + 1, an edge power or a momentum coefficient of i64::MIN or
/// i64::MAX (also 2^31 + 1), and with 1/P(x) (P identically zero) in a vertex
/// factor: each refused before any evaluation (no panic), with the bound's
/// text; the unchanged value and the bound itself (2^31) are accepted.
#[test]
fn the_bounds_reach_the_vacuum_entry_typed_and_through_text() {
    let vakint = Vakint::new().unwrap();
    let settings = no_form_settings();
    let edges = [
        (1, 1, 2, vec![1, 0]),
        (2, 1, 2, vec![0, 1]),
        (3, 2, 1, vec![1, 1]),
    ]
    .iter()
    .map(|(id, s, t, m)| {
        let mut e = edge(
            *id,
            *s,
            *t,
            (2 * id, 2 * id + 1),
            m,
            LineKind::Scalar,
            1,
            Atom::num(1),
        );
        e.mass_squared = atom("witness::msq");
        e
    })
    .collect();
    let mut base = value(
        vec![Some(1), Some(2)],
        Vec::new(),
        edges,
        (1..=2).map(|i| vertex(i, Atom::num(1))).collect(),
        Vec::new(),
    );
    base.label = "2L sunrise, every line massive".into();
    base.record = Some(super::super::DirectRecord {
        original_cycles: None,
        depth: None,
        weights: Vec::new(),
    });
    let big = i64::MAX as usize + 1;
    let renumber_vertex = |v: &mut FactorizedNumerator, id: usize| {
        let old = v.vertices[1].id;
        v.vertices[1].id = id;
        for e in &mut v.edges {
            for end in [&mut e.source, &mut e.target] {
                if *end == super::super::Endpoint::Vertex(old) {
                    *end = super::super::Endpoint::Vertex(id);
                }
            }
        }
    };
    type Change = Box<dyn Fn(&mut FactorizedNumerator)>;
    let cases: Vec<(&str, Change, &str)> = vec![
        (
            "depth usize::MAX",
            Box::new(|v| v.record.as_mut().unwrap().depth = Some(usize::MAX)),
            "the depth",
        ),
        (
            "edge id i64::MAX + 1",
            Box::new(move |v| v.edges[2].id = big),
            "the edge and vertex ids",
        ),
        (
            "vertex id i64::MAX + 1",
            Box::new(move |v| renumber_vertex(v, big)),
            "the edge and vertex ids",
        ),
        (
            "edge power i64::MIN",
            Box::new(|v| v.edges[0].power = i64::MIN),
            "power",
        ),
        (
            "edge power i64::MAX",
            Box::new(|v| v.edges[0].power = i64::MAX),
            "power",
        ),
        (
            "edge power 2^31 + 1",
            Box::new(|v| v.edges[0].power = (1 << 31) + 1),
            "power",
        ),
        (
            "momentum coefficient i64::MIN",
            Box::new(|v| v.edges[0].momentum[0] = i64::MIN),
            "momentum coefficient",
        ),
        (
            "momentum coefficient i64::MAX",
            Box::new(|v| v.edges[0].momentum[0] = i64::MAX),
            "momentum coefficient",
        ),
        (
            "an exponent beyond 64 in a vertex factor",
            Box::new(|v| v.vertices[0].factor = atom("1/((1+UFO::x)^1000000000-UFO::y)")),
            "exceeds 64 in modulus",
        ),
        (
            "a caller-written vakint::slash in a vertex factor",
            Box::new(|v| v.vertices[0].factor = atom("vakint::slash(1,1)")),
            "slash",
        ),
        (
            "1/P(x) in a vertex factor",
            Box::new(|v| {
                v.vertices[0].factor = atom("1/((witness::x-3)^2-witness::x^2+6*witness::x-9)")
            }),
            "a denominator vanishes identically",
        ),
    ];
    for (what, change, text) in &cases {
        let mut v = base.clone();
        change(&mut v);
        // typed (within seconds: the shape rules run first)
        let started = std::time::Instant::now();
        let typed = evaluate_factorized_vacuum(&vakint, &settings, &v, &Atom::num(1));
        assert!(
            started.elapsed().as_secs() < 10,
            "{what}: {:?}",
            started.elapsed()
        );
        match typed {
            Err(e) => assert!(e.to_string().contains(text), "{what} (typed): {e}"),
            Ok(_) => panic!("{what} (typed): accepted"),
        }
        // through the text form: refused when written, read or evaluated
        let refused = match v.to_text() {
            Err(_) => true,
            Ok(written) => match FactorizedNumerator::from_text(&written) {
                Err(_) => true,
                Ok(read) => {
                    evaluate_factorized_vacuum(&vakint, &settings, &read, &Atom::num(1)).is_err()
                }
            },
        };
        assert!(refused, "{what} (text): accepted");
    }
    // the bound itself is accepted (the value is evaluated)
    let mut v = base.clone();
    v.edges[0].power = 1 << 31;
    assert!(v.check_marked_bounds(true).is_ok());
    evaluate_factorized_vacuum(&vakint, &settings, &base, &Atom::num(1)).unwrap();
}
