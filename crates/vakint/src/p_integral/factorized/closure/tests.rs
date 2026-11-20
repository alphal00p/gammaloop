//! The closure kind against the CONTRACTED reference chain, exactly in ε, on
//! hand-built R* remainders at two, three and four loops: root remainders
//! closed by the metric projector (K1), children with open Lorentz ports and
//! independent externals (K2), a pre-expanded sector (form B); every refusal;
//! the v3 text form.
//!
//! The reference is independent of the entry's route: idenso contracts the
//! WHOLE numerator (the same production chain), the Taylor step is
//! symbolica's `series` of the full integrand (scaled momenta and the mass
//! substitution applied to the contracted scalar, every denominator a symbol
//! `a_e + B λ + C λ²`), each denominator pattern goes through the PUBLIC
//! `Vakint::tensor_reduce`, is cut by hand and closed by the ORDINARY
//! `MassiveClosure::evaluate` (the entry: graded blocks, binomial den series,
//! memoised coefficient reduction, `RemainderBoundary` direct lowering).
//! RustRed only, FORM at a path that does not exist.

use std::collections::BTreeMap;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use super::super::vacuum::map_to_vakint_vacuum;
use super::super::{
    ClosureSpec, Endpoint, FactorizedEdge, FactorizedError, FactorizedKind, FactorizedNumerator,
    FactorizedVertex, LineKind, Prefactor, Projector, TaylorOrders, TaylorSpec,
    contract_with_idenso,
};
use super::*;
use crate::p_integral::MassiveClosure;
use crate::utils::vakint_macros::vk_parse;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode,
};

const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const METRIC: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
const MOMENTUM: &str = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const EDGE: &str = "gammalooprs::{spenso::index,spenso::index-label:e}::edge";
const DIM: &str = "gammalooprs::{}::dim";
const R: &str = "gammalooprs::{scalar}::rstar_reference_mass_squared";
const MAUX: &str = "gammalooprs::{scalar}::rstar_mass_squared";

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

fn line(id: usize, s: usize, t: usize, momentum: &[i64], mass: &str) -> FactorizedEdge {
    FactorizedEdge {
        id,
        source: Endpoint::Vertex(s),
        target: Endpoint::Vertex(t),
        source_hedge: Some(2 * id),
        target_hedge: Some(2 * id + 1),
        particle: "x".into(),
        line: LineKind::Scalar,
        power: 1,
        mass_squared: atom(mass),
        momentum: momentum.to_vec(),
        factor: Atom::num(1),
    }
}

/// A leg `id` at vertex `v` (incoming when `into`) with momentum `momentum`;
/// its half-edge is `100 + id`.
fn leg(id: usize, v: usize, into: bool, momentum: &[i64]) -> FactorizedEdge {
    let (source, target, sh, th) = if into {
        (
            Endpoint::External,
            Endpoint::Vertex(v),
            None,
            Some(100 + id),
        )
    } else {
        (
            Endpoint::Vertex(v),
            Endpoint::External,
            Some(100 + id),
            None,
        )
    };
    FactorizedEdge {
        id,
        source,
        target,
        source_hedge: sh,
        target_hedge: th,
        particle: "x".into(),
        line: LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: momentum.to_vec(),
        factor: Atom::num(1),
    }
}

fn vertex(id: usize, factor: &str) -> FactorizedVertex {
    FactorizedVertex {
        id,
        rule: "test".into(),
        factor: if factor.is_empty() {
            Atom::num(1)
        } else {
            atom(factor)
        },
    }
}

fn taylor(degree: u32, scaled_edges: Vec<usize>, compensate: bool) -> Option<TaylorSpec> {
    Some(TaylorSpec {
        degree,
        scaled_edges,
        compensate_mass: compensate,
        orders: TaylorOrders::All,
    })
}

#[allow(clippy::too_many_arguments)]
fn record(
    label: &str,
    loops: usize,
    external_edges: Vec<Option<usize>>,
    mut edges: Vec<FactorizedEdge>,
    vertices: Vec<FactorizedVertex>,
    ports: Vec<Atom>,
    projector: Option<Projector>,
    closing: usize,
    taylor: Option<TaylorSpec>,
    independent: bool,
) -> FactorizedNumerator {
    edges.sort_by_key(|e| e.id);
    let power = edges.iter().find(|e| e.id == closing).unwrap().power;
    let mut loop_edges = vec![None; loops];
    for edge in &edges {
        if let Some(i) = (0..loops).find(|i| {
            edge.momentum
                .iter()
                .enumerate()
                .all(|(j, c)| *c == i64::from(j == *i))
        }) {
            loop_edges[i] = Some(edge.id);
        }
    }
    let auxiliary = if taylor.is_some() {
        atom(MAUX)
    } else {
        atom(R)
    };
    FactorizedNumerator {
        label: label.into(),
        kind: FactorizedKind::Closure(ClosureSpec {
            edge: closing,
            mass_squared: atom(R),
            auxiliary_mass_squared: auxiliary,
            power,
            independent_externals: independent,
            taylor,
            softs: Vec::new(),
            original_cycles: Some((0, 0)),
            depth: None,
            weights: Vec::new(),
        }),
        dimension: FactorizedNumerator::default_dimension(),
        loop_edges,
        external_edges,
        edges,
        vertices,
        ports,
        projector,
        prefactor: Prefactor {
            scalars: atom(&format!("(-1+{DIM})^(-1)*gs^2")),
            ..Prefactor::default()
        },
        colour_branch: None,
        insertions: Vec::new(),
        record: None,
    }
}

fn root_projector() -> Option<Projector> {
    Some(Projector {
        tensor: atom(&format!("{METRIC}({},{})", hedge(120), hedge(121))),
        p2_power: -1,
    })
}

/// L2 root remainder: the sunrise with the closing line 3 (k2, massive),
/// q through line 1 (k1 + q); vertex 1 carries `q1^mu (q1.q3)`, vertex 2
/// `q2^nu`, closed by `g(mu, nu) / p^2`.
fn l2_root(degree: u32, compensate: bool) -> FactorizedNumerator {
    record(
        &format!("L2 root w={degree} compensate={compensate}"),
        2,
        vec![Some(20)],
        vec![
            line(1, 1, 2, &[1, 0, 1], "0"),
            line(2, 1, 2, &[-1, 1, 0], "0"),
            line(3, 2, 1, &[0, 1, 0], R),
            leg(20, 1, true, &[0, 0, 1]),
            leg(21, 2, false, &[0, 0, 1]),
        ],
        vec![
            vertex(
                1,
                &format!(
                    "{}*{}*{}",
                    q(1, &hedge(120)),
                    q(1, &index(1, 7)),
                    q(3, &index(1, 7))
                ),
            ),
            vertex(2, &q(2, &hedge(121))),
        ],
        vec![atom(&hedge(120)), atom(&hedge(121))],
        root_projector(),
        3,
        taylor(degree, vec![20, 21], compensate),
        false,
    )
}

/// L2 child: the same graph with two INDEPENDENT crown momenta (p1 into vertex
/// 1, p2 into vertex 2, line 1 = k1 + p1), the two Lorentz ports open.
fn l2_child(degree: u32) -> FactorizedNumerator {
    record(
        &format!("L2 child w={degree}"),
        2,
        vec![Some(20), Some(21)],
        vec![
            line(1, 1, 2, &[1, 0, 1, 0], "0"),
            line(2, 1, 2, &[-1, 1, 0, 0], "0"),
            line(3, 2, 1, &[0, 1, 0, 0], R),
            leg(20, 1, true, &[0, 0, 1, 0]),
            leg(21, 2, true, &[0, 0, 0, 1]),
        ],
        vec![
            vertex(1, &q(1, &hedge(120))),
            vertex(
                2,
                &format!(
                    "{}*({}*{}+{}*{})",
                    q(2, &hedge(121)),
                    q(21, &index(2, 7)),
                    q(1, &index(2, 7)),
                    q(3, &index(2, 7)),
                    q(1, &index(2, 7))
                ),
            ),
        ],
        vec![atom(&hedge(120)), atom(&hedge(121))],
        None,
        3,
        taylor(degree, vec![20, 21], true),
        true,
    )
}

/// L3 root remainder: the tetrahedron, closing line 3 (k3), q through line 1.
fn l3_root(degree: u32) -> FactorizedNumerator {
    record(
        &format!("L3 root w={degree}"),
        3,
        vec![Some(20)],
        vec![
            line(1, 1, 2, &[1, 0, 0, 1], "0"),
            line(2, 2, 3, &[0, 1, 0, 0], "0"),
            line(3, 3, 1, &[0, 0, 1, 0], R),
            line(4, 1, 4, &[-1, 0, 1, 0], "0"),
            line(5, 2, 4, &[1, -1, 0, 0], "0"),
            line(6, 3, 4, &[0, 1, -1, 0], "0"),
            leg(20, 1, true, &[0, 0, 0, 1]),
            leg(21, 2, false, &[0, 0, 0, 1]),
        ],
        vec![
            vertex(
                1,
                &format!(
                    "{}*{}*{}",
                    q(1, &hedge(120)),
                    q(3, &index(1, 7)),
                    q(4, &index(1, 7))
                ),
            ),
            vertex(2, &q(2, &hedge(121))),
            vertex(3, ""),
            vertex(4, ""),
        ],
        vec![atom(&hedge(120)), atom(&hedge(121))],
        root_projector(),
        3,
        taylor(degree, vec![20, 21], true),
        false,
    )
}

/// L4 root remainder: the tetrahedron with line 2 split into a bubble
/// (2: k2 - k4, 7: k4), closing line 3 (k3).
fn l4_root(degree: u32) -> FactorizedNumerator {
    record(
        &format!("L4 root w={degree}"),
        4,
        vec![Some(20)],
        vec![
            line(1, 1, 2, &[1, 0, 0, 0, 1], "0"),
            line(2, 2, 3, &[0, 1, 0, -1, 0], "0"),
            line(7, 2, 3, &[0, 0, 0, 1, 0], "0"),
            line(3, 3, 1, &[0, 0, 1, 0, 0], R),
            line(4, 1, 4, &[-1, 0, 1, 0, 0], "0"),
            line(5, 2, 4, &[1, -1, 0, 0, 0], "0"),
            line(6, 3, 4, &[0, 1, -1, 0, 0], "0"),
            leg(20, 1, true, &[0, 0, 0, 0, 1]),
            leg(21, 2, false, &[0, 0, 0, 0, 1]),
        ],
        vec![
            vertex(1, &q(1, &hedge(120))),
            vertex(
                2,
                &format!(
                    "{}*{}*{}",
                    q(7, &hedge(121)),
                    q(5, &index(2, 7)),
                    q(2, &index(2, 7))
                ),
            ),
            vertex(3, ""),
            vertex(4, ""),
        ],
        vec![atom(&hedge(120)), atom(&hedge(121))],
        root_projector(),
        3,
        taylor(degree, vec![20, 21], true),
        false,
    )
}

/// A pre-expanded sector (form B) of the L2 root: lambda^0 denominators (line
/// 1 at k1), the den-side `(l1.q)` factor written with Q(20), line 1 at power 2.
fn l2_sector() -> FactorizedNumerator {
    let mut value = record(
        "L2 sector (l1.q) / l1^4",
        2,
        vec![Some(20)],
        vec![
            line(1, 1, 2, &[1, 0, 0], "0"),
            line(2, 1, 2, &[-1, 1, 0], "0"),
            line(3, 2, 1, &[0, 1, 0], R),
            leg(20, 1, true, &[0, 0, 1]),
            leg(21, 2, false, &[0, 0, 1]),
        ],
        vec![
            vertex(
                1,
                &format!(
                    "{}*{}*{}*{}*{}",
                    q(1, &hedge(120)),
                    q(1, &index(1, 7)),
                    q(3, &index(1, 7)),
                    q(1, &index(1, 8)),
                    q(20, &index(1, 8))
                ),
            ),
            vertex(
                2,
                &format!(
                    "{}*{}*{}",
                    q(2, &hedge(121)),
                    q(20, &index(2, 7)),
                    q(2, &index(2, 7))
                ),
            ),
        ],
        vec![atom(&hedge(120)), atom(&hedge(121))],
        root_projector(),
        3,
        None,
        false,
    );
    value.edges[0].power = 2;
    value
}

fn settings() -> VakintSettings {
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

/// The external-tensor split of a reduced numerator (as the vacuum entry's).
fn external_split(numerator: &Atom) -> Vec<(Atom, Atom)> {
    let mut external = std::collections::BTreeSet::new();
    numerator.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            let external_dot = symbol == S.dot
                && function
                    .iter()
                    .all(|v| matches!(v, AtomView::Fun(f) if f.get_symbol() == S.p));
            if symbol == S.g || (symbol == S.p && function.get_nargs() == 2) || external_dot {
                external.insert(part.to_owned());
                return false;
            }
        }
        true
    });
    let external = external.into_iter().collect::<Vec<_>>();
    if external.is_empty() {
        vec![(Atom::num(1), numerator.clone())]
    } else {
        numerator.coefficient_list::<i8>(&external)
    }
}

/// The CONTRACTED reference chain (module documentation). `routing`: the hand
/// cut (vacuum k(i) -> cut momenta / head p(1)); `closing`: the closing edge.
fn reference(value: &FactorizedNumerator, routing: &[(&str, &str)]) -> Atom {
    let vakint = Vakint::new().unwrap();
    let settings = settings();
    let FactorizedKind::Closure(spec) = &value.kind else {
        unreachable!()
    };
    let loops = value.loops();
    let lambda = Atom::var(symbolica::symbol!("reference_lambda"));
    let epsilon = Atom::var(crate::utils::vakint_macros::vk_symbol!(
        &settings.epsilon_symbol
    ));
    let in_epsilon = Atom::num(4) - Atom::num(2) * &epsilon;
    // the whole numerator, contracted by idenso
    let mut network = value.factor_product();
    if let Some(projector) = &value.projector {
        network *= &projector.tensor;
    }
    let contracted = contract_with_idenso(&network, false).unwrap();
    let mut scalar = map_to_vakint_vacuum(value, &contracted)
        .unwrap()
        .replace(Atom::var(value.dimension))
        .with(in_epsilon.clone())
        .expand();
    let (w, compensate) = match &spec.taylor {
        Some(t) => (t.degree as usize, t.compensate_mass),
        None => (0, false),
    };
    let closing_mass = spec.auxiliary_mass_squared.clone();
    // the Taylor step on the contracted scalar: p -> lambda p, R -> ..
    if spec.taylor.is_some() {
        scalar = scalar.replace_map(|part, _, out| {
            let AtomView::Fun(f) = part else { return };
            if f.get_symbol() == S.dot {
                let n = f
                    .iter()
                    .filter(|a| matches!(a, AtomView::Fun(g) if g.get_symbol() == S.p))
                    .count();
                if n > 0 {
                    **out = lambda.pow(Atom::num(n as i64)) * part.to_owned();
                }
            } else if f.get_symbol() == S.p && f.get_nargs() == 2 {
                **out = &lambda * part.to_owned();
            }
        });
        let image = if compensate {
            (Atom::num(1) - lambda.pow(Atom::num(2))) * &closing_mass
        } else {
            closing_mass.clone()
        };
        scalar = scalar.replace(spec.mass_squared.clone()).with(image);
    }
    // every denominator a symbol a_e + B lambda + C lambda^2
    let vector = |head, c: &[i64]| {
        let mut out = Atom::Zero;
        for (i, x) in c.iter().enumerate() {
            if *x != 0 {
                out += Atom::num(*x) * function!(head, Atom::num(i as i64 + 1));
            }
        }
        out
    };
    let a = |id: usize| Atom::var(symbolica::symbol!(format!("reference_den_{id}")));
    let mut integrand = scalar;
    let internal = value
        .edges
        .iter()
        .filter(|e| !e.is_external())
        .collect::<Vec<_>>();
    for edge in &internal {
        let l = vector(S.k, &edge.momentum[..loops]);
        let r = vector(S.p, &edge.momentum[loops..]);
        let mut den = a(edge.id);
        if spec.taylor.is_some() {
            if !r.is_zero() {
                den += Atom::num(2) * &lambda * S.dot.call((&l, &r))
                    + lambda.pow(Atom::num(2)) * S.dot.call((&r, &r));
            }
            if edge.id == spec.edge && compensate {
                den += lambda.pow(Atom::num(2)) * &closing_mass;
            }
        }
        integrand *= den.pow(Atom::num(-edge.power));
    }
    let integrand = if let Some(taylor) = &spec.taylor {
        let series = integrand
            .series(symbolica::symbol!("reference_lambda"), 0, w)
            .unwrap();
        let mut kept = Atom::Zero;
        let orders = taylor.orders;
        for (power, coefficient) in series.terms() {
            let n = power.numerator().to_i64().unwrap();
            assert!(power.is_integer() && n >= 0);
            if (n as usize) <= w && orders.keeps(n as u32) {
                kept += coefficient;
            }
        }
        kept
    } else {
        integrand
    };
    // patterns: a_e -> 1/b_e, then the (positive) powers of the b_e
    let b_of = |id: usize| Atom::var(symbolica::symbol!(format!("reference_inv_{id}")));
    let mut integrand = integrand;
    for edge in &internal {
        integrand = integrand
            .replace(a(edge.id))
            .with(b_of(edge.id).pow(Atom::num(-1)));
    }
    let integrand = integrand.expand();
    let inverses = internal.iter().map(|e| b_of(e.id)).collect::<Vec<_>>();
    let by_pattern = integrand.coefficient_list::<i8>(&inverses);
    let prefactor = value
        .prefactor
        .total()
        .replace(Atom::var(value.dimension))
        .with(in_epsilon.clone());
    let routing = routing
        .iter()
        .map(|(k, image)| symbolica::id::Replacement::new(atom(k).to_pattern(), atom(image)))
        .collect::<Vec<_>>();
    let mut total = Atom::Zero;
    for (monomial, numerator) in by_pattern {
        let mut powers = BTreeMap::new();
        for edge in &internal {
            let mut power = 0;
            monomial.as_view().visitor(&mut |part| {
                if part == b_of(edge.id).as_view() {
                    power = 1;
                }
                if let AtomView::Pow(p) = part {
                    let (base, exp) = p.get_base_exp();
                    if base == b_of(edge.id).as_view() {
                        power = i64::try_from(exp).unwrap();
                        return false;
                    }
                }
                true
            });
            assert!(power > 0, "{monomial}: edge {} not a denominator", edge.id);
            powers.insert(edge.id, power);
        }
        let mut topology = Atom::num(1);
        let mut cut = Atom::num(1);
        for edge in &internal {
            let (Endpoint::Vertex(s), Endpoint::Vertex(t)) = (edge.source, edge.target) else {
                unreachable!()
            };
            let mass = if edge.id == spec.edge {
                closing_mass.clone()
            } else {
                Atom::Zero
            };
            let l = vector(S.k, &edge.momentum[..loops]);
            topology *= S.prop.call((
                edge.id as i64,
                S.edge.call((s as i64, t as i64)),
                l.clone(),
                mass,
                Atom::num(powers[&edge.id]),
            ));
            if edge.id != spec.edge {
                cut *= S.prop.call((
                    edge.id as i64,
                    S.edge.call((s as i64, t as i64)),
                    l.replace_multiple(&routing).expand(),
                    Atom::Zero,
                    Atom::num(powers[&edge.id]),
                ));
            }
        }
        let input = &prefactor * &numerator * S.topo.call(topology);
        let reduced = vakint.tensor_reduce(&settings, input.as_view()).unwrap();
        let expression = crate::VakintExpression::try_from(reduced).unwrap();
        let closure = MassiveClosure {
            momentum: S.p.call(1),
            mass_squared: closing_mass.clone(),
            power: powers[&spec.edge],
        };
        for term in &expression.0 {
            for (tensor, coefficient) in external_split(&term.numerator.expand()) {
                assert!(
                    !coefficient.contains_symbol(S.p),
                    "{}: external momentum left in {coefficient}",
                    value.label
                );
                let hand = (coefficient.replace_multiple(&routing).expand()
                    * S.topo.call(cut.clone()))
                .expand();
                let closed = closure
                    .evaluate(&vakint, &settings, hand.as_view())
                    .unwrap();
                total += tensor * closed;
            }
        }
    }
    if let Some(projector) = &value.projector
        && projector.p2_power != 0
    {
        total *= S
            .dot
            .call((S.p.call(1), S.p.call(1)))
            .pow(Atom::num(projector.p2_power));
    }
    total.expand()
}

fn check(value: &FactorizedNumerator, routing: &[(&str, &str)]) -> ClosureRecordValue {
    value
        .validate()
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &settings(),
        std::slice::from_ref(value),
        &[value.prefactor.total()],
    )
    .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    let started = std::time::Instant::now();
    let expected = reference(value, routing);
    let reference_seconds = started.elapsed().as_secs_f64();
    let record = got.records[0].clone();
    assert!(!expected.is_zero(), "{}: zero reference", value.label);
    assert!(
        (&record.vakint_value - &expected).expand().is_zero(),
        "{}:\nentry     {}\nreference {}",
        value.label,
        record.vakint_value,
        expected
    );
    eprintln!(
        "CLOSURE_KIND label=\"{}\" patterns={} blocks={} reduced={} contract_s={:.4} split_s={:.4} evaluate_s={:.4} lower_s={:.4} reduce_s={:.4} reference_s={:.4}",
        value.label,
        record.patterns,
        record.blocks,
        record.reduced_monomials,
        record.contract_seconds,
        record.split_seconds,
        record.evaluate_seconds,
        record.lower_seconds,
        record.reduce_seconds,
        reference_seconds
    );
    record
}

const L2_CUT: [(&str, &str); 2] = [("k(1)", "k(1)"), ("k(2)", "p(1)")];
const L3_CUT: [(&str, &str); 3] = [("k(1)", "k(1)"), ("k(2)", "k(2)"), ("k(3)", "p(1)")];
const L4_CUT: [(&str, &str); 4] = [
    ("k(1)", "k(1)"),
    ("k(2)", "k(2)"),
    ("k(3)", "p(1)"),
    ("k(4)", "k(3)"),
];

#[test]
fn l2_root_remainder_equals_the_contracted_reference() {
    for (degree, compensate) in [(0, true), (1, true), (2, true), (2, false)] {
        let record = check(&l2_root(degree, compensate), &L2_CUT);
        if degree == 2 {
            assert!(record.patterns > 1, "the Taylor step raised no power");
        }
    }
}

#[test]
fn l2_root_orders_filters() {
    for orders in [TaylorOrders::DropZero, TaylorOrders::DropZeroOdd] {
        let mut value = l2_root(2, false);
        let FactorizedKind::Closure(spec) = &mut value.kind else {
            unreachable!()
        };
        spec.taylor.as_mut().unwrap().orders = orders;
        value.label = format!("L2 root orders={}", orders.as_str());
        check(&value, &L2_CUT);
    }
}

/// K3 shape: an IR hard part at FIXED mass (R -> M_aux, no compensation),
/// independent crown momenta, open Lorentz ports, all orders.
#[test]
fn l2_hard_part_at_fixed_mass_equals_the_contracted_reference() {
    for degree in [0, 1, 2] {
        let mut value = l2_child(degree);
        spec(&mut value).taylor.as_mut().unwrap().compensate_mass = false;
        value.label = format!("L2 hard part fixed mass w={degree}");
        check(&value, &L2_CUT);
    }
}

#[test]
fn l2_child_with_open_ports_equals_the_contracted_reference() {
    for degree in [0, 1, 2] {
        check(&l2_child(degree), &L2_CUT);
    }
}

/// A remainder whose numerator vanishes at the Taylor degree: exactly zero.
#[test]
fn a_vanishing_remainder_is_exactly_zero() {
    let mut value = l2_child(0);
    value.vertices[1].factor = atom(&format!(
        "{}*{}*{}",
        q(2, &hedge(121)),
        q(21, &index(2, 7)),
        q(1, &index(2, 7))
    ));
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &settings(),
        std::slice::from_ref(&value),
        &[Atom::num(1)],
    )
    .unwrap();
    assert!(got.value.is_zero());
    assert!(reference(&value, &L2_CUT).is_zero());
}

#[test]
fn l2_pre_expanded_sector_equals_the_contracted_reference() {
    check(&l2_sector(), &L2_CUT);
}

#[test]
fn l3_root_remainder_equals_the_contracted_reference() {
    for degree in [0, 2] {
        check(&l3_root(degree), &L3_CUT);
    }
}

#[test]
fn l4_root_remainder_equals_the_contracted_reference() {
    check(&l4_root(2), &L4_CUT);
}

/// All records in one call: the sum of the single calls.
#[test]
fn several_records_in_one_call() {
    let vakint = Vakint::new().unwrap();
    let records = [l2_root(2, true), l2_child(1), l2_sector()];
    let prefactors = records
        .iter()
        .map(|r| r.prefactor.total())
        .collect::<Vec<_>>();
    let all = evaluate_factorized_closure(&vakint, &settings(), &records, &prefactors).unwrap();
    let mut sum = Atom::Zero;
    for (record, prefactor) in records.iter().zip(&prefactors) {
        sum += evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(record),
            std::slice::from_ref(prefactor),
        )
        .unwrap()
        .value;
    }
    assert!((&all.value - &sum).expand().is_zero());
    assert_eq!(all.records.len(), 3);
}

fn refused(value: &FactorizedNumerator, what: &str) {
    let vakint = Vakint::new().unwrap();
    match evaluate_factorized_closure(
        &vakint,
        &settings(),
        std::slice::from_ref(value),
        &[Atom::num(1)],
    ) {
        Err(FactorizedError::Closure { detail }) => {
            assert!(detail.contains(what), "{}: {detail}", value.label);
            eprintln!("CLOSURE_REFUSAL {what}: {detail}");
        }
        other => panic!(
            "{}: expected a closure refusal ({what}), got {other:?}",
            value.label
        ),
    }
}

fn spec(value: &mut FactorizedNumerator) -> &mut ClosureSpec {
    let FactorizedKind::Closure(spec) = &mut value.kind else {
        unreachable!()
    };
    spec
}

#[test]
fn the_closure_kind_refuses_what_it_does_not_support() {
    // several massive lines
    let mut v = l2_root(2, true);
    v.edges[1].mass_squared = atom(R);
    refused(&v, "several massive lines");
    // a zero-loop record (OQ-25: the term without remaining loop)
    let mut v = record(
        "zero loops",
        0,
        vec![Some(20)],
        vec![
            line(3, 2, 1, &[1], R),
            leg(20, 2, true, &[1]),
            leg(21, 1, false, &[1]),
        ],
        vec![vertex(1, ""), vertex(2, "")],
        Vec::new(),
        None,
        3,
        taylor(0, vec![20, 21], false),
        false,
    );
    v.label = "zero loops".into();
    refused(&v, "at least one loop");
    // the closing edge massless / wrong mass / wrong power
    let mut v = l2_root(2, true);
    spec(&mut v).edge = 1;
    refused(&v, "mass");
    let mut v = l2_root(2, true);
    spec(&mut v).power = 2;
    refused(&v, "power");
    // Taylor degree above the limit
    let mut v = l2_root(2, true);
    spec(&mut v).taylor.as_mut().unwrap().degree = 3;
    refused(&v, "Taylor degree 3");
    // a scaled-edge set that does not scale every external momentum
    let mut v = l2_child(1);
    spec(&mut v).taylor.as_mut().unwrap().scaled_edges = vec![20];
    refused(&v, "must be scaled");
    // an internal edge among the scaled edges
    let mut v = l2_root(1, true);
    spec(&mut v).taylor.as_mut().unwrap().scaled_edges = vec![1, 20];
    refused(&v, "not an external leg");
    // a momentum in the projector with a Taylor step
    let mut v = l2_root(1, true);
    v.projector.as_mut().unwrap().tensor =
        atom(&format!("{}*{}", q(20, &hedge(120)), q(20, &hedge(121))));
    refused(&v, "momentum-free");
    // a pre-expanded sector whose internal edge still carries q
    let mut v = l2_root(1, true);
    spec(&mut v).taylor = None;
    spec(&mut v).auxiliary_mass_squared = atom(R);
    refused(&v, "carries an external momentum");
    // weights: with original_cycles the prefactor carries none
    let mut v = l2_root(1, true);
    v.prefactor.loop_sign = -1;
    refused(&v, "loop_sign=1");
    let mut v = l2_root(1, true);
    v.prefactor.scalars = atom(&format!("gs^2*{R}"));
    refused(&v, "reference mass");
    // an open bispinor port
    let mut v = l2_child(1);
    v.ports.push(atom(&format!(
        "spenso::{{spenso::representation,spenso::self_dual}}::bis(4,{HEDGE}(150,1))"
    )));
    match v.validate() {
        Err(FactorizedError::Closure { .. }) | Err(_) => {}
        Ok(()) => panic!("an open bispinor port validated"),
    }
    // not a closure record
    let mut v = l2_root(1, true);
    v.kind = FactorizedKind::Direct;
    refused(&v, "not a closure record");
    // an insertion on a vertex that is neither a contracted (ct:) nor a
    // surrogate (ins:) vertex is refused
    let mut v = l2_root(1, true);
    v.insertions.push(super::super::Insertion {
        vertex: 1,
        hedges: Vec::new(),
        cycles: (0, 0),
        coefficient: Atom::num(1),
        factor: Atom::num(1),
    });
    assert!(matches!(
        v.validate(),
        Err(FactorizedError::Insertion { .. })
    ));
    // the adapter and the vacuum entry keep refusing the kind
    let vakint = Vakint::new().unwrap();
    let v = l2_root(1, true);
    assert!(matches!(
        super::super::evaluate_factorized(&vakint, &settings(), &v, &Atom::num(1)),
        Err(FactorizedError::Reserved { .. })
    ));
    assert!(matches!(
        super::super::evaluate_factorized_vacuum(&vakint, &settings(), &v, &Atom::num(1)),
        Err(FactorizedError::Reserved { .. })
    ));
}

#[test]
fn the_v3_text_form_round_trips_every_fixture() {
    let mut child = l2_child(2);
    spec(&mut child).softs = vec![
        super::super::SoftLeg {
            external: 1,
            edge: 20,
            root: true,
        },
        super::super::SoftLeg {
            external: 2,
            edge: 21,
            root: false,
        },
    ];
    spec(&mut child).depth = Some(4);
    spec(&mut child).weights = vec![("scalars".into(), true, child.prefactor.scalars.clone())];
    for value in [l2_root(2, true), child, l2_sector(), l3_root(2), l4_root(2)] {
        value.validate().unwrap();
        let text = value.to_text().unwrap();
        assert!(text.starts_with("schema\tpint6.factorized-numerator.v3\n"));
        let back = FactorizedNumerator::from_text(&text)
            .unwrap_or_else(|e| panic!("{}: {e}\n{text}", value.label));
        assert_eq!(back, value, "{}", value.label);
        assert_eq!(back.to_text().unwrap(), text);
    }
    let text = l2_root(2, true).to_text().unwrap();
    // the beta's cut lines are read and ignored
    let with_cut = text.replace(
        "softs\t0\n",
        "routing\t2\nroute\t1\timage=1,0\nroute\t2\timage=0,1\ncut_topology\tvakint::topo(1)\nsofts\t0\n",
    );
    let lines = with_cut.lines().count() - 1;
    let with_cut = with_cut
        .lines()
        .map(|l| {
            if l.starts_with("lines\t") {
                format!("lines\t{lines}")
            } else {
                l.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    assert_eq!(
        FactorizedNumerator::from_text(&with_cut).unwrap(),
        l2_root(2, true)
    );
    // a wrong open-slot list, kind direct in v3, kind closure in v1
    let wrong = text.replace(
        "open\t0\n",
        &format!("open\t1\nopen_slot\t{}\n", hedge(120)),
    );
    assert!(FactorizedNumerator::parse_text(&wrong).is_err());
    // kind direct is a v3 kind too; with a closure line that is not `none`
    // the record is refused
    let direct = text.replace("kind\tclosure", "kind\tdirect");
    assert!(matches!(
        FactorizedNumerator::parse_text(&direct),
        Err(FactorizedError::TextSyntax { .. })
    ));
    let v1 = text.replace(
        "pint6.factorized-numerator.v3",
        "pint6.factorized-numerator.v1",
    );
    assert!(matches!(
        FactorizedNumerator::parse_text(&v1),
        Err(FactorizedError::Reserved { .. })
    ));
}

/// The beta thread's background-field B -> B DOT graphs (real Feynman rules:
/// gluon, ghost and quark lines, Dirac traces) as R* ROOT remainders: the
/// graph's last loop edge is the closing line (reference mass R), the Taylor
/// step of degree 2 with mass compensation, the metric root projector.
#[cfg(feature = "feynkit-ingress")]
fn dot_root_remainders(
    picks: &[(&str, usize)],
) -> Vec<(FactorizedNumerator, Vec<(String, String)>)> {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .unwrap_or_else(|_| "/common/ben/beta/results/prepared".into());
    let root = std::path::PathBuf::from(root);
    if !root.join("L1").join("model.json").is_file() {
        eprintln!("SKIP: no DOT data under {}", root.display());
        return Vec::new();
    }
    let mut out = Vec::new();
    for (level, n) in picks {
        let directory = root.join(level);
        let model = crate::refpath::load_model(&directory).unwrap();
        let mut files = std::fs::read_dir(&directory)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("diagram-") && n.ends_with(".dot"))
            })
            .collect::<Vec<_>>();
        files.sort_by_key(|p| (std::fs::metadata(p).unwrap().len(), p.clone()));
        out.extend(dot_files_as_roots(
            level,
            &model,
            files.into_iter().take(*n),
        ));
    }
    out
}

/// Named DOT graphs as root remainders ([`dot_root_remainders`]).
#[cfg(feature = "feynkit-ingress")]
fn dot_named_root_remainders(
    names: &[(&str, &str)],
) -> Vec<(FactorizedNumerator, Vec<(String, String)>)> {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .unwrap_or_else(|_| "/common/ben/beta/results/prepared".into());
    let root = std::path::PathBuf::from(root);
    let mut out = Vec::new();
    for (level, name) in names {
        let directory = root.join(level);
        let path = directory.join(name);
        if !path.is_file() {
            eprintln!("SKIP: {}", path.display());
            continue;
        }
        let model = crate::refpath::load_model(&directory).unwrap();
        out.extend(dot_files_as_roots(level, &model, [path]));
    }
    out
}

#[cfg(feature = "feynkit-ingress")]
fn dot_files_as_roots(
    level: &str,
    model: &std::sync::Arc<feynkit_model::Model>,
    files: impl IntoIterator<Item = std::path::PathBuf>,
) -> Vec<(FactorizedNumerator, Vec<(String, String)>)> {
    use super::super::FromFeynkitOptions;
    let model = std::sync::Arc::clone(model);
    let mut out = Vec::new();
    {
        for path in files {
            let text = std::fs::read_to_string(&path).unwrap();
            let diagram =
                feynkit_graph::FeynmanDiagram::from_dot(std::sync::Arc::clone(&model), &text)
                    .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            // the generator DOTs (prepared/L1..L4) carry the unevaluated weight;
            // the regrouped five-loop DOTs an evaluated one with the loop sign
            // inside (the weight is replaced by 1 below either way)
            let values =
                FactorizedNumerator::from_feynkit(&diagram, &model, &FromFeynkitOptions::default())
                    .or_else(|_| {
                        let options = FromFeynkitOptions {
                            loop_sign_in_overall: true,
                            ..FromFeynkitOptions::default()
                        };
                        FactorizedNumerator::from_feynkit_dot(&text, &model, &options)
                            .map(|(_, values)| values)
                    })
                    .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            let mut value = values.into_iter().next().expect("a colour branch");
            let loops = value.loops();
            let (i, closing) = value
                .loop_edges
                .iter()
                .enumerate()
                .rev()
                .find_map(|(i, e)| e.map(|e| (i, e)))
                .expect("a loop edge");
            let cycles = value.closed_loops().unwrap();
            let fermion = cycles
                .iter()
                .filter(|c| c.line == LineKind::Fermion)
                .count() as u32;
            let ghost = cycles.len() as u32 - fermion;
            let legs = value
                .edges
                .iter()
                .filter(|e| e.is_external())
                .map(|e| e.id)
                .collect::<Vec<_>>();
            for edge in value.edges.iter_mut() {
                if edge.id == closing {
                    edge.mass_squared = atom(R);
                }
            }
            value.prefactor.loop_sign = 1;
            value.prefactor.nf_power = 0;
            value.prefactor.overall = Atom::num(1);
            if let Some(projector) = value.projector.as_mut() {
                projector.p2_power = -1;
            }
            value.label = format!(
                "{level} {} closing edge {closing}",
                path.file_name().unwrap().to_string_lossy()
            );
            value.kind = FactorizedKind::Closure(ClosureSpec {
                edge: closing,
                mass_squared: atom(R),
                auxiliary_mass_squared: atom(MAUX),
                power: value.edge(closing).unwrap().power,
                independent_externals: false,
                taylor: taylor(2, legs, true),
                softs: Vec::new(),
                original_cycles: Some((fermion, ghost)),
                depth: None,
                weights: Vec::new(),
            });
            // the cut of the closing loop momentum k(i+1): k(i+1) -> p(1),
            // the later loop momenta shift down by one
            let routing = (0..loops)
                .map(|m| {
                    let image = match m.cmp(&i) {
                        std::cmp::Ordering::Less => format!("k({})", m + 1),
                        std::cmp::Ordering::Equal => "p(1)".to_string(),
                        std::cmp::Ordering::Greater => format!("k({})", m),
                    };
                    (format!("k({})", m + 1), image)
                })
                .collect::<Vec<_>>();
            out.push((value, routing));
        }
    }
    out
}

#[cfg(feature = "feynkit-ingress")]
#[test]
fn beta_dot_graphs_as_root_remainders_equal_the_contracted_reference() {
    for (value, routing) in dot_root_remainders(&[("L2", 3), ("L3", 2), ("L4", 1)]) {
        let routing = routing
            .iter()
            .map(|(a, b)| (a.as_str(), b.as_str()))
            .collect::<Vec<_>>();
        check(&value, &routing);
        let text = value.to_text().unwrap();
        assert_eq!(FactorizedNumerator::from_text(&text).unwrap(), value);
    }
}

/// Single-core timing, entry (uncontracted record) vs the contracted reference
/// chain, alternating, median of `VAKINT_CLOSURE_TIMING_REPS` (default 3).
/// Run explicitly (`--run-ignored only`), on one CPU.
#[cfg(feature = "feynkit-ingress")]
#[test]
#[ignore]
fn closure_kind_timing() {
    let reps = std::env::var("VAKINT_CLOSURE_TIMING_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3usize);
    let mut cases = vec![
        (
            l2_root(2, true),
            L2_CUT
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect::<Vec<_>>(),
        ),
        (
            l2_child(2),
            L2_CUT
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect(),
        ),
        (
            l3_root(2),
            L3_CUT
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect(),
        ),
        (
            l4_root(2),
            L4_CUT
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect(),
        ),
    ];
    cases.extend(dot_root_remainders(&[("L2", 3), ("L3", 2), ("L4", 1)]));
    // the beta's quick-turnaround L4 set (their Q7: 40, 42, 49, 100)
    cases.extend(dot_named_root_remainders(&[
        ("L4", "diagram-000040.dot"),
        ("L4", "diagram-000042.dot"),
        ("L4", "diagram-000049.dot"),
        ("L4", "diagram-000100.dot"),
    ]));
    // VAKINT_CLOSURE_TIMING_DOTS="L5/diagram-000051.dot,..." (under
    // VAKINT_TEST_FACTORIZED_DOT_ROOT) replaces the cases
    if let Ok(list) = std::env::var("VAKINT_CLOSURE_TIMING_DOTS") {
        let named = list
            .split(',')
            .filter_map(|item| item.split_once('/'))
            .collect::<Vec<_>>();
        cases = dot_named_root_remainders(&named);
    }
    let with_reference = std::env::var("VAKINT_CLOSURE_TIMING_NO_REFERENCE").is_err();
    let peak_kb = || {
        std::fs::read_to_string("/proc/self/status")
            .ok()
            .and_then(|s| {
                s.lines()
                    .find_map(|l| l.strip_prefix("VmHWM:"))
                    .and_then(|v| v.trim().trim_end_matches(" kB").trim().parse::<u64>().ok())
            })
            .unwrap_or(0)
    };
    let vakint = Vakint::new().unwrap();
    let median = |mut v: Vec<f64>| {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    let only = std::env::var("VAKINT_CLOSURE_TIMING_ONLY").ok();
    for (value, routing) in &cases {
        if only
            .as_ref()
            .is_some_and(|o| !value.label.contains(o.as_str()))
        {
            continue;
        }
        let routing = routing
            .iter()
            .map(|(a, b)| (a.as_str(), b.as_str()))
            .collect::<Vec<_>>();
        let (mut entry, mut reference_times) = (Vec::new(), Vec::new());
        let (mut contract, mut lower, mut reduce) = (Vec::new(), Vec::new(), Vec::new());
        let mut rustred = (0, 0, 0, 0, 0);
        for _ in 0..reps {
            let started = std::time::Instant::now();
            let got = evaluate_factorized_closure(
                &vakint,
                &settings(),
                std::slice::from_ref(value),
                &[value.prefactor.total()],
            )
            .unwrap();
            entry.push(started.elapsed().as_secs_f64());
            contract.push(got.records[0].contract_seconds);
            rustred = (
                got.records[0].patterns,
                got.records[0].blocks,
                got.records[0].key_plans,
                got.records[0].lowered_entries,
                got.records[0].distinct_monomials,
            );
            lower.push(got.records[0].lower_seconds);
            reduce.push(got.records[0].reduce_seconds);
            if with_reference {
                let started = std::time::Instant::now();
                let expected = reference(value, &routing);
                reference_times.push(started.elapsed().as_secs_f64());
                assert!(
                    (&got.records[0].vakint_value - &expected)
                        .expand()
                        .is_zero()
                );
            } else {
                reference_times.push(f64::NAN);
            }
        }
        eprintln!(
            "CLOSURE_TIMING label=\"{}\" peak_rss_kb={} patterns_blocks_keyplans_lowered_monomials={rustred:?} reps={reps} entry_median_s={:.4} reference_median_s={:.4} entry_contract_s={:.4} entry_lower_s={:.4} entry_reduce_s={:.4} entry_all={:?} reference_all={:?}",
            value.label,
            peak_kb(),
            median(entry.clone()),
            median(reference_times.clone()),
            median(contract),
            median(lower),
            median(reduce),
            entry,
            reference_times
        );
    }
}

/// The caller's settings of the beta's R* sites (gammalooprs
/// `VakintSettings::true_settings`: ε = `gammalooprs::ε`, μ_r² =
/// `gammalooprs::μᵣ²`, MSbar), RustRed only; the record's `depth` decides
/// the number of terms.
fn rstar_site_settings() -> VakintSettings {
    VakintSettings {
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        ..settings()
    }
}

fn fixture_value(text: &str) -> Atom {
    let line = text
        .lines()
        .find_map(|l| l.strip_prefix("value\t"))
        .expect("a value line");
    atom(line)
}

/// THE BETA'S REAL K3 RECORD: the gluonic IR hard part h_s of
/// B_self_energy_L3_036 (correction 0, term 0; fixed mass, degree 1, three
/// independent crown momenta, three open Lorentz ports, couplings and i
/// inside the factors), against the value their contracted path returns at
/// that call (`finish_tensor`), exactly.
#[test]
fn the_beta_k3_hs_record_equals_its_contracted_value() {
    let record = include_str!("../fixtures/wh-k3-hs.record.txt");
    let expected = fixture_value(include_str!("../fixtures/wh-k3-hs.value.txt"));
    let value = FactorizedNumerator::from_text(record).unwrap_or_else(|e| panic!("{e}"));
    value.validate().unwrap();
    assert_eq!(
        value.to_text().unwrap(),
        record,
        "byte-identical round trip"
    );
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&value),
        &[value.prefactor.total()],
    )
    .unwrap();
    eprintln!("K3_HS entry {}\nK3_HS beta  {}", got.value, expected);
    assert_equal_through(
        &got.value,
        &expected,
        valid_through(include_str!("../fixtures/wh-k3-hs.value.txt")),
        "k3-hs",
    );
}

/// The order in ε through which a fixture's value is valid (`valid_through
/// eps^n`).
fn valid_through(text: &str) -> i64 {
    let line = text
        .lines()
        .find_map(|l| l.strip_prefix("valid_through\t"))
        .expect("a valid_through line");
    let order = line
        .trim_start_matches("eps^")
        .split_whitespace()
        .next()
        .unwrap();
    order.parse().unwrap()
}

/// Exact equality coefficient by coefficient of the Laurent series in ε
/// (`gammalooprs::ε`) through order `order`, AND NO FURTHER (the fixtures'
/// identity rule: higher orders are not valid in their value).
fn assert_equal_through(got: &Atom, expected: &Atom, order: i64, what: &str) {
    let epsilon = symbolica::symbol!("gammalooprs::ε");
    let difference = (got - expected).expand();
    let mut checked = 0;
    if !difference.is_zero() {
        let series = difference
            .series(epsilon, 0, order as u32 as usize)
            .unwrap();
        for (power, coefficient) in series.terms() {
            let n = power.numerator().to_i64().unwrap();
            if n <= order {
                assert!(
                    coefficient.expand().is_zero(),
                    "{what}: coefficient of eps^{n} differs: {coefficient}\nentry {got}\nbeta  {expected}"
                );
                checked += 1;
            }
        }
    }
    eprintln!("{what}: equal through eps^{order} ({checked} nonzero-difference orders inspected)");
}

/// THE BETA'S REAL K3 DIRAC RECORD: the IR hard part of B_self_energy_L4_4966
/// with its open quark line closed by `gamma(bis, bis, nu)` (one fresh
/// Lorentz index left open, the 1/Tr[1] normalisation applied by the caller,
/// not in the record), fixed mass, degree 1, two independent crown momenta:
/// against their contracted value, exactly through eps^0 (its valid_through:
/// a (4 - 2 eps) trace factor multiplies part of their series, so eps^1 is
/// incomplete there) and NO FURTHER.
#[test]
fn the_beta_k3_dirac_record_equals_its_contracted_value() {
    let record = include_str!("../fixtures/wh-k3-dirac.record.txt");
    let value_text = include_str!("../fixtures/wh-k3-dirac.value.txt");
    let expected = fixture_value(value_text);
    let value = FactorizedNumerator::from_text(record).unwrap_or_else(|e| panic!("{e}"));
    value.validate().unwrap();
    assert_eq!(
        value.to_text().unwrap(),
        record,
        "byte-identical round trip"
    );
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&value),
        &[value.prefactor.total()],
    )
    .unwrap();
    assert_equal_through(&got.value, &expected, valid_through(value_text), "k3-dirac");
}

fn real_record(record: &str, value_text: &str, what: &str) -> ClosureRecordValue {
    let expected = fixture_value(value_text);
    let value = FactorizedNumerator::from_text(record).unwrap_or_else(|e| panic!("{what}: {e}"));
    value.validate().unwrap_or_else(|e| panic!("{what}: {e}"));
    assert_eq!(
        value.to_text().unwrap(),
        record,
        "{what}: byte-identical round trip"
    );
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&value),
        &[value.prefactor.total()],
    )
    .unwrap_or_else(|e| panic!("{what}: {e}"));
    assert_equal_through(&got.value, &expected, valid_through(value_text), what);
    // the caller's notation: no vakint-notation vector or scalar product survives
    let mut vakint_notation = None;
    got.value.as_view().visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && (f.get_symbol() == S.dot || f.get_symbol() == S.p || f.get_symbol() == S.k)
        {
            vakint_notation = Some(part.to_owned());
        }
        vakint_notation.is_none()
    });
    assert!(
        vakint_notation.is_none(),
        "{what}: vakint notation {vakint_notation:?} in the returned value"
    );
    got.records[0].clone()
}

/// THE BETA'S REAL K1 RECORD: the root remainder of B_self_energy_L2_22
/// (closed by the metric projector, p2_power -1, w = 2 compensated, ghost
/// loop): equal to `evaluate_boundary`'s value at that site through its
/// valid_through, in the caller's notation (p^2 = spenso::dot(Q(0),Q(0))).
#[test]
fn the_beta_k1_root_record_equals_its_contracted_value() {
    real_record(
        include_str!("../fixtures/wh-k1-root.record.txt"),
        include_str!("../fixtures/wh-k1-root.value.txt"),
        "k1-root",
    );
}

/// THE BETA'S REAL K2 RECORD: the isotropic closed UV child {2,3,5} of
/// B_self_energy_L2_22 (three independent crown momenta, w = 0, the channel
/// markers uv_isotropic_channel_1..3 carried through the projector).
#[test]
fn the_beta_k2_child_record_equals_its_contracted_value() {
    real_record(
        include_str!("../fixtures/wh-k2-child.record.txt"),
        include_str!("../fixtures/wh-k2-child.value.txt"),
        "k2-child",
    );
}

/// A unimodular change of the loop basis `k(i) -> k(i) + c p(j)`: the edge
/// momenta are data, the factors keep their edge momenta `Q(e)`.
fn shift_loop(value: &FactorizedNumerator, i: usize, j: usize, c: i64) -> FactorizedNumerator {
    let mut out = value.clone();
    let loops = out.loops();
    for edge in out.edges.iter_mut() {
        edge.momentum[loops + j] += c * edge.momentum[i];
    }
    out.loop_edges[i] = None;
    out.label = format!(
        "{} rerouted k({}) -> k({}) + {c} p({})",
        value.label,
        i + 1,
        i + 1,
        j + 1
    );
    out
}

/// C4: the closing line carrying the external momentum (r_c != 0, the
/// production root routing can do this). The same remainder rerouted so that
/// q flows through the closing edge gives the SAME value, and equals the
/// contracted reference chain (which expands the closing line's l.q terms);
/// also on the beta's real K1 record.
#[test]
fn a_closing_line_carrying_the_external_momentum_gives_the_same_value() {
    let vakint = Vakint::new().unwrap();
    let evaluate = |value: &FactorizedNumerator| {
        value
            .validate()
            .unwrap_or_else(|e| panic!("{}: {e}", value.label));
        evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(value),
            &[value.prefactor.total()],
        )
        .unwrap_or_else(|e| panic!("{}: {e}", value.label))
        .value
    };
    for (degree, compensate) in [(1, true), (2, true), (2, false)] {
        let value = l2_root(degree, compensate);
        // the closing edge 3 is k(2): k(2) -> k(2) + p(1) sends q through it
        let rerouted = shift_loop(&value, 1, 0, 1);
        assert!(
            rerouted
                .edges
                .iter()
                .any(|e| e.id == 3 && e.momentum[2] != 0)
        );
        let a = evaluate(&value);
        let b = evaluate(&rerouted);
        assert!(
            (&a - &b).expand().is_zero(),
            "{}:\n{a}\n{b}",
            rerouted.label
        );
        check(&rerouted, &L2_CUT);
    }
    let value = l2_child(2);
    let rerouted = shift_loop(&value, 1, 1, -1);
    assert!(
        (&evaluate(&value) - &evaluate(&rerouted))
            .expand()
            .is_zero()
    );
    check(&rerouted, &L2_CUT);
    // the beta's K1 record, its closing edge 2 = k(1) rerouted onto k(1) + p(1)
    let text = include_str!("../fixtures/wh-k1-root.value.txt");
    let record =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k1-root.record.txt")).unwrap();
    let rerouted = shift_loop(&record, 0, 0, 1);
    assert_equal_through(
        &evaluate(&rerouted),
        &fixture_value(text),
        valid_through(text),
        "k1-root rerouted",
    );
}

/// An open-port record keeps every order: the beta's k3-hs hard part (three
/// open Lorentz ports, omega 1) is entirely of ODD order 1 in its crown
/// momenta: at degree 0 the entry gives exactly 0, at degree 1 the beta's
/// nonzero value (the boundary path's odd-order drop is only for records
/// without an open slot).
#[test]
fn an_open_port_record_keeps_its_odd_order() {
    let vakint = Vakint::new().unwrap();
    let record =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k3-hs.record.txt")).unwrap();
    let at = |degree: u32| {
        let mut value = record.clone();
        spec(&mut value).taylor.as_mut().unwrap().degree = degree;
        evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(&value),
            &[value.prefactor.total()],
        )
        .unwrap()
        .value
    };
    assert!(at(0).is_zero(), "order 0 of k3-hs is not zero");
    let text = include_str!("../fixtures/wh-k3-hs.value.txt");
    let odd = at(1);
    assert!(!odd.is_zero());
    assert_equal_through(
        &odd,
        &fixture_value(text),
        valid_through(text),
        "k3-hs order 1",
    );
}

/// THE BETA'S VERIFIED C4 FIXTURE: the K1 root of B_self_energy_L2_21 under
/// their production root routing (num-est), whose closing edge 2 carries
/// r_c = +p(1), and its twin under the default routing (r_c = 0): each equals
/// its value.txt, and the two values are identical.
#[test]
fn the_beta_k1_root_with_q_through_the_closing_line_equals_its_twin() {
    let rc = real_record(
        include_str!("../fixtures/wh-k1-root-rc.record.txt"),
        include_str!("../fixtures/wh-k1-root-rc.value.txt"),
        "k1-root-rc (r_c = +p(1))",
    );
    // the control's value file has no valid_through line: the same remainder,
    // compared through the r_c != 0 fixture's (eps^-1)
    let order = valid_through(include_str!("../fixtures/wh-k1-root-rc.value.txt"));
    let control =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k1-root-rc0.record.txt"))
            .unwrap();
    control.validate().unwrap();
    assert!(
        control.edge(2).unwrap().momentum[2] == 0,
        "the control's closing edge carries no q"
    );
    let vakint = Vakint::new().unwrap();
    let rc0 = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&control),
        &[control.prefactor.total()],
    )
    .unwrap()
    .records[0]
        .clone();
    assert_equal_through(
        &rc0.value,
        &fixture_value(include_str!("../fixtures/wh-k1-root-rc0.value.txt")),
        order,
        "k1-root-rc control (r_c = 0)",
    );
    let record =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k1-root-rc.record.txt"))
            .unwrap();
    assert!(
        record.edge(2).unwrap().momentum[2] != 0,
        "the fixture's closing edge carries q"
    );
    assert!(
        (&rc.value - &rc0.value).expand().is_zero(),
        "r_c != 0 {}\nr_c = 0 {}",
        rc.value,
        rc0.value
    );
}

/// The value does not depend on the caller's `tensor_reduction_method` or
/// `use_dot_product_notation` (format v1.4 section 7): k1-root and k3-hs with
/// each flipped give the same value.
#[test]
fn the_value_does_not_depend_on_the_reduction_method_or_the_notation() {
    let vakint = Vakint::new().unwrap();
    for (record, what) in [
        (include_str!("../fixtures/wh-k1-root.record.txt"), "k1-root"),
        (include_str!("../fixtures/wh-k3-hs.record.txt"), "k3-hs"),
    ] {
        let value = FactorizedNumerator::from_text(record).unwrap();
        let run = |method: crate::TensorReductionMethod, dots: bool| {
            let mut settings = rstar_site_settings();
            settings.tensor_reduction_method = method;
            settings.use_dot_product_notation = dots;
            evaluate_factorized_closure(
                &vakint,
                &settings,
                std::slice::from_ref(&value),
                &[value.prefactor.total()],
            )
            .unwrap()
            .value
        };
        let reference = run(crate::TensorReductionMethod::FeynKit, true);
        for (method, dots) in [
            (crate::TensorReductionMethod::FeynKit, false),
            (crate::TensorReductionMethod::AlphaLoop, true),
            (crate::TensorReductionMethod::AlphaLoop, false),
        ] {
            let other = run(method, dots);
            assert!(
                (&other - &reference).expand().is_zero(),
                "{what}: {method:?} dots={dots}: {other}\nvs {reference}"
            );
        }
    }
}

/// zeta(n) is a legal scalar inside a factor (the beta's K5 coefficient
/// vertices carry symbolica::zeta(3)): accepted, and it multiplies the value.
#[test]
fn zeta_values_are_legal_scalars_in_factors() {
    let record = include_str!("../fixtures/wh-k3-hs.record.txt");
    let text = include_str!("../fixtures/wh-k3-hs.value.txt");
    let mut value = FactorizedNumerator::from_text(record).unwrap();
    let zeta = atom("symbolica::{}::zeta(3)");
    value.vertices[0].factor = &value.vertices[0].factor * &zeta;
    value.validate().unwrap();
    let back = FactorizedNumerator::from_text(&value.to_text().unwrap()).unwrap();
    assert_eq!(back, value);
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&value),
        &[value.prefactor.total()],
    )
    .unwrap();
    assert_equal_through(
        &got.value,
        &(fixture_value(text) * &zeta),
        valid_through(text),
        "k3-hs x zeta(3)",
    );
}

/// Both projector forms work: K1 emits p2_power = -1 (the 1/p^2 applied by
/// the entry), the K5 surrogate root p2_power = 0 (applied by the caller):
/// k1-root at p2_power 0, times 1/p^2 in the caller's notation, is its value.
#[test]
fn p2_power_zero_and_minus_one_projectors_agree() {
    let text = include_str!("../fixtures/wh-k1-root.value.txt");
    let mut value =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k1-root.record.txt")).unwrap();
    value.projector.as_mut().unwrap().p2_power = 0;
    value.validate().unwrap();
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&value),
        &[value.prefactor.total()],
    )
    .unwrap();
    let q = format!("{MOMENTUM}(0,{MINK}({DIM}))");
    let p2 = atom(&format!(
        "spenso::{{symmetric,linear,scalar}}::dot({q},{q})"
    ));
    assert_equal_through(
        &(got.value / p2),
        &fixture_value(text),
        valid_through(text),
        "k1-root at p2_power 0, / p^2",
    );
}

/// Contraction probe (ignored; run explicitly): the entry on the v3 record in
/// `VAKINT_CLOSURE_PROBE_RECORD`, and, when `VAKINT_CLOSURE_PROBE_ROUTED`
/// names a file whose first line `routed_numerator<TAB><atom>` (or a bare atom)
/// is the producer's routed product, the SAME idenso chain on that product;
/// both timed, `VAKINT_CLOSURE_TIMING_REPS` repetitions alternating.
#[test]
#[ignore]
fn closure_contraction_probe() {
    let Ok(path) = std::env::var("VAKINT_CLOSURE_PROBE_RECORD") else {
        eprintln!("SKIP: VAKINT_CLOSURE_PROBE_RECORD unset");
        return;
    };
    let reps = std::env::var("VAKINT_CLOSURE_TIMING_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3usize);
    let value = FactorizedNumerator::from_text(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let theirs = std::env::var("VAKINT_CLOSURE_PROBE_ROUTED").ok().map(|p| {
        let text = std::fs::read_to_string(p).unwrap();
        let line = text
            .lines()
            .find_map(|l| l.strip_prefix("routed_numerator\t"))
            .unwrap_or_else(|| text.lines().next().unwrap());
        atom(line)
    });
    let vakint = Vakint::new().unwrap();
    if let Some(routed) = theirs.as_ref()
        && std::env::var("VAKINT_CLOSURE_PROBE_THEIRS_FIRST").is_ok()
    {
        let started = std::time::Instant::now();
        let contracted =
            super::super::oracle::contract_with_idenso_unexpanded(routed, false).unwrap();
        eprintln!(
            "CLOSURE_PROBE_FIRST their_product_contract_s={:.4} terms={}",
            started.elapsed().as_secs_f64(),
            crate::census::summands(contracted.as_view())
        );
    }
    for _ in 0..reps {
        let got = evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(&value),
            &[value.prefactor.total()],
        )
        .unwrap();
        let r = &got.records[0];
        let their_seconds = theirs.as_ref().map(|routed| {
            let started = std::time::Instant::now();
            let contracted =
                super::super::oracle::contract_with_idenso_unexpanded(routed, false).unwrap();
            (
                started.elapsed().as_secs_f64(),
                crate::census::summands(contracted.as_view()),
                contracted.as_view().get_byte_size(),
                routed.as_view().get_byte_size(),
            )
        });
        eprintln!(
            "CLOSURE_PROBE label=\"{}\" entry_contract_s={:.4} entry_lower_s={:.4} entry_reduce_s={:.4} patterns={} key_plans={} lowered={} their_product_contract_s_terms_bytes_inbytes={their_seconds:?}",
            value.label,
            r.contract_seconds,
            r.lower_seconds,
            r.reduce_seconds,
            r.patterns,
            r.key_plans,
            r.lowered_entries
        );
    }
}

/// The beta's k1-quick set (every root remainder of the L4 graphs 40, 42, 49,
/// 100 that reaches their boundary path, as v3 records with the value and the
/// RR_CONTRACT seconds of that call; 30 of the 35 with child insertions):
/// each record through the entry, its value equal to value.txt through
/// valid_through, and its contraction seconds printed beside theirs. Data: `VAKINT_TEST_WH_K1_QUICK` (default the beta's directory);
/// SKIP without it.
#[test]
fn the_beta_k1_quick_records_equal_their_values() {
    let root = std::env::var("VAKINT_TEST_WH_K1_QUICK")
        .unwrap_or_else(|_| "/common/ben/beta/results/consolidation/wh-fixtures/k1-quick".into());
    let root = std::path::PathBuf::from(root);
    if !root.is_dir() {
        eprintln!("SKIP: no k1-quick data under {}", root.display());
        return;
    }
    let mut records = Vec::new();
    for graph in std::fs::read_dir(&root).unwrap().map(|e| e.unwrap().path()) {
        if let Ok(entries) = std::fs::read_dir(&graph) {
            for record in entries.map(|e| e.unwrap().path()) {
                if record.join("input.txt").is_file() && record.join("value.txt").is_file() {
                    records.push(record);
                }
            }
        }
    }
    records.sort();
    let only = std::env::var("VAKINT_TEST_WH_K1_ONLY").ok();
    let vakint = Vakint::new().unwrap();
    let (mut mine, mut theirs, mut checked) = (0.0, 0.0, 0);
    for record in &records {
        if only
            .as_ref()
            .is_some_and(|o| !record.display().to_string().ends_with(o.as_str()))
        {
            continue;
        }
        let input = std::fs::read_to_string(record.join("input.txt")).unwrap();
        let value_text = std::fs::read_to_string(record.join("value.txt")).unwrap();
        let value = match FactorizedNumerator::from_text(&input) {
            Ok(value) => value,
            Err(e) => panic!("{}: {e}", record.display()),
        };
        let meta = std::fs::read_to_string(record.join("meta.json")).unwrap_or_default();
        let their_seconds = meta
            .split("\"rr_contract_seconds\": ")
            .nth(1)
            .and_then(|t| t.split([',', '\n']).next())
            .and_then(|t| t.trim().parse::<f64>().ok())
            .unwrap_or(f64::NAN);
        let got = evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(&value),
            &[value.prefactor.total()],
        )
        .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        let what = record.display().to_string();
        assert_equal_through(
            &got.value,
            &fixture_value(&value_text),
            valid_through(&value_text),
            &what,
        );
        let r = &got.records[0];
        eprintln!(
            "K1_QUICK {what} insertions={} entry_contract_s={:.4} their_rr_contract_s={their_seconds:.4} ratio={:.2} split_s={:.4} evaluate_s={:.4} peak_rss_kb={} patterns={} blocks={} key_plans={} lowered={} lower_s={:.4} reduce_s={:.4}",
            value.insertions.len(),
            r.contract_seconds,
            r.contract_seconds / their_seconds,
            r.split_seconds,
            r.evaluate_seconds,
            std::fs::read_to_string("/proc/self/status")
                .ok()
                .and_then(|s| s
                    .lines()
                    .find_map(|l| l.strip_prefix("VmHWM:"))
                    .map(|v| v.trim().to_string()))
                .unwrap_or_default(),
            r.patterns,
            r.blocks,
            r.key_plans,
            r.lowered_entries,
            r.lower_seconds,
            r.reduce_seconds
        );
        mine += r.contract_seconds;
        theirs += their_seconds;
        checked += 1;
    }
    eprintln!(
        "K1_QUICK SUM records={checked} entry_contract_s={mine:.4} their_rr_contract_s={theirs:.4}"
    );
}

/// THE BETA'S REAL RECORDS WITH INSERTIONS (closure kind + insertion
/// vertices), through `real_record` (exact through `valid_through`, caller
/// notation, byte-identical round trip), with the summed pole order printed.
fn real_insertion_record(record: &str, value_text: &str, what: &str) {
    let value = FactorizedNumerator::from_text(record).unwrap();
    let epsilon = crate::utils::vakint_macros::vk_symbol!("gammalooprs::ε");
    let started = std::time::Instant::now();
    real_record(record, value_text, what);
    eprintln!(
        "{what}: insertions {} pole order {} seconds {:.3}",
        value.insertions.len(),
        value.insertion_pole_order(epsilon).unwrap(),
        started.elapsed().as_secs_f64()
    );
}

#[test]
fn the_beta_k2_child_with_an_insertion_equals_its_contracted_value() {
    real_insertion_record(
        include_str!("../fixtures/wh-k2-child-ins.record.txt"),
        include_str!("../fixtures/wh-k2-child-ins.value.txt"),
        "k2-child-ins",
    );
}

#[test]
fn the_beta_k5_surrogate_root_equals_its_contracted_value() {
    real_insertion_record(
        include_str!("../fixtures/wh-k5-surrogate.record.txt"),
        include_str!("../fixtures/wh-k5-surrogate.value.txt"),
        "k5-surrogate",
    );
}

#[test]
fn the_beta_k4_quotient_equals_its_contracted_value() {
    real_insertion_record(
        include_str!("../fixtures/wh-k4-quotient.record.txt"),
        include_str!("../fixtures/wh-k4-quotient.value.txt"),
        "k4-quotient",
    );
}

/// A child insertion of momentum degree 2 in a record with Taylor degree
/// w = 2 goes through the same routing, lambda scaling and truncation as a
/// vertex factor: the record with the polynomial on a contracted (ct:)
/// vertex equals the record with it as that vertex's own factor.
#[test]
fn a_degree_two_insertion_is_routed_and_truncated_like_a_factor() {
    let mut base = l2_root(2, true);
    // the port of vertex 2 moves to its leg; the vertex keeps a degree-2
    // scalar polynomial of its boundary momenta
    let compact = format!("{MINK}({DIM})");
    base.edges.iter_mut().find(|e| e.id == 21).unwrap().factor = atom(&q(2, &hedge(121)));
    let polynomial = atom(&format!(
        "({DIM}-2)/3*spenso::dot({},{})+spenso::dot({},{})",
        q(1, &compact),
        q(2, &compact),
        q(3, &compact),
        q(3, &compact)
    ));
    base.vertices.iter_mut().find(|v| v.id == 2).unwrap().factor = polynomial.clone();
    base.validate().unwrap();
    let mut with = base.clone();
    let hedges = with
        .local_hedges(super::super::Location::Vertex(2))
        .into_iter()
        .collect::<Vec<_>>();
    let v = with.vertices.iter_mut().find(|v| v.id == 2).unwrap();
    v.factor = Atom::num(1);
    v.rule = "ct:child".into();
    with.insertions.push(super::super::Insertion {
        vertex: 2,
        hedges,
        cycles: (0, 0),
        coefficient: Atom::num(1),
        factor: polynomial,
    });
    with.validate().unwrap();
    let vakint = Vakint::new().unwrap();
    let settings = rstar_site_settings();
    let run = |v: &FactorizedNumerator| {
        evaluate_factorized_closure(
            &vakint,
            &settings,
            std::slice::from_ref(v),
            &[v.prefactor.total()],
        )
        .unwrap()
        .value
    };
    let (a, b) = (run(&with), run(&base));
    assert!(!b.is_zero());
    assert!((&a - &b).expand().is_zero(), "insertion {a}\nvertex {b}");
}

/// The beta's k5-stress set (L5 graph 3793, closed surrogate root, five
/// records with four insertion vertices each; r0051 has five loops and needs
/// the 1-4 loop package of `VAKINT_TEST_P_LIBRARY_1_4`): each record equal to
/// its value through valid_through, with the entry's seconds per stage.
/// Data: `VAKINT_TEST_WH_K5_STRESS` (default the beta's directory); SKIP
/// without it.
#[test]
fn the_beta_k5_stress_records_equal_their_values() {
    let root = std::env::var("VAKINT_TEST_WH_K5_STRESS")
        .unwrap_or_else(|_| "/common/ben/beta/results/consolidation/wh-fixtures/k5-stress".into());
    let root = std::path::PathBuf::from(root);
    if !root.is_dir() {
        eprintln!("SKIP: no k5-stress data under {}", root.display());
        return;
    }
    let mut vakint = Vakint::new().unwrap();
    let package = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4");
    if let Some(directory) = &package {
        let library = rustred_app::propagator::load_p_library(
            rustred_app::propagator::PLibrarySource::Directory(directory.into()),
            &rustred_app::propagator::PLibraryOptions {
                verify_descent: false,
            },
        )
        .unwrap();
        vakint.inject_p_library(library).unwrap();
    }
    let mut records = std::fs::read_dir(&root)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.join("input.txt").is_file())
        .collect::<Vec<_>>();
    records.sort();
    let epsilon = crate::utils::vakint_macros::vk_symbol!("gammalooprs::ε");
    for record in &records {
        let input = std::fs::read_to_string(record.join("input.txt")).unwrap();
        let value_text = std::fs::read_to_string(record.join("value.txt")).unwrap();
        let value = FactorizedNumerator::from_text(&input)
            .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        if value.loops() > 4 && package.is_none() {
            eprintln!(
                "K5_STRESS {} SKIP: five loops need VAKINT_TEST_P_LIBRARY_1_4",
                record.display()
            );
            continue;
        }
        let started = std::time::Instant::now();
        let got = evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(&value),
            &[value.prefactor.total()],
        )
        .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        let total = started.elapsed().as_secs_f64();
        let what = record.display().to_string();
        assert_equal_through(
            &got.value,
            &fixture_value(&value_text),
            valid_through(&value_text),
            &what,
        );
        let r = &got.records[0];
        eprintln!(
            "K5_STRESS {what} loops={} insertions={} pole_order={} total_s={total:.3} contract_s={:.3} split_s={:.3} lower_s={:.3} reduce_s={:.3} evaluate_s={:.3} patterns={}",
            value.loops(),
            value.insertions.len(),
            value.insertion_pole_order(epsilon).unwrap(),
            r.contract_seconds,
            r.split_seconds,
            r.lower_seconds,
            r.reduce_seconds,
            r.evaluate_seconds,
            r.patterns
        );
    }
}

/// Two contracted children with a pole each (ε^-1 and 3 ε^-1) on the L2 root
/// remainder: the entry raises the integral's depth by the SUM of their pole
/// orders (2). The contracted route: the record without the insertions (their
/// factors left on the vertices) evaluated two terms deeper, times the product
/// of the coefficients, truncated at eps^(depth - L - 1), as the beta's
/// integrate.rs does. With the maximum (1) the last kept order would differ.
#[test]
fn two_poled_insertions_raise_the_depth_by_their_sum() {
    let eps = "gammalooprs::ε";
    let mut base = l2_root(1, true);
    base.vertices.iter_mut().find(|v| v.id == 1).unwrap().factor = atom(&q(1, &hedge(120)));
    let mut with = base.clone();
    let coefficients = [atom(&format!("{eps}^-1")), atom(&format!("3*{eps}^-1+1/2"))];
    for (vertex, coefficient) in [1usize, 2].into_iter().zip(&coefficients) {
        let hedges = with
            .local_hedges(super::super::Location::Vertex(vertex))
            .into_iter()
            .collect::<Vec<_>>();
        let v = with.vertices.iter_mut().find(|v| v.id == vertex).unwrap();
        let factor = std::mem::replace(&mut v.factor, Atom::num(1));
        v.rule = format!("ct:child{vertex}");
        with.insertions.push(super::super::Insertion {
            vertex,
            hedges,
            cycles: (0, 0),
            coefficient: coefficient.clone(),
            factor,
        });
    }
    with.validate().unwrap();
    let epsilon = crate::utils::vakint_macros::vk_symbol!(eps);
    assert_eq!(with.insertion_pole_order(epsilon).unwrap(), 2);
    let vakint = Vakint::new().unwrap();
    let settings = rstar_site_settings();
    let terms = settings.number_of_terms_in_epsilon_expansion;
    let got = evaluate_factorized_closure(
        &vakint,
        &settings,
        std::slice::from_ref(&with),
        &[with.prefactor.total()],
    )
    .unwrap();
    let mut deeper = settings.clone();
    deeper.number_of_terms_in_epsilon_expansion = terms + 2;
    let bare = evaluate_factorized_closure(
        &vakint,
        &deeper,
        std::slice::from_ref(&base),
        &[base.prefactor.total()],
    )
    .unwrap();
    let target = terms - base.loops() as i64 - 1;
    let reference = (&coefficients[0] * &coefficients[1] * &bare.value)
        .series(
            epsilon,
            Atom::Zero.as_view(),
            symbolica::domains::rational::Rational::from(target),
        )
        .unwrap()
        .to_atom()
        .expand();
    assert!(!reference.is_zero());
    assert!(
        (&got.value - &reference).expand().is_zero(),
        "insertions {}\ncontracted {}",
        got.value,
        reference
    );
    eprintln!("TWO_POLED identical through eps^{target}: {}", got.value);
}

/// Weights with insertions: carried cycles + the insertions' cycles may be
/// FEWER than `original_cycles` (a scope misses a loop of the diagram), never
/// more.
#[test]
fn insertion_cycles_count_towards_original_cycles_at_most() {
    let record = include_str!("../fixtures/wh-k2-child-ins.record.txt");
    let value = FactorizedNumerator::from_text(record).unwrap();
    let mut fewer = value.clone();
    if let FactorizedKind::Closure(spec) = &mut fewer.kind {
        spec.original_cycles = Some((1, 1));
    }
    fewer.insertions[0].cycles = (1, 0);
    fewer.validate().unwrap();
    let mut more = fewer.clone();
    more.insertions[0].cycles = (2, 0);
    assert!(
        more.validate().is_err(),
        "the insertions' cycles exceed original_cycles"
    );
}

/// The beta's k1-l5 set: K1 records of every boundary root remainder of a
/// FIVE-LOOP diagram (most with child insertions), each against value.txt
/// through valid_through, its seconds printed beside the RR_CONTRACT of
/// meta.json. Heavy (five-loop reductions): runs only when
/// `VAKINT_TEST_WH_K1_L5` names a graph directory (e.g. .../k1-l5/L5-60), with
/// the 1-4 loop package of `VAKINT_TEST_P_LIBRARY_1_4`.
#[test]
fn the_beta_k1_l5_records_equal_their_values() {
    let Some(root) = std::env::var_os("VAKINT_TEST_WH_K1_L5") else {
        eprintln!("SKIP: VAKINT_TEST_WH_K1_L5 not set");
        return;
    };
    let Some(package) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") else {
        eprintln!("SKIP: five-loop records need VAKINT_TEST_P_LIBRARY_1_4");
        return;
    };
    let mut vakint = Vakint::new().unwrap();
    let library = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(package.into()),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap();
    vakint.inject_p_library(library).unwrap();
    let mut records = std::fs::read_dir(std::path::PathBuf::from(root))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.join("input.txt").is_file() && p.join("value.txt").is_file())
        .collect::<Vec<_>>();
    records.sort();
    let (mut mine, mut theirs) = (0.0, 0.0);
    for record in &records {
        let input = std::fs::read_to_string(record.join("input.txt")).unwrap();
        let value_text = std::fs::read_to_string(record.join("value.txt")).unwrap();
        let value = FactorizedNumerator::from_text(&input)
            .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        let meta = std::fs::read_to_string(record.join("meta.json")).unwrap_or_default();
        let their_seconds = meta
            .split("\"rr_contract_seconds\": ")
            .nth(1)
            .and_then(|t| t.split([',', '\n', '}']).next())
            .and_then(|t| t.trim().parse::<f64>().ok())
            .unwrap_or(f64::NAN);
        let started = std::time::Instant::now();
        let got = evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(&value),
            &[value.prefactor.total()],
        )
        .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        let total = started.elapsed().as_secs_f64();
        let what = record.display().to_string();
        assert_equal_through(
            &got.value,
            &fixture_value(&value_text),
            valid_through(&value_text),
            &what,
        );
        let r = &got.records[0];
        eprintln!(
            "K1_L5 {what} loops={} insertions={} total_s={total:.3} contract_s={:.3} their_rr_contract_s={their_seconds:.3} lower_s={:.3} reduce_s={:.3} patterns={}",
            value.loops(),
            value.insertions.len(),
            r.contract_seconds,
            r.lower_seconds,
            r.reduce_seconds,
            r.patterns
        );
        mine += r.contract_seconds;
        theirs += their_seconds;
    }
    eprintln!(
        "K1_L5 SUM records={} entry_contract_s={mine:.3} their_rr_contract_s={theirs:.3}",
        records.len()
    );
}
