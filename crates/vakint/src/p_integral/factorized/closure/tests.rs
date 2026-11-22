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

/// THE shared fixture helper (load or fail, compare through valid_through
/// with a signed bound and a nonzero reference, count): one file for the unit
/// and the integration tests.
#[path = "../../../../tests/support/wh_fixtures.rs"]
mod wh_fixtures;

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
            scalars: atom(&format!("(-1+{DIM})^(-1)*UFO::gs^2")),
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

/// The vakint vacuum notation of an idenso-contracted network, by a
/// substitution written HERE (not `map_to_vakint_vacuum`): a scalar product
/// of two compact momenta `Q(a, mink(D))` -> `dot(v_a, v_b)`, a momentum on
/// port i of `value.ports` (record order) -> `v_e` with the raw index
/// 9001 + i, a metric of two ports -> `g(9001 + i, 9001 + j)`; `v_e` from
/// the edge's momentum vector (`k(l)` loops, `p(j)` externals). Any other
/// surviving momentum or slot is a failure.
fn independent_vacuum_map(value: &FactorizedNumerator, contracted: &Atom) -> Atom {
    let head_of = |text: &str| {
        let a = atom(text);
        let AtomView::Fun(f) = a.as_view() else {
            panic!("{a}")
        };
        f.get_symbol()
    };
    let q_head = head_of(&q(1, &hedge(1)));
    let metric = head_of(&format!("{METRIC}({},{})", hedge(1), hedge(2)));
    let ports = value.ports.clone();
    let loops = value.loops();
    let vector = |id: i64, index: Option<i64>| -> Atom {
        let edge = value
            .edges
            .iter()
            .find(|e| e.id as i64 == id)
            .unwrap_or_else(|| panic!("unknown edge {id}"));
        let mut out = Atom::Zero;
        for (i, c) in edge.momentum.iter().enumerate() {
            if *c == 0 {
                continue;
            }
            let (h, n) = if i < loops {
                (S.k, i as i64 + 1)
            } else {
                (S.p, (i - loops) as i64 + 1)
            };
            out += Atom::num(*c)
                * match index {
                    Some(x) => function!(h, Atom::num(n), Atom::num(x)),
                    None => function!(h, Atom::num(n)),
                };
        }
        out
    };
    // Q(e, rep): (e, rep)
    let momentum = |v: AtomView<'_>| -> Option<(i64, Atom)> {
        let AtomView::Fun(f) = v else { return None };
        if f.get_symbol() != q_head || f.get_nargs() != 2 {
            return None;
        }
        let args = f.iter().collect::<Vec<_>>();
        Some((i64::try_from(args[0]).ok()?, args[1].to_owned()))
    };
    let compact = |rep: &Atom| matches!(rep.as_view(), AtomView::Fun(f) if f.get_nargs() == 1);
    let port = |slot: &Atom| {
        ports
            .iter()
            .position(|p| p == slot)
            .map(|i| 9001 + i as i64)
    };
    let mapped = contracted.replace_map(|part, _, out| {
        let AtomView::Fun(f) = part else { return };
        let args = f.iter().collect::<Vec<_>>();
        let is_product = f.get_symbol() == metric || f.get_symbol().get_name() == "spenso::dot";
        if is_product && args.len() == 2 {
            if let (Some((a, ra)), Some((b, rb))) = (momentum(args[0]), momentum(args[1]))
                && compact(&ra)
                && compact(&rb)
            {
                **out = S.dot.call((vector(a, None), vector(b, None)));
                return;
            }
            if f.get_symbol() == metric
                && let (Some(i), Some(j)) = (port(&args[0].to_owned()), port(&args[1].to_owned()))
            {
                **out = S.g.call((Atom::num(i), Atom::num(j)));
            }
            return;
        }
        if let Some((e, slot)) = momentum(part)
            && let Some(i) = port(&slot)
        {
            **out = vector(e, Some(i));
        }
    });
    let mut leftover = None;
    mapped.as_view().visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && (f.get_symbol() == q_head || f.get_symbol() == metric)
        {
            leftover = Some(part.to_owned());
        }
        leftover.is_none()
    });
    assert!(leftover.is_none(), "left after the mapping: {leftover:?}");
    mapped
}

/// THE INDEPENDENT REFERENCE: the contracted reference chain
/// sharing NONE of `TaylorOrders::keeps`, `map_to_vakint_vacuum` and
/// `Prefactor::total` with the entry (orders written out, the mapping
/// [`independent_vacuum_map`], the prefactor ARGUMENT multiplied in by plain
/// multiplication AFTER the evaluation; exact for an ε-free argument).
fn independent_reference(
    value: &FactorizedNumerator,
    routing: &[(&str, &str)],
    argument: &Atom,
) -> Atom {
    independent_reference_with(value, routing, argument, settings())
}

/// [`independent_reference`] with the given settings (another depth).
fn independent_reference_with(
    value: &FactorizedNumerator,
    routing: &[(&str, &str)],
    argument: &Atom,
    settings: VakintSettings,
) -> Atom {
    let vakint = Vakint::new().unwrap();
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
    let mut scalar = independent_vacuum_map(value, &contracted)
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
        // the kept orders written out (not TaylorOrders::keeps)
        let kept_orders: Vec<i64> = match taylor.orders {
            TaylorOrders::All => (0..=w as i64).collect(),
            TaylorOrders::DropZero => (1..=w as i64).collect(),
            TaylorOrders::DropZeroOdd => (2..=w as i64).filter(|n| n % 2 == 0).collect(),
        };
        for (power, coefficient) in series.terms() {
            let n = power.numerator().to_i64().unwrap();
            assert!(power.is_integer() && n >= 0);
            if kept_orders.contains(&n) {
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
    // NO prefactor inside: the argument multiplies the result afterwards
    let prefactor = Atom::num(1);
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
    let argument = argument
        .replace(Atom::var(value.dimension))
        .with(in_epsilon.clone());
    (total * argument).expand()
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

/// The Taylor order filters other than `orders=all` (the top-mass drop
/// modes) are refused, by validation and by the entry, until a verified
/// record exists; a record carrying one still parses. `orders=all` on the
/// same record evaluates and equals the reference.
#[test]
fn l2_root_orders_filters_are_refused() {
    let vakint = Vakint::new().unwrap();
    for orders in [TaylorOrders::DropZero, TaylorOrders::DropZeroOdd] {
        let mut value = l2_root(2, false);
        let FactorizedKind::Closure(spec) = &mut value.kind else {
            unreachable!()
        };
        spec.taylor.as_mut().unwrap().orders = orders;
        value.label = format!("L2 root orders={}", orders.as_str());
        let prefix = format!("orders={} is not supported: ", orders.as_str());
        let refused = |result: std::result::Result<(), FactorizedError>, what: &str| match result {
            Err(FactorizedError::Closure { detail }) => {
                assert!(detail.starts_with(&prefix), "{what}: {detail}")
            }
            other => panic!("{what}: {other:?}"),
        };
        refused(value.validate(), "validate");
        refused(
            evaluate_factorized_closure(
                &vakint,
                &settings(),
                std::slice::from_ref(&value),
                &[value.prefactor.total()],
            )
            .map(|_| ()),
            "entry",
        );
        // the text form keeps the field: parse_text reads it, from_text
        // (which validates) refuses it (the diagnostic writer's lines with the
        // record schema; its own schema line is refused by every reader)
        let text = value.to_text_unvalidated_for_diagnostics().unwrap();
        assert!(FactorizedNumerator::parse_text(&text).is_err());
        let text = text.replacen(
            crate::p_integral::factorized::DIAGNOSTIC_SCHEMA_SUFFIX,
            "",
            1,
        );
        assert!(text.contains(&format!("orders={}", orders.as_str())));
        let parsed = FactorizedNumerator::parse_text(&text).unwrap();
        refused(parsed.validate(), "parsed");
        refused(
            FactorizedNumerator::from_text(&text).map(|_| ()),
            "from_text",
        );
    }
    check(&l2_root(2, false), &L2_CUT);
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
    v.prefactor.scalars = atom(&format!("UFO::gs^2*{R}"));
    refused(&v, "reference mass");
    // a bispinor port in no factor: refused by the port rule (the open-slot
    // rule itself: `refusal_open_spinor_slot`)
    let mut v = l2_child(1);
    v.ports.push(atom(&format!(
        "spenso::{{spenso::representation,spenso::self_dual}}::bis(4,{HEDGE}(150,1))"
    )));
    match v.validate() {
        Err(FactorizedError::Ports { detail }) => assert!(
            detail.starts_with("port ") && detail.contains("does not occur in the factors"),
            "{detail}"
        ),
        other => panic!("a bispinor port in no factor: {other:?}"),
    }
    // not a closure record
    let mut v = l2_root(1, true);
    v.kind = FactorizedKind::Direct;
    refused(&v, "not a closure record");
    // an insertion on a vertex that is neither a contracted (ct:) nor a
    // surrogate (ins:) vertex is refused BY THAT RULE: vertex 2's polynomial
    // (degree 1) moved into an insertion with the vertex's own half-edges and
    // factor 1, the rule name left as it is; with the rule `ct:` the same
    // record is accepted
    let mut v = l2_root(1, true);
    let hedges = v
        .local_hedges(super::super::Location::Vertex(2))
        .into_iter()
        .collect::<Vec<_>>();
    let vertex = v.vertices.iter_mut().find(|x| x.id == 2).unwrap();
    let factor = std::mem::replace(&mut vertex.factor, Atom::num(1));
    v.insertions.push(super::super::Insertion {
        vertex: 2,
        hedges,
        cycles: (0, 0),
        coefficient: Atom::num(1),
        factor,
    });
    match v.validate() {
        Err(FactorizedError::Insertion { vertex: 2, detail }) => {
            assert!(detail.starts_with("the vertex rule "), "{detail}")
        }
        other => panic!("neither ct: nor ins:: {other:?}"),
    }
    v.vertices.iter_mut().find(|x| x.id == 2).unwrap().rule = "ct:child".into();
    v.validate().unwrap();
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
    // a wrong open-slot list (one line more, the `lines` count raised with
    // it: only the open list is wrong): refused by the open-slot rule; the
    // same record with its true (empty) list reads
    assert!(FactorizedNumerator::parse_text(&text).is_ok());
    let wrong = text
        .replace(
            "open\t0\n",
            &format!("open\t1\nopen_slot\t{}\n", hedge(120)),
        )
        .lines()
        .map(|l| match l.strip_prefix("lines\t") {
            Some(n) => format!("lines\t{}", n.parse::<usize>().unwrap() + 1),
            None => l.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    assert_eq!(wrong.lines().count(), text.lines().count() + 1);
    match FactorizedNumerator::parse_text(&wrong) {
        Err(FactorizedError::TextSyntax { detail, .. }) => assert!(
            detail.starts_with("open slots ")
                && detail.contains("differ from the ports the projector leaves open"),
            "{detail}"
        ),
        other => panic!("a wrong open-slot list: {other:?}"),
    }
    // kind direct in v3, kind closure in v1
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
    let root = dot_root();
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
        assert!(
            files.len() >= *n,
            "{}: {} DOT graphs, {n} requested",
            directory.display(),
            files.len()
        );
        out.extend(dot_files_as_roots(
            level,
            &model,
            files.into_iter().take(*n),
        ));
    }
    out
}

/// The R* producer's prepared DOT graphs (`<root>/L<n>/diagram-*.dot` and
/// `model.json`): `VAKINT_TEST_FACTORIZED_DOT_ROOT`, REQUIRED (the tests that
/// read it are `#[ignore]`; run, a missing root is a failure).
#[cfg(feature = "feynkit-ingress")]
fn dot_root() -> std::path::PathBuf {
    let root = wh_fixtures::required_env(
        "VAKINT_TEST_FACTORIZED_DOT_ROOT",
        "the DOT suites read the producer's prepared DOT graphs",
    );
    assert!(
        root.join("L1").join("model.json").is_file(),
        "no DOT data under {}",
        root.display()
    );
    root
}

/// Named DOT graphs as root remainders ([`dot_root_remainders`]).
#[cfg(feature = "feynkit-ingress")]
fn dot_named_root_remainders(
    names: &[(&str, &str)],
) -> Vec<(FactorizedNumerator, Vec<(String, String)>)> {
    let root = dot_root();
    let mut out = Vec::new();
    for (level, name) in names {
        let directory = root.join(level);
        let path = directory.join(name);
        assert!(path.is_file(), "missing DOT graph {}", path.display());
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

/// External data (`VAKINT_TEST_FACTORIZED_DOT_ROOT`): `#[ignore]`; run, the
/// six graphs (L2 x3, L3 x2, L4 x1) must all be read and compared.
#[cfg(feature = "feynkit-ingress")]
#[test]
#[ignore = "external data: VAKINT_TEST_FACTORIZED_DOT_ROOT"]
fn beta_dot_graphs_as_root_remainders_equal_the_contracted_reference() {
    let cases = dot_root_remainders(&[("L2", 3), ("L3", 2), ("L4", 1)]);
    let mut compared = 0;
    for (value, routing) in cases {
        let routing = routing
            .iter()
            .map(|(a, b)| (a.as_str(), b.as_str()))
            .collect::<Vec<_>>();
        check(&value, &routing);
        let text = value.to_text().unwrap();
        assert_eq!(FactorizedNumerator::from_text(&text).unwrap(), value);
        compared += 1;
    }
    wh_fixtures::assert_compared("DOT_ROOTS", compared, Some(6));
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
    if !with_reference {
        eprintln!(
            "CLOSURE_TIMING: VAKINT_CLOSURE_TIMING_NO_REFERENCE set: NO VALUE IS COMPARED (timing only)"
        );
    }
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
    atom(wh_fixtures::value_line(text))
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
/// eps^n`, signed).
fn valid_through(text: &str) -> i64 {
    wh_fixtures::valid_through(text)
}

/// A fixture's value file against `got` in `gammalooprs::ε`: a DECLARED zero
/// ([`wh_fixtures::is_declared_zero`]) must be zero through valid_through,
/// every other value equal through valid_through ([`assert_equal_through`]).
fn assert_fixture(got: &Atom, value_text: &str, what: &str) {
    let order = valid_through(value_text);
    if wh_fixtures::is_declared_zero(value_text) {
        let epsilon = symbolica::symbol!("gammalooprs::ε");
        wh_fixtures::assert_declared_zero(got, epsilon, order, what);
        eprintln!("{what}: declared zero, zero through eps^{order}");
    } else {
        assert_equal_through(got, &fixture_value(value_text), order, what);
    }
}

/// [`wh_fixtures::assert_equal_through`] in `gammalooprs::ε`: exact through
/// `order` (SIGNED) and NO FURTHER, the reference nonzero through `order`.
fn assert_equal_through(got: &Atom, expected: &Atom, order: i64, what: &str) {
    let epsilon = symbolica::symbol!("gammalooprs::ε");
    let compared = wh_fixtures::assert_equal_through(got, expected, epsilon, order, what);
    eprintln!("{what}: equal through eps^{order} ({compared} orders compared)");
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
    let path = wh_fixtures::required_env("VAKINT_CLOSURE_PROBE_RECORD", "the record");
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
/// valid_through, and its contraction seconds printed beside theirs. Data:
/// `$VAKINT_TEST_WH_FIXTURES/k1-quick` (`#[ignore]`; run, missing data is a
/// failure and all 35 records must be compared, unless
/// `VAKINT_TEST_WH_K1_ONLY` selects some by suffix).
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn the_beta_k1_quick_records_equal_their_values() {
    let records = wh_fixtures::record_dirs(&wh_fixtures::fixture_dir("k1-quick"));
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
        assert_fixture(&got.value, &value_text, &what);
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
    wh_fixtures::assert_compared("K1_QUICK", checked, only.is_none().then_some(35));
}

/// A child insertion of momentum degree 2 in a record with Taylor degree
/// w = 2 goes through the same routing, lambda scaling and truncation as a
/// vertex factor: the entry on the record with the polynomial on a
/// contracted (ct:) vertex equals the INDEPENDENT reference of the record
/// with the polynomial multiplied into that vertex's own factor by hand.
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
    let argument = atom("UFO::gs^2*witness::x");
    let got = evaluate_factorized_closure(
        &vakint,
        &settings(),
        std::slice::from_ref(&with),
        std::slice::from_ref(&argument),
    )
    .unwrap();
    let expected = independent_reference(&base, &L2_CUT, &argument);
    assert!(!expected.is_zero());
    assert!(
        (&got.records[0].vakint_value - &expected)
            .expand()
            .is_zero(),
        "insertion {}\nindependent (by hand) {expected}",
        got.records[0].vakint_value
    );
}

/// An insertion factor that depends on the dimension outside its slots
/// through integer powers of sums is a Laurent factor, and its VALUE is right
/// through the expansion: on a contracted (ct:) vertex of the L2 root,
/// - `1/((D-1)(D-2))` (a projector normalisation, pole order 0): the entry
///   equals the INDEPENDENT reference of the record with the factor on the
///   vertex by hand, at the settings' depth, where the O(ε) terms of
///   1/((3-2ε)(2-2ε)) multiply the integral's poles (the value with the
///   factor at ε = 0, 1/6, differs);
/// - `1/(D-4)` (pole order 1): the entry raises the depth by one; it equals
///   1/(-2ε) times the independent reference of the record without it one
///   term deeper, truncated at ε^(depth - L - 1).
#[test]
fn dimension_dependent_insertion_factors_equal_the_independent_reference() {
    let compact = format!("{MINK}({DIM})");
    let polynomial = atom(&format!(
        "spenso::dot({},{})+spenso::dot({},{})",
        q(1, &compact),
        q(2, &compact),
        q(3, &compact),
        q(3, &compact)
    ));
    let record = |vertex_factor: Atom| {
        let mut base = l2_root(2, true);
        base.edges.iter_mut().find(|e| e.id == 21).unwrap().factor = atom(&q(2, &hedge(121)));
        base.vertices.iter_mut().find(|v| v.id == 2).unwrap().factor = vertex_factor;
        base
    };
    let inserted = |factor: Atom| {
        let mut with = record(Atom::num(1));
        let hedges = with
            .local_hedges(super::super::Location::Vertex(2))
            .into_iter()
            .collect::<Vec<_>>();
        let v = with.vertices.iter_mut().find(|v| v.id == 2).unwrap();
        v.rule = "ct:child".into();
        with.insertions.push(super::super::Insertion {
            vertex: 2,
            hedges,
            cycles: (0, 0),
            coefficient: Atom::num(1),
            factor,
        });
        with.validate().unwrap();
        with
    };
    let vakint = Vakint::new().unwrap();
    let argument = atom("UFO::gs^2*witness::x");
    let entry = |with: &FactorizedNumerator| {
        evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(with),
            std::slice::from_ref(&argument),
        )
        .unwrap()
        .records[0]
            .vakint_value
            .clone()
    };
    let epsilon = crate::utils::vakint_macros::vk_symbol!(&settings().epsilon_symbol);
    let in_epsilon = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    // pole order 0
    let normalisation = atom(&format!("1/(({DIM}-1)*({DIM}-2))"));
    let with = inserted(&polynomial * &normalisation);
    assert_eq!(with.insertion_pole_order(epsilon).unwrap(), 0);
    let got = entry(&with);
    let expected = independent_reference(&record(&polynomial * &normalisation), &L2_CUT, &argument);
    assert!(!expected.is_zero());
    assert!(
        (&got - &expected).expand().is_zero(),
        "1/((D-1)(D-2)): entry {got}\nindependent {expected}"
    );
    let at_zero = independent_reference(&record(&polynomial / Atom::num(6)), &L2_CUT, &argument);
    assert!(
        !(&expected - &at_zero).expand().is_zero(),
        "the O(eps) terms of the normalisation must matter at this depth"
    );
    // pole order 1
    let pole = atom(&format!("1/({DIM}-4)"));
    let with = inserted(&polynomial * &pole);
    assert_eq!(with.insertion_pole_order(epsilon).unwrap(), 1);
    let got = entry(&with);
    let mut deeper = settings();
    deeper.number_of_terms_in_epsilon_expansion += 1;
    let bare = independent_reference_with(&record(polynomial.clone()), &L2_CUT, &argument, deeper);
    let target = settings().number_of_terms_in_epsilon_expansion - with.loops() as i64 - 1;
    let pole_in_epsilon = pole
        .replace(Atom::var(with.dimension))
        .with(in_epsilon.clone());
    let expected = (&pole_in_epsilon * &bare)
        .series(
            epsilon,
            Atom::Zero.as_view(),
            symbolica::domains::rational::Rational::from(target),
        )
        .unwrap()
        .to_atom()
        .expand();
    assert!(!expected.is_zero());
    assert!(
        (&got - &expected).expand().is_zero(),
        "1/(D-4): entry {got}\nindependent, one term deeper {expected}"
    );
}

/// A negative power of a sum whose terms cancel at ε = 0 keeps its EXACT pole
/// order through the closure entry: ε^-55·((1-ε)^-1 - 1 - ε)^-5 (degree
/// count 60, pole order 65) is refused by the entry before any reduction
/// (`Insertion`, the total above 64); the twin ((1-ε)^-1 - 1 - ε)^-1 (count
/// 1, pole order 2: the depth is raised by 2) equals the coefficient times the INDEPENDENT
/// reference of the record without it, two terms deeper, truncated at
/// ε^(depth - L - 1).
#[test]
fn a_cancelling_coefficient_raises_the_depth_by_its_exact_pole_order() {
    let eps = settings().epsilon_symbol.clone();
    let witness =
        closure_with_insertion(&format!("{eps}^-55*((1-{eps})^-1-1-{eps})^-5"), &settings())
            .unwrap_err();
    match &witness {
        FactorizedError::Insertion { vertex: 1, detail } => assert!(
            detail.starts_with(
                "an ε power out of range: the insertions' pole orders add up to 65, above 64"
            ),
            "{detail}"
        ),
        other => panic!("expected the Insertion refusal at vertex 1, got {other:?}"),
    }
    let mut base = l2_root(1, true);
    base.vertices.iter_mut().find(|v| v.id == 1).unwrap().factor = atom(&q(1, &hedge(120)));
    let coefficient = atom(&format!("((1-{eps})^-1-1-{eps})^-1"));
    let got = closure_with_insertion(&coefficient.to_canonical_string(), &settings()).unwrap();
    let epsilon = crate::utils::vakint_macros::vk_symbol!(&eps);
    let mut deeper = settings();
    deeper.number_of_terms_in_epsilon_expansion += 2;
    let bare = independent_reference_with(&base, &L2_CUT, &base.prefactor.total(), deeper);
    let target = settings().number_of_terms_in_epsilon_expansion - base.loops() as i64 - 1;
    let expected = (&coefficient * &bare)
        .series(
            epsilon,
            Atom::Zero.as_view(),
            symbolica::domains::rational::Rational::from(target),
        )
        .unwrap()
        .to_atom()
        .expand();
    assert!(!expected.is_zero());
    assert!(
        (&got.records[0].vakint_value - &expected)
            .expand()
            .is_zero(),
        "entry {}\nindependent, two terms deeper {expected}",
        got.records[0].vakint_value
    );
}

/// The beta's k5-stress set (L5 graph 3793, closed surrogate root, five
/// records with four insertion vertices each; r0051 has five loops and needs
/// the 1-4 loop package of `VAKINT_TEST_P_LIBRARY_1_4`): each record equal to
/// its value through valid_through, with the entry's seconds per stage.
/// Data: `$VAKINT_TEST_WH_FIXTURES/k5-stress` (`#[ignore]`; run, missing data
/// is a failure, all 5 records are compared, r0051 requires
/// `VAKINT_TEST_P_LIBRARY_1_4`).
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES (+ VAKINT_TEST_P_LIBRARY_1_4)"]
fn the_beta_k5_stress_records_equal_their_values() {
    let records = wh_fixtures::record_dirs(&wh_fixtures::fixture_dir("k5-stress"));
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
    let epsilon = crate::utils::vakint_macros::vk_symbol!("gammalooprs::ε");
    let mut compared = 0;
    for record in &records {
        let input = std::fs::read_to_string(record.join("input.txt")).unwrap();
        let value_text = std::fs::read_to_string(record.join("value.txt")).unwrap();
        let value = FactorizedNumerator::from_text(&input)
            .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        assert!(
            value.loops() <= 4 || package.is_some(),
            "{}: five loops need VAKINT_TEST_P_LIBRARY_1_4",
            record.display()
        );
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
        assert_fixture(&got.value, &value_text, &what);
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
        compared += 1;
    }
    wh_fixtures::assert_compared("K5_STRESS", compared, Some(5));
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
    // the base record (no insertion) against the independent reference
    check_independent(&base, &L2_CUT, &atom("UFO::gs^2*witness::x"));
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

/// The closure entry is the only caller of `insertion_pole_order` (its depth
/// raise). The insertions whose Laurent series in ε would fail there are
/// refused first, by the entry's validation, `Insertion`: exp(1/ε) (an
/// essential singularity, not a Laurent polynomial) and ε^(-10^30) (an
/// integer power beyond i64, which symbolica does not expand: an ε power out
/// of range). The pole order's `Evaluation` is not reached through the entry
/// by any input found.
#[test]
fn a_non_laurent_insertion_is_refused_by_the_closure_entry_before_its_pole_order() {
    let eps = "gammalooprs::ε";
    let refused = closure_insertion_error(&format!("exp(1/{eps})"));
    match &refused {
        FactorizedError::Insertion { vertex: 1, detail } => assert!(
            detail.starts_with("the dependence on ε is not a Laurent polynomial"),
            "{detail}"
        ),
        other => panic!("expected the Insertion refusal at vertex 1, got {other:?}"),
    }
    let huge = closure_insertion_error(&format!("{eps}^(-10^30)"));
    match &huge {
        FactorizedError::Insertion { vertex: 1, detail } => {
            assert!(detail.contains("exceeds 64 in modulus"), "{detail}")
        }
        other => panic!("expected the Insertion refusal at vertex 1, got {other:?}"),
    }
}

/// The closure entry's error for the L2 root with vertex 1 as a contracted
/// child of coefficient `coefficient` (its polynomial as the factor).
fn closure_insertion_error(coefficient: &str) -> FactorizedError {
    closure_with_insertion(coefficient, &rstar_site_settings()).unwrap_err()
}

/// The closure entry on the L2 root with vertex 1 as a contracted child of
/// coefficient `coefficient` (its polynomial as the factor), `settings`.
fn closure_with_insertion(
    coefficient: &str,
    settings: &VakintSettings,
) -> std::result::Result<ClosureValue, FactorizedError> {
    let mut with = l2_root(1, true);
    with.vertices.iter_mut().find(|v| v.id == 1).unwrap().factor = atom(&q(1, &hedge(120)));
    let hedges = with
        .local_hedges(super::super::Location::Vertex(1))
        .into_iter()
        .collect::<Vec<_>>();
    let v = with.vertices.iter_mut().find(|v| v.id == 1).unwrap();
    let factor = std::mem::replace(&mut v.factor, Atom::num(1));
    v.rule = "ct:child1".into();
    with.insertions.push(super::super::Insertion {
        vertex: 1,
        hedges,
        cycles: (0, 0),
        coefficient: atom(coefficient),
        factor,
    });
    let vakint = Vakint::new().unwrap();
    evaluate_factorized_closure(
        &vakint,
        settings,
        std::slice::from_ref(&with),
        &[with.prefactor.total()],
    )
}

/// Validation recognises ε by its name; with an epsilon setting of another
/// name (`user::ep`) the closure entry applies the pole order's rules with it:
/// log(user::ep) refused (`Insertion`: `log` is a function head that is not
/// allowed, whatever its argument); the twin with a pole (user::ep^-1)
/// evaluates.
#[test]
fn the_closure_entry_applies_the_laurent_rule_to_an_epsilon_setting_of_another_name() {
    let mut user = rstar_site_settings();
    user.epsilon_symbol = "user::ep".into();
    match closure_with_insertion("log(user::ep)", &user).unwrap_err() {
        FactorizedError::Insertion { vertex: 1, detail } => assert!(
            detail.starts_with(
                "a function head that is not allowed in insertion 1 coefficient: symbolica::log"
            ),
            "{detail}"
        ),
        other => panic!("expected the Insertion refusal at vertex 1, got {other:?}"),
    }
    closure_with_insertion("user::ep^-1", &user).expect("user::ep^-1 is a Laurent coefficient");
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
    match more.validate() {
        Err(FactorizedError::Closure { detail }) => assert!(
            detail.starts_with("the carried graph and its insertions have "),
            "{detail}"
        ),
        other => panic!("the insertions' cycles exceed original_cycles: {other:?}"),
    }
}

/// The beta's k1-l5 set: K1 records of every boundary root remainder of a
/// FIVE-LOOP diagram (most with child insertions), each against value.txt
/// through valid_through, its seconds printed beside the RR_CONTRACT of
/// meta.json. Heavy (five-loop reductions). Data:
/// `$VAKINT_TEST_WH_FIXTURES/k1-l5/<graph>` and the 1-4 loop package of
/// `VAKINT_TEST_P_LIBRARY_1_4` (`#[ignore]`; run, missing data is a failure);
/// `VAKINT_TEST_WH_K1_L5_GRAPHS` (comma list, default all three) selects
/// graphs, `VAKINT_TEST_WH_K1_L5_INSERTIONS_ONLY` the records with insertions;
/// the count compared is asserted for the selection (L5-60 15 / 14, L5-51
/// 6 / 5, L5-34342 1 / 0).
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES + VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k1_l5_records_equal_their_values() {
    // the graphs (all three by default) and their record counts: all
    // records, or (VAKINT_TEST_WH_K1_L5_INSERTIONS_ONLY) those with insertions
    const GRAPHS: [(&str, usize, usize); 3] =
        [("L5-60", 15, 14), ("L5-51", 6, 5), ("L5-34342", 1, 0)];
    let insertions_only = std::env::var_os("VAKINT_TEST_WH_K1_L5_INSERTIONS_ONLY").is_some();
    let selected = std::env::var("VAKINT_TEST_WH_K1_L5_GRAPHS").ok();
    let graphs = GRAPHS
        .iter()
        .filter(|(g, _, _)| {
            selected
                .as_ref()
                .is_none_or(|s| s.split(',').any(|x| x == *g))
        })
        .collect::<Vec<_>>();
    assert!(
        !graphs.is_empty(),
        "VAKINT_TEST_WH_K1_L5_GRAPHS selects no graph"
    );
    let expected = graphs
        .iter()
        .map(|(_, all, with)| if insertions_only { *with } else { *all })
        .sum::<usize>();
    let package = wh_fixtures::required_env(
        "VAKINT_TEST_P_LIBRARY_1_4",
        "five-loop records need the 1-4 P library",
    );
    let root = wh_fixtures::fixture_dir("k1-l5");
    let records = graphs
        .iter()
        .flat_map(|(g, _, _)| wh_fixtures::record_dirs(&root.join(g)))
        .collect::<Vec<_>>();
    let mut vakint = Vakint::new().unwrap();
    let library = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(package),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap();
    vakint.inject_p_library(library).unwrap();
    let (mut mine, mut theirs, mut compared) = (0.0, 0.0, 0);
    for record in &records {
        let input = std::fs::read_to_string(record.join("input.txt")).unwrap();
        let value_text = std::fs::read_to_string(record.join("value.txt")).unwrap();
        let value = FactorizedNumerator::from_text(&input)
            .unwrap_or_else(|e| panic!("{}: {e}", record.display()));
        if value.insertions.is_empty() && insertions_only {
            eprintln!(
                "K1_L5 {} SKIP: no insertion (VAKINT_TEST_WH_K1_L5_INSERTIONS_ONLY)",
                record.display()
            );
            continue;
        }
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
        assert_fixture(&got.value, &value_text, &what);
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
        compared += 1;
    }
    eprintln!(
        "K1_L5 SUM records={compared} entry_contract_s={mine:.3} their_rr_contract_s={theirs:.3}"
    );
    wh_fixtures::assert_compared("K1_L5", compared, Some(expected));
}

/// REGRESSION of the placeholder collision: a record naming one of the
/// entry's internal symbols (`vakint::closure_tensor_<n>`, the Taylor
/// parameter) as a legal scalar used to change the value silently (the
/// placeholder substitution erased or replaced it); it is now refused, in the
/// prefactor and inside a factor, closed and with open ports.
#[test]
fn internal_symbol_names_are_refused_in_a_record() {
    let vakint = Vakint::new().unwrap();
    for (mut value, name, inside) in [
        (l2_root(2, true), "vakint::closure_tensor_0", false),
        (l2_root(2, true), "vakint::closure_tensor_0", true),
        (l2_child(2), "vakint::closure_tensor_1", false),
        (l2_child(2), "vakint::closure_tensor_1", true),
        (l2_root(2, true), "vakint::closure_lambda", true),
    ] {
        let named = atom(name);
        if inside {
            value.vertices[0].factor = &value.vertices[0].factor * &named;
        } else {
            value.prefactor.scalars = &value.prefactor.scalars * &named;
        }
        let refusal = evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(&value),
            &[value.prefactor.total()],
        );
        match refusal {
            Err(FactorizedError::Graph { detail }) => {
                assert!(
                    detail.starts_with("reserved symbol ") && detail.contains(name),
                    "{name}: {detail}"
                );
                eprintln!("RESERVED {name} inside_factor={inside}: {detail}");
            }
            other => panic!("{name} inside_factor={inside}: not refused: {other:?}"),
        }
        // the renamed record still evaluates
        let mut renamed = value.clone();
        let innocent = atom("witness::collision_scalar");
        for factor in renamed
            .vertices
            .iter_mut()
            .map(|v| &mut v.factor)
            .chain([&mut renamed.prefactor.scalars])
        {
            *factor = factor.replace(named.clone()).with(innocent.clone());
        }
        evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(&renamed),
            &[renamed.prefactor.total()],
        )
        .unwrap();
    }
}

/// Every internal name of the sweep, as a variable and (for the vakint
/// notation's heads) as a function head.
fn reserved_names() -> Vec<&'static str> {
    vec![
        "vakint::closure_tensor_0",
        "vakint::closure_tensor_7",
        "vakint::closure_lambda",
        "vakint::slash(-1,0)",
        "vakint::nf",
        "vakint::p(1)",
        "vakint::k(1)",
        "vakint::dot(vakint::p(1),vakint::p(1))",
        "vakint::g(9001,9002)",
        "vakint::g(1,2)",
        "vakint::dot(vakint::p(1),UFO::x)",
        "UFO::coupling(vakint::p(1))",
        "vakint::topo(1)",
        "vakint::prop(1,vakint::edge(0,1),vakint::k(1),0,1)",
        "vakint::mursq",
        "vakint::graph_entry_branch_0",
        "vakint::logMassiveClosure",
        "vakint::any_name_at_all",
    ]
}

fn assert_reserved(value: &FactorizedNumerator, name: &str, field: &str) {
    match value.validate() {
        Err(FactorizedError::Graph { detail }) if detail.starts_with("reserved symbol ") => {
            assert!(
                detail.contains(field),
                "{name}: {detail} (expected field {field})"
            );
        }
        other => panic!("{name} in {field}: not refused as reserved: {other:?}"),
    }
}

/// One refusal per reserved name in a closure record's numerator factors
/// (edge factor and vertex factor).
#[test]
fn reserved_symbols_are_refused_in_factors() {
    for name in reserved_names() {
        let symbol = atom(name);
        let mut v = l2_root(2, true);
        v.edges[1].factor = &v.edges[1].factor * &symbol;
        assert_reserved(&v, name, "edge 2 factor");
        let mut v = l2_root(2, true);
        v.vertices[1].factor = &v.vertices[1].factor * &symbol;
        assert_reserved(&v, name, "vertex 2 factor");
    }
}

/// ... in the prefactor fields (overall, scalars, colour).
#[test]
fn reserved_symbols_are_refused_in_the_prefactor() {
    for name in reserved_names() {
        let symbol = atom(name);
        for field in ["prefactor_overall", "prefactor_scalars", "prefactor_colour"] {
            let mut v = l2_root(2, true);
            let target = match field {
                "prefactor_overall" => &mut v.prefactor.overall,
                "prefactor_scalars" => &mut v.prefactor.scalars,
                _ => &mut v.prefactor.colour,
            };
            *target = &*target * &symbol;
            assert_reserved(&v, name, field);
        }
    }
}

/// ... in the projector, the ports, the masses (an edge's, the closure
/// line's, the auxiliary one), the weights and the dimension.
#[test]
fn reserved_symbols_are_refused_in_projector_ports_masses_weights() {
    for name in reserved_names() {
        let symbol = atom(name);
        let mut v = l2_root(2, true);
        let projector = v.projector.as_mut().unwrap();
        projector.tensor = &projector.tensor * &symbol;
        assert_reserved(&v, name, "projector");
        let mut v = l2_child(2);
        v.ports.push(function!(
            symbolica::symbol!(names::MINK),
            Atom::var(v.dimension),
            symbol.clone()
        ));
        if name.starts_with("vakint::slash(") {
            // the reader's own dummy slot form mink(D, vakint::slash(..)) is
            // exempt from the namespace: a forged one is refused by
            // the index rules (here: a port that occurs in no factor)
            match v.validate() {
                Err(FactorizedError::Ports { detail }) => {
                    assert!(
                        detail.ends_with("does not occur in the factors"),
                        "{detail}"
                    )
                }
                other => panic!("forged slash port: {other:?}"),
            }
        } else {
            assert_reserved(&v, name, "ports");
        }
        let mut v = l2_root(2, true);
        v.edges[0].mass_squared = symbol.clone();
        assert_reserved(&v, name, "edge 1 mass");
        let mut v = l2_root(2, true);
        spec(&mut v).auxiliary_mass_squared = symbol.clone();
        assert_reserved(&v, name, "taylor auxiliary mass");
        let mut v = l2_root(2, true);
        spec(&mut v).mass_squared = symbol.clone();
        assert_reserved(&v, name, "closure mass");
        let mut v = l2_root(2, true);
        spec(&mut v).weights = vec![("scalars".into(), false, symbol.clone())];
        assert_reserved(&v, name, "weight scalars");
    }
    let mut v = l2_root(2, true);
    v.dimension = symbolica::symbol!("vakint::closure_lambda");
    assert_reserved(&v, "vakint::closure_lambda", "dimension");
}

/// ... in an insertion's coefficient and factor (a record with insertions is
/// a lifted kind: closure and direct records alike).
#[test]
fn reserved_symbols_are_refused_in_insertions() {
    let record = include_str!("../fixtures/wh-k2-child-ins.record.txt");
    let base = FactorizedNumerator::from_text(record).unwrap();
    base.validate().unwrap();
    assert!(!base.insertions.is_empty());
    let vertex = base.insertions[0].vertex;
    for name in reserved_names() {
        let symbol = atom(name);
        let mut v = base.clone();
        v.insertions[0].coefficient = &v.insertions[0].coefficient * &symbol;
        assert_reserved(&v, name, &format!("insertion {vertex} coefficient"));
        let mut v = base.clone();
        v.insertions[0].factor = &v.insertions[0].factor * &symbol;
        assert_reserved(&v, name, &format!("insertion {vertex} factor"));
    }
}

/// A v3 DIRECT record is a lifted kind too; a v1 direct record (the
/// pre-existing adapter's input) keeps its earlier behaviour.
#[test]
fn reserved_symbols_in_direct_records() {
    let record = include_str!("../fixtures/wh-k5-direct.record.txt");
    let base = FactorizedNumerator::from_text(record).unwrap();
    base.validate().unwrap();
    assert!(base.record.is_some() || !base.insertions.is_empty());
    let mut v = base.clone();
    v.prefactor.scalars = &v.prefactor.scalars * atom("vakint::closure_tensor_0");
    assert_reserved(&v, "vakint::closure_tensor_0", "prefactor_scalars");
    // v1, no insertion: not a lifted kind, the namespace is not checked
    let mut v1 = base.clone();
    v1.record = None;
    v1.insertions.clear();
    v1.prefactor.scalars = &v1.prefactor.scalars * atom("vakint::closure_tensor_0");
    assert!(!matches!(
        v1.validate(),
        Err(FactorizedError::Graph { detail }) if detail.starts_with("reserved symbol")
    ));
}

/// vakint's default epsilon is the one vakint symbol a record may carry (a
/// default-settings caller writes its insertion coefficients in it).
#[test]
fn the_default_epsilon_is_not_reserved() {
    let mut v = l2_root(2, true);
    v.prefactor.scalars = &v.prefactor.scalars * atom("1+vakint::ε");
    v.validate().unwrap();
}

/// The comparator FAILS FAST on a wrong oracle: k1-root's value perturbed at
/// its lowest valid order (eps^-2, valid through eps^-1) is reported as
/// different within seconds (an unsigned cast of the negative bound used to
/// request ~4e9 orders instead).
#[test]
fn a_wrong_oracle_fails_fast() {
    let text = include_str!("../fixtures/wh-k1-root.value.txt");
    let expected = fixture_value(text);
    let wrong = &expected + atom("UFO::gs_minus^2*gammalooprs::ε^(-2)");
    let order = valid_through(text);
    assert!(order < 0, "the fixture is valid through a negative order");
    let started = std::time::Instant::now();
    let outcome =
        std::panic::catch_unwind(|| assert_equal_through(&expected, &wrong, order, "wrong oracle"));
    let seconds = started.elapsed().as_secs_f64();
    assert!(outcome.is_err(), "the wrong oracle was not reported");
    assert!(
        seconds < 10.0,
        "the wrong oracle took {seconds} s to report"
    );
    eprintln!("WRONG_ORACLE reported in {seconds:.3} s");
}

/// The comparator refuses a ZERO reference (a fixture whose value parsed to
/// zero, or truncated above `valid_through`, cannot pass against a zero
/// entry), and compares a nonzero one order by order.
#[test]
fn the_comparator_needs_a_nonzero_reference() {
    let text = include_str!("../fixtures/wh-k1-root.value.txt");
    let order = valid_through(text);
    let zero = std::panic::catch_unwind(|| {
        assert_equal_through(&Atom::Zero, &Atom::Zero, order, "zero reference")
    });
    assert!(zero.is_err(), "a zero reference passed");
    // a reference that is nonzero only ABOVE valid_through is zero for the test
    let above = atom(&format!("UFO::x*gammalooprs::ε^({})", order + 1));
    let truncated = std::panic::catch_unwind(|| {
        assert_equal_through(&above, &above, order, "reference above valid_through")
    });
    assert!(
        truncated.is_err(),
        "a reference zero through valid_through passed"
    );
    let expected = fixture_value(text);
    assert_equal_through(&expected, &expected, order, "the fixture against itself");
}

// ---------------------------------------------------------------------------
// The reserved namespace and the caller-symbol rules: priority (1) tests.
// Each refusal test asserts the VARIANT and the DETAIL PREFIX on a minimal
// input (a valid record with one defect) and that the same input without the
// defect is accepted.
// ---------------------------------------------------------------------------

/// `Err(Graph)` whose detail starts with `prefix`.
fn assert_graph_prefix(result: std::result::Result<(), FactorizedError>, prefix: &str, what: &str) {
    match result {
        Err(FactorizedError::Graph { detail }) => {
            assert!(
                detail.starts_with(prefix),
                "{what}: {detail}\n(expected prefix {prefix:?})"
            )
        }
        other => panic!("{what}: expected Graph {prefix:?}, got {other:?}"),
    }
}

/// The full reserved-symbol refusal: a RECORD field names the rule for the
/// reader's default namespace, an ARGUMENT (the prefactor argument, a
/// setting) the reserved namespace.
fn reserved_prefix(name: &str, field: &str) -> String {
    if field == "prefactor argument" || field.ends_with("setting") {
        format!(
            "reserved symbol {name} in {field}: the vakint:: namespace is reserved for vakint's internal symbols"
        )
    } else {
        format!(
            "reserved symbol {name} in {field}: the vakint:: namespace is reserved for vakint's internal symbols; a name written without a namespace is read into vakint::"
        )
    }
}

/// vakint's default epsilon symbol (the default settings' epsilon).
fn vakint_epsilon() -> symbolica::atom::Symbol {
    crate::utils::vakint_macros::vk_symbol!(&VakintSettings::default().epsilon_symbol)
}

/// F1: `vakint::closure_lambda` AND `vakint::closure_tensor_0` in
/// EVERY atom of a closure record, one case per field, with the exact message.
#[test]
fn closure_lambda_and_closure_tensor_0_are_refused_in_every_field() {
    for name in ["vakint::closure_lambda", "vakint::closure_tensor_0"] {
        let x = atom(name);
        let mut cases: Vec<(FactorizedNumerator, String)> = Vec::new();
        let mut v = l2_root(2, true);
        v.prefactor.overall = &v.prefactor.overall * &x;
        cases.push((v, "prefactor_overall".into()));
        let mut v = l2_root(2, true);
        v.prefactor.scalars = &v.prefactor.scalars * &x;
        cases.push((v, "prefactor_scalars".into()));
        let mut v = l2_root(2, true);
        v.prefactor.colour = &v.prefactor.colour * &x;
        cases.push((v, "prefactor_colour".into()));
        let mut v = l2_root(2, true);
        v.edges[0].mass_squared = x.clone();
        cases.push((v, "edge 1 mass".into()));
        let mut v = l2_root(2, true);
        spec(&mut v).auxiliary_mass_squared = x.clone();
        cases.push((v, "taylor auxiliary mass".into()));
        let mut v = l2_root(2, true);
        spec(&mut v).mass_squared = x.clone();
        cases.push((v, "closure mass".into()));
        let mut v = l2_root(2, true);
        v.edges[0].factor = &v.edges[0].factor * &x;
        cases.push((v, "edge 1 factor".into()));
        let mut v = l2_root(2, true);
        v.vertices[0].factor = &v.vertices[0].factor * &x;
        cases.push((v, "vertex 1 factor".into()));
        let mut v = l2_root(2, true);
        let projector = v.projector.as_mut().unwrap();
        projector.tensor = &projector.tensor * &x;
        cases.push((v, "projector".into()));
        let mut v = l2_child(2);
        v.ports[0] = function!(
            symbolica::symbol!(names::MINK),
            Atom::var(v.dimension),
            x.clone()
        );
        cases.push((v, "ports".into()));
        let mut v = l2_root(2, true);
        spec(&mut v).weights = vec![("scalars".into(), false, x.clone())];
        cases.push((v, "weight scalars".into()));
        let mut v = l2_root(2, true);
        v.dimension = symbolica::symbol!(name);
        cases.push((v, "dimension".into()));
        let base =
            FactorizedNumerator::from_text(include_str!("../fixtures/wh-k2-child-ins.record.txt"))
                .unwrap();
        let vertex = base.insertions[0].vertex;
        let mut v = base.clone();
        v.insertions[0].coefficient = &v.insertions[0].coefficient * &x;
        cases.push((v, format!("insertion {vertex} coefficient")));
        let mut v = base.clone();
        v.insertions[0].factor = &v.insertions[0].factor * &x;
        cases.push((v, format!("insertion {vertex} factor")));
        for (v, field) in cases {
            assert_graph_prefix(v.validate(), &reserved_prefix(name, &field), &field);
        }
        // without the defect: accepted
        l2_root(2, true).validate().unwrap();
        l2_child(2).validate().unwrap();
        base.validate().unwrap();
    }
}

/// F1 for the separately supplied prefactor ARGUMENT of every
/// entry: closure (l2_root), evaluate_factorized, the ports entry and
/// the vacuum entry (a v3 direct record): refused with the exact message;
/// the record with a clean argument passes the argument checks.
#[test]
fn the_prefactor_argument_is_checked_at_every_entry() {
    let vakint = Vakint::new().unwrap();
    let field = "prefactor argument";
    let direct =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k5-direct.record.txt"))
            .unwrap();
    assert!(super::super::is_marked(&direct));
    for name in [
        "vakint::closure_lambda",
        "vakint::closure_tensor_0",
        "vakint::dot(vakint::p(1),vakint::p(1))",
        "UFO::coupling(vakint::p(1))",
    ] {
        let x = atom(name);
        let reported = match x.as_view() {
            AtomView::Var(v) => v.get_symbol().get_name().to_string(),
            _ if name.starts_with("UFO::") => "vakint::p".to_string(),
            _ => "vakint::dot".to_string(),
        };
        let prefix = reserved_prefix(&reported, field);
        let root = l2_root(2, true);
        let closure = evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(&root),
            &[root.prefactor.total() * &x],
        )
        .map(|_| ());
        assert_graph_prefix(closure, &prefix, "closure entry");
        let argument = direct.prefactor.total() * &x;
        let entries: [(&str, std::result::Result<(), FactorizedError>); 3] = [
            (
                "evaluate_factorized",
                crate::p_integral::factorized::evaluate_factorized(
                    &vakint,
                    &settings(),
                    &direct,
                    &argument,
                )
                .map(|_| ()),
            ),
            (
                "evaluate_factorized_ports",
                crate::p_integral::factorized::evaluate_factorized_ports(
                    &vakint,
                    &settings(),
                    &direct,
                    &argument,
                )
                .map(|_| ()),
            ),
            (
                "evaluate_factorized_vacuum",
                crate::p_integral::factorized::evaluate_factorized_vacuum(
                    &vakint,
                    &settings(),
                    &direct,
                    &argument,
                )
                .map(|_| ()),
            ),
        ];
        for (entry, result) in entries {
            assert_graph_prefix(result, &prefix, entry);
        }
    }
    // without the defect the argument checks pass (marked records)
    super::super::check_prefactor_argument(&direct, &direct.prefactor.total(), vakint_epsilon())
        .unwrap();
    let root = l2_root(2, true);
    super::super::check_prefactor_argument(&root, &root.prefactor.total(), vakint_epsilon())
        .unwrap();
    // an UNMARKED value keeps the base behaviour: no argument check
    let mut v1 = direct.clone();
    v1.record = None;
    v1.insertions.clear();
    assert!(!super::super::is_marked(&v1));
    super::super::check_prefactor_argument(
        &v1,
        &atom("vakint::closure_tensor_0"),
        vakint_epsilon(),
    )
    .unwrap();
}

/// every item of the prefactor argument's checks, on the closure
/// entry's l2_root (variant AND detail prefix; the clean argument accepted).
#[test]
fn the_prefactor_argument_gets_the_record_prefactor_checks() {
    let root = l2_root(2, true);
    let clean = root.prefactor.total();
    let check = |x: &str| {
        super::super::check_prefactor_argument(&root, &(&clean * atom(x)), vakint_epsilon())
    };
    let not_scalar = |x: &str, prefix: &str| match check(x) {
        Err(FactorizedError::NotScalar { location, detail }) => {
            assert_eq!(location, "prefactor argument", "{x}");
            assert!(detail.starts_with(prefix), "{x}: {detail}");
        }
        other => panic!("{x}: expected NotScalar {prefix:?}, got {other:?}"),
    };
    // the reserved namespace and wildcards: Graph
    assert_graph_prefix(
        check("vakint::closure_lambda"),
        &reserved_prefix("vakint::closure_lambda", "prefactor argument"),
        "lambda",
    );
    assert_graph_prefix(
        check("UFO::x_"),
        "a symbol ending in _ cannot be used as prefactor argument: UFO::x_ is a pattern wildcard",
        "wildcard",
    );
    // loop momenta and indices (open colour indices included): NotScalar
    not_scalar(&q(1, &hedge(120)), "must be free of indices and momenta");
    not_scalar(
        &format!("{METRIC}({},{})", hedge(120), hedge(121)),
        "must be free of indices and momenta",
    );
    not_scalar(
        "spenso::t(spenso::coad(8,gammalooprs::a),spenso::cof(3,gammalooprs::i),spenso::dind(spenso::cof(3,gammalooprs::j)))",
        "must be free of indices and momenta",
    );
    // inexact numbers
    match check("0.5") {
        Err(FactorizedError::InexactNumber { location, .. }) => {
            assert_eq!(location, "prefactor argument")
        }
        other => panic!("0.5: {other:?}"),
    }
    // the flavour symbol (nf_power = 0), mUVexp, the reference mass R
    not_scalar(
        "beta::nf",
        "the flavour symbol beta::nf with the record's nf_power=0",
    );
    not_scalar(
        "gammalooprs::mUVexp",
        "gammalooprs::mUVexp (the producer maps it to mUV)",
    );
    not_scalar(R, "the reference mass ");
    // accepted: the clean argument, vakint's epsilon, a colour invariant and
    // the flavour factor of a record with nf_power > 0 (the auxiliary mass
    // may occur only as the taylor line's mass)
    for x in [
        "1",
        "vakint::ε",
        "spenso::cas(2,spenso::coad(8))",
        "witness::x",
    ] {
        check(x).unwrap_or_else(|e| panic!("{x}: {e}"));
    }
    let mut with_nf = root.clone();
    with_nf.prefactor.nf_power = 1;
    super::super::check_prefactor_argument(&with_nf, &with_nf.prefactor.total(), vakint_epsilon())
        .unwrap();
    super::super::check_prefactor_argument(
        &with_nf,
        &(&clean * atom("beta::nf")),
        vakint_epsilon(),
    )
    .unwrap();
}

/// THE DELTA WITNESS as acceptance tests on the REAL k1-root
/// record: with `vakint::closure_tensor_0` in prefactor_scalars or in the
/// prefactor argument it is REFUSED; with the control name
/// `user::closure_tensor_0` in the argument the value is A x (A = the
/// fixture's value) through valid_through, and with it in prefactor_scalars
/// too.
#[test]
fn the_k1_root_collision_witness_is_refused_and_its_control_is_exact() {
    let record = include_str!("../fixtures/wh-k1-root.record.txt");
    let value_text = include_str!("../fixtures/wh-k1-root.value.txt");
    let value = FactorizedNumerator::from_text(record).unwrap();
    let expected = fixture_value(value_text);
    let order = valid_through(value_text);
    let vakint = Vakint::new().unwrap();
    let run = |v: &FactorizedNumerator, argument: Atom| {
        evaluate_factorized_closure(
            &vakint,
            &rstar_site_settings(),
            std::slice::from_ref(v),
            &[argument],
        )
    };
    let internal = atom("vakint::closure_tensor_0");
    let control = atom("user::closure_tensor_0");
    let mut named = value.clone();
    named.prefactor.scalars = &named.prefactor.scalars * &internal;
    assert_graph_prefix(
        run(&named, named.prefactor.total()).map(|_| ()),
        &reserved_prefix("vakint::closure_tensor_0", "prefactor_scalars"),
        "k1-root scalars",
    );
    assert_graph_prefix(
        run(&value, value.prefactor.total() * &internal).map(|_| ()),
        &reserved_prefix("vakint::closure_tensor_0", "prefactor argument"),
        "k1-root argument",
    );
    let got = run(&value, value.prefactor.total() * &control).unwrap();
    assert_equal_through(
        &got.value,
        &(&expected * &control),
        order,
        "k1-root x user:: argument",
    );
    let mut renamed = value.clone();
    renamed.prefactor.scalars = &renamed.prefactor.scalars * &control;
    let got = run(&renamed, renamed.prefactor.total()).unwrap();
    assert_equal_through(
        &got.value,
        &(&expected * &control),
        order,
        "k1-root x user:: scalars",
    );
}

/// a symbol whose name ends in `_` is a symbolica WILDCARD when an
/// atom is used as a replacement pattern (witness: `replace(R_)` rewrites a
/// sum that does not contain R_), so it is refused in a marked record and as
/// the settings' epsilon; an unmarked value keeps the base behaviour.
#[test]
fn wildcard_names_are_refused() {
    let wild = atom("UFO::R_");
    let target = atom("UFO::a+UFO::b");
    let replaced = target.replace(wild.clone()).with(Atom::num(1));
    eprintln!("WILDCARD witness: (a+b).replace(R_).with(1) = {replaced}");
    assert_ne!(replaced, target, "R_ did not act as a wildcard");
    // the reference mass (closing edge and closure line): the first atom met
    // is the closing edge's mass
    let mut v = l2_root(2, true);
    v.edges[2].mass_squared = wild.clone();
    spec(&mut v).mass_squared = wild.clone();
    assert_graph_prefix(
        v.validate(),
        "a symbol ending in _ cannot be used as edge 3 mass: UFO::R_ is a pattern wildcard",
        "R_",
    );
    let mut v = l2_root(2, true);
    spec(&mut v).auxiliary_mass_squared = wild.clone();
    assert_graph_prefix(
        v.validate(),
        "a symbol ending in _ cannot be used as taylor auxiliary mass",
        "M_aux",
    );
    let mut v = l2_root(2, true);
    v.dimension = symbolica::symbol!("UFO::D_");
    assert_graph_prefix(
        v.validate(),
        "a symbol ending in _ cannot be used as dimension: UFO::D_ is a pattern wildcard",
        "D_",
    );
    let mut v = l2_root(2, true);
    v.vertices[1].factor = &v.vertices[1].factor * atom("UFO::g_");
    assert_graph_prefix(
        v.validate(),
        "a symbol ending in _ cannot be used as vertex 2 factor",
        "factor",
    );
    l2_root(2, true).validate().unwrap();
    // the settings' epsilon
    let vakint = Vakint::new().unwrap();
    let root = l2_root(0, true);
    for (epsilon, prefix) in [
        (
            "UFO::eps_",
            "a symbol ending in _ cannot be used as the epsilon_symbol setting",
        ),
        (
            "closure_tensor_0",
            "reserved symbol vakint::closure_tensor_0 in the epsilon_symbol setting",
        ),
    ] {
        let settings = VakintSettings {
            epsilon_symbol: epsilon.into(),
            ..settings()
        };
        assert_graph_prefix(
            evaluate_factorized_closure(
                &vakint,
                &settings,
                std::slice::from_ref(&root),
                &[root.prefactor.total()],
            )
            .map(|_| ()),
            prefix,
            epsilon,
        );
    }
    // an unmarked (v1 direct) value is not checked: the base behaviour
    let mut v1 =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k5-direct.record.txt"))
            .unwrap();
    v1.record = None;
    v1.insertions.clear();
    v1.prefactor.scalars = atom("UFO::x_");
    assert!(!matches!(
        v1.validate(),
        Err(FactorizedError::Graph { detail }) if detail.starts_with("a symbol ending in _")
    ));
}

/// Rename edge `old` to `new` everywhere in a record (ids, loop / external
/// edges, the closing edge, every `Q(old, ..)`).
fn renumber(value: &FactorizedNumerator, old: usize, new: usize) -> FactorizedNumerator {
    let mut v = value.clone();
    let q_atom = atom(&q(1, &hedge(1)));
    let AtomView::Fun(q_fun) = q_atom.as_view() else {
        unreachable!()
    };
    let q_symbol = q_fun.get_symbol();
    let rename = |a: &Atom| {
        a.replace_map(|part, _, out| {
            if let AtomView::Fun(f) = part
                && f.get_symbol() == q_symbol
                && f.get_nargs() == 2
            {
                let args = f.iter().collect::<Vec<_>>();
                if i64::try_from(args[0]).ok() == Some(old as i64) {
                    **out = function!(q_symbol, Atom::num(new as i64), args[1].to_owned());
                }
            }
        })
    };
    for e in &mut v.edges {
        if e.id == old {
            e.id = new;
        }
        e.factor = rename(&e.factor);
    }
    for vx in &mut v.vertices {
        vx.factor = rename(&vx.factor);
    }
    if let Some(p) = &mut v.projector {
        p.tensor = rename(&p.tensor);
    }
    for slot in v.loop_edges.iter_mut().chain(v.external_edges.iter_mut()) {
        if *slot == Some(old) {
            *slot = Some(new);
        }
    }
    if let FactorizedKind::Closure(spec) = &mut v.kind {
        if spec.edge == old {
            spec.edge = new;
        }
        if let Some(taylor) = &mut spec.taylor {
            for e in &mut taylor.scaled_edges {
                if *e == old {
                    *e = new;
                }
            }
        }
    }
    v.edges.sort_by_key(|e| e.id);
    v
}

/// edge ids are written as `Q(id, ..)` (i64) and a loop without a loop
/// edge gets the fresh id max + 1 + i: an id near usize::MAX (which used to
/// wrap the fresh id onto an existing edge) is refused, an id whose fresh ids
/// stay at or below i64::MAX is accepted.
#[test]
fn edge_ids_beyond_i64_are_refused() {
    // l2_root with loop 0 without a loop edge: line 1 (k1 + q) carries loop 0
    let mut base = l2_root(0, true);
    base.loop_edges[0] = None;
    base.validate().unwrap();
    let fresh = base.loop_edges.iter().filter(|e| e.is_none()).count();
    assert_eq!(fresh, 1);
    let prefix = |id: usize| {
        format!(
            "edge id {id}: the ids (with fresh loop-vector ids up to {id} + 1 above them) must not exceed i64::MAX"
        )
    };
    for id in [usize::MAX, i64::MAX as usize] {
        match renumber(&base, 21, id).validate() {
            Err(FactorizedError::Closure { detail }) => {
                assert!(detail.starts_with(&prefix(id)), "{id}: {detail}")
            }
            other => panic!("{id}: {other:?}"),
        }
    }
    renumber(&base, 21, i64::MAX as usize - 1)
        .validate()
        .unwrap();
}

/// a marked value's edge power is at most 2^31 in absolute
/// value (powers are summed when lines are merged in series): a power of
/// i64::MAX / 2 and -(2^31 + 1) are refused by validate ("edge 1: power ..
/// out of range"), 2^31 and 3 are accepted with a Taylor step of degree 2. With
/// that bound and a degree of at most 2 the Taylor coefficients cannot
/// overflow (the bound function's own limit at degree 2 is about 3.04e9 >
/// 2^31); the bound function stays checked on its own.
#[test]
fn powers_that_overflow_the_taylor_coefficients_are_refused() {
    let huge = i64::MAX / 2;
    for power in [huge, -(1i64 << 31) - 1] {
        let mut v = l2_root(2, true);
        v.edges[0].power = power;
        match v.validate() {
            Err(FactorizedError::Closure { detail }) => assert!(
                detail.starts_with(&format!(
                    "edge 1: power {power} out of range (at most 2^31 in absolute value)"
                )),
                "{detail}"
            ),
            other => panic!("{power}: {other:?}"),
        }
    }
    for power in [1i64 << 31, 3] {
        let mut v = l2_root(2, true);
        v.edges[0].power = power;
        v.validate().unwrap_or_else(|e| panic!("{power}: {e}"));
    }
    // the bound function itself
    assert!(super::taylor_coefficient_bound(3_037_000_499, 2).is_some());
    assert!(super::taylor_coefficient_bound(i64::MAX, 1).is_none());
    assert!(super::taylor_coefficient_bound(huge, 2).is_none());
}

/// a closure-kind value's VERTEX ids are i64 too (they are written as
/// edge(source, target)): a vertex renumbered to usize::MAX or i64::MAX + 1 is
/// refused by validate ("id <n>: the edge and vertex ids"), i64::MAX is
/// accepted (no fresh ids are counted on vertices).
#[test]
fn closure_vertex_ids_beyond_i64_are_refused() {
    let base = l2_root(0, true);
    base.validate().unwrap();
    let renumbered = |id: usize| {
        let mut v = base.clone();
        // the last vertex (ids stay sorted)
        let last = v.vertices.len() - 1;
        let old = v.vertices[last].id;
        v.vertices[last].id = id;
        for e in &mut v.edges {
            for end in [&mut e.source, &mut e.target] {
                if *end == Endpoint::Vertex(old) {
                    *end = Endpoint::Vertex(id);
                }
            }
        }
        v
    };
    for id in [usize::MAX, i64::MAX as usize + 1] {
        match renumbered(id).validate() {
            Err(FactorizedError::Closure { detail }) => assert!(
                detail.starts_with(&format!("id {id}: the edge and vertex ids")),
                "{id}: {detail}"
            ),
            other => panic!("{id}: {other:?}"),
        }
    }
    renumbered(i64::MAX as usize)
        .validate()
        .unwrap_or_else(|e| panic!("i64::MAX: {e}"));
}

/// a degree read off an exponent saturates (a degree is only compared
/// with w <= 2) instead of truncating `x^(2^32)` to order 0.
#[test]
fn degree_exponents_saturate() {
    let exponent = |text: &str| {
        let pow = atom(&format!("UFO::x^({text})"));
        let AtomView::Pow(p) = pow.as_view() else {
            panic!("{pow}")
        };
        super::degree_exponent(p.get_base_exp().1)
    };
    assert_eq!(exponent("2"), 2);
    assert_eq!(exponent("4294967296"), u32::MAX);
    assert_eq!(exponent("4294967298"), u32::MAX);
    assert_eq!(exponent("100000000000000000000000"), u32::MAX);
    assert_eq!(exponent("-2"), 0);
}

/// the insertions' cycles are summed with checked arithmetic (two
/// insertions of u32::MAX and 1 fermion cycles used to wrap to 0).
#[test]
fn insertion_cycles_do_not_wrap() {
    let base =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    let mut a = base.insertions[0].clone();
    let mut b = a.clone();
    a.cycles = (u32::MAX, 0);
    b.cycles = (1, 0);
    assert_eq!(super::insertion_cycles(&[a.clone(), b.clone()]), None);
    b.cycles = (0, 1);
    assert_eq!(super::insertion_cycles(&[a, b]), Some((u32::MAX, 1)));
}

/// Raw port indices `PORT + i` (9001, 9002, ..) cannot meet a
/// caller's integer index. (a) A caller's OPEN slot cannot be an integer: a
/// port is a half-edge slot (`IndexNotLocal` otherwise, asserted here). (b) A
/// caller's integer DUMMY is contracted away by idenso before the mapping;
/// the mapping goes BY POSITION of the open slot, and the way back rewrites
/// only vakint's `g` / `p` heads, which a marked record cannot carry.
/// Witness of (b): l2_child with its dummy renamed to the integer 9001 (the
/// first raw port index, while port 0 maps to 9001) gives the same value.
#[test]
fn integer_indices_in_the_port_range_do_not_collide() {
    let vakint = Vakint::new().unwrap();
    let original = l2_child(1);
    let slot = |n: i64| atom(&format!("{MINK}({DIM},{n})"));
    // (a) an integer port is refused
    let mut port = original.clone();
    for v in &mut port.vertices {
        v.factor = v.factor.replace(atom(&hedge(120))).with(slot(9002));
    }
    port.ports[0] = slot(9002);
    assert!(
        matches!(port.validate(), Err(FactorizedError::IndexNotLocal { .. })),
        "{:?}",
        port.validate()
    );
    // (b) an integer dummy at the raw port base
    let mut dummy = original.clone();
    for v in &mut dummy.vertices {
        v.factor = v.factor.replace(atom(&index(2, 7))).with(slot(9001));
    }
    assert_ne!(dummy, original);
    dummy.validate().unwrap();
    let run = |v: &FactorizedNumerator| {
        evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(v),
            &[v.prefactor.total()],
        )
        .unwrap()
        .value
    };
    let expected = run(&original);
    let got = run(&dummy);
    assert!(!expected.is_zero());
    assert!(
        (&got - &expected).expand().is_zero(),
        "integer dummy 9001\n got      {got}\n expected {expected}"
    );
}

/// a FORGED `vakint::slash` dummy: as a bare symbol it is reserved; in
/// the reader's slot form it is refused by the index rules (it can only meet
/// a reader dummy), never evaluated.
#[test]
fn a_forged_slash_dummy_is_refused() {
    let mut v = l2_root(2, true);
    v.vertices[1].factor = &v.vertices[1].factor * atom("vakint::slash(-1,0)");
    assert_graph_prefix(
        v.validate(),
        &reserved_prefix("vakint::slash", "vertex 2 factor"),
        "bare slash",
    );
    let slash = format!("{MINK}({DIM},vakint::slash(2,0))");
    let mut v = l2_root(2, true);
    v.vertices[1].factor = &v.vertices[1].factor
        * atom(&format!("{}*{}", q(2, &slash), q(1, &slash)))
        * atom(&format!("{}*{}", q(3, &slash), q(1, &slash)));
    assert!(
        matches!(
            v.validate(),
            Err(FactorizedError::IndexOverContracted { .. })
                | Err(FactorizedError::DummyIndexShared { .. })
        ),
        "{:?}",
        v.validate()
    );
}

/// The entry against the INDEPENDENT reference, exactly, one record per
/// kind: root (projector, Taylor step, compensated), child (open ports,
/// independent externals), pre-expanded sector (no Taylor step), and the
/// order filters; the argument is ε-free (so that multiplying it after the
/// evaluation is exact).
fn check_independent(value: &FactorizedNumerator, routing: &[(&str, &str)], argument: &Atom) {
    value
        .validate()
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    let vakint = Vakint::new().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &settings(),
        std::slice::from_ref(value),
        std::slice::from_ref(argument),
    )
    .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    let expected = independent_reference(value, routing, argument);
    assert!(!expected.is_zero(), "{}: zero reference", value.label);
    assert!(
        (&got.records[0].vakint_value - &expected)
            .expand()
            .is_zero(),
        "{}: entry vs the independent reference\nentry     {}\nreference {}",
        value.label,
        got.records[0].vakint_value,
        expected
    );
}

#[test]
fn every_kind_equals_the_independent_reference() {
    let argument = atom("UFO::gs^2*witness::x");
    check_independent(&l2_root(2, true), &L2_CUT, &argument);
    check_independent(&l2_root(1, false), &L2_CUT, &argument);
    check_independent(&l2_child(2), &L2_CUT, &argument);
    check_independent(&l2_sector(), &L2_CUT, &argument);
    // the order filters BELOW the validation (the public path refuses
    // the non-`all` filters on the stacked branch; the retained evaluation
    // code stays covered by the independent reference)
    for orders in [TaylorOrders::DropZero, TaylorOrders::DropZeroOdd] {
        let mut v = l2_root(2, true);
        spec(&mut v).taylor.as_mut().unwrap().orders = orders;
        check_independent_below_validation(&v, &L2_CUT, &argument);
    }
}

/// [`check_independent`] through the evaluation BELOW the validation
/// (`evaluate_validated`), for record shapes the validation refuses.
fn check_independent_below_validation(
    value: &FactorizedNumerator,
    routing: &[(&str, &str)],
    argument: &Atom,
) {
    let vakint = Vakint::new().unwrap();
    let FactorizedKind::Closure(spec) = &value.kind else {
        unreachable!()
    };
    let got = super::evaluate_validated(&vakint, &settings(), value, spec, argument)
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    let expected = independent_reference(value, routing, argument);
    assert!(!expected.is_zero(), "{}: zero reference", value.label);
    assert!(
        (&got.vakint_value - &expected).expand().is_zero(),
        "{}: below validation vs the independent reference\nentry     {}\nreference {}",
        value.label,
        got.vakint_value,
        expected
    );
}

/// THE COLLISION WITNESSES against the independent reference: the control
/// names `user::closure_tensor_<n>` (what a caller may legally write) in the
/// prefactor argument and inside a factor, closed (root) and with open ports
/// (child): the entry equals the reference exactly (the internal names are
/// refused: `internal_symbol_names_are_refused_in_a_record`).
#[test]
fn the_collision_controls_equal_the_independent_reference() {
    for (value, name) in [
        (l2_root(2, true), "user::closure_tensor_0"),
        (l2_child(2), "user::closure_tensor_1"),
    ] {
        let x = atom(name);
        check_independent(&value, &L2_CUT, &(atom("UFO::gs^2") * &x));
        let mut inside = value.clone();
        inside.vertices[0].factor = &inside.vertices[0].factor * &x;
        check_independent(&inside, &L2_CUT, &atom("UFO::gs^2"));
    }
}

// ---------------------------------------------------------------------------
// REFUSAL TESTS THAT CAN FAIL: `validate` of a MINIMAL input (a valid
// record with one defect) gives the variant AND the detail prefix of the
// REFUSALS list, and the same input without the defect is ACCEPTED.
// Priority (2): the refusals the R* emitter can reach.
// ---------------------------------------------------------------------------

/// `bad.validate()` = `Closure { detail }` with `detail` starting with
/// `prefix`; `good.validate()` = Ok.
fn assert_closure_refusal(bad: &FactorizedNumerator, good: &FactorizedNumerator, prefix: &str) {
    good.validate()
        .unwrap_or_else(|e| panic!("the twin without the defect is refused: {e}"));
    match bad.validate() {
        Err(FactorizedError::Closure { detail }) => assert!(
            detail.starts_with(prefix),
            "{}: {detail}\n(expected prefix {prefix:?})",
            bad.label
        ),
        other => panic!("{}: expected Closure {prefix:?}, got {other:?}", bad.label),
    }
}

/// A closure record with no loop (an inserted counterterm alone).
fn zero_loop_record() -> FactorizedNumerator {
    record(
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
    )
}

#[test]
fn refusal_zero_loops() {
    assert_closure_refusal(
        &zero_loop_record(),
        &l2_root(2, true),
        "a closure record has at least one loop",
    );
}

#[test]
fn refusal_taylor_degree_above_two() {
    let good = l2_root(2, true);
    let mut bad = good.clone();
    spec(&mut bad).taylor.as_mut().unwrap().degree = 3;
    assert_closure_refusal(&bad, &good, "Taylor degree 3 above 2");
}

/// The order filter that keeps nothing (by a UNIT test of the function
/// the rule lives in, `orders_keep_any`, so that it stays covered when the
/// non-`all` filters are refused earlier in validate_closure): drop0-odd at
/// degree 1 keeps nothing, at degree 2 order 2; drop0 at degree 0 nothing,
/// at 1 order 1; all always order 0.
#[test]
fn refusal_orders_keep_nothing() {
    let taylor = |degree: u32, orders: TaylorOrders| TaylorSpec {
        degree,
        scaled_edges: vec![20],
        compensate_mass: true,
        orders,
    };
    for (degree, orders, keeps) in [
        (1, TaylorOrders::DropZeroOdd, false),
        (2, TaylorOrders::DropZeroOdd, true),
        (0, TaylorOrders::DropZero, false),
        (1, TaylorOrders::DropZero, true),
        (0, TaylorOrders::All, true),
    ] {
        assert_eq!(
            super::orders_keep_any(&taylor(degree, orders)),
            keeps,
            "degree {degree} orders={}",
            orders.as_str()
        );
    }
}

/// An open spinor chain: a projector leaving two bispinor slots open (on
/// the pre-expanded sector: without a Taylor step the projector may carry
/// the leg momentum); the twin closes the root with the metric alone.
#[test]
fn refusal_open_spinor_slot() {
    let good = l2_sector();
    let mut bad = good.clone();
    let bis = |k: usize| {
        format!("spenso::{{spenso::representation,spenso::self_dual}}::bis(4,{HEDGE}(120,{k}))")
    };
    let projector = bad.projector.as_mut().unwrap();
    projector.tensor = &projector.tensor
        * atom(&format!(
            "spenso::{{linear,spenso::tensor}}::gamma({},{},{MINK}({DIM},{HEDGE}(120,4)))",
            bis(2),
            bis(3)
        ))
        * atom(&q(20, &format!("{MINK}({DIM},{HEDGE}(120,4))")));
    match bad.validate() {
        Err(FactorizedError::Closure { detail }) => assert!(
            detail.starts_with("open slot ") && detail.contains("is not a Lorentz slot"),
            "{detail}"
        ),
        other => panic!("expected the open-slot refusal, got {other:?}"),
    }
    good.validate().unwrap();
}

/// "every line is massive": a one-loop record whose only internal line is
/// the closing line (the back end would be a single-scale family, not a
/// closure); refused by validate (R-15: not later at the evaluation). The
/// twin: the same loop with a massless line beside the closing one.
#[test]
fn refusal_every_line_is_massive() {
    let bad = record(
        "one loop, the closing line alone",
        1,
        vec![Some(20)],
        vec![
            line(3, 1, 1, &[1, 0], R),
            leg(20, 1, true, &[0, 1]),
            leg(21, 1, false, &[0, 1]),
        ],
        vec![vertex(1, "")],
        Vec::new(),
        None,
        3,
        taylor(0, vec![20, 21], false),
        false,
    );
    let good = record(
        "one loop, closing line and a massless line",
        1,
        vec![Some(20)],
        vec![
            line(1, 1, 2, &[1, 1], "0"),
            line(3, 2, 1, &[1, 0], R),
            leg(20, 1, true, &[0, 1]),
            leg(21, 2, false, &[0, 1]),
        ],
        vec![vertex(1, ""), vertex(2, "")],
        Vec::new(),
        None,
        3,
        taylor(0, vec![20, 21], false),
        false,
    );
    assert_closure_refusal(&bad, &good, "every line is massive: not a closure");
}

/// A compact slash `gamma(a, b, Q(x, mink(D)))` in the projector of a record
/// with a Taylor step (the reader's rule R4, `Ports`); without the Taylor
/// step the same projector is written out.
#[test]
fn refusal_compact_slash_under_a_taylor_step() {
    let mut bad = l2_root(1, true);
    let slash = atom(&format!(
        "spenso::{{linear,spenso::tensor}}::gamma(spenso::{{spenso::representation,spenso::self_dual}}::bis(4,{HEDGE}(120,2)),spenso::{{spenso::representation,spenso::self_dual}}::bis(4,{HEDGE}(121,2)),{MOMENTUM}(20,{MINK}({DIM})))"
    ));
    let projector = bad.projector.as_mut().unwrap();
    projector.tensor = &projector.tensor * &slash;
    let mut good = bad.clone();
    spec(&mut good).taylor = None;
    good.expand_compact_slashes().unwrap();
    assert_ne!(
        good.projector, bad.projector,
        "the twin's slash was not written out"
    );
    match bad.expand_compact_slashes() {
        Err(FactorizedError::Ports { detail }) => assert!(
            detail.starts_with(
                "a compact slash gamma(a, b, Q(x, mink(D))) in the projector of a record with a Taylor step"
            ),
            "{detail}"
        ),
        other => panic!("expected the compact-slash refusal, got {other:?}"),
    }
}

// Priority (3): every other designed refusal of the closure kind, the same
// rule (variant + detail prefix on a minimal input, the twin accepted).

/// `bad.validate()` = `Closure` whose detail starts with `prefix` and
/// contains `contains`; `good.validate()` = Ok.
fn assert_closure_refusal_with(
    bad: &FactorizedNumerator,
    good: &FactorizedNumerator,
    prefix: &str,
    contains: &str,
) {
    good.validate()
        .unwrap_or_else(|e| panic!("the twin without the defect is refused: {e}"));
    match bad.validate() {
        Err(FactorizedError::Closure { detail }) => assert!(
            detail.starts_with(prefix) && detail.contains(contains),
            "{}: {detail}\n(expected prefix {prefix:?} containing {contains:?})",
            bad.label
        ),
        other => panic!("{}: expected Closure {prefix:?}, got {other:?}", bad.label),
    }
}

/// The record with one change.
fn with(
    base: &FactorizedNumerator,
    change: impl FnOnce(&mut FactorizedNumerator),
) -> FactorizedNumerator {
    let mut v = base.clone();
    change(&mut v);
    v
}

#[test]
fn refusal_original_cycles_weights() {
    let good = l2_root(1, true);
    let bad = with(&good, |v| v.prefactor.loop_sign = -1);
    assert_closure_refusal(
        &bad,
        &good,
        "with original_cycles the prefactor has loop_sign=1, nf_power=0, overall 1",
    );
}

#[test]
fn refusal_more_cycles_than_original_cycles() {
    let value =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    let good = with(&value, |v| {
        spec(v).original_cycles = Some((1, 1));
        v.insertions[0].cycles = (1, 0);
    });
    let bad = with(&good, |v| v.insertions[0].cycles = (2, 0));
    assert_closure_refusal(&bad, &good, "the carried graph and its insertions have ");
}

/// Without original_cycles the record's own loop sign is checked (variant
/// `LoopSign`, the contract's).
#[test]
fn refusal_loop_sign_without_original_cycles() {
    let good = with(&l2_root(1, true), |v| spec(v).original_cycles = None);
    good.validate().unwrap();
    let bad = with(&good, |v| v.prefactor.loop_sign = -1);
    assert!(
        matches!(
            bad.validate(),
            Err(FactorizedError::LoopSign {
                found: -1,
                loops: 0
            })
        ),
        "{:?}",
        bad.validate()
    );
}

#[test]
fn refusal_weight_applied_by_vakint_differs() {
    let base = l2_root(1, true);
    let scalars = base.prefactor.scalars.clone();
    let good = with(&base, |v| {
        spec(v).weights = vec![("scalars".into(), true, scalars.clone())]
    });
    let bad = with(&base, |v| {
        spec(v).weights = vec![("scalars".into(), true, &scalars * atom("UFO::two"))]
    });
    assert_closure_refusal(&bad, &good, "weight scalars applied_by=vakint ");
}

#[test]
fn refusal_flavour_symbol_in_the_record() {
    let good = l2_root(1, true);
    let bad = with(&good, |v| {
        v.prefactor.scalars = &v.prefactor.scalars * atom("beta::nf")
    });
    assert_closure_refusal(
        &bad,
        &good,
        "a flavour symbol in the record (n_f stays with the caller)",
    );
}

#[test]
fn refusal_muvexp_in_the_record() {
    let good = l2_root(1, true);
    let bad = with(&good, |v| {
        v.vertices[0].factor = &v.vertices[0].factor * atom("gammalooprs::mUVexp")
    });
    assert_closure_refusal(
        &bad,
        &good,
        "mUVexp in the record (the producer maps it to mUV)",
    );
}

#[test]
fn refusal_closing_edge_missing_external_massless_power() {
    let good = l2_root(2, true);
    assert_closure_refusal(
        &with(&good, |v| spec(v).edge = 99),
        &good,
        "closing edge 99 does not exist",
    );
    assert_closure_refusal(
        &with(&good, |v| spec(v).edge = 20),
        &good,
        "closing edge 20 is an external leg",
    );
    assert_closure_refusal(
        &with(&good, |v| spec(v).edge = 1),
        &good,
        "closing edge 1: mass 0 is zero or differs from the closure line's",
    );
    assert_closure_refusal(
        &with(&good, |v| spec(v).power = 2),
        &good,
        "closing edge 3: power 1 differs from the closure line's 2",
    );
}

#[test]
fn refusal_several_massive_lines() {
    let good = l2_root(2, true);
    let bad = with(&good, |v| v.edges[1].mass_squared = atom(R));
    assert_closure_refusal(&bad, &good, "several massive lines ([2, 3])");
}

#[test]
fn refusal_auxiliary_mass_zero() {
    let good = l2_root(2, true);
    let bad = with(&good, |v| spec(v).auxiliary_mass_squared = Atom::Zero);
    assert_closure_refusal(&bad, &good, "the auxiliary mass is zero");
}

#[test]
fn refusal_reference_mass_in_the_prefactor() {
    let good = l2_root(1, true);
    let bad = with(&good, |v| {
        v.prefactor.scalars = &v.prefactor.scalars * atom(R)
    });
    assert_closure_refusal(
        &bad,
        &good,
        &format!(
            "the reference mass {} in the prefactor (it is substituted in the integrand only)",
            atom(R)
        ),
    );
}

#[test]
fn refusal_momentum_in_the_projector() {
    let good = l2_root(1, true);
    let bad = with(&good, |v| {
        v.projector.as_mut().unwrap().tensor =
            atom(&format!("{}*{}", q(20, &hedge(120)), q(20, &hedge(121))))
    });
    assert_closure_refusal(
        &bad,
        &good,
        "the projector contains Q(20, ..): with a Taylor step",
    );
    // without a Taylor step, an internal edge's momentum: the contract's
    // rule refuses a LOOP momentum first (the closure's "the projector
    // contains the momentum of internal edge" is reachable only for an
    // internal edge without loop momentum: a listed limit, by reasoning)
    let good = l2_sector();
    let bad = with(&good, |v| {
        let p = v.projector.as_mut().unwrap();
        p.tensor = &p.tensor * atom(&format!("{}*{}", q(1, &index(9, 1)), q(1, &index(9, 1))))
    });
    good.validate().unwrap();
    match bad.validate() {
        Err(FactorizedError::UnsupportedTensor { location, tensor }) => {
            assert_eq!(location, "projector");
            assert!(
                tensor.ends_with("(loop momentum in the projector)"),
                "{tensor}"
            );
        }
        other => panic!("{other:?}"),
    }
}

/// p2_power on a record without an external momentum (a vacuum record).
#[test]
fn refusal_p2_power_without_an_external_momentum() {
    let vacuum = |p2_power: i64| {
        record(
            "vacuum sunrise",
            2,
            Vec::new(),
            vec![
                line(1, 1, 2, &[1, 0], "0"),
                line(2, 1, 2, &[-1, 1], "0"),
                line(3, 2, 1, &[0, 1], R),
            ],
            vec![vertex(1, ""), vertex(2, "")],
            Vec::new(),
            Some(Projector {
                tensor: Atom::num(1),
                p2_power,
            }),
            3,
            None,
            false,
        )
    };
    let mut good = vacuum(0);
    spec(&mut good).auxiliary_mass_squared = atom(R);
    let mut bad = vacuum(-1);
    spec(&mut bad).auxiliary_mass_squared = atom(R);
    assert_closure_refusal(&bad, &good, "p2_power without an external momentum");
}

#[test]
fn refusal_soft_line_mismatch() {
    let soft = |edge: usize| super::super::SoftLeg {
        external: 1,
        edge,
        root: false,
    };
    let good = with(&l2_root(1, true), |v| spec(v).softs = vec![soft(20)]);
    let bad = with(&good, |v| spec(v).softs = vec![soft(21)]);
    assert_closure_refusal(&bad, &good, "soft line 1 edge 21: not the leg of p(1)");
}

#[test]
fn refusal_pre_expanded_sector_with_an_auxiliary_mass() {
    let good = l2_sector();
    let bad = with(&good, |v| spec(v).auxiliary_mass_squared = atom(MAUX));
    assert_closure_refusal(
        &bad,
        &good,
        "a pre-expanded sector (taylor none) has no auxiliary mass",
    );
}

#[test]
fn refusal_pre_expanded_sector_with_an_external_momentum_inside() {
    let good = l2_sector();
    let bad = with(&l2_root(1, true), |v| {
        spec(v).taylor = None;
        spec(v).auxiliary_mass_squared = atom(R);
    });
    assert_closure_refusal(
        &bad,
        &good,
        "taylor none (lambda^0 denominators): internal edge 1 carries an external momentum",
    );
}

#[test]
fn refusal_scaled_edges() {
    let good = l2_root(1, true);
    assert_closure_refusal(
        &with(&good, |v| {
            spec(v).taylor.as_mut().unwrap().scaled_edges = vec![1, 20]
        }),
        &good,
        "scaled_edges: 1 is not an external leg",
    );
    assert_closure_refusal(
        &with(&good, |v| {
            spec(v).taylor.as_mut().unwrap().scaled_edges = vec![20, 20]
        }),
        &good,
        "scaled_edges: 20 twice",
    );
    let good = l2_child(1);
    assert_closure_refusal_with(
        &with(&good, |v| {
            spec(v).taylor.as_mut().unwrap().scaled_edges = vec![20]
        }),
        &good,
        "scaled_edges=[20] scale p",
        "must be scaled",
    );
}

#[test]
fn refusal_constant_vector_under_a_taylor_step() {
    let good = l2_root(1, true);
    let bad = with(&good, |v| {
        v.vertices[0].factor =
            &v.vertices[0].factor * atom(&format!("{}*{}", q(99, &index(1, 9)), q(1, &index(1, 9))))
    });
    assert_closure_refusal(
        &bad,
        &good,
        "constant external vectors [99] (edges outside the graph) with a Taylor step",
    );
}

#[test]
fn refusal_reference_mass_not_a_symbol() {
    let good = l2_root(1, true);
    let twice = atom(&format!("2*{R}"));
    let bad = with(&good, |v| {
        v.edges[2].mass_squared = twice.clone();
        spec(v).mass_squared = twice.clone();
    });
    assert_closure_refusal(
        &bad,
        &good,
        &format!("the reference mass {twice} is not a symbol (it is substituted everywhere)"),
    );
}

/// The entry's own refusals of its arguments (before any record is read).
#[test]
fn refusal_entry_arguments() {
    let vakint = Vakint::new().unwrap();
    let v = l2_root(0, true);
    let run = |records: &[FactorizedNumerator], prefactors: &[Atom]| {
        evaluate_factorized_closure(&vakint, &settings(), records, prefactors).map(|_| ())
    };
    run(std::slice::from_ref(&v), &[v.prefactor.total()]).unwrap();
    match run(std::slice::from_ref(&v), &[]) {
        Err(FactorizedError::Closure { detail }) => {
            assert!(detail.starts_with("1 records and 0 prefactors"), "{detail}")
        }
        other => panic!("{other:?}"),
    }
    let direct = with(&v, |v| v.kind = FactorizedKind::Direct);
    match run(std::slice::from_ref(&direct), &[v.prefactor.total()]) {
        Err(FactorizedError::Closure { detail }) => assert!(
            detail.starts_with(&format!(
                "{}: not a closure record (kind direct",
                direct.label
            )),
            "{detail}"
        ),
        other => panic!("{other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Delta review of 4f3f54110.
// ---------------------------------------------------------------------------

/// a BARE name in a marked record text is read into `vakint::` and
/// refused with the advice to qualify it.
#[test]
fn a_bare_name_in_a_marked_record_is_refused_with_advice() {
    let record = include_str!("../fixtures/wh-k1-root.record.txt");
    let line = record
        .lines()
        .find(|l| l.starts_with("prefactor_scalars\t"))
        .unwrap();
    let bare = record.replace(line, &format!("{line}*gs"));
    assert_ne!(bare, record);
    match FactorizedNumerator::from_text(&bare) {
        Err(FactorizedError::Graph { detail }) => assert_eq!(
            detail,
            "reserved symbol vakint::gs in prefactor_scalars: the vakint:: namespace is reserved for vakint's internal symbols; a name written without a namespace is read into vakint::"
        ),
        other => panic!("{other:?}"),
    }
    // qualified: accepted
    let qualified = record.replace(line, &format!("{line}*UFO::gs"));
    FactorizedNumerator::from_text(&qualified).unwrap();
}

/// an UNMARKED v1 record with bare names parses, validates and
/// evaluates as on the base (bare names are read into `vakint::` there too):
/// its value is the qualified record's value with the symbol renamed.
#[test]
fn an_unmarked_v1_record_with_bare_names_is_unchanged() {
    let qualified_text = include_str!("../fixtures/l1-gluon-bubble.txt");
    assert!(qualified_text.contains("UFO::{}::gs_minus"));
    let bare_text = qualified_text.replace("UFO::{}::gs_minus", "gs_minus");
    let bare = FactorizedNumerator::from_text(&bare_text).unwrap();
    assert!(!super::super::is_marked(&bare));
    bare.validate().unwrap();
    let qualified = FactorizedNumerator::from_text(qualified_text).unwrap();
    let vakint = Vakint::new().unwrap();
    let run = |v: &FactorizedNumerator| {
        super::super::evaluate_factorized(&vakint, &settings(), v, &v.prefactor.total())
            .unwrap()
            .value
    };
    let got = run(&bare);
    let expected = run(&qualified)
        .replace(atom("UFO::gs_minus"))
        .with(atom("vakint::gs_minus"));
    assert!(!expected.is_zero());
    assert!(
        (&got - &expected).expand().is_zero(),
        "bare {got}\nqualified renamed {expected}"
    );
}

/// `vakint::ε` in a marked record is exempt ONLY when it is the
/// settings' epsilon; with another epsilon it is refused (record and
/// argument), with vakint's own it is accepted.
#[test]
fn vakint_epsilon_is_exempt_only_as_the_settings_epsilon() {
    let vakint = Vakint::new().unwrap();
    let mut v = l2_root(1, true);
    v.vertices[0].factor = &v.vertices[0].factor * atom("1+vakint::ε");
    v.validate().unwrap();
    let run = |settings: &VakintSettings| {
        evaluate_factorized_closure(
            &vakint,
            settings,
            std::slice::from_ref(&v),
            &[v.prefactor.total()],
        )
        .map(|_| ())
    };
    assert_graph_prefix(
        run(&rstar_site_settings()),
        "reserved symbol vakint::ε in vertex 1 factor: the settings' epsilon is gammalooprs::ε",
        "record, other epsilon",
    );
    run(&settings()).unwrap();
    // the argument
    let root = l2_root(1, true);
    let argument = root.prefactor.total() * atom("1+vakint::ε");
    let other = crate::utils::vakint_macros::vk_symbol!("gammalooprs::ε");
    assert_graph_prefix(
        super::super::check_prefactor_argument(&root, &argument, other),
        &reserved_prefix("vakint::ε", "prefactor argument"),
        "argument, other epsilon",
    );
    super::super::check_prefactor_argument(&root, &argument, vakint_epsilon()).unwrap();
}

/// a FORMAT defect found by the lines' back end inside validate (a
/// closing line whose loop momentum has no unit coefficient) is a `Closure`
/// refusal, not `Evaluation`.
#[test]
fn a_closing_line_without_a_unit_coefficient_is_a_closure_refusal() {
    let good = l2_root(2, true);
    let mut bad = good.clone();
    for e in &mut bad.edges {
        match e.id {
            2 => e.momentum = vec![-1, 2, 0],
            3 => e.momentum = vec![0, 2, 0],
            _ => {}
        }
    }
    bad.loop_edges = vec![None, None];
    bad.label = "closing line 2 k2".into();
    assert_closure_refusal(
        &bad,
        &good,
        "the closure's vacuum lines are not one closing line among massless ones: ",
    );
}

/// the fresh loop-vector id of loop i is max + 1 + i: the bound uses
/// the LARGEST index without a loop edge (here loop 1 of the sector, whose
/// loop 0 has its edge), not their count.
#[test]
fn the_fresh_id_bound_uses_the_top_loop_index() {
    let mut base = l2_sector();
    assert!(base.loop_edges[0].is_some());
    base.loop_edges[1] = None;
    base.validate().unwrap();
    let id = i64::MAX as usize - 1;
    match renumber(&base, 21, id).validate() {
        Err(FactorizedError::Closure { detail }) => assert!(
            detail.starts_with(&format!(
                "edge id {id}: the ids (with fresh loop-vector ids up to {id} + 2 above them) must not exceed i64::MAX"
            )),
            "{detail}"
        ),
        other => panic!("{other:?}"),
    }
    renumber(&base, 21, i64::MAX as usize - 2)
        .validate()
        .unwrap();
}

/// without original_cycles the insertions' cycles enter the loop-sign
/// parity in exact arithmetic (a u32 sum wrapped at 2^32): the refusal names
/// the true count.
#[test]
fn insertion_cycles_without_original_cycles_are_counted_exactly() {
    let mut v =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    spec(&mut v).original_cycles = None;
    v.insertions[0].cycles = (u32::MAX, 1);
    let carried = v.closed_loops().unwrap().len();
    let n = carried + u32::MAX as usize + 1;
    let sign = if n.is_multiple_of(2) { 1 } else { -1 };
    v.prefactor.loop_sign = -sign;
    match v.validate() {
        Err(FactorizedError::LoopSign { found, loops }) => {
            assert_eq!((found, loops), (-sign, n))
        }
        other => panic!("{other:?}"),
    }
}

/// The one-loop bubble: the closing line and one massless line (its cut has
/// no loop left: the ORDINARY closure path).
fn one_loop_bubble() -> FactorizedNumerator {
    record(
        "one loop, closing line and a massless line",
        1,
        vec![Some(20)],
        vec![
            line(1, 1, 2, &[1, 1], "0"),
            line(3, 2, 1, &[1, 0], R),
            leg(20, 1, true, &[0, 1]),
            leg(21, 2, false, &[0, 1]),
        ],
        vec![vertex(1, ""), vertex(2, "")],
        Vec::new(),
        None,
        3,
        taylor(0, vec![20, 21], false),
        false,
    )
}

/// The L2 root with two poled contracted children (the truncation path).
fn two_poled_root() -> FactorizedNumerator {
    let eps = "gammalooprs::ε";
    let mut with = l2_root(1, true);
    with.vertices.iter_mut().find(|v| v.id == 1).unwrap().factor = atom(&q(1, &hedge(120)));
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
    with
}

/// each of the four internal-failure CALL SITES of the evaluation:
/// the test-only fault injection makes the WRAPPED call itself return its
/// native error (a `VakintError` from the boundary, the lowering, the
/// ordinary closure; a `SeriesError` from the truncation), and the site's
/// own `map_err` must report `Evaluation` with the stage's prefix; the same
/// record without the injection evaluates (the site is on its path).
#[test]
fn every_internal_failure_site_reports_evaluation() {
    let vakint = Vakint::new().unwrap();
    for (stage, value, settings) in [
        ("boundary", l2_root(2, true), settings()),
        ("lowering and reduction", l2_root(2, true), settings()),
        ("ordinary closure", one_loop_bubble(), settings()),
        (
            "truncation of the insertion product",
            two_poled_root(),
            rstar_site_settings(),
        ),
    ] {
        let run = || {
            evaluate_factorized_closure(
                &vakint,
                &settings,
                std::slice::from_ref(&value),
                &[value.prefactor.total()],
            )
            .map(|_| ())
        };
        run().unwrap_or_else(|e| panic!("{stage}: without injection: {e}"));
        super::INJECT_INTERNAL_FAILURE.with(|c| c.set(Some(stage)));
        let injected = run();
        super::INJECT_INTERNAL_FAILURE.with(|c| c.set(None));
        // the wrapped call itself failed with its NATIVE error; the site's
        // own map_err made it Evaluation with the stage's prefix
        match injected {
            Err(FactorizedError::Evaluation { detail }) => {
                assert!(
                    detail.starts_with(&format!("closure {stage}: ")),
                    "{stage}: {detail}"
                );
                assert!(
                    detail.contains("injected failure") || detail.contains("depth"),
                    "{stage}: not the injected native error: {detail}"
                );
            }
            other => panic!("{stage}: {other:?}"),
        }
    }
}

/// `vakint_input` (public) called DIRECTLY with a marked
/// value checks the prefactor argument too; unmarked, it does not.
#[test]
fn vakint_input_checks_a_marked_value_s_argument() {
    let direct =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k5-direct.record.txt"))
            .unwrap();
    let argument = direct.prefactor.total() * atom("vakint::k(1)");
    // the record's depth is the number of terms
    let settings = VakintSettings {
        number_of_terms_in_epsilon_expansion: direct.record.as_ref().unwrap().depth.unwrap() as i64,
        ..settings()
    };
    assert_graph_prefix(
        super::super::vakint_input(&settings, &direct, &argument).map(|_| ()),
        &reserved_prefix("vakint::k", "prefactor argument"),
        "vakint_input",
    );
    let _ = super::super::vakint_input(&settings, &direct, &direct.prefactor.total()).unwrap();
}

// ---------------------------------------------------------------------------
// names Symbolica's parser resolves to BUILT-IN symbols are
// not in vakint's namespace and are accepted in a marked record. Symbolica
// (rev 394f6e1, `State::BUILTIN_NAMES_AND_ALIASES`): arg, coeff, exp, log,
// sin, cos, sqrt, conj, abs, if, der, 𝑒, 𝜋, ‖, opt, alt, and the aliases
// euler_e (-> 𝑒) and pi (-> 𝜋), all resolved to `symbolica::<name>`; `𝑖` is
// not a symbol at all (the parser reads the complex number 1i).
// ---------------------------------------------------------------------------

/// The parser's resolution, as the rule relies on it.
#[test]
fn bare_builtins_resolve_outside_vakint() {
    let head = |text: &str| -> String {
        let a = atom(text);
        match a.as_view() {
            AtomView::Fun(f) => f.get_symbol().get_name().to_string(),
            AtomView::Var(v) => v.get_symbol().get_name().to_string(),
            other => format!("not a symbol: {other:?}"),
        }
    };
    // Symbolica's own table of the names its parser resolves to built-ins
    for name in [
        "arg", "coeff", "exp", "log", "sin", "cos", "sqrt", "conj", "abs", "if", "der", "𝑒", "𝜋",
        "‖", "opt", "alt", "euler_e", "pi",
    ] {
        assert!(symbolica::state::State::is_builtin(name), "{name}");
    }
    assert!(!symbolica::state::State::is_builtin("gs"));
    for name in ["log", "sin", "cos", "abs", "conj", "arg"] {
        assert_eq!(
            head(&format!("{name}(UFO::x)")),
            format!("symbolica::{name}")
        );
    }
    // exp(x) is normalised by the parser to 𝑒^x (the built-in constant)
    let e = atom("exp(UFO::x)");
    assert_eq!(e, atom("𝑒^UFO::x"));
    for (text, name) in [
        ("𝜋", "symbolica::𝜋"),
        ("pi", "symbolica::𝜋"),
        ("𝑒", "symbolica::𝑒"),
    ] {
        assert_eq!(head(text), name, "{text}");
    }
    assert!(
        matches!(atom("𝑖").as_view(), AtomView::Num(_)),
        "𝑖 is a number"
    );
    // and a bare non-built-in name is read into vakint::
    assert_eq!(head("gs"), "vakint::gs");
    // the walk accepts every built-in
    for text in [
        "𝑖*𝜋*𝑒*pi",
        "log(UFO::x)*exp(UFO::x)*sqrt(UFO::x)*sin(UFO::x)*cos(UFO::x)*abs(UFO::x)*conj(UFO::x)*arg(UFO::x)",
    ] {
        assert_eq!(
            super::super::first_forbidden_symbol(atom(text).as_view(), &[]),
            None,
            "{text}"
        );
    }
}

/// 𝑖, 𝜋 and zeta(3) in a FACTOR, in an INSERTION coefficient and in the
/// PREFACTOR ARGUMENT of a marked record: accepted, and the value is the
/// record's value times the constant; a built-in FUNCTION of a caller scalar
/// (log, cos) is a function head that is not allowed: refused in the factor
/// and in the argument.
#[test]
fn builtins_are_accepted_in_a_marked_record_and_leave_the_value_unchanged() {
    let vakint = Vakint::new().unwrap();
    let base = l2_root(1, true);
    let run = |v: &FactorizedNumerator, argument: &Atom| {
        v.validate().unwrap_or_else(|e| panic!("{}: {e}", v.label));
        evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(v),
            std::slice::from_ref(argument),
        )
        .unwrap_or_else(|e| panic!("{}: {e}", v.label))
        .vakint_value
    };
    let reference = run(&base, &base.prefactor.total());
    assert!(!reference.is_zero());
    for (function, head) in [
        ("log(UFO::x)", "symbolica::log"),
        ("symbolica::zeta(3)*cos(UFO::x)", "symbolica::cos"),
    ] {
        let c = atom(function);
        let mut in_factor = base.clone();
        in_factor.vertices[1].factor = &in_factor.vertices[1].factor * &c;
        let id = in_factor.vertices[1].id;
        match in_factor.validate() {
            Err(FactorizedError::Graph { detail }) => assert!(
                detail.starts_with(&format!(
                    "a function head that is not allowed in vertex {id} factor: {head}"
                )),
                "{detail}"
            ),
            other => panic!("factor {function}: {other:?}"),
        }
        match evaluate_factorized_closure(
            &vakint,
            &settings(),
            std::slice::from_ref(&base),
            &[base.prefactor.total() * &c],
        ) {
            Err(FactorizedError::Graph { detail }) => assert!(
                detail.starts_with(&format!(
                    "a function head that is not allowed in prefactor argument: {head}"
                )),
                "{detail}"
            ),
            other => panic!("argument {function}: {other:?}"),
        }
    }
    for constant in ["𝑖*𝜋", "𝜋^-2", "symbolica::zeta(3)"] {
        let c = atom(constant);
        let mut in_factor = base.clone();
        in_factor.vertices[1].factor = &in_factor.vertices[1].factor * &c;
        let got = run(&in_factor, &in_factor.prefactor.total());
        assert!(
            (&got - &reference * &c).expand().is_zero(),
            "factor {constant}"
        );
        let got = run(&base, &(base.prefactor.total() * &c));
        assert!(
            (&got - &reference * &c).expand().is_zero(),
            "argument {constant}"
        );
    }
    // an insertion coefficient (the producer's k2-child-ins already carries 𝑖)
    let record = include_str!("../fixtures/wh-k2-child-ins.record.txt");
    let value_text = include_str!("../fixtures/wh-k2-child-ins.value.txt");
    assert!(record.contains('𝑖'));
    let mut v = FactorizedNumerator::from_text(record).unwrap();
    let c = atom("𝑖*𝜋");
    v.insertions[0].coefficient = &v.insertions[0].coefficient * &c;
    v.validate().unwrap();
    let got = evaluate_factorized_closure(
        &vakint,
        &rstar_site_settings(),
        std::slice::from_ref(&v),
        &[v.prefactor.total()],
    )
    .unwrap();
    assert_equal_through(
        &got.value,
        &(fixture_value(value_text) * &c),
        valid_through(value_text),
        "k2-child-ins x 𝑖𝜋 in the insertion coefficient",
    );
}

/// The namespace rule against a set of producer record files: every v3
/// record is marked and passes the walk; returns (records, 𝑖 and 𝜋 in the
/// text).
fn namespace_rule_on(files: &[std::path::PathBuf]) -> (usize, usize, usize) {
    let (mut records, mut i, mut pi) = (0, 0, 0);
    for path in files {
        let text = wh_fixtures::read_fixture(path);
        if !text.starts_with("schema\tpint6.factorized-numerator.v3") {
            continue;
        }
        let value = FactorizedNumerator::parse_text(&text)
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert!(
            super::super::is_marked(&value),
            "{}: not marked",
            path.display()
        );
        value
            .check_reserved_symbols()
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        records += 1;
        i += text.matches('𝑖').count();
        pi += text.matches('𝜋').count();
    }
    (records, i, pi)
}

/// Every record file under `dir` (recursively): `*record*.txt`, `input.txt`.
fn record_files(dir: &std::path::Path) -> Vec<std::path::PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        for entry in std::fs::read_dir(&d).unwrap() {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_string_lossy().to_string();
            if path.is_dir() {
                if name != "runs" {
                    stack.push(path);
                }
            } else if name == "input.txt" || (name.contains("record") && name.ends_with(".txt")) {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

/// The vendored producer records (this crate's fixtures): all 10 v3 records
/// pass the namespace rule (they carry 𝑖 and 𝜋 bare).
#[test]
fn the_vendored_records_pass_the_namespace_rule() {
    let dir =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/p_integral/factorized/fixtures");
    let (records, i, pi) = namespace_rule_on(&record_files(&dir));
    eprintln!("NAMESPACE_RULE vendored records={records} 𝑖={i} 𝜋={pi}");
    assert_eq!(records, 10);
    assert!(i > 0 && pi > 0);
}

/// EVERY producer record under `$VAKINT_TEST_WH_FIXTURES` (all 75 v3 records
/// of the wh-fixtures root, `runs/` excluded) passes the namespace rule.
#[test]
#[ignore = "external data: VAKINT_TEST_WH_FIXTURES"]
fn every_producer_record_passes_the_namespace_rule() {
    let (records, i, pi) = namespace_rule_on(&record_files(&wh_fixtures::fixture_root()));
    eprintln!("NAMESPACE_RULE producer records={records} 𝑖={i} 𝜋={pi}");
    wh_fixtures::assert_compared("NAMESPACE_RULE", records, Some(75));
}

/// `contract_to_vakint_dots` (public, takes a value) walks a MARKED value
/// like the other input functions; an unmarked value is not walked.
#[test]
fn contract_to_vakint_dots_walks_a_marked_value() {
    let direct =
        FactorizedNumerator::from_text(include_str!("../fixtures/wh-k5-direct.record.txt"))
            .unwrap();
    assert!(super::super::is_marked(&direct));
    let _ = super::super::evaluate::contract_to_vakint_dots(&direct).unwrap();
    let mut bad = direct.clone();
    bad.vertices[0].factor =
        &bad.vertices[0].factor * atom("vakint::dot(vakint::p(1),vakint::p(1))");
    assert_graph_prefix(
        super::super::evaluate::contract_to_vakint_dots(&bad).map(|_| ()),
        &reserved_prefix(
            "vakint::dot",
            &format!("vertex {} factor", bad.vertices[0].id),
        ),
        "contract_to_vakint_dots",
    );
    let mut unmarked = bad.clone();
    unmarked.record = None;
    unmarked.insertions.clear();
    assert!(!super::super::is_marked(&unmarked));
    assert!(!matches!(
        super::super::evaluate::contract_to_vakint_dots(&unmarked),
        Err(FactorizedError::Graph { ref detail }) if detail.starts_with("reserved symbol")
    ));
}

/// A graph for the carried-cycle count: `lines` (id, source, target,
/// source half-edge, target half-edge, line kind; source / target `None` =
/// external), vertices 1..=`vertices`, and `insertions` (vertex, boundary
/// half-edges).
type CycleLine = (
    usize,
    Option<usize>,
    Option<usize>,
    Option<usize>,
    Option<usize>,
    LineKind,
);

fn cycle_graph(
    vertices: usize,
    lines: &[CycleLine],
    insertions: &[(usize, Vec<usize>)],
) -> FactorizedNumerator {
    let end = |v: Option<usize>| v.map_or(Endpoint::External, Endpoint::Vertex);
    FactorizedNumerator {
        label: "cycle graph".into(),
        kind: FactorizedKind::Direct,
        dimension: FactorizedNumerator::default_dimension(),
        loop_edges: Vec::new(),
        external_edges: Vec::new(),
        edges: lines
            .iter()
            .map(|(id, s, t, hs, ht, line)| FactorizedEdge {
                id: *id,
                source: end(*s),
                target: end(*t),
                source_hedge: *hs,
                target_hedge: *ht,
                particle: "x".into(),
                line: *line,
                power: 1,
                mass_squared: Atom::Zero,
                momentum: Vec::new(),
                factor: Atom::num(1),
            })
            .collect(),
        vertices: (1..=vertices)
            .map(|id| FactorizedVertex {
                id,
                rule: if insertions.iter().any(|(v, _)| *v == id) {
                    "ins:x".into()
                } else {
                    "v".into()
                },
                factor: Atom::num(1),
            })
            .collect(),
        ports: Vec::new(),
        projector: None,
        prefactor: Prefactor::default(),
        colour_branch: None,
        insertions: insertions
            .iter()
            .map(|(vertex, hedges)| super::super::Insertion {
                vertex: *vertex,
                hedges: hedges.clone(),
                cycles: (0, 0),
                coefficient: Atom::num(1),
                factor: Atom::num(1),
            })
            .collect(),
        record: None,
    }
}

/// THE CARRIED-CYCLE COUNT (R2: an insertion's span is ONE node): a
/// fermion cycle closed through an insertion (1 -> 2, the span {2, 3}, 3 ->
/// 1); a ghost cycle through TWO insertions (1 -> 2, {2, 3}, 3 -> 4, {4, 5},
/// 5 -> 1); a ghost cycle through a span of THREE vertices (1 -> 2, {2, 3,
/// 4}, 4 -> 1); each counts one cycle, and none without the insertions
/// (the open paths). An open crown ghost line (external -> 1 -> 2 ->
/// external) is no cycle.
#[test]
fn the_carried_cycle_count_identifies_each_insertion_span() {
    use LineKind::{Fermion, Ghost};
    // half-edges: edge e has source half-edge 10e and target half-edge 10e + 1
    let h = |e: usize| (Some(10 * e), Some(10 * e + 1));
    let line = |e: usize, s: usize, t: usize, kind: LineKind| {
        let (hs, ht) = h(e);
        (e, Some(s), Some(t), hs, ht, kind)
    };
    let fermion = [line(1, 1, 2, Fermion), line(2, 3, 1, Fermion)];
    // the insertion at vertex 2 with the boundary half-edge 20 (edge 2 at 3)
    assert_eq!(
        carried_cycles(&cycle_graph(3, &fermion, &[(2, vec![11, 20])])),
        (1, 0)
    );
    assert_eq!(carried_cycles(&cycle_graph(3, &fermion, &[])), (0, 0));
    let two = [
        line(1, 1, 2, Ghost),
        line(2, 3, 4, Ghost),
        line(3, 5, 1, Ghost),
    ];
    let both = [(2, vec![11, 20]), (4, vec![21, 30])];
    assert_eq!(carried_cycles(&cycle_graph(5, &two, &both)), (0, 1));
    assert_eq!(carried_cycles(&cycle_graph(5, &two, &both[..1])), (0, 0));
    // a span of three vertices: the insertion at 2 with half-edges at 3 and 4
    let wide = [
        line(1, 1, 2, Ghost),
        line(2, 4, 1, Ghost),
        line(3, 3, 5, LineKind::Scalar),
    ];
    assert_eq!(
        carried_cycles(&cycle_graph(5, &wide, &[(2, vec![11, 30, 20])])),
        (0, 1)
    );
    assert_eq!(carried_cycles(&cycle_graph(5, &wide, &[])), (0, 0));
    // an open crown ghost line
    let open = [
        (1, None, Some(1), None, Some(11), Ghost),
        line(2, 1, 2, Ghost),
        (3, Some(2), None, Some(30), None, Ghost),
    ];
    assert_eq!(carried_cycles(&cycle_graph(2, &open, &[])), (0, 0));
}
