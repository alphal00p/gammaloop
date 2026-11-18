//! The adapter from a [`FactorizedNumerator`] to a value: idenso contracts,
//! vakint evaluates (helper VH, numerator assignment, lane rstar-interface).
//!
//! 1. CONTRACT with idenso only: `projector.tensor · Π factors` through
//!    [`contract_with_idenso`] (the R* production chain: cook with
//!    reversible encoding and representation payloads, `infer_indexed`,
//!    `simplify_algebra{gamma: canonical(), epsilon: true}`, `to_dots`,
//!    uncook). The result is a closed scalar in EDGE momentum dot products
//!    `g(Q(e, mink(D)), Q(f, mink(D)))` (idenso's dot notation).
//! 2. MAP each edge dot to vakint's `dot` of the edges' momenta in vakint's
//!    p-integral notation: the value's loop momenta become `k(1..L)` in order
//!    and its external momentum `orientation · p(1)`, `p(1)` flowing INTO the
//!    graph through the external leg of smallest edge id: exactly the
//!    convention of vakint's FeynKit ingress (`VakintExpression::from_diagram`
//!    for a propagator graph), so that a DOT graph reaches vakint here in the
//!    same notation as through that ingress (lane refpath's chain).
//! 3. EVALUATE with vakint's existing p path: `prefactor ×
//!    dot(p(1), p(1))^projector.p2_power × scalar × topo(Π prop(id, edge(s,t),
//!    momentum, 0, power))` with `D → 4 − 2ε`, through `Vakint::evaluate`
//!    (the direct lowering of plain requests, its zero-sector pruning switch,
//!    `fold_requests`, one RustRed reduction per folded channel, the p
//!    finaliser). The folded channel inputs are available, as for every
//!    evaluation, through the existing gate observers
//!    ([`crate::rustred_evaluation::p_gate::observe`]).
//!
//! No numerator algebra happens here beyond idenso's; the mapping of step 2 is
//! linear substitution. `prefactor` is the caller's scalar (colour, couplings,
//! graph weight, projection normalisation: e.g. `value.prefactor.total()`
//! times a normalisation), multiplied once.

use std::time::Instant;

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};

use super::{Endpoint, FactorizedError, FactorizedKind, FactorizedNumerator, Result, names};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintSettings};

/// A value of [`evaluate_factorized`] and what it was made of.
#[derive(Clone, Debug)]
pub struct FactorizedValue {
    /// The Laurent expansion vakint returns (ε of the settings).
    pub value: Atom,
    /// The integrand handed to `Vakint::evaluate` (numerator × topology, in
    /// ε), for dumps and checks.
    pub input: Atom,
    /// Wall seconds of the idenso contraction and of the vakint evaluation.
    pub contract_seconds: f64,
    pub evaluate_seconds: f64,
}

fn evaluation(detail: impl Into<String>) -> FactorizedError {
    FactorizedError::Evaluation {
        detail: detail.into(),
    }
}

/// The vakint momentum of every edge (module documentation, step 2), and
/// `p(1)`.
fn vakint_momenta(value: &FactorizedNumerator) -> Result<(Vec<(usize, Atom)>, Atom)> {
    let loops = value.loops();
    if value.externals() != 1 {
        return Err(evaluation(format!(
            "a two-point value with one independent external momentum, not {}",
            value.externals()
        )));
    }
    let head = S.p.call(1);
    let first = value
        .edges
        .iter()
        .filter(|edge| edge.is_external())
        .min_by_key(|edge| edge.id)
        .ok_or_else(|| evaluation("a value without external legs"))?;
    if first.momentum[..loops].iter().any(|c| *c != 0) {
        return Err(evaluation("an external leg carries a loop momentum"));
    }
    let along = first.momentum[loops];
    let into = if first.source == Endpoint::External {
        along
    } else {
        -along
    };
    if into != 1 && into != -1 {
        return Err(evaluation("the first external leg does not carry ±p"));
    }
    let mut out = Vec::with_capacity(value.edges.len());
    for edge in &value.edges {
        let mut momentum = Atom::Zero;
        for (i, c) in edge.momentum[..loops].iter().enumerate() {
            if *c != 0 {
                momentum += Atom::num(*c) * S.k.call(i as i64 + 1);
            }
        }
        let c = edge.momentum[loops] * into;
        if c != 0 {
            momentum += Atom::num(c) * &head;
        }
        out.push((edge.id, momentum));
    }
    Ok((out, head))
}

/// The edge id of idenso's bare momentum `Q(e, mink(D))`.
fn bare_momentum(view: AtomView<'_>) -> Option<i64> {
    let AtomView::Fun(function) = view else {
        return None;
    };
    if function.get_symbol().get_name() != names::MOMENTUM || function.get_nargs() != 2 {
        return None;
    }
    let mut args = function.iter();
    let edge = i64::try_from(args.next()?).ok()?;
    let AtomView::Fun(rep) = args.next()? else {
        return None;
    };
    (rep.get_symbol().get_name() == names::MINK && rep.get_nargs() == 1).then_some(edge)
}

/// The closed scalar of `projector.tensor · Π factors` contracted by idenso
/// and mapped to vakint's dot products of the edge momenta (steps 1 and 2),
/// in the value's dimension symbol; and the seconds of the contraction.
pub fn contract_to_vakint_dots(value: &FactorizedNumerator) -> Result<(Atom, f64)> {
    let (momenta, _) = vakint_momenta(value)?;
    let started = Instant::now();
    let contracted = value.contract_naively()?;
    let seconds = started.elapsed().as_secs_f64();
    let lookup = |id: i64| {
        momenta
            .iter()
            .find(|(edge, _)| *edge as i64 == id)
            .map(|(_, momentum)| momentum.clone())
    };
    let unknown = std::cell::Cell::new(None);
    let mapped = contracted.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        // idenso writes a contracted pair as `g(Q(e, mink(D)), Q(f, mink(D)))`
        // (`AlgebraContraction::Fully`) or `dot(..)` of the same arguments.
        let symbol = function.get_symbol();
        let name = symbol.get_name();
        if (name != names::METRIC && name != "spenso::dot") || function.get_nargs() != 2 {
            return;
        }
        let args = function.iter().collect::<Vec<_>>();
        let (Some(a), Some(b)) = (bare_momentum(args[0]), bare_momentum(args[1])) else {
            return;
        };
        match (lookup(a), lookup(b)) {
            (Some(left), Some(right)) => **out = S.dot.call((&left, &right)),
            _ => unknown.set(Some((a, b))),
        }
    });
    if let Some((a, b)) = unknown.get() {
        return Err(evaluation(format!("a dot of unknown edges {a}, {b}")));
    }
    // Open ports (no projector; D6, helper VH): an edge momentum carrying a
    // port index, `Q(e, mink(D, i))`, becomes the edge's vakint momentum with
    // that index (`k(j, i)`, `±p(1, i)`), a metric between two ports
    // `g(mink(D, a), mink(D, b))` vakint's `g(a, b)`; vakint's tensor prepass
    // reduces them after its denominator match. With a projector every index
    // is contracted and this changes nothing.
    let mapped = if value.projector.is_none() && !value.ports.is_empty() {
        let port_index = |view: AtomView<'_>| match view {
            AtomView::Fun(rep) if rep.get_symbol().get_name() == names::MINK && rep.get_nargs() == 2 => {
                rep.iter().nth(1).map(|index| index.to_owned())
            }
            _ => None,
        };
        let indexed = |momentum: &Atom, index: &Atom| {
            momentum.replace_map(|term, _, out| {
                if let AtomView::Fun(vector) = term
                    && (vector.get_symbol() == S.k || vector.get_symbol() == S.p)
                    && vector.get_nargs() == 1
                {
                    **out = FunctionBuilder::new(vector.get_symbol())
                        .add_arg(vector.iter().next().unwrap())
                        .add_arg(index)
                        .finish();
                }
            })
        };
        let unknown = std::cell::Cell::new(None);
        let open = mapped.replace_map(|part, _, out| {
            let AtomView::Fun(function) = part else {
                return;
            };
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            let args = function.iter().collect::<Vec<_>>();
            if args.len() != 2 {
                return;
            }
            if name == names::METRIC
                && let (Some(a), Some(b)) = (port_index(args[0]), port_index(args[1]))
            {
                **out = FunctionBuilder::new(S.g).add_arg(&a).add_arg(&b).finish();
            } else if name == names::MOMENTUM
                && let (Ok(edge), Some(index)) = (i64::try_from(args[0]), port_index(args[1]))
            {
                match lookup(edge) {
                    Some(momentum) => **out = indexed(&momentum, &index),
                    None => unknown.set(Some(edge)),
                }
            }
        });
        if let Some(edge) = unknown.get() {
            return Err(evaluation(format!("an open index on the unknown edge {edge}")));
        }
        open
    } else {
        mapped
    };
    // Nothing of the tensor notation may survive.
    let mut leftover = None;
    mapped.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            if name == names::MOMENTUM
                || name == "spenso::dot"
                || name == names::MINK
                || name == names::BIS
                || name == names::METRIC
                || name == names::GAMMA
            {
                leftover = Some(part.to_owned());
                return false;
            }
        }
        leftover.is_none()
    });
    if let Some(part) = leftover {
        return Err(evaluation(format!(
            "the contraction left tensor notation: {part}"
        )));
    }
    Ok((mapped, seconds))
}

/// The integrand vakint receives for `value` (module documentation): the
/// prefactor, `dot(p(1), p(1))^p2_power`, the contracted scalar and the
/// topology, with `D → 4 − 2ε` (ε of the settings); and the contraction's
/// seconds.
pub fn vakint_input(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<(Atom, f64)> {
    if let FactorizedKind::Closure(_) = value.kind {
        return Err(FactorizedError::Reserved {
            feature: "the closure kind (R* remainder with a closing line)".into(),
        });
    }
    value.validate()?;
    let (momenta, head) = vakint_momenta(value)?;
    let (scalar, seconds) = contract_to_vakint_dots(value)?;
    let mut topology = Atom::num(1);
    for edge in value.edges.iter().filter(|edge| !edge.is_external()) {
        let (Endpoint::Vertex(source), Endpoint::Vertex(target)) = (edge.source, edge.target)
        else {
            unreachable!("internal edges have two vertices (validated)");
        };
        if !edge.mass_squared.is_zero() {
            return Err(evaluation("a massive line (p-integrals are massless)"));
        }
        let momentum = &momenta
            .iter()
            .find(|(id, _)| *id == edge.id)
            .expect("every edge has a momentum")
            .1;
        topology *= S.prop.call((
            edge.id as i64,
            S.edge.call((source as i64, target as i64)),
            momentum,
            Atom::Zero,
            Atom::num(edge.power),
        ));
    }
    let p2_power = value.projector.as_ref().map_or(0, |p| p.p2_power);
    let square = S.dot.call((&head, &head));
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * epsilon;
    let input = (prefactor * square.pow(Atom::num(p2_power)) * scalar * S.topo.call(topology))
        .replace(Atom::var(value.dimension))
        .with(in_epsilon);
    Ok((input, seconds))
}

/// Evaluate `value` (a DIRECT-kind value: one graph, one colour branch) times
/// `prefactor` with vakint's existing p path (module documentation). Wrap the
/// call in [`crate::rustred_evaluation::p_gate::observe`] to also get the
/// folded channel inputs RustRed received.
pub fn evaluate_factorized(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<FactorizedValue> {
    let (input, contract_seconds) = vakint_input(settings, value, prefactor)?;
    let loops = value.loops();
    if value
        .edges
        .iter()
        .filter(|edge| !edge.is_external())
        .all(|edge| edge.momentum[loops] == 0)
    {
        // No internal line carries the external momentum: a scaleless
        // massless vacuum integral, an exact zero in dimensional
        // regularisation (refpath's `scaleless_contact`).
        return Ok(FactorizedValue {
            value: Atom::Zero,
            input,
            contract_seconds,
            evaluate_seconds: 0.0,
        });
    }
    let started = Instant::now();
    let result = vakint
        .evaluate(settings, input.as_view())
        .map_err(|e| evaluation(e.to_string()))?;
    Ok(FactorizedValue {
        value: result,
        input,
        contract_seconds,
        evaluate_seconds: started.elapsed().as_secs_f64(),
    })
}
