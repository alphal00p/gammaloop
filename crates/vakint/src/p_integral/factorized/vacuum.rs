//! The uncontracted-numerator input mode for SINGLE-SCALE VACUUM integrals
//! (helper VH, numerator assignment, lane rstar-kinds): the vacuum sibling of
//! [`super::evaluate_factorized`] (which takes massless propagator graphs).
//!
//! A [`FactorizedNumerator`] whose internal edges carry loop momenta only
//! (no external momentum flows through a denominator) and massive lines of
//! one mass scale; its numerator is the product of its spenso factors (gamma
//! chains, metrics, momenta), optionally closed by the caller's projector.
//! OPEN INDICES are the declared `ports` (`mink(D, ..)` slots of the factor
//! product); EXTERNAL VECTORS are the independent external momenta `p(j)`,
//! written in the factors as `Q(x, ..)` of the external leg `x` with momentum
//! `p(j)`. The same declaration serves the propagator case (lane d6) and the
//! R* hard parts (D5 section 14).
//!
//! 1. CONTRACT with idenso: `projector · PROD factors` through
//!    [`super::contract_with_idenso`] (cook, `infer_indexed`,
//!    `simplify_algebra{gamma canonical, epsilon}`, `to_dots`, uncook): every
//!    trace, chain and contraction among the factors; what is left is the
//!    scope vakint's vacuum evaluation always supported: scalar products,
//!    momenta with an open index and metrics between open indices.
//! 2. MAP to vakint's vacuum notation (linear substitution): an edge dot
//!    becomes `dot` of the edges' momenta in `k(1..L)` and `p(j)`, an edge
//!    momentum with open port i becomes the same combination with the raw
//!    index `PORT + i` (`k(l, PORT + i)`, `p(j, PORT + i)`), a metric between
//!    ports `g(PORT + i, PORT + j)`; anything else is an error.
//! 3. EVALUATE with vakint's EXISTING vacuum path (`Vakint::evaluate`): its
//!    vacuum tensor reduction (method of the settings: FeynKit by default,
//!    feynkit `TensorReducer`) projects the loop-momentum tensors onto tensors
//!    of purely external objects, which ride along as constant coefficients;
//!    matching, the backend of the evaluation order (RustRed first by default:
//!    one to three loops from its vendored artifacts, four with its external
//!    data, five with an application package; FORM-based backends only where
//!    the order names them) and the finaliser.
//! 4. The value keeps its OPEN indices; it is returned in vakint notation and
//!    mapped back to the caller's notation (`g(slot, slot)`, `Q(x, slot)`,
//!    `spenso::dot(Q(x, mink(D)), Q(y, mink(D)))`).

use std::time::Instant;

use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use super::{Endpoint, FactorizedError, FactorizedKind, FactorizedNumerator, Result, names};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintSettings};

/// First raw index of the open ports in vakint notation (`PORT + i`).
pub const PORT: i64 = 9001;

fn error(detail: impl Into<String>) -> FactorizedError {
    FactorizedError::Evaluation {
        detail: format!("vacuum: {}", detail.into()),
    }
}

/// A value of [`evaluate_factorized_vacuum`] and what it was made of.
#[derive(Clone, Debug)]
pub struct VacuumValue {
    /// The Laurent expansion with the open ports and external vectors in the
    /// caller's notation.
    pub value: Atom,
    /// The same in vakint notation (raw port indices `PORT + i`, `p(j)`).
    pub vakint_value: Atom,
    /// The integrand handed to `Vakint::evaluate` (in epsilon).
    pub input: Atom,
    /// Wall seconds of the idenso contraction and of the vakint evaluation.
    pub contract_seconds: f64,
    pub evaluate_seconds: f64,
}

/// The vakint momentum `SUM c_i k(i) + SUM d_j p(j)` of every edge, with an
/// optional raw index on every vector.
fn momentum(value: &FactorizedNumerator, id: i64, index: Option<i64>) -> Option<Atom> {
    let edge = value.edges.iter().find(|edge| edge.id as i64 == id)?;
    let loops = value.loops();
    let mut out = Atom::Zero;
    for (i, c) in edge.momentum.iter().enumerate() {
        if *c == 0 {
            continue;
        }
        let (head, n) = if i < loops {
            (S.k, i as i64 + 1)
        } else {
            (S.p, (i - loops) as i64 + 1)
        };
        let vector = match index {
            Some(index) => function!(head, Atom::num(n), Atom::num(index)),
            None => function!(head, Atom::num(n)),
        };
        out += Atom::num(*c) * vector;
    }
    Some(out)
}

/// `(edge id, its mink argument)` of `Q(e, mink(D, ..))` / `Q(e, mink(D))`.
fn momentum_parts(view: AtomView<'_>) -> Option<(i64, AtomView<'_>)> {
    let AtomView::Fun(function) = view else {
        return None;
    };
    if function.get_symbol().get_name() != names::MOMENTUM || function.get_nargs() != 2 {
        return None;
    }
    let mut args = function.iter();
    let edge = i64::try_from(args.next()?).ok()?;
    Some((edge, args.next()?))
}

fn is_compact(rep: AtomView<'_>) -> bool {
    matches!(rep, AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK && f.get_nargs() == 1)
}

/// The checks of the vacuum mode (the contract's graph checks apply where they
/// do not assume a propagator graph).
fn validate_vacuum(value: &FactorizedNumerator) -> Result<()> {
    if value.kind != FactorizedKind::Direct {
        return Err(FactorizedError::Reserved {
            feature: "the closure kind in the vacuum mode".into(),
        });
    }
    if !value.insertions.is_empty() {
        return Err(FactorizedError::Reserved {
            feature: "insertions".into(),
        });
    }
    let loops = value.loops();
    if loops == 0 {
        return Err(error("a vacuum value has at least one loop"));
    }
    let mut mass = None::<&Atom>;
    for edge in &value.edges {
        if edge.momentum.len() != loops + value.externals() {
            return Err(error(format!(
                "edge {}: momentum of the wrong length",
                edge.id
            )));
        }
        if edge.is_external() {
            continue;
        }
        if edge.momentum[loops..].iter().any(|c| *c != 0) {
            return Err(error(format!(
                "edge {} carries an external momentum (not a vacuum integral)",
                edge.id
            )));
        }
        if !edge.mass_squared.is_zero() {
            match mass {
                None => mass = Some(&edge.mass_squared),
                Some(m) if *m == edge.mass_squared => {}
                Some(m) => {
                    return Err(error(format!(
                        "two mass scales ({m} and {}): single-scale vacua only",
                        edge.mass_squared
                    )));
                }
            }
        }
    }
    if mass.is_none() {
        return Err(error(
            "no massive line: a massless vacuum integral is scaleless",
        ));
    }
    for port in &value.ports {
        if !matches!(port.as_view(), AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK && f.get_nargs() == 2)
        {
            return Err(error(format!(
                "port {port}: only Lorentz ports mink(D, index) are in the vacuum scope"
            )));
        }
    }
    Ok(())
}

/// The closed-except-ports scalar of `projector · PROD factors` contracted by
/// idenso and mapped to vakint's vacuum notation (steps 1 and 2), in the
/// value's dimension symbol; and the seconds of the contraction.
pub fn contract_to_vakint_vacuum(value: &FactorizedNumerator) -> Result<(Atom, f64)> {
    validate_vacuum(value)?;
    let started = Instant::now();
    let contracted = value.contract_naively()?;
    let seconds = started.elapsed().as_secs_f64();
    let port_index = |slot: AtomView<'_>| {
        value
            .ports
            .iter()
            .position(|port| port.as_view() == slot)
            .map(|i| PORT + i as i64)
    };
    let failure = std::cell::Cell::new(None::<String>);
    let mapped = contracted.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let symbol = function.get_symbol();
        let name = symbol.get_name();
        // a scalar product of two edge momenta
        if (name == names::METRIC || name == "spenso::dot") && function.get_nargs() == 2 {
            let args = function.iter().collect::<Vec<_>>();
            if let (Some((a, ra)), Some((b, rb))) =
                (momentum_parts(args[0]), momentum_parts(args[1]))
                && is_compact(ra)
                && is_compact(rb)
            {
                match (momentum(value, a, None), momentum(value, b, None)) {
                    (Some(x), Some(y)) => **out = S.dot.call((&x, &y)),
                    _ => failure.set(Some(format!("a dot of unknown edges {a}, {b}"))),
                }
                return;
            }
            // a metric between two open ports
            if name == names::METRIC
                && let (Some(i), Some(j)) = (port_index(args[0]), port_index(args[1]))
            {
                **out = S.g.call((Atom::num(i), Atom::num(j)));
                return;
            }
            return;
        }
        // an edge momentum carrying an open port
        if let Some((edge, slot)) = momentum_parts(part)
            && let Some(index) = port_index(slot)
        {
            match momentum(value, edge, Some(index)) {
                Some(vector) => **out = vector,
                None => failure.set(Some(format!("unknown edge {edge}"))),
            }
        }
    });
    if let Some(detail) = failure.take() {
        return Err(error(detail));
    }
    let mut leftover = None;
    mapped.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            if [
                names::MOMENTUM,
                "spenso::dot",
                names::MINK,
                names::BIS,
                names::METRIC,
                names::GAMMA,
            ]
            .contains(&name)
            {
                leftover = Some(part.to_owned());
            }
        }
        leftover.is_none()
    });
    if let Some(part) = leftover {
        return Err(error(format!(
            "outside the vacuum scope (metrics and momenta with open indices) after idenso: {part}"
        )));
    }
    Ok((mapped, seconds))
}

/// One propagator line of the vacuum topology handed to vakint.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumLine {
    pub id: usize,
    pub source: usize,
    pub target: usize,
    /// Over the loop momenta only.
    pub momentum: Vec<i64>,
    pub mass_squared: Atom,
    pub power: i64,
}

/// The lines of the vacuum graph: its internal edges with every chain of
/// lines in SERIES merged (a vertex with exactly two internal line ends whose
/// lines carry the same momentum through it and the same mass is removed, the
/// two lines become one with the summed power). This is graph normalisation
/// (vakint's matching takes one line per distinct propagator), as the merged
/// lines of D5 section 1.3; external legs carry no momentum through a vacuum
/// and do not count.
pub fn vacuum_lines(value: &FactorizedNumerator) -> Vec<VacuumLine> {
    let loops = value.loops();
    let mut lines = value
        .edges
        .iter()
        .filter(|edge| !edge.is_external())
        .map(|edge| {
            let (Endpoint::Vertex(source), Endpoint::Vertex(target)) = (edge.source, edge.target)
            else {
                unreachable!("internal edges have two vertices")
            };
            VacuumLine {
                id: edge.id,
                source,
                target,
                momentum: edge.momentum[..loops].to_vec(),
                mass_squared: edge.mass_squared.clone(),
                power: edge.power,
            }
        })
        .collect::<Vec<_>>();
    loop {
        let mut merged = false;
        let mut vertices = lines
            .iter()
            .flat_map(|line| [line.source, line.target])
            .collect::<Vec<_>>();
        vertices.sort_unstable();
        vertices.dedup();
        for v in vertices {
            let ends = lines
                .iter()
                .enumerate()
                .filter(|(_, l)| l.source == v || l.target == v)
                .collect::<Vec<_>>();
            if ends.len() != 2 || ends.iter().any(|(_, l)| l.source == l.target) {
                continue;
            }
            let (ia, a) = (ends[0].0, ends[0].1.clone());
            let (ib, b) = (ends[1].0, ends[1].1.clone());
            if a.mass_squared != b.mass_squared {
                continue;
            }
            // a into v, b out of v (reversing a line flips its momentum)
            let (a_from, a_in) = if a.target == v {
                (a.source, a.momentum.clone())
            } else {
                (a.target, a.momentum.iter().map(|c| -c).collect())
            };
            let (b_to, b_out) = if b.source == v {
                (b.target, b.momentum.clone())
            } else {
                (b.source, b.momentum.iter().map(|c| -c).collect())
            };
            if a_in != b_out {
                continue;
            }
            let line = VacuumLine {
                id: a.id.min(b.id),
                source: a_from,
                target: b_to,
                momentum: a_in,
                mass_squared: a.mass_squared.clone(),
                power: a.power + b.power,
            };
            let (first, second) = (ia.min(ib), ia.max(ib));
            lines.remove(second);
            lines[first] = line;
            merged = true;
            break;
        }
        if !merged {
            break;
        }
    }
    lines.sort_by_key(|line| line.id);
    lines
}

/// The integrand vakint receives (module documentation): `prefactor ×
/// scalar × topo(PROD prop(id, edge(s,t), momentum, mass², power))`, `D → 4 −
/// 2ε`; and the contraction's seconds.
pub fn vakint_vacuum_input(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<(Atom, f64)> {
    let (scalar, seconds) = contract_to_vakint_vacuum(value)?;
    let mut topology = Atom::num(1);
    for line in vacuum_lines(value) {
        let mut momentum = Atom::Zero;
        for (i, c) in line.momentum.iter().enumerate() {
            if *c != 0 {
                momentum += Atom::num(*c) * S.k.call(i as i64 + 1);
            }
        }
        topology *= S.prop.call((
            line.id as i64,
            S.edge.call((line.source as i64, line.target as i64)),
            momentum,
            line.mass_squared.clone(),
            Atom::num(line.power),
        ));
    }
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * epsilon;
    let input = (prefactor * scalar * S.topo.call(topology))
        .replace(Atom::var(value.dimension))
        .with(in_epsilon);
    Ok((input, seconds))
}

/// The vakint-notation value mapped back to the caller's notation: raw port
/// indices to the port slots, `p(j)` to `Q(x_j, ..)` of the external leg of
/// `p(j)`; an external vector without a leg stays `p(j)`.
pub fn to_caller_notation(value: &FactorizedNumerator, vakint_value: &Atom) -> Atom {
    let q = symbolica::symbol!(names::MOMENTUM);
    let rep = value.ports.first().and_then(|port| match port.as_view() {
        AtomView::Fun(f) => Some(function!(
            f.get_symbol(),
            f.iter().next().unwrap().to_owned()
        )),
        _ => None,
    });
    // a raw index, or a representation slot `mink(D', raw index)`
    let port = |index: AtomView<'_>| {
        let raw = match index {
            AtomView::Fun(f) if f.get_nargs() == 2 => {
                f.iter().nth(1).and_then(|i| i64::try_from(i).ok())
            }
            _ => i64::try_from(index).ok(),
        };
        raw.and_then(|i| usize::try_from(i - PORT).ok())
            .and_then(|i| value.ports.get(i).cloned())
    };
    let leg = |j: AtomView<'_>| {
        i64::try_from(j)
            .ok()
            .and_then(|j| usize::try_from(j - 1).ok())
            .and_then(|j| value.external_edges.get(j).copied().flatten())
    };
    vakint_value.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let args = function.iter().collect::<Vec<_>>();
        if function.get_symbol() == S.g && args.len() == 2 {
            if let (Some(a), Some(b)) = (port(args[0]), port(args[1])) {
                **out = function!(ETS.metric, a, b);
            }
        } else if function.get_symbol() == S.p && args.len() == 2 {
            if let (Some(x), Some(slot)) = (leg(args[0]), port(args[1])) {
                **out = function!(q, Atom::num(x as i64), slot);
            }
        } else if function.get_symbol() == S.dot && args.len() == 2 {
            let vector = |v: AtomView<'_>| match v {
                AtomView::Fun(f) if f.get_symbol() == S.p && f.get_nargs() == 1 => {
                    leg(f.iter().next().unwrap())
                }
                _ => None,
            };
            if let (Some(x), Some(y), Some(rep)) = (vector(args[0]), vector(args[1]), rep.as_ref())
            {
                **out = function!(
                    SPENSO_TAG.dot,
                    function!(q, Atom::num(x as i64), rep),
                    function!(q, Atom::num(y as i64), rep)
                );
            }
        }
    })
}

/// Evaluate a single-scale VACUUM value (module documentation) times
/// `prefactor` with vakint's existing vacuum evaluation; open ports stay open.
pub fn evaluate_factorized_vacuum(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<VacuumValue> {
    let (input, contract_seconds) = vakint_vacuum_input(settings, value, prefactor)?;
    let started = Instant::now();
    let vakint_value = vakint
        .evaluate(settings, input.as_view())
        .map_err(|e| error(e.to_string()))?;
    let evaluate_seconds = started.elapsed().as_secs_f64();
    Ok(VacuumValue {
        value: to_caller_notation(value, &vakint_value),
        vakint_value,
        input,
        contract_seconds,
        evaluate_seconds,
    })
}

#[cfg(test)]
mod tests;
