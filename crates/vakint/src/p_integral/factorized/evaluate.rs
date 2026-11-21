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
//!    graph through the external leg of smallest edge id (among the legs
//!    carrying it: soft legs are skipped): exactly the
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
//! 4. OPEN SLOTS ([`FactorizedNumerator::open_slots`]: the ports no
//!    projector closes, and slots of the projector that are not ports; e.g. a
//!    propagator-type R* remainder whose external indices the caller keeps):
//!    the edge momenta and metrics idenso leaves open carry the CALLER's
//!    indices (`mink(D, a)` -> `a`) into vakint's `k(i, a)`, `p(1, a)`,
//!    `g(a, b)`; vakint's tensor prepass (feynkit-tensor's reducer with `p(1)`
//!    external) reduces them. [`evaluate_factorized`] returns `Σ c_S T_S`,
//!    `T_S` products of vakint's `g(a, b)` and `p(1, a)` in the orientation of
//!    step 2; [`evaluate_factorized_ports`] returns the same value in the
//!    CALLER's notation (`spenso::g(slot, slot)`, `Q(x, slot)` of the
//!    external leg `x` of the value's own `p(1)`, `spenso::dot(Q(x, mink(D)),
//!    Q(x, mink(D)))`). A soft external leg (all-zero momentum vector) may
//!    carry ports. For a MARKED value (and on the ports entry, which marks
//!    what it is given) at most [`MAX_OPEN_PORTS`] open Lorentz slots, no
//!    open spinor slot, no integer index; a value without new features keeps
//!    the base's behaviour on the old entries (no such limit).
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
    if !value.has_new_features() {
        return vakint_momenta_v1(value);
    }
    vakint_momenta_oriented(value).map(|(momenta, head, _)| (momenta, head))
}

/// [`vakint_momenta`] and the orientation: vakint's `p(1)` is `into ×` the
/// value's own `p(1)` (the last component of its momentum
/// vectors). vakint's `p(1)` flows INTO the graph through the first external
/// leg (by edge id) that carries the external momentum; a leg whose momentum
/// vector is all zero is a SOFT leg (e.g. the zero-momentum legs of an R*
/// child made propagator-type) and is skipped.
/// Every edge's vakint momentum, vakint's `p(1)` and its orientation.
type OrientedMomenta = (Vec<(usize, Atom)>, Atom, Orientation);

fn vakint_momenta_oriented(value: &FactorizedNumerator) -> Result<OrientedMomenta> {
    let loops = value.loops();
    if value.externals() != 1 {
        return Err(evaluation(format!(
            "a two-point value with one independent external momentum, not {}",
            value.externals()
        )));
    }
    let head = S.p.call(1);
    let legs = value.edges.iter().filter(|edge| edge.is_external());
    let first_leg = legs
        .clone()
        .min_by_key(|edge| edge.id)
        .ok_or_else(|| evaluation("a value without external legs"))?;
    // Unchanged when the first leg by id carries momentum; a soft first leg
    // (all components zero) was refused before and now defers to the first
    // leg that carries the external momentum.
    let soft = |edge: &&super::FactorizedEdge| edge.momentum.iter().all(|c| *c == 0);
    let first = if soft(&first_leg) {
        legs.filter(|edge| !soft(edge))
            .min_by_key(|edge| edge.id)
            .ok_or_else(|| evaluation("every external leg is soft (zero momentum)"))?
    } else {
        first_leg
    };
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
    Ok((
        out,
        head,
        Orientation {
            into,
            leg: first.id,
        },
    ))
}

/// How vakint's `p(1)` relates to the value's own: `p_vakint(1) = into ×
/// p_value(1)`, flowing in through the external leg `leg`.
#[derive(Clone, Copy, Debug)]
struct Orientation {
    into: i64,
    leg: usize,
}

/// The most open Lorentz ports a MARKED propagator-type value may keep (the
/// unmarked values keep the base's behaviour: no limit): the transverse Gram
/// matrix of the tensor reduction with one external momentum is regular at
/// `d = 4` through rank 7 and singular at rank 8 (a spurious `1/(d-4)`, N/notes
/// D6 section 3); ranks above 6 are not validated. A CONSERVATIVE guard, not
/// the exact condition: the exact one is the open LOOP-tensor rank (an open
/// slot carried by an external momentum adds no loop rank).
pub const MAX_OPEN_PORTS: usize = 6;

/// The open slots ([`FactorizedNumerator::open_slots`]) the p-integral route
/// can return for a MARKED value: Lorentz slots `mink(D, a)` with a
/// non-integer index `a`, at most [`MAX_OPEN_PORTS`]. A closed value passes.
/// Applied in the marked arm of [`vakint_input`] (so in
/// [`evaluate_factorized`] and [`evaluate_factorized_ports`] for marked
/// values, as at its introduction); an unmarked value runs the base's body,
/// which has none of these rules.
fn check_open_ports(value: &FactorizedNumerator) -> Result<()> {
    let open = value.open_slots()?;
    if open.is_empty() {
        return Ok(());
    }
    let ports = |detail: String| FactorizedError::Ports { detail };
    if open.len() > MAX_OPEN_PORTS {
        return Err(ports(format!(
            "{} open slots; a propagator-type value keeps at most {MAX_OPEN_PORTS} open \
             Lorentz slots (close the others with a projector)",
            open.len()
        )));
    }
    for port in &open {
        let Some(slot) = super::as_slot(port.as_view()) else {
            return Err(ports(format!("{port} is not a slot")));
        };
        if slot.rep != names::MINK {
            return Err(ports(format!(
                "open port {port} is not a Lorentz slot: the p-integral route returns open \
                 Lorentz indices only (close the spinor line first)"
            )));
        }
        if i64::try_from(slot.index.as_view()).is_ok() {
            return Err(ports(format!(
                "open port {port} has an integer index: vakint would read a metric g(i, j) of \
                 integer indices as a vector; name the index (hedge(h, k) or a symbol)"
            )));
        }
    }
    Ok(())
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
    if !value.has_new_features() {
        return contract_to_vakint_dots_v1(value);
    }
    // a marked value: the whole validation (the reserved names, the bounds,
    // the contract, the atom rules, the defined-and-finite fields), as at the
    // entries; an unmarked value is not walked (the base behaviour)
    if super::is_marked(value) {
        value.validate()?;
    }
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
    // Open slots (D6, helper VH; `open_slots`: the ports without a projector,
    // or what a projector leaves open): an edge momentum carrying an open
    // index, `Q(e, mink(D, i))`, becomes the edge's vakint momentum with that
    // index (`k(j, i)`, `±p(1, i)`), a metric between two open slots
    // `g(mink(D, a), mink(D, b))` vakint's `g(a, b)`; vakint's tensor prepass
    // reduces them after its denominator match. A closed value (every index
    // contracted) is unchanged by this.
    let open_slots = if value.projector.is_none() {
        !value.ports.is_empty()
    } else {
        !value.open_slots()?.is_empty()
    };
    let mapped = if open_slots {
        let port_index = |view: AtomView<'_>| match view {
            AtomView::Fun(rep)
                if rep.get_symbol().get_name() == names::MINK && rep.get_nargs() == 2 =>
            {
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
            return Err(evaluation(format!(
                "an open index on the unknown edge {edge}"
            )));
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
    if !value.has_new_features() {
        return vakint_input_v1(settings, value, prefactor);
    }
    if let FactorizedKind::Closure(_) = value.kind {
        return Err(FactorizedError::Reserved {
            feature: "the closure kind (R* remainder with a closing line)".into(),
        });
    }
    // a marked value called here directly: the entry check too
    if super::is_marked(value) {
        // the settings' rule before any computing rule
        super::check_master_modes(settings)?;
    }
    super::check_marked_entry(
        value,
        prefactor,
        vk_symbol!(&settings.epsilon_symbol),
        Some(settings.number_of_terms_in_epsilon_expansion),
    )?;
    super::check_entry_names(value, prefactor, settings)?;
    value.validate()?;
    check_open_ports(value)?;
    value.check_insertion_epsilon(vk_symbol!(&settings.epsilon_symbol).get_name())?;
    vakint_input_unchecked(settings, value, prefactor)
}

/// [`vakint_input`] without the contract checks (the reference of the tests:
/// a numerator contracted into one vertex factor by hand, as a caller's
/// contracted path holds it).
pub(crate) fn vakint_input_unchecked(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<(Atom, f64)> {
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

/// Whether `value` is one of the record kinds of the R* interface (a record
/// with insertions, or a v3 direct record): its value must not depend on the
/// caller's tensor-reduction method or notation (record format, section 7).
/// Values accepted before (no insertion, v1) keep the caller's settings.
pub(crate) fn is_new_kind(value: &FactorizedNumerator) -> bool {
    value.has_new_features()
}

/// The caller's settings with the tensor reduction and the notation pinned to
/// what the R* sites use (FeynKit, dot notation), on a clone.
pub(crate) fn pinned_settings(settings: &VakintSettings) -> VakintSettings {
    let mut pinned = settings.clone();
    pinned.tensor_reduction_method = crate::TensorReductionMethod::FeynKit;
    pinned.use_dot_product_notation = true;
    pinned
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
    if !value.has_new_features() {
        return evaluate_factorized_v1(vakint, settings, value, prefactor);
    }
    // a marked value: the separately supplied prefactor gets the checks of
    // the record's prefactor fields, the settings' epsilon the caller-symbol
    // rule (no-ops for an unmarked value: the base behaviour; the closure
    // kind keeps its earlier refusal of this entry, `Reserved`)
    if !matches!(value.kind, super::FactorizedKind::Closure(_)) {
        if super::is_marked(value) {
            // the settings' rule before any computing rule
            super::check_master_modes(settings)?;
        }
        super::check_marked_entry(
            value,
            prefactor,
            vk_symbol!(&settings.epsilon_symbol),
            Some(settings.number_of_terms_in_epsilon_expansion),
        )?;
        super::check_entry_names(value, prefactor, settings)?;
        super::check_data_atoms_catalogue(value, vakint)?;
    }
    let pinned;
    let settings = if is_new_kind(value) {
        pinned = pinned_settings(settings);
        &pinned
    } else {
        settings
    };
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
    super::backstop::check_result(value, prefactor, settings, &result, true)?;
    Ok(FactorizedValue {
        value: result,
        input,
        contract_seconds,
        evaluate_seconds: started.elapsed().as_secs_f64(),
    })
}

/// A value of [`evaluate_factorized_ports`]: the open tensor in the caller's
/// notation.
#[derive(Clone, Debug)]
pub struct PortsValue {
    /// The Laurent expansion in the CALLER's notation (module documentation,
    /// step 4): `spenso::g(slot, slot)`, `Q(x, slot)`, `spenso::dot(Q(x,
    /// mink(D)), Q(x, mink(D)))`, `Q(x)` the external leg of the value's own
    /// `p(1)`.
    pub value: Atom,
    /// The same value in vakint notation (what [`evaluate_factorized`]
    /// returns).
    pub vakint_value: Atom,
    /// The open slots ([`FactorizedNumerator::open_slots`]), sorted.
    pub open: Vec<Atom>,
    pub input: Atom,
    pub contract_seconds: f64,
    pub evaluate_seconds: f64,
    /// The `vakint::timing` counters of this call on the calling thread:
    /// contraction / lowering / reduction / finalisation
    /// ([`crate::timing::Delta::stage`]), thread CPU and sampled RSS when the
    /// resources are enabled.
    pub timing: crate::timing::Delta,
}

/// [`evaluate_factorized`] with the result mapped to the caller's notation
/// ([`to_caller_notation_p`]): the entry for a propagator-type value whose
/// open slots the caller keeps (no projector, or a projector that closes only
/// some ports). A closed value comes back as [`evaluate_factorized`]'s scalar
/// with `dot(p(1), p(1))` written `spenso::dot(Q(x, mink(D)), Q(x, mink(D)))`.
/// A value without new features is evaluated here as a MARKED one (the entry
/// marks a local clone with an empty v3 direct record), so the record kinds'
/// rules, the open-slot rules and the reserved-namespace and bare-name
/// refusals apply to whatever goes through this entry.
pub fn evaluate_factorized_ports(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<PortsValue> {
    // a ports record is a new kind whatever it carries: a value without new
    // features takes the record kinds' rules here (an empty v3 provenance:
    // the weights stay the v1 rule, `original_cycles` is none)
    let settings = &pinned_settings(settings);
    let marked;
    let value = if value.has_new_features() {
        value
    } else {
        marked = FactorizedNumerator {
            record: Some(super::DirectRecord {
                original_cycles: None,
                depth: None,
                weights: Vec::new(),
            }),
            ..value.clone()
        };
        &marked
    };
    value.validate()?;
    check_open_ports(value)?;
    let counted = crate::timing::Snapshot::now();
    let started = Instant::now();
    let evaluated = evaluate_factorized(vakint, settings, value, prefactor)?;
    let open = value.open_slots()?;
    let caller = to_caller_notation_p(value, &evaluated.value)?;
    super::backstop::check_result(value, prefactor, settings, &caller, false)?;
    let timing = counted.elapsed();
    if crate::timing::record_log_enabled() {
        eprintln!(
            "{}",
            crate::timing::record_line(
                "direct",
                &value.label,
                value.loops(),
                Some(settings.number_of_terms_in_epsilon_expansion),
                None,
                &timing,
                started.elapsed().as_secs_f64(),
                evaluated.input.as_view().get_byte_size(),
            )
        );
    }
    Ok(PortsValue {
        value: caller,
        vakint_value: evaluated.value,
        open,
        input: evaluated.input,
        contract_seconds: evaluated.contract_seconds,
        evaluate_seconds: evaluated.evaluate_seconds,
        timing,
    })
}

/// A value of [`evaluate_factorized`] (vakint notation) in the caller's
/// notation: `g(a, b)` of two open indices -> `spenso::g(slot_a, slot_b)`;
/// `p(1, a)` -> `σ Q(x, slot_a)` and `dot(p(1), p(1))` ->
/// `spenso::dot(Q(x, mink(D)), Q(x, mink(D)))`, with `x` the edge named for
/// `p(1)` (`external 1 edge=x`), else the leg vakint's `p(1)` flows through,
/// and `σ = ±1` such that the result is in the VALUE's own `p(1)` (the last
/// component of its momentum vectors), whatever leg vakint chose. Fails
/// closed if any of vakint's `g`, `p`, `k` survives.
pub fn to_caller_notation_p(value: &FactorizedNumerator, vakint_value: &Atom) -> Result<Atom> {
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
    let (_, _, orientation) = vakint_momenta_oriented(value)?;
    let loops = value.loops();
    let x = value.external_edges[0].unwrap_or(orientation.leg);
    let leg = value
        .edge(x)
        .ok_or_else(|| evaluation(format!("the external leg {x} of p(1) is not an edge")))?;
    let xi = leg.momentum[loops];
    if (xi != 1 && xi != -1) || leg.momentum[..loops].iter().any(|c| *c != 0) {
        return Err(evaluation(format!(
            "the leg {x} of p(1) does not carry ±p(1)"
        )));
    }
    // p_vakint(1) = into p_value(1) = into xi Q(x).
    let sigma = Atom::num(orientation.into * xi);
    let q = symbolica::symbol!(names::MOMENTUM);
    let open = value.open_slots()?;
    let slot_of = |index: AtomView<'_>| {
        open.iter()
            .find(|slot| super::as_slot(slot.as_view()).is_some_and(|s| s.index.as_view() == index))
            .cloned()
    };
    let mink = symbolica::function!(symbolica::symbol!(names::MINK), Atom::var(value.dimension));
    let head = |view: AtomView<'_>| match view {
        AtomView::Fun(f) if f.get_symbol() == S.p && f.get_nargs() == 1 => {
            f.iter().next().and_then(|id| i64::try_from(id).ok()) == Some(1)
        }
        _ => false,
    };
    let mapped = vakint_value.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let args = function.iter().collect::<Vec<_>>();
        if function.get_symbol() == S.g && args.len() == 2 {
            if let (Some(a), Some(b)) = (slot_of(args[0]), slot_of(args[1])) {
                **out = symbolica::function!(ETS.metric, a, b);
            }
        } else if function.get_symbol() == S.p && args.len() == 2 {
            if i64::try_from(args[0]).ok() == Some(1)
                && let Some(slot) = slot_of(args[1])
            {
                **out = &sigma * symbolica::function!(q, Atom::num(x as i64), slot);
            }
        } else if function.get_symbol() == S.dot
            && args.len() == 2
            && head(args[0])
            && head(args[1])
        {
            let vector = symbolica::function!(q, Atom::num(x as i64), &mink);
            **out = symbolica::function!(SPENSO_TAG.dot, &vector, &vector);
        }
    });
    let mut leftover = None;
    mapped.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            if symbol == S.g || symbol == S.p || symbol == S.k {
                leftover = Some(part.to_owned());
                return false;
            }
        }
        leftover.is_none()
    });
    if let Some(part) = leftover {
        return Err(evaluation(format!(
            "vakint notation left after the mapping to the caller's notation: {part}"
        )));
    }
    Ok(mapped)
}

// The base's old entries, verbatim, for values without new features.

/// The base's `vakint_momenta`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn vakint_momenta_v1(value: &FactorizedNumerator) -> Result<(Vec<(usize, Atom)>, Atom)> {
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

/// The base's `contract_to_vakint_dots`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn contract_to_vakint_dots_v1(value: &FactorizedNumerator) -> Result<(Atom, f64)> {
    let (momenta, _) = vakint_momenta_v1(value)?;
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
            AtomView::Fun(rep)
                if rep.get_symbol().get_name() == names::MINK && rep.get_nargs() == 2 =>
            {
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
            return Err(evaluation(format!(
                "an open index on the unknown edge {edge}"
            )));
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

/// The base's `vakint_input`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn vakint_input_v1(
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
    let (momenta, head) = vakint_momenta_v1(value)?;
    let (scalar, seconds) = contract_to_vakint_dots_v1(value)?;
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

/// The base's `evaluate_factorized`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn evaluate_factorized_v1(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<FactorizedValue> {
    let (input, contract_seconds) = vakint_input_v1(settings, value, prefactor)?;
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
