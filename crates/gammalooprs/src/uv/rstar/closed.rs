//! `GAMMALOOP_IR_HARD=closed|check`: IR hard moments through polarized
//! metric closures instead of an open vacuum tensor projection.
//!
//! The hard Taylor polynomial `T_omega H` of an IR correction is a sum of
//! homogeneous sectors, one per multidegree `m` in the soft (boundary)
//! momenta `Q_e`. For a parity-even hard part whose open slots are Minkowski
//! slots only, the vacuum average of a sector is
//!
//! ```text
//!   <T_m H>(Q) = (1/m!) Q_{a_1}^{alpha_1} .. Q_{a_n}^{alpha_n} <D_{alpha mu}>,
//!   D = d^n T_m H / dQ_{a_1}^{alpha_1} .. dQ_{a_n}^{alpha_n}   (Q-free),
//! ```
//!
//! and the vacuum average `<D>` of the Q-free tensor `D` with open slots
//! `{mu} u {alpha}` is isotropic: `<D> = sum_P c_P P` over the metric pairings
//! `P` of those slots, `c = G^-1 (P' . D)` with the pairing Gram matrix `G`.
//! Every `P' . D` is a closed scalar: the hard loop momenta occur only in
//! scalar products, so no tensor reduction is needed. Re-substituting the
//! momenta gives the insertion of each `c_P` (`g(alpha,mu) -> Q_a(mu)`,
//! `g(alpha,beta) -> Q_a.Q_b`). The `(coefficient, insertion)` pairs are
//! the same as those of `project_hard_vacuum` + `hard_moments`: the pairing
//! projection is the vacuum projection, and the coefficient of each formal
//! covariant monomial is unique. The `d` of `G^-1` stays in the coefficient,
//! outside the quotient `K`, exactly as the projector's `1/d` does today.
//!
//! Supported: open slots all Minkowski, total rank `r = slots + n <= 4` per
//! sector (one, one or three pairings), no gamma5/projector/Levi-Civita
//! symbol. Anything else declines with a reason and the caller uses today's
//! path.

use std::sync::LazyLock;

use eyre::{Result, eyre};
use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use linnet::half_edge::involution::EdgeIndex;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibraryRep, Minkowski, RepName},
        slot::ParseableAind,
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    symbol,
};

use super::taylor::TaylorExpansion;
use crate::{numerator::aind::Aind, utils::{GS, W_}, uv::approx::integrated::Integrated};
use spenso::structure::{partial::PartialStructureExt, slot::IsAbstractSlot};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClosedMode {
    Off,
    On,
    Check,
}

/// `GAMMALOOP_IR_HARD`: unset/`off` (today), `closed`, or `check` (both
/// paths, compared at a random rational point; today's moments are used).
pub fn closed_mode() -> ClosedMode {
    static MODE: LazyLock<ClosedMode> = LazyLock::new(|| {
        match std::env::var("GAMMALOOP_IR_HARD").ok().as_deref().map(str::trim) {
            Some("closed") => ClosedMode::On,
            Some("check") => ClosedMode::Check,
            _ => ClosedMode::Off,
        }
    });
    *MODE
}

/// Counts of one closed evaluation, for the `IR_CLOSED` line.
#[derive(Debug, Default)]
pub struct ClosedStats {
    pub slots: usize,
    pub sectors: usize,
    pub odd_sectors: usize,
    pub closures: usize,
    pub taylor_seconds: f64,
    pub polarize_seconds: f64,
    pub close_seconds: f64,
    pub taylor_bytes: usize,
    pub closed_bytes: usize,
    /// Closures idenso certified a closed scalar (the producer fact of
    /// `simplify_view_certified`), and those that needed the slot scan.
    pub certified: usize,
    pub uncertified: usize,
    /// Moments whose coefficient expands to zero, dropped before the quotient.
    pub dropped_zero: usize,
}

/// A slot of a sector tensor: an open slot of the hard part, or the
/// polarization slot of the `k`-th derivative (momentum of soft edge `edge`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Leg {
    Open(usize),
    Alpha { index: usize, edge: usize },
}

fn cooking() -> CookSettings {
    CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true)
}

fn bare_minkowski() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

fn alpha_slot(index: usize) -> Atom {
    Minkowski {}
        .new_rep(GS.dim)
        .to_symbolic([Aind::Normal(9_000_000 + index).to_atom()])
}

fn soft_vector(edge: usize, slot: Option<&Atom>) -> Atom {
    match slot {
        Some(slot) => function!(GS.emr_mom, edge, slot),
        None => function!(GS.emr_mom, edge, bare_minkowski()),
    }
}

/// `Q(edge, ..)` with the given soft edge, and its last argument.
fn soft_momentum<'a>(view: AtomView<'a>, edge: usize) -> Option<AtomView<'a>> {
    let AtomView::Fun(function) = view else {
        return None;
    };
    if function.get_symbol() != GS.emr_mom || function.get_nargs() < 2 {
        return None;
    }
    (usize::try_from(function.get(0)).ok() == Some(edge)).then(|| function.get(function.get_nargs() - 1))
}

/// A slot-less vector `V(.., mink(D))` with its last argument replaced by `slot`.
fn indexed(vector: AtomView<'_>, slot: &Atom) -> Option<Atom> {
    let AtomView::Fun(function) = vector else {
        return None;
    };
    let nargs = function.get_nargs();
    if nargs < 2 || function.get(nargs - 1) != bare_minkowski().as_view() {
        return None;
    }
    let mut builder = FunctionBuilder::new(function.get_symbol());
    for argument in function.iter().take(nargs - 1) {
        builder = builder.add_arg(argument);
    }
    Some(builder.add_arg(slot.as_view()).finish())
}

/// The derivative of a polynomial in `Q_edge` along the polarization slot
/// `slot`: every primitive occurrence `X(Q)` becomes `X(Q) + t dX`, and the
/// coefficient of `t` is read off (product rule). Errors when `Q_edge` occurs
/// in any other position (for example inside a denominator).
fn derivative(expression: &Atom, edge: usize, slot: &Atom) -> std::result::Result<Atom, String> {
    let t = symbol!("gammalooprs::ir_closed_t");
    let mut failure = None;
    let substituted = expression.replace_map(|view, _, out| {
        if failure.is_some() {
            return;
        }
        let AtomView::Fun(function) = view else {
            return;
        };
        if let Some(last) = soft_momentum(view, edge) {
            if last == bare_minkowski().as_view() {
                failure = Some("bare_soft_vector".to_owned());
                return;
            }
            // Q(edge, mink(D, i)) -> Q + t g(mink(D, i), alpha)
            **out = view.to_owned() + Atom::var(t) * function!(ETS.metric, last, slot);
            return;
        }
        if function.get_symbol() == SPENSO_TAG.dot && function.get_nargs() == 2 {
            let (a, b) = (function.get(0), function.get(1));
            let a_soft = soft_momentum(a, edge).is_some();
            let b_soft = soft_momentum(b, edge).is_some();
            let increment = match (a_soft, b_soft) {
                (true, true) => indexed(a, slot).map(|vector| Atom::num(2) * vector),
                (true, false) => indexed(b, slot),
                (false, true) => indexed(a, slot),
                (false, false) => return,
            };
            match increment {
                Some(increment) => **out = view.to_owned() + Atom::var(t) * increment,
                None => failure = Some(format!("dot_shape:{view}")),
            }
            return;
        }
        if function.get_symbol() == GS.den {
            let mut found = false;
            view.visitor(&mut |inner| {
                found |= soft_momentum(inner, edge).is_some();
                !found
            });
            if found {
                failure = Some("soft_momentum_in_den".to_owned());
            }
            // Leave den untouched (it is Q-free after the Taylor step).
            **out = view.to_owned();
        }
    });
    if let Some(failure) = failure {
        return Err(failure);
    }
    let derived = substituted.derivative(t).replace(t).with(Atom::Zero);
    // Every Q_edge of a homogeneous sector is consumed by its own derivative;
    // any remaining one marks a non-polynomial position.
    Ok(derived)
}

fn contains_soft(expression: AtomView<'_>, edge: usize) -> bool {
    let mut found = false;
    expression.visitor(&mut |inner| {
        found |= soft_momentum(inner, edge).is_some();
        !found
    });
    found
}

/// The perfect matchings of `legs` (all, in a fixed order); `r <= 4`.
/// `None` for any other rank: the caller declines (fails closed). A rank > 4
/// sector is never truncated to fewer pairings.
fn pairings(legs: &[Leg]) -> Option<Vec<Vec<(Leg, Leg)>>> {
    match legs.len() {
        0 => Some(vec![vec![]]),
        2 => Some(vec![vec![(legs[0], legs[1])]]),
        4 => Some(vec![
            vec![(legs[0], legs[1]), (legs[2], legs[3])],
            vec![(legs[0], legs[2]), (legs[1], legs[3])],
            vec![(legs[0], legs[3]), (legs[1], legs[2])],
        ]),
        _ => None,
    }
}

/// The metric product of a pairing, in the sector tensor's slots.
fn pairing_tensor(pairing: &[(Leg, Leg)], open: &[Atom], alphas: &[Atom]) -> Atom {
    let slot = |leg: Leg| match leg {
        Leg::Open(i) => open[i].clone(),
        Leg::Alpha { index, .. } => alphas[index].clone(),
    };
    pairing.iter().fold(Atom::num(1), |product, (a, b)| {
        product * function!(ETS.metric, slot(*a), slot(*b))
    })
}

/// The covariant insertion of a pairing: polarization slots carry their
/// soft momenta again.
fn pairing_insertion(pairing: &[(Leg, Leg)], open: &[Atom]) -> Atom {
    pairing.iter().fold(Atom::num(1), |product, pair| {
        let factor = match *pair {
            (Leg::Open(i), Leg::Open(j)) => function!(ETS.metric, &open[i], &open[j]),
            (Leg::Open(i), Leg::Alpha { edge, .. }) | (Leg::Alpha { edge, .. }, Leg::Open(i)) => {
                soft_vector(edge, Some(&open[i]))
            }
            (Leg::Alpha { edge: a, .. }, Leg::Alpha { edge: b, .. }) => {
                function!(SPENSO_TAG.dot, soft_vector(a, None), soft_vector(b, None))
            }
        };
        product * factor
    })
}

/// The open slots of `expression` when all are explicit Minkowski slots.
fn open_minkowski_slots(expression: &Atom) -> std::result::Result<Vec<Atom>, String> {
    let cooking = cooking();
    let cooked = cooking
        .try_cook(expression.as_view())
        .map_err(|error| format!("cook:{error}"))?;
    let tensor = SymbolicTensor::infer(cooked).map_err(|error| format!("infer:{error}"))?;
    let minkowski: LibraryRep = Minkowski {}.into();
    let mut slots = Vec::new();
    for slot in tensor.structure().logical_slots() {
        let spenso::structure::partial::PartialIndex::Explicit(index) = slot.aind else {
            return Err("implicit_slot".into());
        };
        if slot.rep().rep != minkowski {
            return Err(format!("slot_rep:{}", slot.rep().rep));
        }
        let atom = slot.rep().slot::<AbstractIndex, _>(index).to_atom();
        slots.push(cooking.uncook(atom.as_view()));
    }
    Ok(slots)
}

/// The `(coefficient, insertion)` moments of the hard routed product
/// `routed` (open slots kept, soft momenta `Q(e)` for `e` in `soft_edges`)
/// at Taylor order `degree`, or a decline reason. `simplify` is the hard
/// path's contraction (`RStarIntegrand::simplify_hard_taylor`).
pub(super) fn closed_moments(
    routed: &Atom,
    soft_edges: &[EdgeIndex],
    degree: usize,
    simplify: &dyn Fn(&Atom) -> Result<(Atom, bool)>,
    stats: &mut ClosedStats,
) -> Result<std::result::Result<Vec<(Atom, Atom)>, String>> {
    let parity_odd = [
        idenso::dirac::AGS.gamma5,
        idenso::dirac::AGS.projp,
        idenso::dirac::AGS.projm,
        *idenso::epsilon::EPSILON_SYMBOL,
    ];
    if parity_odd.iter().any(|symbol| routed.contains_symbol(*symbol)) {
        return Ok(Err("parity_odd".into()));
    }
    let open = match open_minkowski_slots(routed) {
        Ok(open) => open,
        Err(reason) => return Ok(Err(reason)),
    };
    stats.slots = open.len();
    if open.len() > 4 {
        return Ok(Err(format!("slots:{}", open.len())));
    }
    let edges = soft_edges.iter().map(|edge| usize::from(*edge)).collect::<Vec<_>>();

    // Sectors: the multisets of soft momenta of size n <= degree. An odd
    // total rank (slots + n) averages to zero for a parity-even hard part.
    let mut sectors = Vec::new();
    for n in 0..=degree {
        let mut current = Vec::new();
        multisets(edges.len(), n, 0, &mut current, &mut |indices| {
            sectors.push(indices.iter().map(|i| edges[*i]).collect::<Vec<_>>());
        });
    }
    if sectors
        .iter()
        .any(|sector| (open.len() + sector.len()) % 2 == 0 && open.len() + sector.len() > 4)
    {
        return Ok(Err(format!("rank>4 slots={}", open.len())));
    }

    let dimension = Atom::var(GS.dim);
    let lambda = symbol!("rstar_hard_taylor"; Scalar);
    let mu = symbol!("gammalooprs::ir_closed_mu"; Scalar);
    let scale = |expression: &Atom, edges: &[usize], parameter: Symbol| {
        let replacements = edges
            .iter()
            .map(|edge| {
                let momentum = function!(GS.emr_mom, *edge, W_.x___);
                Replacement::new(momentum.to_pattern(), Atom::var(parameter) * momentum)
            })
            .collect::<Vec<_>>();
        expression.replace_multiple(&replacements)
    };
    let mut moments: Vec<(Atom, Atom)> = Vec::new();
    for multidegree in sectors {
        let n = multidegree.len();
        stats.sectors += 1;
        if (open.len() + n) % 2 == 1 {
            stats.odd_sectors += 1;
            continue;
        }
        // The homogeneous sector: the other soft momenta at zero, then the
        // exact power of each sector momentum (the series' own coefficients).
        let started = std::time::Instant::now();
        let others = edges
            .iter()
            .copied()
            .filter(|edge| !multidegree.contains(edge))
            .collect::<Vec<_>>();
        let mut sector = routed.clone();
        if !others.is_empty() {
            sector = TaylorExpansion::new(mu, 0).orders(&scale(&sector, &others, mu))?.swap_remove(0);
        }
        let mut counts = std::collections::BTreeMap::new();
        for edge in &multidegree {
            *counts.entry(*edge).or_insert(0usize) += 1;
        }
        for (edge, count) in &counts {
            if sector.is_zero() {
                break;
            }
            sector = TaylorExpansion::new(lambda, *count)
                .orders(&scale(&sector, &[*edge], lambda))?
                .swap_remove(*count);
        }
        // Slash momenta inside Dirac chains become explicit vectors, as for
        // project_hard_vacuum's rank count.
        let sector = Integrated::materialize_tensor(&sector)?;
        stats.taylor_bytes += sector.as_view().get_byte_size();
        stats.taylor_seconds += started.elapsed().as_secs_f64();
        if sector.is_zero() {
            continue;
        }
        // Polarize: one derivative per soft momentum of the multidegree.
        let polarize_started = std::time::Instant::now();
        let alphas = (0..n).map(alpha_slot).collect::<Vec<_>>();
        let mut tensor = sector;
        for (index, edge) in multidegree.iter().enumerate() {
            tensor = match derivative(&tensor, *edge, &alphas[index]) {
                Ok(tensor) => tensor,
                Err(reason) => return Ok(Err(reason)),
            };
        }
        if edges.iter().any(|edge| contains_soft(tensor.as_view(), *edge)) {
            return Ok(Err("soft_momentum_left".into()));
        }
        // 1/m! for the multiplicities of the multidegree.
        let factorial = counts
            .values()
            .map(|count| (1..=*count as i64).product::<i64>())
            .product::<i64>();
        stats.polarize_seconds += polarize_started.elapsed().as_secs_f64();

        // Close with the metric pairings of all slots (channel markers).
        let close_started = std::time::Instant::now();
        let legs = (0..open.len())
            .map(Leg::Open)
            .chain(
                multidegree
                    .iter()
                    .enumerate()
                    .map(|(index, edge)| Leg::Alpha { index, edge: *edge }),
            )
            .collect::<Vec<_>>();
        // Fails closed: every rank other than 0, 2, 4 declines the whole term.
        let Some(matchings) = pairings(&legs) else {
            return Ok(Err(format!("rank:{}", legs.len())));
        };
        debug_assert_eq!(matchings.len(), [1, 0, 1, 0, 3][legs.len()]);
        let channels = (0..matchings.len())
            .map(|j| symbol!(format!("gammalooprs::ir_closed_tau{j}").as_str()))
            .collect::<Vec<_>>();
        let closure = matchings
            .iter()
            .zip(&channels)
            .fold(Atom::Zero, |sum, (pairing, channel)| {
                sum + Atom::var(*channel) * pairing_tensor(pairing, &open, &alphas)
            });
        let (closed, certified) = simplify(&(tensor * closure))?;
        // The closure must leave a closed scalar. Without idenso's certificate,
        // scan for open slots: a leaked spinor index (an open trace or chain)
        // or any other open slot declines the term.
        if certified {
            stats.certified += 1;
        } else {
            stats.uncertified += 1;
            match open_slots(&closed) {
                Ok(0) => {}
                Ok(count) => return Ok(Err(format!("open_after_closure:{count}"))),
                Err(reason) => return Ok(Err(format!("open_after_closure:{reason}"))),
            }
        }
        stats.closures += matchings.len();
        stats.closed_bytes += closed.as_view().get_byte_size();
        // The channel markers occur linearly, possibly inside unexpanded
        // products: collect in them only.
        let channel_atoms = channels.iter().map(|channel| Atom::var(*channel)).collect::<Vec<_>>();
        let mut projections = vec![Atom::Zero; channels.len()];
        for (monomial, coefficient) in closed.coefficient_list::<u8>(&channel_atoms) {
            let Some(position) = channel_atoms.iter().position(|channel| *channel == monomial) else {
                if coefficient.is_zero() {
                    continue;
                }
                return Ok(Err(format!("channel_monomial:{monomial}")));
            };
            projections[position] = &projections[position] + coefficient;
        }
        // c = G^-1 (P' . D): r = 0: 1; r = 2: 1/d; r = 4: the isotropic duals.
        let coefficients = match matchings.len() {
            1 if legs.is_empty() => vec![projections[0].clone()],
            1 => vec![&projections[0] / &dimension],
            3 => {
                let sum = projections.iter().fold(Atom::Zero, |sum, p| sum + p);
                let normalization = Atom::one() / (&dimension * (&dimension - 1) * (&dimension + 2));
                projections
                    .iter()
                    .map(|p| ((&dimension + 2) * p - &sum) * &normalization)
                    .collect()
            }
            _ => return Ok(Err(format!("pairings:{}", matchings.len()))),
        };
        stats.close_seconds += close_started.elapsed().as_secs_f64();
        for (pairing, coefficient) in matchings.iter().zip(coefficients) {
            if coefficient.is_zero() {
                continue;
            }
            // Zero coefficients that are not literally zero (equal channels):
            // a nonzero value at the random point proves nonzero; only a zero
            // value pays for the exact expansion.
            if random_point(&coefficient).is_zero() && coefficient.expand().is_zero() {
                stats.dropped_zero += 1;
                continue;
            }
            let coefficient =
                Integrated::dimensionally_regularized(&(coefficient / Atom::num(factorial)));
            let insertion = pairing_insertion(pairing, &open);
            match moments.iter_mut().find(|(_, existing)| *existing == insertion) {
                Some((sum, _)) => *sum = &*sum + coefficient,
                None => moments.push((coefficient, insertion)),
            }
        }
    }
    // The sum over sectors can cancel as well (equal insertions summed).
    let before = moments.len();
    moments.retain(|(coefficient, _)| {
        !(random_point(coefficient).is_zero() && coefficient.expand().is_zero())
    });
    stats.dropped_zero += before - moments.len();
    Ok(Ok(moments))
}

/// The number of open slots of a (closed) expression.
fn open_slots(expression: &Atom) -> std::result::Result<usize, String> {
    let cooking = cooking();
    let cooked = cooking
        .try_cook(expression.as_view())
        .map_err(|error| format!("cook:{error}"))?;
    let tensor = SymbolicTensor::infer(cooked).map_err(|error| format!("infer:{error}"))?;
    Ok(tensor.structure().logical_slots().len())
}

/// A fixed rational for every variable and function of `expression` (by its
/// canonical string): one exact random point. Equal expressions give equal
/// numbers; it is a probabilistic equality test for check mode.
pub(super) fn random_point(expression: &Atom) -> Atom {
    use std::hash::{Hash, Hasher};
    let value = |key: String, variable: bool| -> Atom {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        key.hash(&mut hasher);
        let n = (hasher.finish() % 1_000_003) as i64 + 2;
        if variable {
            Atom::num(n) / Atom::num(97) + Atom::num(1) / Atom::num(3)
        } else {
            Atom::num(n)
        }
    };
    expression
        .replace_map(|node, _, out| match node {
            AtomView::Fun(_) => **out = value(node.to_canonical_string(), false),
            AtomView::Var(_) => **out = value(node.to_canonical_string(), true),
            _ => {}
        })
        .expand()
}

/// Every non-decreasing index sequence of length `n` from `start..count`.
fn multisets(
    count: usize,
    n: usize,
    start: usize,
    current: &mut Vec<usize>,
    emit: &mut dyn FnMut(&[usize]),
) {
    if current.len() == n {
        emit(current);
        return;
    }
    for index in start..count {
        current.push(index);
        multisets(count, n, index, current, emit);
        current.pop();
    }
}

#[allow(dead_code)]
fn unused() -> Result<()> {
    Err(eyre!("unused"))
}

#[cfg(test)]
mod tests {
    use super::super::{
        expression::RStarIntegrand, moments::hard_moments, project::project_hard_vacuum,
    };
    use super::*;
    use crate::initialisation::test_initialise;
    use idenso::representations::Bispinor;

    fn slot(hedge: u16) -> Atom {
        Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(hedge, 0).to_atom()])
    }
    fn spinor(hedge: u16) -> Atom {
        Bispinor {}.new_rep(4).to_symbolic([Aind::Hedge(hedge, 1).to_atom()])
    }
    fn vector(edge: usize) -> Atom {
        soft_vector(edge, None)
    }
    fn indexed(edge: usize, slot: &Atom) -> Atom {
        soft_vector(edge, Some(slot))
    }
    fn dot(a: &Atom, b: &Atom) -> Atom {
        function!(SPENSO_TAG.dot, a, b)
    }
    fn gamma(left: &Atom, right: &Atom, argument: &Atom) -> Atom {
        function!(idenso::dirac::AGS.gamma, left, right, argument)
    }
    fn edges(list: &[usize]) -> Vec<EdgeIndex> {
        list.iter().map(|edge| EdgeIndex::from(*edge)).collect()
    }

    /// Today's path: Taylor, contraction, vacuum projection, moments.
    fn today(expression: &Atom, soft: &[usize], hard: &[usize], degree: usize) -> Result<Vec<(Atom, Atom)>> {
        let lambda = symbol!("rstar_hard_taylor"; Scalar);
        let scaled = RStarIntegrand::scale_edges(expression, &edges(soft), lambda);
        let polynomial = TaylorExpansion::new(lambda, degree).evaluate(&scaled)?;
        let simplified = RStarIntegrand::simplify(&polynomial)?;
        let projected = project_hard_vacuum(&simplified, &edges(hard))?;
        hard_moments(&projected, &edges(soft))
    }

    fn closed(
        expression: &Atom,
        soft: &[usize],
        degree: usize,
    ) -> Result<(std::result::Result<Vec<(Atom, Atom)>, String>, ClosedStats)> {
        let mut stats = ClosedStats::default();
        let moments = closed_moments(
            expression,
            &edges(soft),
            degree,
            &|expression| RStarIntegrand::simplify_view_certified(expression.as_view()),
            &mut stats,
        )?;
        Ok((moments, stats))
    }

    fn value(rows: &[(Atom, Atom)]) -> Atom {
        rows.iter().fold(Atom::Zero, |sum, (coefficient, insertion)| {
            sum + random_point(coefficient) * random_point(insertion)
        })
    }

    /// The closed moments equal today's at a random rational point, and
    /// neither is zero.
    fn assert_same(expression: &Atom, soft: &[usize], hard: &[usize], degree: usize) -> Result<ClosedStats> {
        let reference = today(expression, soft, hard, degree)?;
        let (moments, stats) = closed(expression, soft, degree)?;
        let moments = moments.map_err(|reason| eyre!("declined: {reason}"))?;
        assert!(!value(&reference).is_zero(), "the reference is zero: {reference:?}");
        assert!(
            (value(&moments) - value(&reference)).is_zero(),
            "closed {moments:?}\ntoday {reference:?}"
        );
        Ok(stats)
    }

    /// Mixed Q_a.Q_b (b != a), a slot vector Q_a(mu), and the multidegree
    /// normalization prod_a m_a! = 1 (not |m|! = 2) of the {a,b} sector.
    #[test]
    fn mixed_soft_momenta_rank_four() -> Result<()> {
        test_initialise()?;
        let (mu, nu) = (slot(10), slot(11));
        let k = vector(0);
        let expression = indexed(0, &mu) * indexed(0, &nu) * dot(&k, &vector(4)) * dot(&k, &vector(5))
            + dot(&vector(4), &vector(5)) * dot(&k, &k) * function!(ETS.metric, &mu, &nu)
            + indexed(4, &mu) * dot(&k, &vector(5)) * indexed(0, &nu)
            + dot(&k, &k) * function!(ETS.metric, &mu, &nu);
        let stats = assert_same(&expression, &[4, 5], &[0], 2)?;
        assert!(stats.closures > 0);
        Ok(())
    }

    /// The {a,a} sector: two derivatives give 2 (k.Q)^2 -> 1/2! restores it.
    #[test]
    fn repeated_soft_momentum_normalization() -> Result<()> {
        test_initialise()?;
        let (mu, nu) = (slot(10), slot(11));
        let k = vector(0);
        let q = vector(4);
        let expression = dot(&k, &q).pow(2) * indexed(0, &mu) * indexed(0, &nu)
            + dot(&q, &q) * dot(&k, &k) * function!(ETS.metric, &mu, &nu)
            + indexed(4, &mu) * indexed(4, &nu) * dot(&k, &k);
        assert_same(&expression, &[4], &[0], 2)?;
        Ok(())
    }

    /// A slashed soft momentum inside a closed Dirac trace (gamma.Q_a ->
    /// gamma^alpha after polarization), with a hard slash and open mu, nu.
    #[test]
    fn slashed_soft_momentum_in_a_trace() -> Result<()> {
        test_initialise()?;
        let (mu, nu) = (slot(10), slot(11));
        let (a, b, c, e) = (spinor(20), spinor(21), spinor(22), spinor(23));
        let k = vector(0);
        let trace = gamma(&a, &b, &mu) * gamma(&b, &c, &vector(4)) * gamma(&c, &e, &k) * gamma(&e, &a, &nu);
        let expression = trace * dot(&k, &vector(5));
        let stats = assert_same(&expression, &[4, 5], &[0], 2)?;
        assert_eq!(stats.slots, 2);
        Ok(())
    }

    /// A rank-6 sector declines (fails closed), it is never truncated.
    #[test]
    fn rank_above_four_declines() -> Result<()> {
        test_initialise()?;
        let slots = [slot(10), slot(11), slot(12), slot(13)];
        let k = vector(0);
        let expression = slots.iter().fold(dot(&k, &vector(4)).pow(2), |product, slot| product * indexed(0, slot));
        let (moments, _) = closed(&expression, &[4], 2)?;
        let reason = moments.expect_err("a rank-6 sector must decline");
        assert!(reason.starts_with("rank"), "{reason}");
        Ok(())
    }

    /// Open spinor slots (a quark leg) and parity-odd content decline.
    #[test]
    fn spinor_legs_and_parity_odd_decline() -> Result<()> {
        test_initialise()?;
        let mu = slot(10);
        let k = vector(0);
        let open_chain = gamma(&spinor(20), &spinor(21), &mu) * dot(&k, &vector(4));
        let (moments, _) = closed(&open_chain, &[4], 1)?;
        assert!(moments.expect_err("open spinor slots").starts_with("slot_rep"));
        let odd = dot(&k, &k) * function!(ETS.metric, &mu, &slot(11))
            * Atom::var(idenso::dirac::AGS.gamma5);
        let (moments, _) = closed(&odd, &[4], 0)?;
        assert_eq!(moments.expect_err("parity odd"), "parity_odd");
        Ok(())
    }
}
