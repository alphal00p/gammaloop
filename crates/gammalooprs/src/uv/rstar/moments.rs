//! Scalar hard coefficients and complete covariant quotient insertions.

use eyre::{Result, WrapErr, bail, ensure};
use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use linnet::half_edge::involution::EdgeIndex;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing::TensorCollectFilter,
    structure::{TensorStructure, partial::PartialStructure},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::utils::GS;

fn is_momentum(value: AtomView<'_>) -> bool {
    matches!(value, AtomView::Fun(momentum)
        if [GS.emr_mom, GS.loop_mom, GS.external_mom].contains(&momentum.get_symbol()))
}

fn is_soft_momentum(value: AtomView<'_>, soft_edges: &[EdgeIndex]) -> bool {
    let AtomView::Fun(momentum) = value else {
        return false;
    };
    momentum.get_symbol() == GS.emr_mom
        && momentum.get_nargs() >= 2
        && usize::try_from(momentum.get(0))
            .is_ok_and(|edge| soft_edges.contains(&EdgeIndex::from(edge)))
}

fn contains_soft_momentum(expression: AtomView<'_>, soft_edges: &[EdgeIndex]) -> bool {
    let mut found = false;
    expression.visitor(&mut |value| {
        found |= is_soft_momentum(value, soft_edges);
        !found
    });
    found
}

/// Vacuum angular projection leaves hard momenta only in scalar invariants.
fn is_hard_invariant(value: AtomView<'_>, soft_edges: &[EdgeIndex]) -> bool {
    matches!(value, AtomView::Fun(dot)
        if [SPENSO_TAG.dot, ETS.metric].contains(&dot.get_symbol())
            && dot.get_nargs() == 2
            && dot.iter().all(|vector| is_momentum(vector) && !is_soft_momentum(vector, soft_edges)))
}

/// Momentum of a dimensionless Euclidean self-square raised to -m epsilon.
pub(super) fn analytic_line_edge(value: AtomView<'_>) -> Option<EdgeIndex> {
    let AtomView::Pow(power) = value else {
        return None;
    };
    if !power.get_exp().contains_symbol(GS.dim_epsilon) {
        return None;
    }
    let count = -power.get_exp().to_owned() / Atom::var(GS.dim_epsilon);
    if i64::try_from(count.as_view()).ok()? <= 0 {
        return None;
    }
    let square = -power.get_base().to_owned() * Atom::var(GS.mu_r_sq);
    let AtomView::Fun(dot) = square.as_view() else {
        return None;
    };
    if dot.get_symbol() != SPENSO_TAG.dot || dot.get_nargs() != 2 || dot.get(0) != dot.get(1) {
        return None;
    }
    let AtomView::Fun(q) = dot.get(0) else {
        return None;
    };
    if q.get_symbol() != GS.emr_mom || q.get_nargs() != 2 {
        return None;
    }
    usize::try_from(q.get(0)).ok().map(EdgeIndex::from)
}

/// The covariant collector can put a noninteger power of a closed dot in
/// the tensor monomial. Its hard line weight belongs outside the quotient K.
fn separate_hard_line_weights(insertion: Atom, soft: &[EdgeIndex]) -> (Atom, Atom) {
    if !insertion.contains_symbol(GS.dim_epsilon) || !insertion.contains_symbol(GS.emr_mom) {
        return (Atom::one(), insertion);
    }
    let is_hard_weight = |part| analytic_line_edge(part).is_some_and(|edge| !soft.contains(&edge));
    if is_hard_weight(insertion.as_view()) {
        return (insertion, Atom::one());
    }
    if let AtomView::Mul(product) = insertion.as_view()
        && product.iter().any(is_hard_weight)
    {
        let mut weight = Atom::one();
        let mut soft_insertion = Atom::one();
        for factor in product.iter() {
            if is_hard_weight(factor) {
                weight *= factor.to_owned();
            } else {
                soft_insertion *= factor.to_owned();
            }
        }
        return (weight, soft_insertion);
    }
    (Atom::one(), insertion)
}

/// Return `(scalar hard coefficient, complete covariant insertion)` pairs.
///
/// Hard loops must first undergo vacuum tensor projection. Their scalar
/// invariants, denominators and Laurent coefficients remain outside the next
/// quotient's pole projection; all soft vectors and attachment tensors enter
/// its local vertex. This gives bare hard terms and already integrated local
/// counterterms the same covariant basis. Splitting unprojected hard components
/// instead would move dimension-dependent traces across the quotient's `K`.
pub(super) fn hard_moments(
    expression: &Atom,
    soft_edges: &[EdgeIndex],
) -> Result<Vec<(Atom, Atom)>> {
    if expression.is_zero() {
        return Ok(Vec::new());
    }

    let mut unprojected = None;
    let mut soft_denominator = false;
    expression.visitor(&mut |value| {
        if unprojected.is_some() || soft_denominator {
            return false;
        }
        if matches!(value, AtomView::Fun(function) if function.get_symbol() == GS.den) {
            soft_denominator = contains_soft_momentum(value, soft_edges);
            return false;
        }
        if is_hard_invariant(value, soft_edges) {
            return false;
        }
        if is_momentum(value) && !is_soft_momentum(value, soft_edges) {
            unprojected = Some(value.to_owned());
            return false;
        }
        true
    });
    ensure!(
        !soft_denominator,
        "a soft R* momentum remains in a hard Taylor denominator"
    );
    if let Some(unprojected) = unprojected {
        // A vanishing hard part can keep a closed hard dummy between two sum
        // factors, which contraction cannot close without distributing them.
        // Literal zero children upstream should leave no such term.
        if expression.expand().is_zero() {
            tracing::warn!(
                uv = true,
                rstar = true,
                "an R* hard part vanishes only after expansion; it has no hard moments"
            );
            return Ok(Vec::new());
        }
        bail!(
            "R* hard moments require vacuum tensor projection before quotient UV maps: {unprojected}"
        );
    }

    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let tensor: SymbolicTensor<PartialStructure> =
        SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view())?)?.into();
    let rows = tensor
        .coefficient_list_with(|value| {
            !is_hard_invariant(value, soft_edges)
                && !matches!(value, AtomView::Fun(function) if function.get_symbol() == GS.den)
                && TensorCollectFilter::<0>::Tensors.matches(value)
        })
        .wrap_err("could not separate scalar hard R* coefficients from covariant insertions")?;
    rows.into_iter()
        .map(|(insertion, coefficient)| {
            ensure!(
                coefficient.structure().canonical().is_scalar(),
                "an open tensor remains in an R* hard coefficient: {}",
                coefficient.expression()
            );
            let insertion = cooking.uncook(insertion.expression().as_view());
            let (hard_weight, insertion) = separate_hard_line_weights(insertion, soft_edges);
            let coefficient = hard_weight * cooking.uncook(coefficient.expression().as_view());
            ensure!(
                !contains_soft_momentum(coefficient.as_view(), soft_edges),
                "a soft R* momentum remains in the hard coefficient: {coefficient}"
            );
            let mut has_hard_momentum = false;
            insertion.visitor(&mut |value| {
                has_hard_momentum |= is_momentum(value) && !is_soft_momentum(value, soft_edges);
                !has_hard_momentum
            });
            ensure!(
                !has_hard_momentum,
                "a hard R* momentum remains in the soft insertion: {insertion}"
            );
            Ok((coefficient, insertion))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{initialisation::test_initialise, numerator::aind::Aind};
    use idenso::{representations::ColorAdjoint, tensor::ContractSettings};
    use spenso::structure::{
        representation::{Minkowski, RepName},
        slot::ParseableAind,
    };
    use symbolica::{function, parse};

    fn normalized(expression: &Atom) -> Result<Atom> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view())?)?
            .contract(ContractSettings::default())?
            .to_dots()?
            .canonize()?;
        Ok(cooking.uncook(tensor.expression().as_view()).expand())
    }

    fn reconstruct(rows: &[(Atom, Atom)]) -> Atom {
        rows.iter()
            .fold(Atom::Zero, |sum, (coefficient, insertion)| {
                sum + coefficient * insertion
            })
    }

    #[test]
    fn projected_moments_put_complete_attachment_and_colour_tensors_inside() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::Hedge(0, 0).to_atom()]);
        let nu = mink.to_symbolic([Aind::Hedge(1, 0).to_atom()]);
        let alpha = mink.to_symbolic([Aind::Hedge(2, 0).to_atom()]);
        let beta = mink.to_symbolic([Aind::Hedge(3, 0).to_atom()]);
        let colour = ColorAdjoint {}.new_rep(8);
        let a = colour.to_symbolic([Aind::Hedge(2, 1).to_atom()]);
        let b = colour.to_symbolic([Aind::Hedge(3, 1).to_atom()]);
        let attachment = function!(ETS.metric, &alpha, &beta) * function!(ETS.metric, a, b);
        let hard = function!(GS.emr_mom, 0, mink.to_symbolic([]));
        let scalar = function!(SPENSO_TAG.dot, &hard, &hard) / Atom::var(GS.dim_epsilon);
        let insertion = &attachment * function!(GS.emr_mom, 4, &mu) * function!(GS.emr_mom, 5, &nu);
        let source = &scalar * &insertion;
        let rows = hard_moments(&source, &[EdgeIndex::from(4), EdgeIndex::from(5)])?;
        assert_eq!(rows.len(), 1);
        assert_eq!(normalized(&rows[0].0)?, normalized(&scalar)?);
        assert_eq!(normalized(&rows[0].1)?, normalized(&insertion)?);
        assert_eq!(normalized(&reconstruct(&rows))?, normalized(&source)?);
        Ok(())
    }

    #[test]
    fn analytic_weights_on_hard_and_soft_edges_stay_on_their_own_side() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let weight = |edge| {
            let q = function!(GS.emr_mom, edge, &mink);
            (-function!(SPENSO_TAG.dot, &q, &q) / Atom::var(GS.mu_r_sq))
                .pow(-Atom::var(GS.dim_epsilon))
        };
        let hard = weight(3);
        let soft = weight(4);
        let rows = hard_moments(&(&hard * &soft), &[EdgeIndex::from(4)])?;
        assert_eq!(rows, vec![(hard, soft)]);
        Ok(())
    }

    #[test]
    fn local_metric_counterterms_use_a_covariant_insertion_without_soft_momenta() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::Hedge(0, 0).to_atom()]);
        let nu = mink.to_symbolic([Aind::Hedge(1, 0).to_atom()]);
        let metric = function!(ETS.metric, mu, nu);
        let coefficient = Atom::one() / Atom::var(GS.dim_epsilon);
        let rows = hard_moments(&(&coefficient * &metric), &[])?;
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].0, coefficient);
        assert_eq!(normalized(&rows[0].1)?, normalized(&metric)?);
        Ok(())
    }

    #[test]
    fn mixed_hard_moments_are_projected_before_collecting_the_soft_square() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let hard = function!(GS.emr_mom, 0, &mink);
        let soft = function!(GS.emr_mom, 4, &mink);
        let mixed = function!(SPENSO_TAG.dot, &hard, &soft).pow(2);
        assert!(hard_moments(&mixed, &[EdgeIndex::from(4)]).is_err());
        let projected = super::super::project::project_hard_vacuum(&mixed, &[EdgeIndex::from(0)])?;
        let rows = hard_moments(&projected, &[EdgeIndex::from(4)])?;
        assert_eq!(rows.len(), 1);
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let coefficient = function!(SPENSO_TAG.dot, &hard, &hard) / dimension;
        let soft_square = function!(SPENSO_TAG.dot, &soft, &soft);
        assert_eq!(normalized(&rows[0].0)?, normalized(&coefficient)?);
        assert_eq!(normalized(&rows[0].1)?, normalized(&soft_square)?);
        assert_eq!(normalized(&reconstruct(&rows))?, normalized(&projected)?);
        Ok(())
    }

    #[test]
    fn scalar_coefficients_and_zero_need_no_covariant_insertion() -> Result<()> {
        test_initialise()?;
        let coefficient = parse!("a/(eps*(q2-m2)^2)");
        assert_eq!(
            hard_moments(&coefficient, &[EdgeIndex::from(4)])?,
            vec![(coefficient, Atom::one())]
        );
        assert!(hard_moments(&Atom::Zero, &[EdgeIndex::from(4)])?.is_empty());
        Ok(())
    }

    /// The shape of B_self_energy_L5_105327's failing hard part: a child
    /// counterterm `X/epsilon + (-X)/epsilon` with differently distributed
    /// sums, whose hard vectors close only against a second sum factor.
    #[test]
    fn hard_part_vanishing_after_expansion_has_no_moments() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let index = |aind: Aind| mink.to_symbolic([aind.to_atom()]);
        let (dummy, port, open) = (
            index(Aind::Dummy(0)),
            index(Aind::Hedge(3, 1)),
            index(Aind::Hedge(0, 1)),
        );
        let hard = |slot: &Atom| function!(GS.emr_mom, 2, slot);
        let vector = hard(&mink.to_symbolic([]));
        let hard_square = function!(SPENSO_TAG.dot, &vector, &vector);
        let metric = function!(ETS.metric, &dummy, &port);
        let tensor = hard(&dummy) * hard(&port);
        let trace = &hard_square * &metric;
        let epsilon = Atom::var(GS.dim_epsilon);
        let child = (&tensor - &trace) / &epsilon + (&trace - &tensor) / &epsilon;
        assert!(!child.is_zero() && child.expand().is_zero());
        let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let projector = &metric - &dimension * &metric;
        let soft = function!(GS.emr_mom, 4, &open);
        let soft_edges = [EdgeIndex::from(4)];
        assert!(hard_moments(&(&projector * &child * &soft), &soft_edges)?.is_empty());
        // A nonzero hard part of the same shape still needs its projection.
        let nonzero = &projector * (&tensor - &trace) / &epsilon * &soft;
        assert!(hard_moments(&nonzero, &soft_edges).is_err());
        Ok(())
    }
}
