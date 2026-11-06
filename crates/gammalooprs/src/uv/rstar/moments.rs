//! Scalar hard coefficients and complete covariant quotient insertions.

use eyre::{Result, WrapErr, ensure};
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
    ensure!(
        unprojected.is_none(),
        "R* hard moments require vacuum tensor projection before quotient UV maps: {}",
        unprojected.unwrap_or(Atom::Zero)
    );

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
            let coefficient = cooking.uncook(coefficient.expression().as_view());
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
}
