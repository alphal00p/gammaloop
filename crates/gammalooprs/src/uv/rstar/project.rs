//! Vacuum projection of hard Taylor tensors before quotient counterterms.

use eyre::{Result, WrapErr};
use feynkit_tensor::TensorReducer;
use idenso::{CookMode, CookSettings};
use linnet::half_edge::involution::EdgeIndex;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{atom::Atom, function};

use crate::{utils::GS, uv::approx::integrated::Integrated};

/// Average hard loop directions while retaining their scalar integrals.
///
/// The resulting complete local covariants can be inserted in the next UV
/// counterterm. For example, `(k.q)^2` becomes `k^2 q^2/d`: the scalar hard
/// factor is `k^2/d`, and the insertion is the whole polynomial `q^2`.
/// Projecting only after a quotient pole operation would instead expose its
/// open-index counterterm and give a different evanescent contribution.
pub(super) fn project_hard_vacuum(expression: &Atom, hard_edges: &[EdgeIndex]) -> Result<Atom> {
    if hard_edges.is_empty() || expression.is_zero() {
        return Ok(expression.clone());
    }
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
    let mut reducer = TensorReducer::new(Atom::var(GS.dim));
    for edge in hard_edges {
        reducer = reducer.with_integrated_vector(function!(
            GS.emr_mom,
            usize::from(*edge),
            &minkowski
        ))?;
    }
    let source = reducer.admit(cooking.try_cook(expression.as_view())?)?;
    let reduced = reducer
        .reduce(&source)
        .wrap_err("could not project an R* hard Taylor coefficient onto vacuum covariants")?
        .into_expression()?;
    let result = cooking.uncook(reduced.expression().as_view());
    // Rational projector coefficients must participate in subsequent Laurent
    // expansions. The old R helper preserves d in Lorentz slot metadata.
    Ok(Integrated::dimensionally_regularized(&result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{initialisation::test_initialise, numerator::aind::Aind};
    use idenso::representations::ColorAdjoint;
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::slot::ParseableAind,
    };
    use symbolica::{atom::AtomCore, symbol};

    fn momentum(edge: usize) -> Atom {
        function!(
            GS.emr_mom,
            edge,
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        )
    }

    fn dot(left: &Atom, right: &Atom) -> Atom {
        function!(SPENSO_TAG.dot, left, right)
    }

    #[test]
    fn quadratic_and_quartic_hard_moments_are_covariant() -> Result<()> {
        test_initialise()?;
        let hard = momentum(0);
        let soft = momentum(4);
        let regulated = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let mixed = dot(&hard, &soft);
        let invariant = dot(&hard, &hard) * dot(&soft, &soft);
        let quadratic = project_hard_vacuum(&mixed.pow(2), &[EdgeIndex::from(0)])?;
        assert!((quadratic - &invariant / &regulated).together().is_zero());
        let quartic = project_hard_vacuum(&mixed.pow(4), &[EdgeIndex::from(0)])?;
        let expected = Atom::num(3) * invariant.pow(2) / (&regulated * (&regulated + Atom::num(2)));
        assert!((quartic - expected).together().is_zero());
        let local = dot(&soft, &soft);
        assert_eq!(project_hard_vacuum(&local, &[EdgeIndex::from(0)])?, local);
        Ok(())
    }

    #[test]
    fn projection_preserves_propagators_and_attachment_indices() -> Result<()> {
        test_initialise()?;
        let hard = momentum(0);
        let soft = momentum(4);
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Aind::Hedge(7, 0).to_atom()]);
        let nu = minkowski.to_symbolic([Aind::Hedge(8, 0).to_atom()]);
        let colour = ColorAdjoint {}.new_rep(8);
        let a = colour.to_symbolic([Aind::Hedge(7, 1).to_atom()]);
        let b = colour.to_symbolic([Aind::Hedge(8, 1).to_atom()]);
        let attachment = function!(ETS.metric, &mu, &nu) * function!(ETS.metric, a, b);
        let mass = Atom::var(symbol!("rstar_projection_mass_squared"; Scalar));
        let denominator = function!(
            GS.den,
            0,
            function!(GS.emr_mom, 0),
            mass.sqrt(),
            dot(&hard, &hard) - &mass
        );
        let source = &attachment * dot(&hard, &soft).pow(2) / denominator.pow(3);
        let reduced = project_hard_vacuum(&source, &[EdgeIndex::from(0)])?;
        let regulated = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let expected =
            attachment * dot(&hard, &hard) * dot(&soft, &soft) / (regulated * denominator.pow(3));
        assert!((reduced - expected).together().is_zero());
        Ok(())
    }

    #[test]
    fn separate_hard_loop_carriers_share_the_vacuum_projection() -> Result<()> {
        test_initialise()?;
        let left = momentum(0);
        let right = momentum(1);
        let soft = momentum(4);
        let source = dot(&left, &soft) * dot(&right, &soft);
        let regulated = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let expected = dot(&left, &right) * dot(&soft, &soft) / regulated;
        let reduced = project_hard_vacuum(&source, &[EdgeIndex::from(0), EdgeIndex::from(1)])?;
        assert!((reduced - expected).together().is_zero());
        Ok(())
    }
}
