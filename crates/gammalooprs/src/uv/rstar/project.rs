//! Vacuum projection of hard Taylor tensors before quotient counterterms.

use eyre::{Result, WrapErr};
use feynkit_tensor::TensorReducer;
use idenso::{
    CookMode, CookSettings,
    dirac::GammaSimplifySettings,
    shorthands::UndoShorthands,
    tensor::{AlgebraSettings, SymbolicNetExt, SymbolicNetParse},
};
use linnet::half_edge::involution::EdgeIndex;
use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::representation::{Minkowski, RepName},
};
use symbolica::{atom::Atom, function};

use crate::{numerator::aind::Aind, utils::GS, uv::approx::integrated::Integrated};

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
    // The hard simplifier can collect slash momenta inside Dirac chains.
    // Expose them before counting tensor rank: k_mu slash(k) has rank two,
    // although only one vector is visible outside its compact gamma tensor.
    // Use the same native shorthand expansion as the vacuum-integral bridge.
    let explicit = expression.undo_trace::<Aind>()?.undo_chain::<Aind>()?;
    let explicit = explicit
        .parse_to_symbolic_net::<Aind>(&ParseSettings {
            shorthand_parsing: ShorthandParsing::Expand {
                schoonschip: SchoonschipExpansionMode {
                    inner_products: false,
                    expand_inside_chains: true,
                    expand_schoonship: true,
                },
                trace: true,
                chain: true,
            },
            ..Default::default()
        })?
        .simple_execute::<()>()?;
    let source = reducer.admit(cooking.try_cook(explicit.as_view())?)?;
    let reduced = reducer
        .reduce(&source)
        .wrap_err("could not project an R* hard Taylor coefficient onto vacuum covariants")?
        .into_expression()?;
    // Finish scalar and Clifford contractions before separating the hard
    // coefficient from its quotient insertion. In particular, gamma_mu
    // gamma^mu must cancel the projector's 1/d here; leaving it inside the
    // quotient would put that d on the wrong side of its Laurent projection.
    let reduced = reduced
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        })?
        .to_dots()?;
    let result = cooking.uncook(reduced.expression().as_view());
    // Rational projector coefficients must participate in subsequent Laurent
    // expansions. The old R helper preserves d in Lorentz slot metadata.
    Ok(Integrated::dimensionally_regularized(&result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{initialisation::test_initialise, numerator::aind::Aind};
    use idenso::representations::{Bispinor, ColorAdjoint};
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
    fn hard_clifford_contractions_precede_coefficient_separation() -> Result<()> {
        test_initialise()?;
        let hard_left = momentum(0);
        let hard_right = momentum(1);
        let soft = momentum(2);
        let spin = Bispinor {}.new_rep(4);
        let left = spin.to_symbolic([Aind::Hedge(8, 0).to_atom()]);
        let right = spin.to_symbolic([Aind::Hedge(9, 0).to_atom()]);
        let gamma = |argument: &Atom| {
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                argument
            )
        };
        let slash = spenso::chain!(&left, &right, gamma(&soft));
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let scalar = dot(&hard_left, &hard_right);
        for (source, coefficient) in [
            (
                spenso::chain!(
                    &left,
                    &right,
                    gamma(&hard_left),
                    gamma(&hard_right),
                    gamma(&soft)
                ),
                Atom::one(),
            ),
            (
                spenso::chain!(
                    &left,
                    &right,
                    gamma(&hard_left),
                    gamma(&soft),
                    gamma(&hard_right)
                ),
                Atom::num(2) / &dimension - 1,
            ),
        ] {
            let actual = project_hard_vacuum(&source, &[EdgeIndex::from(0), EdgeIndex::from(1)])?;
            let expected = &scalar * coefficient * &slash;
            assert!(
                (&actual - &expected).together().is_zero(),
                "actual={actual}, expected={expected}"
            );
        }
        Ok(())
    }

    #[test]
    fn hard_projection_counts_slash_momenta_inside_dirac_chains() -> Result<()> {
        test_initialise()?;
        let hard = momentum(0);
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::Hedge(7, 0).to_atom()]);
        let spin = Bispinor {}.new_rep(4);
        let left = spin.to_symbolic([Aind::Hedge(8, 0).to_atom()]);
        let right = spin.to_symbolic([Aind::Hedge(9, 0).to_atom()]);
        let colour = ColorAdjoint {}.new_rep(8);
        let a = colour.to_symbolic([Aind::Hedge(8, 1).to_atom()]);
        let b = colour.to_symbolic([Aind::Hedge(9, 1).to_atom()]);
        let attachment = function!(ETS.metric, a, b);
        let slash = function!(idenso::dirac::AGS.gamma, &left, &right, &hard);
        let chain = function!(
            SPENSO_TAG.chain,
            &left,
            &right,
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                &hard
            )
        );
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let expected = &attachment
            * dot(&hard, &hard)
            * function!(idenso::dirac::AGS.gamma, &left, &right, &mu)
            / &dimension;
        for compact in [slash, chain] {
            let source = &attachment * function!(GS.emr_mom, 0, &mu) * compact;
            let projected = project_hard_vacuum(&source, &[EdgeIndex::from(0)])?;
            assert!(!projected.is_zero(), "an even hard tensor was discarded");
            super::super::tests::assert_same_counterterm(&projected, &expected)?;
        }
        Ok(())
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

    #[test]
    fn factored_hard_vector_sums_are_scalar_invariants_before_moment_separation() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::Hedge(20, 0).to_atom()]);
        let nu = mink.to_symbolic([Aind::Hedge(21, 0).to_atom()]);
        let vector = |edge: usize, slot: &Atom| function!(GS.emr_mom, edge, slot);
        let source = (vector(0, &mu) + vector(1, &mu))
            * (vector(0, &nu) + vector(1, &nu))
            * vector(4, &mu)
            * vector(4, &nu);
        let projected = project_hard_vacuum(&source, &[EdgeIndex::from(0), EdgeIndex::from(1)])?;
        let moments = super::super::moments::hard_moments(&projected, &[EdgeIndex::from(4)])?;
        let reconstructed = moments
            .iter()
            .fold(Atom::Zero, |sum, (coefficient, insertion)| {
                sum + coefficient * insertion
            });
        let left = momentum(0);
        let right = momentum(1);
        let soft = momentum(4);
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let expected =
            (dot(&left, &left) + Atom::num(2) * dot(&left, &right) + dot(&right, &right))
                * dot(&soft, &soft)
                / dimension;
        super::super::tests::assert_same_counterterm(&reconstructed, &expected)
    }
}
