//! Native auxiliary angular projector for the quadratic gluon two-point term.
//! The direction of q is averaged at fixed q²; q is NOT a sixth physical loop.

use super::*;
use feynkit_tensor::TensorReducer;
use rustred::algebra::Coefficient;
use spenso::structure::representation::{BaseRepName, Minkowski};
use symbolica::{
    atom::Symbol,
    prelude::{Q, Z},
};

pub(super) fn algebra() -> AlgebraSettings {
    AlgebraSettings {
        contract: AlgebraContraction::Dots,
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        ..AlgebraSettings::hep()
    }
}

pub(super) fn scalar(expression: Atom) -> Result<Atom> {
    let tensor = SymbolicTensor::infer(expression)?.simplify_algebra(&algebra())?;
    if !tensor.is_scalar() {
        return Err("native gamma/color contraction left open tensor ports".into());
    }
    Ok(tensor.into_expression())
}

fn angular(expression: &Atom, q: Symbol, d: &Atom, guards: &mut BTreeSet<Atom>) -> Result<Atom> {
    if expression.is_zero() {
        return Ok(Atom::Zero);
    }
    let rep = Minkowski::selfless_symbol().call(d);
    let reduction = TensorReducer::new(d.clone())
        .with_integrated_vector(q.call(rep))
        .reduce(expression.as_view())?;
    if reduction.max_rank() > 4 {
        return Err("the quadratic transverse projector exceeded auxiliary q rank four".into());
    }
    for term in reduction.terms() {
        let coefficient: Coefficient = term
            .coefficient()
            .try_to_rational_polynomial(&Q, &Z, None)?;
        guards.insert(coefficient.denominator.to_expression());
    }
    let radial = Kinematics::in_dimension(d)?
        .with_momenta([Atom::var(q)])?
        .scalar_product(&Atom::var(q), &Atom::var(q))?;
    // For a homogeneous rank-two/four source, division by q²/q⁴ in the
    // transverse projector is exactly evaluation of the remaining radial
    // invariant at one AFTER angular contraction. No physical limit is taken
    // on an uncontracted tensor, and no extra loop normalization is introduced.
    let result = scalar(reduction.into_expression())?
        .replace(radial)
        .with(Atom::one())
        .together();
    if result.contains_symbol(q) {
        return Err("auxiliary angular projection retained external q dependence".into());
    }
    Ok(result)
}

pub(super) fn quadratic_transverse(
    metric_part: &Atom,
    longitudinal_part: &Atom,
    q: Symbol,
    d: &Atom,
    adjoint_dimension: i64,
    guards: &mut BTreeSet<Atom>,
) -> Result<Atom> {
    guards.insert(d.clone());
    guards.insert(d - Atom::one());
    guards.insert(d + Atom::num(2));
    let metric = angular(metric_part, q, d, guards)?;
    let longitudinal = angular(longitudinal_part, q, d, guards)?;
    Ok(((metric - longitudinal) / (Atom::num(adjoint_dimension) * (d - Atom::one()))).together())
}

/// Independent analytic controls are rerun by the preparation command.
pub(super) fn controls(d: &Atom) -> Result<()> {
    let q = spenso::vector_symbol!("qcd_ladder::q");
    let k = spenso::vector_symbol!("qcd_ladder_control::k");
    let p = spenso::vector_symbol!("qcd_ladder_control::p");
    let mu = Minkowski::selfless_symbol().call((d, Atom::var(symbol!("qcd_ladder_control::mu"))));
    let nu = Minkowski::selfless_symbol().call((d, Atom::var(symbol!("qcd_ladder_control::nu"))));
    let kin =
        Kinematics::in_dimension(d)?.with_momenta([Atom::var(q), Atom::var(k), Atom::var(p)])?;
    let dot = |a, b| kin.scalar_product(&Atom::var(a), &Atom::var(b));
    let a = Atom::var(symbol!("qcd_ladder_control::A"));
    let b = Atom::var(symbol!("qcd_ladder_control::B"));
    let tensor = &a * dot(q, q)? * ETS.metric(&mu, &nu) + &b * q.call(&mu) * q.call(&nu);
    let metric = scalar(&tensor * ETS.metric(&mu, &nu))?;
    let longitudinal = scalar(&tensor * q.call(&mu) * q.call(&nu))?;
    let mut guards = BTreeSet::new();
    let value = quadratic_transverse(&metric, &longitudinal, q, d, 1, &mut guards)?;
    if !(value - a).expand().together().is_zero() {
        return Err("native projector control did not retain A and remove B".into());
    }

    // This genuinely uses native rank-two AND rank-four q moments with other
    // vectors as spectators, and matches the local derivative projector.
    let anisotropic = k.call(&mu) * k.call(&nu) * dot(q, p)?.pow(2);
    let metric = scalar(&anisotropic * ETS.metric(&mu, &nu))?;
    let longitudinal = scalar(&anisotropic * q.call(&mu) * q.call(&nu))?;
    let value = quadratic_transverse(&metric, &longitudinal, q, d, 1, &mut guards)?;
    let expected = ((d + Atom::one()) * dot(k, k)? * dot(p, p)? - Atom::num(2) * dot(k, p)?.pow(2))
        / (d * (d - Atom::one()) * (d + Atom::num(2)));
    if !(value - expected).expand().together().is_zero() {
        return Err(
            "native auxiliary angular average disagrees with the local quadratic projector".into(),
        );
    }
    Ok(())
}
