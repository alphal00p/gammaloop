//! Native auxiliary angular projector for the quadratic gluon two-point term.
//! The direction of q is averaged at fixed q²; q is NOT a sixth physical loop.

use super::*;
use feynkit_tensor::TensorReducer;
use idenso::shorthands::UndoShorthands;
use rustred::algebra::Coefficient;
use spenso::structure::representation::{BaseRepName, Minkowski};
use spenso::structure::{abstract_index::AbstractIndex, slot::DummyAind};
use symbolica::{
    atom::Symbol,
    prelude::{Q, Z},
};

pub(super) fn scalar(expression: Atom) -> Result<Atom> {
    let reduced = SymbolicTensor::infer(expression)?.simplify_algebra(&AlgebraSettings {
        gamma: Some(Default::default()),
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::Dots,
        ..Default::default()
    })?;
    if reduced.reduction_status() != ReductionStatus::Complete {
        return Err(format!(
            "native tensor contraction is incomplete: {:?}",
            reduced.reduction_status()
        )
        .into());
    }
    let expression = reduced.expression().clone();
    if !expression
        .list_dangling::<spenso::structure::abstract_index::AbstractIndex>()?
        .is_empty()
    {
        return Err("native gamma/color contraction left open tensor ports".into());
    }
    Ok(expression)
}

/// Materialize each scalar-product copy with an independent contracted index.
/// A network power of an already contracted scalar retains that scalar's dummy
/// scope; the indexed reducer instead needs explicit copies of that scope.
fn indexed_dot_polynomial(expression: AtomView<'_>) -> Result<Atom> {
    match expression {
        AtomView::Add(sum) => sum.iter().try_fold(Atom::Zero, |sum, term| {
            Ok(sum + indexed_dot_polynomial(term)?)
        }),
        AtomView::Mul(product) => product.iter().try_fold(Atom::one(), |product, factor| {
            Ok(product * indexed_dot_polynomial(factor)?)
        }),
        AtomView::Pow(power) if power.get_base().contains_symbol(SPENSO_TAG.dot) => {
            let count = usize::try_from(i64::try_from(power.get_exp())?)?;
            (0..count).try_fold(Atom::one(), |product, _| {
                Ok(product * indexed_dot_polynomial(power.get_base())?)
            })
        }
        AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.dot => Ok(expression
            .to_owned()
            .undo_dots::<AbstractIndex>()?
            .canonize::<AbstractIndex>(|_| AbstractIndex::new_dummy())?),
        _ => Ok(expression.to_owned()),
    }
}

fn angular(expression: &Atom, q: Symbol, d: &Atom, guards: &mut BTreeSet<Atom>) -> Result<Atom> {
    if expression.is_zero() {
        return Ok(Atom::Zero);
    }
    let rep = Minkowski::selfless_symbol().call(d);
    // Retain the validated indexed boundary while the compact native path is
    // checked independently against these controls and the complete source.
    let indexed = indexed_dot_polynomial(expression.as_view())?;
    let reduction = TensorReducer::new(d.clone())
        .with_output_term_limit(ANGULAR_OUTPUT_TERM_LIMIT)
        .with_integrated_vector(q.call(rep))?
        .reduce(&SymbolicTensor::infer(indexed)?)?;
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
    let result = scalar(reduction.into_expression()?.expression().clone())?
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

#[cfg(test)]
#[test]
fn local_projector_preserves_exact_angular_dirac_and_color_controls() -> Result<()> {
    controls(&dimension()?)?;
    let expression = Atom::parse(
        "spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(d,mu))*spenso::gamma(spenso::bis(4,b),spenso::bis(4,a),spenso::mink(d,mu))",
        "qcd_ladder_control",
        Default::default(),
    )?;
    let d = Atom::var(symbol!("qcd_ladder_control::d"));
    assert!((scalar(expression)? - Atom::num(4) * d).expand().is_zero());
    let expression = Atom::parse(
        "spenso::t(spenso::coad(8,a),spenso::cof(3,i),spenso::dind(spenso::cof(3,j)))*spenso::t(spenso::coad(8,a),spenso::cof(3,j),spenso::dind(spenso::cof(3,i)))",
        "qcd_ladder_control",
        Default::default(),
    )?;
    assert!((scalar(expression)? - Atom::num(4)).expand().is_zero());
    Ok(())
}

#[cfg(test)]
#[test]
fn scalar_contracts_metric_sums_created_by_six_gamma_trace() -> Result<()> {
    let word = [0, 1, 2, 0, 1, 2];
    let factors = word
        .iter()
        .enumerate()
        .map(|(position, index)| {
            format!(
                "spenso::gamma(spenso::bis(4,s{position}),spenso::bis(4,s{}),spenso::mink(d,mu{index}))",
                (position + 1) % word.len()
            )
        })
        .collect::<Vec<_>>();
    let trace = Atom::parse(factors.join("*"), "qcd_ladder_control", Default::default())?;
    let d = Atom::var(symbol!("qcd_ladder_control::d"));
    // Independent Clifford pairing recurrence gives tr(1)=4 and
    // tr(gamma_mu gamma_nu gamma_rho gamma_mu gamma_nu gamma_rho).
    let expected = Atom::num(4) * &d * (-d.clone().pow(2) + Atom::num(6) * d - Atom::num(4));
    assert!((scalar(trace)? - expected).expand().together().is_zero());
    Ok(())
}

#[cfg(test)]
#[test]
fn scalar_contracts_metric_prerequisites_before_twenty_gamma_trace() -> Result<()> {
    let mut factors = (0..20)
        .map(|position| {
            format!(
                "spenso::gamma(spenso::bis(4,t{position}),spenso::bis(4,t{}),spenso::mink(d,m{position}))",
                (position + 1) % 20
            )
        })
        .collect::<Vec<_>>();
    for pair in 0..10 {
        factors.push(format!(
            "spenso::g(spenso::mink(d,m{}),spenso::mink(d,m{}))",
            2 * pair,
            2 * pair + 1
        ));
    }
    let trace = Atom::parse(factors.join("*"), "qcd_ladder_control", Default::default())?;
    let d = Atom::var(symbol!("qcd_ladder_control::d"));
    // Each metric identifies an adjacent gamma pair, giving d times the
    // identity. Ten pairs therefore leave tr(1)*d^10, with tr(1)=4.
    let expected = Atom::num(4) * d.pow(10);
    assert!((scalar(trace)? - expected).expand().together().is_zero());
    Ok(())
}

#[cfg(test)]
#[test]
fn scalar_preserves_factored_vector_coefficients_before_trace_collection() -> Result<()> {
    let d = Atom::var(symbol!("qcd_ladder_factored_control::D"));
    let p = spenso::vector_symbol!("qcd_ladder_factored_control::p");
    let q = spenso::vector_symbol!("qcd_ladder_factored_control::q");
    let expression = Atom::parse(
        "(p(spenso::mink(D,mu))*q(spenso::mink(D,nu))+p(spenso::mink(D,nu))*q(spenso::mink(D,mu)))*spenso::gamma(spenso::bis(4,s0),spenso::bis(4,s1),spenso::mink(D,mu))*spenso::gamma(spenso::bis(4,s1),spenso::bis(4,s2),spenso::mink(D,a))*spenso::gamma(spenso::bis(4,s2),spenso::bis(4,s3),spenso::mink(D,nu))*spenso::gamma(spenso::bis(4,s3),spenso::bis(4,s0),spenso::mink(D,a))",
        "qcd_ladder_factored_control",
        Default::default(),
    )?;
    let kinematics = Kinematics::in_dimension(&d)?.with_momenta([Atom::var(p), Atom::var(q)])?;
    // Tr(gamma_mu gamma_a gamma_nu gamma_a) = 4(2-D)g_mu_nu.
    let expected = Atom::num(8)
        * (Atom::num(2) - d)
        * kinematics.scalar_product(&Atom::var(p), &Atom::var(q))?;
    assert!(
        (scalar(expression)? - expected)
            .expand()
            .together()
            .is_zero()
    );
    Ok(())
}
