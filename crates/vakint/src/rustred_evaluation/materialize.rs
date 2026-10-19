use rustred::algebra::ExactAlgebraLimits;
use rustred::family::IntegralKey;
use rustred::reduction::Reducer;
use rustred::scalar_numerator::FamilyScalarNumeratorService;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use crate::matad::MATAD;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{MATADOptions, VakintSettings};

use super::matching::MatchedScalarFamily;
use super::{RustRedEvaluationError, RustRedEvaluationOptions};

pub(super) fn evaluate(
    vakint: &crate::Vakint,
    matched: MatchedScalarFamily,
    settings: &VakintSettings,
    numerator: AtomView,
    _options: &RustRedEvaluationOptions,
) -> Result<Atom, RustRedEvaluationError> {
    let dimension = vk_parse!("4-2*ep").expect("Vakint's internal dimension expression parses");
    let rustred_dimension = vk_parse!("rustred::d")
        .expect("RustRed's stable coefficient namespace parses in Symbolica");
    let data = vakint.rustred_data()?;
    let snapshot = data.artifact(matched.family)?;
    let artifact = snapshot.artifact();
    let base_integral = IntegralKey::try_new(matched.powers.iter().copied()).map_err(|error| {
        RustRedEvaluationError::IntegralKey {
            detail: error.to_string(),
        }
    })?;
    let loop_momenta = (1..=matched.loop_count)
        .map(|index| function!(S.k, Atom::num(index as i64)))
        .collect::<Vec<_>>();
    let routed_numerator = super::numerator::route_to_parent(
        numerator,
        &loop_momenta,
        matched.parent_coordinates.as_deref(),
        settings,
        super::numerator::ScalarNumeratorForm::ExpandedPolynomial,
    )?;
    let internal_numerator = routed_numerator
        .replace(Atom::var(vk_symbol!(settings.epsilon_symbol.as_str())).to_pattern())
        .with(vk_parse!("ep").unwrap().to_pattern());
    // Lower coordinates and scalar algebra independently of reduction coverage.
    // The artifact-backed service rejects uncertified roots before we can
    // preserve them as exact residuals below.
    let scalar_service = FamilyScalarNumeratorService::try_new(
        artifact.family(),
        S.dot,
        loop_momenta,
        settings.rustred_scalar_numerator_limits,
    )
    .map_err(|error| RustRedEvaluationError::ScalarNumerator {
        detail: error.to_string(),
    })?;
    let lowering = scalar_service
        .lower(&internal_numerator, &base_integral)
        .map_err(|error| RustRedEvaluationError::ScalarNumerator {
            detail: error.to_string(),
        })?;
    // The ordinary reducer is still request-local. An explicitly attached
    // native cache supplies the same typed homogeneity result without creating
    // a second, unused seeded reducer alongside its configured RAM allowance.
    let mut reducer = None;
    let mut raw_masters = Atom::Zero;
    for lowered in lowering.terms() {
        let admitted = match snapshot.program().validate_entry(lowered.integral()) {
            Ok(()) => true,
            Err(
                rustred::reduction::ReductionError::OutsideCertifiedRootDomain { .. }
                | rustred::reduction::ReductionError::OutsideCertifiedTotalExcessDomain { .. },
            ) => false,
            Err(error) => return Err(super::data::failure(error)),
        };
        let terms = if !admitted {
            log::warn!(
                "RustRed {}L input {:?} is outside the certified root domain; there is no guarantee of a complete reduction to the known masters",
                matched.loop_count,
                lowered.integral().powers()
            );
            vec![(
                lowered.integral().clone(),
                artifact.coefficient_context().one(),
            )]
        } else {
            let decomposition =
                match vakint.cached_closed_homogeneous(matched.family, lowered.integral())? {
                    Some(decomposition) => decomposition,
                    None => {
                        if reducer.is_none() {
                            reducer = Some(Reducer::new(artifact).map_err(super::data::failure)?);
                        }
                        reducer
                            .as_mut()
                            .unwrap()
                            .reduce_with_common_mass_homogeneity(lowered.integral())
                            .map_err(super::data::failure)?
                    }
                };
            decomposition
                .terms()
                .iter()
                .map(|(key, coefficient)| {
                    (key.clone(), coefficient.unit_mass_coefficient().clone())
                })
                .collect()
        };
        let exponent = lowered
            .integral()
            .powers()
            .iter()
            .try_fold(2 * matched.loop_count as i128, |sum, &power| {
                sum.checked_sub(i128::from(power))
            })
            .and_then(|sum| sum.checked_add(i128::from(lowered.common_mass_squared_power())))
            .and_then(|sum| i64::try_from(sum).ok())
            .ok_or_else(|| super::data::failure("lower-loop master mass exponent overflow"))?;
        let scale = matched.mass_squared.clone().pow(Atom::num(exponent));
        for (master, coefficient) in terms {
            let combined = artifact
                .coefficient_context()
                .try_mul(
                    lowered.coefficient(),
                    &coefficient,
                    ExactAlgebraLimits::default(),
                )
                .map_err(super::data::failure)?;
            let exact_coefficient = combined
                .to_expression()
                .replace(rustred_dimension.to_pattern())
                .with(dimension.to_pattern());
            let terminal = data
                .masters
                .symbol(artifact.family_fingerprint(), master.powers());
            raw_masters +=
                exact_coefficient * &scale * terminal * lowered.scalar_spectator().clone();
        }
    }

    let loop_count =
        i64::try_from(matched.loop_count).map_err(|_| RustRedEvaluationError::Reduction {
            detail: "matched loop count does not fit Vakint's normalization exponent".to_owned(),
        })?;
    let raw = MATAD::with_settings(settings.clone())
        .finalize_reduced_masters(
            super::combine_exact(raw_masters),
            loop_count,
            None,
            &matched.mass_squared,
            &MATADOptions {
                expand_masters: false,
                susbstitute_masters: false,
                substitute_hpls: false,
                direct_numerical_substition: false,
            },
        )
        .map_err(super::data::failure)?;
    Ok(raw)
}

#[cfg(test)]
#[path = "materialize_tests.rs"]
mod tests;
