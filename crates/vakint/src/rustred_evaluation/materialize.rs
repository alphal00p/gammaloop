use std::collections::BTreeMap;

use rustred::algebra::ExactAlgebraLimits;
use rustred::family::IntegralKey;
use rustred::reduction::Reducer;
use rustred::scalar_numerator::{FamilyScalarNumeratorService, ScalarNumeratorLowering};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use crate::matad::MATAD;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{MATADOptions, VakintSettings};

use super::experimental::timing::MathSpan;
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
    let lowering = lower_input(&matched, settings, numerator, artifact)?;
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
        // Immutable artifacts and scalar lowering precede this boundary.
        // Request-local native state/dispatch setup is execution work in both
        // engines, alongside IBPs and the old weighted terminal assembly.
        // An attached DFS backend can perform reductions on a cache miss.
        let _math = MathSpan::new();
        let cached = if admitted {
            vakint.cached_closed_homogeneous(matched.family, lowered.integral())?
        } else {
            None
        };
        if admitted && cached.is_none() && reducer.is_none() {
            reducer = Some(Reducer::new(artifact).map_err(super::data::failure)?);
        }
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
            let decomposition = match cached {
                Some(decomposition) => decomposition,
                None => reducer
                    .as_mut()
                    .unwrap()
                    .reduce_with_common_mass_homogeneity(lowered.integral())
                    .map_err(super::data::failure)?,
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
    let raw_masters = {
        let _math = MathSpan::new();
        super::combine_exact(raw_masters)
    };
    let raw = MATAD::with_settings(settings.clone())
        .finalize_reduced_masters(
            raw_masters,
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
    // This request-owned memo would naturally drop on return. Time that
    // cleanup without including the preceding MATAD normalization interval.
    {
        let _math = MathSpan::new();
        drop(reducer);
    }
    Ok(raw)
}

fn lower_input(
    matched: &MatchedScalarFamily,
    settings: &VakintSettings,
    numerator: AtomView,
    artifact: &rustred::foundry::artifact::ClosedArtifact,
) -> Result<ScalarNumeratorLowering, RustRedEvaluationError> {
    let (scalar, mut input) = prepare_lower_input(
        matched,
        settings,
        numerator,
        artifact,
        super::numerator::ScalarNumeratorForm::ExpandedPolynomial,
    )?;
    input.numerator = input
        .numerator
        .replace(Atom::var(vk_symbol!(settings.epsilon_symbol.as_str())).to_pattern())
        .with(vk_parse!("ep").unwrap().to_pattern());
    scalar
        .lower(&input.numerator, &input.target)
        .map_err(RustRedEvaluationError::from_native)
}

fn prepare_lower_input<'family>(
    matched: &MatchedScalarFamily,
    settings: &VakintSettings,
    numerator: AtomView,
    artifact: &'family rustred::foundry::artifact::ClosedArtifact,
    form: super::numerator::ScalarNumeratorForm,
) -> Result<
    (
        FamilyScalarNumeratorService<'family>,
        super::numerator::MatchedScalarInput,
    ),
    RustRedEvaluationError,
> {
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
        form,
    )?;
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
    Ok((
        scalar_service,
        super::numerator::MatchedScalarInput {
            target: base_integral,
            mass_squared: matched.mass_squared.clone(),
            numerator: routed_numerator,
        },
    ))
}

pub(super) fn evaluate_ordered_batch(
    vakint: &crate::Vakint,
    settings: &VakintSettings,
    contributions: &[super::MatchedContribution],
    mut observation: Option<&mut rustred::solver::CandidateRoutedApplicationObserver<'_>>,
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<Atom, RustRedEvaluationError> {
    let data = vakint.rustred_data()?;
    let first = contributions
        .first()
        .ok_or_else(|| super::data::failure("empty ordered family"))?;
    let first = MatchedScalarFamily::try_new(&first.matched)?;
    let snapshot = data.artifact(first.family)?;
    let artifact = snapshot.artifact();
    let dimension = vk_parse!("4-2*ep").unwrap();
    let rustred_dimension = vk_parse!("rustred::d").unwrap();
    let mut channels =
        BTreeMap::<(Atom, Atom, i64), Vec<(IntegralKey, rustred::algebra::Coefficient)>>::new();
    let mut anchors = BTreeMap::<Atom, i64>::new();
    let mut raw = BTreeMap::<Atom, Atom>::new();
    for contribution in contributions {
        let matched = MatchedScalarFamily::try_new(&contribution.matched)?;
        if matched.family != first.family {
            return Err(super::data::failure(
                "ordered lower-loop batch changes family",
            ));
        }
        let (scalar, mut input) = prepare_lower_input(
            &matched,
            settings,
            contribution.numerator.as_view(),
            artifact,
            super::numerator::ScalarNumeratorForm::PreserveProducts,
        )?;
        let epsilon = Atom::var(vk_symbol!(settings.epsilon_symbol.as_str()));
        let physical_dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        // The Atom-only lower-loop API has no public guard field. It still
        // rejects an original pole that vanishes under the dimension relation.
        drop(input.input_conditions_after_mapping(&scalar, |condition| {
            condition
                .replace(vk_parse!("ep").unwrap().to_pattern())
                .with(epsilon.to_pattern())
                .replace(rustred_dimension.to_pattern())
                .with(physical_dimension.to_pattern())
        })?);
        // MATAD finalization reserves this exact internal regulator symbol.
        let numerator = input
            .numerator
            .replace(vk_parse!("ep").unwrap().to_pattern())
            .with(epsilon.to_pattern());
        input.numerator = super::numerator::MatchedScalarInput::numerator_in_family_dimension(
            numerator.as_view(),
            settings,
            &rustred_dimension,
        );
        let lowering = scalar
            .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
            .map_err(RustRedEvaluationError::from_native)?;
        for term in lowering.terms() {
            let exponent = term
                .integral()
                .powers()
                .iter()
                .try_fold(2 * matched.loop_count as i128, |sum, power| {
                    sum.checked_sub(i128::from(*power))
                })
                .and_then(|sum| sum.checked_add(i128::from(term.common_mass_squared_power())))
                .and_then(|sum| i64::try_from(sum).ok())
                .ok_or_else(|| super::data::failure("lower-loop master mass exponent overflow"))?;
            let admitted = match snapshot.program().validate_entry(term.integral()) {
                Ok(()) => true,
                Err(
                    rustred::reduction::ReductionError::OutsideCertifiedRootDomain { .. }
                    | rustred::reduction::ReductionError::OutsideCertifiedTotalExcessDomain {
                        ..
                    },
                ) => false,
                Err(error) => return Err(super::data::failure(error)),
            };
            if !admitted {
                log::warn!(
                    "RustRed {}L input {:?} is outside the certified root domain; there is no guarantee of a complete reduction to the known masters",
                    matched.loop_count,
                    term.integral().powers()
                );
                *raw.entry(matched.mass_squared.clone())
                    .or_insert(Atom::Zero) += term
                    .coefficient()
                    .to_expression()
                    .replace(rustred_dimension.to_pattern())
                    .with(dimension.to_pattern())
                    * term.scalar_spectator()
                    * matched.mass_squared.clone().pow(exponent)
                    * data
                        .masters
                        .symbol(artifact.family_fingerprint(), term.integral().powers());
                continue;
            }
            let anchor = *anchors
                .entry(matched.mass_squared.clone())
                .or_insert(exponent);
            let (spectator, exponent) = super::numerator::MatchedScalarInput::mass_channel(
                term.scalar_spectator(),
                &matched.mass_squared,
                exponent,
                anchor,
                settings
                    .rustred_scalar_numerator_limits
                    .diagnostics_enabled
                    .then_some(
                        settings
                            .rustred_scalar_numerator_limits
                            .max_scalar_product_degree,
                    ),
            );
            let (spectator, coefficient) =
                super::numerator::MatchedScalarInput::channel_coefficient(
                    spectator,
                    term.coefficient(),
                    artifact.coefficient_context(),
                    settings
                        .rustred_application_options
                        .reduction_limits
                        .exact_algebra,
                )?;
            channels
                .entry((matched.mass_squared.clone(), spectator, exponent))
                .or_default()
                .push((term.integral().clone(), coefficient));
        }
    }
    let cohort_total = channels.len();
    for (cohort_index, ((mass, spectator, exponent), terms)) in channels.into_iter().enumerate() {
        if let Some(progress) = progress {
            progress.cohorts(cohort_index, cohort_total);
        }
        let _math = MathSpan::new();
        let reduction = match observation.as_deref_mut() {
            Some(observer) => {
                observer.clear_native_summary();
                let mut envelope = rustred::solver::CandidateRoutedApplicationEnvelope::new(
                    rustred::solver::CandidateRoutedApplicationEvent::ActivityStarted,
                    artifact.family_fingerprint(),
                );
                envelope.caller_phase = Some("entry_admission");
                observer.emit(envelope);
                Reducer::reduce_artifact_linear_combination_ordered_with_observation(
                    artifact,
                    terms,
                    settings.rustred_application_options.reduction_limits,
                    &mut |progress, current| {
                        observer.emit_ordered_progress(
                            artifact.family_fingerprint(),
                            *progress,
                            current,
                        );
                        Ok(())
                    },
                )
            }
            None => Reducer::reduce_artifact_linear_combination_ordered(
                artifact,
                terms,
                settings.rustred_application_options.reduction_limits,
            ),
        }
        .map_err(RustRedEvaluationError::from_native)?;
        let masters = reduction
            .terms
            .into_iter()
            .map(|(master, coefficient)| {
                coefficient
                    .to_expression()
                    .replace(rustred_dimension.to_pattern())
                    .with(dimension.to_pattern())
                    * data
                        .masters
                        .symbol(artifact.family_fingerprint(), master.powers())
            })
            .sum::<Atom>();
        let spectator = spectator
            .replace(rustred_dimension.to_pattern())
            .with(dimension.to_pattern());
        *raw.entry(mass.clone()).or_insert(Atom::Zero) += masters * spectator * mass.pow(exponent);
    }
    if let Some(progress) = progress {
        progress.cohorts(cohort_total, cohort_total);
    }
    let mut value = Atom::Zero;
    for (mass, raw) in raw {
        let raw = {
            let _math = MathSpan::new();
            super::combine_exact(raw)
        };
        let normalized = MATAD::with_settings(settings.clone())
            .finalize_reduced_masters(
                raw,
                first.loop_count as i64,
                None,
                &mass,
                &MATADOptions {
                    expand_masters: false,
                    susbstitute_masters: false,
                    substitute_hpls: false,
                    direct_numerical_substition: false,
                },
            )
            .map_err(super::data::failure)?;
        {
            let _math = MathSpan::new();
            value += normalized;
        }
    }
    Ok(value)
}

#[cfg(test)]
#[path = "materialize_tests.rs"]
mod tests;
