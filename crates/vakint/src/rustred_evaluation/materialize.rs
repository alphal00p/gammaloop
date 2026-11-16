use std::collections::BTreeMap;

use rustred::family::IntegralKey;
use rustred::reduction::Reducer;
use rustred::scalar_numerator::FamilyScalarNumeratorService;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use crate::matad::MATAD;
use crate::stage_timing::StageSpan;
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{MATADOptions, VakintSettings};

use super::RustRedEvaluationError;
use super::experimental::timing::MathSpan;
use super::matching::MatchedScalarFamily;
use super::request_dump::{self, DumpChannel, DumpProgram, DumpRequest, RequestDump};

#[cfg(feature = "weighted-lowering-diagnostics")]
static CAPTURED_GENERIC_SCALAR_INPUTS: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);

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
    let _stage = StageSpan::new("rustred.lower.prepare_input");
    let base_integral = IntegralKey::try_new(matched.powers.iter().copied()).map_err(|error| {
        RustRedEvaluationError::IntegralKey {
            detail: error.to_string(),
        }
    })?;
    let loop_momenta = (1..=matched.loop_count)
        .map(|index| function!(S.k, Atom::num(index as i64)))
        .collect::<Vec<_>>();
    let routed_numerator = {
        let _stage = StageSpan::new("rustred.lower.route_numerator");
        super::numerator::route_to_parent(
            numerator,
            &loop_momenta,
            matched.parent_coordinates.as_deref(),
            settings,
            form,
        )?
    };
    // Lower coordinates and scalar algebra independently of reduction coverage.
    // The artifact-backed service rejects uncertified roots before we can
    // preserve them as exact residuals below.
    let scalar_service = {
        let _stage = StageSpan::new("rustred.lower.scalar_service");
        FamilyScalarNumeratorService::try_new(
            artifact.family(),
            S.dot,
            loop_momenta,
            settings.rustred_scalar_numerator_limits,
        )
        .map_err(|error| RustRedEvaluationError::ScalarNumerator {
            detail: error.to_string(),
        })?
    };
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
    let _stage = StageSpan::new("rustred.lower.ordered_batch");
    let setup = StageSpan::new("rustred.lower.artifact_and_setup");
    let data = vakint.rustred_data()?;
    let first = contributions
        .first()
        .ok_or_else(|| super::data::failure("empty ordered family"))?;
    let first = MatchedScalarFamily::try_new(&first.matched)?;
    let snapshot = data.artifact(first.family)?;
    let artifact = snapshot.artifact();
    let dimension = S.epsilon_dimension.clone();
    let rustred_dimension = S.rustred_dimension.clone();
    let epsilon = Atom::var(vk_symbol!(settings.epsilon_symbol.as_str()));
    let physical_dimension = Atom::num(4) - Atom::num(2) * &epsilon;
    let mut channels =
        BTreeMap::<(Atom, Atom, i64), Vec<(IntegralKey, rustred::algebra::Coefficient)>>::new();
    let mut anchors = BTreeMap::<Atom, i64>::new();
    // Contributions per mass, summed once at the end (see four_loop).
    let mut raw = BTreeMap::<Atom, Vec<Atom>>::new();
    drop(setup);
    for contribution in contributions {
        let _stage = StageSpan::new("rustred.lower.admission_and_channels");
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
        // The Atom-only lower-loop API has no public guard field. It still
        // rejects an original pole that vanishes under the dimension relation.
        {
            let _stage = StageSpan::new("rustred.lower.input_conditions");
            drop(input.input_conditions_after_mapping(&scalar, |condition| {
                condition
                    .replace(Atom::var(S.form_epsilon).to_pattern())
                    .with(epsilon.to_pattern())
                    .replace(rustred_dimension.to_pattern())
                    .with(physical_dimension.to_pattern())
            })?);
        }
        // MATAD finalization reserves this exact internal regulator symbol.
        let dimensions = StageSpan::new("rustred.lower.dimension_mapping");
        // The intermediate numerator is dropped before the old one is
        // replaced, so only the mapped numerator is alive in the lowering.
        input.numerator = {
            let numerator = input
                .numerator
                .replace(Atom::var(S.form_epsilon).to_pattern())
                .with(epsilon.to_pattern());
            super::numerator::MatchedScalarInput::numerator_in_family_dimension(
                numerator.as_view(),
                settings,
                &rustred_dimension,
            )
        };
        drop(dimensions);
        let lowering = {
            let _stage = StageSpan::new("rustred.lower.scalar_lowering");
            let lower =
                || scalar.lower_weighted(&input.numerator, &input.target, &input.mass_squared);
            #[cfg(feature = "weighted-lowering-diagnostics")]
            let result = {
                use rustred::scalar_numerator::lowering_diagnostics;
                let capture_enabled = tracing::enabled!(
                    target: "vakint::coefficient_lane_inputs", tracing::Level::INFO
                );
                if (tracing::enabled!(target: "vakint::coefficient_lanes", tracing::Level::INFO)
                    || capture_enabled)
                    && lowering_diagnostics::current().is_none()
                {
                    let (result, counts, opaque, rejection) =
                        lowering_diagnostics::scoped_with_opaque_rejection(false, lower);
                    tracing::info!(target: "vakint::coefficient_lanes",
                        rational = counts.rational, grouped = counts.grouped,
                        exact = counts.exact, opaque, generic = counts.generic,
                        opaque_reject_reason = ?rejection.map(|r| r.reason),
                        opaque_reject_stage = ?rejection.map(|r| r.stage),
                        opaque_reject_observed = rejection.and_then(|r| r.observed),
                        opaque_reject_limit = rejection.and_then(|r| r.limit),
                        opaque_reject_support_n = rejection.and_then(|r| r.shape).map(|s| s.numerator),
                        opaque_reject_support_d = rejection.and_then(|r| r.shape).map(|s| s.denominator),
                        opaque_reject_degree_n = rejection.and_then(|r| r.shape).map(|s| s.numerator_degree),
                        opaque_reject_degree_d = rejection.and_then(|r| r.shape).map(|s| s.denominator_degree),
                        opaque_reject_den_relation = ?rejection.and_then(|r| r.denominators).map(|d| d.relation),
                        opaque_reject_left_factors = rejection.and_then(|r| r.denominators).map(|d| d.left_factors),
                        opaque_reject_right_factors = rejection.and_then(|r| r.denominators).map(|d| d.right_factors),
                        opaque_reject_shared_factors = rejection.and_then(|r| r.denominators).map(|d| d.shared_factors),
                        "scalar numerator coefficient arithmetic");
                    if capture_enabled
                        && counts.generic != 0
                        && CAPTURED_GENERIC_SCALAR_INPUTS
                            .fetch_add(1, std::sync::atomic::Ordering::Relaxed)
                            < 8
                    {
                        tracing::info!(target: "vakint::coefficient_lane_inputs", {
                            file.numerator = %input.numerator.to_canonical_string(),
                            file.mass = %input.mass_squared.to_canonical_string(),
                            },
                            "Captured generic scalar numerator input");
                    }
                    result
                } else {
                    lower()
                }
            };
            #[cfg(not(feature = "weighted-lowering-diagnostics"))]
            let result = lower();
            result.map_err(RustRedEvaluationError::from_native)?
        };
        let _channels = StageSpan::new("rustred.lower.channel_assembly");
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
                raw.entry(matched.mass_squared.clone()).or_default().push(
                    term.coefficient()
                        .to_expression()
                        .replace(rustred_dimension.to_pattern())
                        .with(dimension.to_pattern())
                        * term.scalar_spectator()
                        * matched.mass_squared.clone().pow(exponent)
                        * data
                            .masters
                            .symbol(artifact.family_fingerprint(), term.integral().powers()),
                );
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
            for (spectator, coefficient) in
                super::numerator::MatchedScalarInput::channel_coefficients(
                    spectator,
                    term.coefficient(),
                    artifact.coefficient_context(),
                    settings
                        .rustred_application_options
                        .reduction_limits
                        .exact_algebra,
                )?
            {
                channels
                    .entry((matched.mass_squared.clone(), spectator, exponent))
                    .or_default()
                    .push((term.integral().clone(), coefficient));
            }
        }
    }
    let step_cache = {
        let _stage = StageSpan::new("rustred.lower.step_cache_setup");
        if channels.is_empty() {
            None
        } else {
            vakint.ordered_step_cache(
                settings,
                &snapshot.ordered_cache_program(
                    first.loop_count,
                    first.family.name(),
                    &data.version,
                ),
            )?
        }
    };
    let cohort_total = channels.len();
    let mut skipped = 0;
    let arithmetic = settings.rustred_ordered_weight_arithmetic();
    for (cohort_index, ((mass, spectator, exponent), terms)) in channels.into_iter().enumerate() {
        let setup = StageSpan::new("rustred.lower.cohort_setup");
        if let Some(progress) = progress {
            progress.cohorts(cohort_index, cohort_total);
        }
        let mut dump = RequestDump::begin(DumpRequest {
            loop_count: first.loop_count,
            program: DumpProgram::ClosedArtifact {
                family_fingerprint: artifact.family_fingerprint(),
                artifact_sha256: snapshot.source_digest(),
                data: &data.loops[&first.loop_count],
            },
            weight_arithmetic: arithmetic,
            master_mode: request_dump::master_mode(settings),
            channel: DumpChannel {
                index: cohort_index,
                total: cohort_total,
                mass_squared: &mass,
                spectator: &spectator,
                exponent,
            },
            limits: settings.rustred_application_options.reduction_limits,
            coefficient_variables: artifact.coefficient_context().parameter_names().to_vec(),
            inputs: &terms,
            caller_conditions: Vec::new(),
        })?;
        if dump.as_ref().is_some_and(RequestDump::only) {
            skipped += 1;
            continue;
        }
        drop(setup);
        let _math = MathSpan::new();
        let reduction = request_dump::observe_native(
            dump.as_mut(),
            observation.as_deref_mut(),
            |observation| match observation {
                Some(observer) => {
                    observer.clear_native_summary();
                    let mut envelope = rustred::solver::CandidateRoutedApplicationEnvelope::new(
                        rustred::solver::CandidateRoutedApplicationEvent::ActivityStarted,
                        artifact.family_fingerprint(),
                    );
                    envelope.caller_phase = Some("entry_admission");
                    observer.emit(envelope);
                    let _stage = StageSpan::new("rustred.lower.native_ibp");
                    Reducer::reduce_artifact_linear_combination_ordered_with_step_cache(
                        artifact,
                        terms,
                        settings.rustred_application_options.reduction_limits,
                        arithmetic,
                        step_cache.as_deref(),
                        &mut |progress, current| {
                            observer.emit_ordered_progress(
                                artifact.family_fingerprint(),
                                *progress,
                                current,
                            );
                            Ok(())
                        },
                        None,
                    )
                }
                None => {
                    let _stage = StageSpan::new("rustred.lower.native_ibp");
                    Reducer::reduce_artifact_linear_combination_ordered_with_step_cache(
                        artifact,
                        terms,
                        settings.rustred_application_options.reduction_limits,
                        arithmetic,
                        step_cache.as_deref(),
                        &mut |_, _| Ok(()),
                        None,
                    )
                }
            },
        )?
        .map_err(RustRedEvaluationError::from_native)?;
        if let Some(dump) = dump {
            let _stage = StageSpan::new("rustred.lower.request_dump");
            dump.finish(
                request_dump::family_terms(artifact.family_fingerprint(), &reduction.terms),
                &reduction.nonzero_conditions,
                Some(reduction.progress),
            )?;
        }
        let masters = {
            let _stage = StageSpan::new("rustred.lower.materialize_coefficients");
            reduction
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
                .sum::<Atom>()
        };
        let _accumulation = StageSpan::new("rustred.lower.accumulate_channel");
        let spectator = spectator
            .replace(rustred_dimension.to_pattern())
            .with(dimension.to_pattern());
        raw.entry(mass.clone())
            .or_default()
            .push(masters * spectator * mass.pow(exponent));
    }
    if let Some(progress) = progress {
        progress.cohorts(cohort_total, cohort_total);
    }
    if skipped != 0 {
        return Err(request_dump::skipped(skipped));
    }
    let mut value = Atom::Zero;
    for (mass, raw) in raw {
        let raw = {
            let _stage = StageSpan::new("rustred.lower.combine_before_normalization");
            let _math = MathSpan::new();
            super::combine_exact(Atom::add_many(raw))
        };
        let normalized = {
            let _stage = StageSpan::new("rustred.lower.matad_normalization");
            MATAD::with_settings(settings.clone())
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
                .map_err(super::data::failure)?
        };
        {
            let _stage = StageSpan::new("rustred.lower.accumulate_mass");
            let _math = MathSpan::new();
            value += normalized;
        }
    }
    Ok(value)
}

#[cfg(test)]
#[path = "materialize_tests.rs"]
mod tests;
