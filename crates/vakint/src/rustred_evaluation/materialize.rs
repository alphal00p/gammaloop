use rustred::algebra::ExactAlgebraLimits;
use rustred::family::IntegralKey;
use rustred::reduction::Reducer;
use rustred::scalar_numerator::{ScalarNumeratorLimits, ScalarNumeratorService};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use crate::matad::MATAD;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{MATADOptions, Vakint, VakintSettings};

use super::artifact::ArtifactFamily;
use super::matching::MatchedScalarFamily;
use super::{RustRedEvaluationError, RustRedEvaluationOptions};

pub(super) fn evaluate(
    matched: MatchedScalarFamily,
    settings: &VakintSettings,
    numerator: AtomView,
    options: &RustRedEvaluationOptions,
) -> Result<Atom, RustRedEvaluationError> {
    let dimension = vk_parse!("4-2*ep").expect("Vakint's internal dimension expression parses");
    let rustred_dimension = vk_parse!("rustred::d")
        .expect("RustRed's stable coefficient namespace parses in Symbolica");
    let artifact = matched.family.artifact()?;
    let base_integral = IntegralKey::try_new(matched.powers.iter().copied()).map_err(|error| {
        RustRedEvaluationError::IntegralKey {
            detail: error.to_string(),
        }
    })?;
    let loop_momenta = (1..=matched.loop_count)
        .map(|index| function!(S.k, Atom::num(index as i64)))
        .collect();
    let internal_numerator = Vakint::convert_to_dot_notation(numerator)
        .replace(Atom::var(vk_symbol!(settings.epsilon_symbol.as_str())).to_pattern())
        .with(vk_parse!("ep").unwrap().to_pattern());
    let scalar_service = ScalarNumeratorService::try_new(
        artifact,
        S.dot,
        loop_momenta,
        ScalarNumeratorLimits::default(),
    )
    .map_err(|error| RustRedEvaluationError::ScalarNumerator {
        detail: error.to_string(),
    })?;
    let lowering = scalar_service
        .lower(&internal_numerator, &base_integral)
        .map_err(|error| RustRedEvaluationError::ScalarNumerator {
            detail: error.to_string(),
        })?;
    let mut reducer =
        Reducer::new(artifact).map_err(|error| RustRedEvaluationError::Reduction {
            detail: error.to_string(),
        })?;
    let mut raw_masters = Atom::Zero;
    for lowered in lowering.terms() {
        let decomposition = reducer
            .reduce_with_common_mass_homogeneity(lowered.integral())
            .map_err(|error| RustRedEvaluationError::Reduction {
                detail: error.to_string(),
            })?;
        for (master, coefficient) in decomposition.terms() {
            let combined = artifact
                .coefficient_context()
                .try_mul(
                    lowered.coefficient(),
                    coefficient.unit_mass_coefficient(),
                    ExactAlgebraLimits::default(),
                )
                .map_err(|error| RustRedEvaluationError::Reduction {
                    detail: error.to_string(),
                })?;
            let exact_coefficient = combined
                .to_expression()
                .replace(rustred_dimension.to_pattern())
                .with(dimension.to_pattern());
            let exponent = coefficient
                .common_mass_squared_power()
                .checked_add(i128::from(lowered.common_mass_squared_power()))
                .ok_or(RustRedEvaluationError::MassExponentAdditionOverflow {
                    reduction: coefficient.common_mass_squared_power(),
                    numerator: lowered.common_mass_squared_power(),
                })?;
            let exponent = i64::try_from(exponent)
                .map_err(|_| RustRedEvaluationError::MassExponentOverflow { exponent })?;
            let scale = matched.mass_squared.clone().pow(Atom::num(exponent));
            let terminal = terminal_master(matched.family, master, &matched.mass_squared)?;
            raw_masters +=
                exact_coefficient * scale * terminal * lowered.scalar_spectator().clone();
        }
    }

    let loop_count =
        i64::try_from(matched.loop_count).map_err(|_| RustRedEvaluationError::Reduction {
            detail: "matched loop count does not fit Vakint's normalization exponent".to_owned(),
        })?;
    let matad_options = MATADOptions {
        expand_masters: options.substitute_masters,
        susbstitute_masters: options.substitute_masters,
        substitute_hpls: options.substitute_masters,
        direct_numerical_substition: options.substitute_masters,
    };
    MATAD::with_settings(settings.clone())
        .finalize_reduced_masters(
            raw_masters,
            loop_count,
            None,
            &matched.mass_squared,
            &matad_options,
        )
        .map_err(|error| RustRedEvaluationError::Reduction {
            detail: format!("Vakint master materialization failed: {error}"),
        })
}

fn terminal_master(
    family: ArtifactFamily,
    master: &IntegralKey,
    mass_squared: &Atom,
) -> Result<Atom, RustRedEvaluationError> {
    // MATAD's exact output canonicalizes Gamma(-1+ep) to Gamma(1+ep).
    // Retaining that canonical basis makes raw backend comparisons purely
    // rational, without teaching the adapter a second Gamma simplifier.
    let tadpole = mass_squared.clone()
        * vk_parse!("Gam(1,1)/(ep*(ep-1))")
            .expect("Vakint's canonical MATAD tadpole master parses");
    match (family, master.powers()) {
        (ArtifactFamily::UnitMassVacuumK1, [1]) => Ok(-tadpole),
        (ArtifactFamily::UnitMassVacuumK3, [1, 1, 1]) => Ok(-mass_squared.clone()
            * vk_parse!("miT111").expect("Vakint's MATAD sunset master parses")),
        (ArtifactFamily::UnitMassVacuumK3, [0, 1, 1]) => Ok(tadpole.pow(Atom::num(2))),
        _ => Err(RustRedEvaluationError::UnsupportedMaster {
            family: family.name(),
            powers: master.powers().to_vec(),
        }),
    }
}
