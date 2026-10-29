//! Select already shipped lower-loop programs by a verified momentum map.
//! RustRed owns guarded application, sum-of-products algebra and recursion.

use std::{collections::BTreeMap, sync::Arc};

use feynkit_graph::{IntegralFamily as MomentumFamily, IntegralMapping};
use feynkit_kinematics::Kinematics;
use rustred::family::{IntegralFamily, IntegralFamilyLimits, IntegralKey, ScalarProductCoordinate};
use rustred::reduction::terminal_normalization::{
    VacuumFactorProvider, VacuumFactorReduction, VacuumFactorTerm, VacuumFactorWork,
    VacuumProductReductionError, reduce_transported_vacuum_factor,
};
use rustred::reduction::{Reducer, ReductionLimits, ReductionStatistics};
use rustred::sector::{Mask, symmetry};
use symbolica::atom::{Atom, AtomCore};
use symbolica::prelude::{Q, Z};
use symbolica::utils::resources::{self, ResourceLimitError, ResourceLimits};

use super::{artifact::shipped_family_for_loop_count, data::RustRedData};
use crate::symbols::S;

/// Only immutable family preparation is retained here. Guarded native actions
/// use fresh request state; the enclosing product-reduction session owns reuse.
#[derive(Default)]
pub(super) struct ShippedFactorProvider {
    data: Option<Arc<RustRedData>>,
    closed_families: BTreeMap<usize, Arc<IntegralFamily>>,
}

impl ShippedFactorProvider {
    pub(super) fn with_data(data: Arc<RustRedData>) -> Self {
        Self {
            data: Some(data),
            closed_families: BTreeMap::new(),
        }
    }
}

fn error(value: impl std::error::Error + 'static) -> FactorProviderError {
    if let Some(resource) = FactorProviderError::find_resource(&value) {
        FactorProviderError::Resource(resource)
    } else {
        FactorProviderError::Invalid(value.to_string())
    }
}

fn share_family(family: &IntegralFamily) -> Result<Arc<IntegralFamily>, FactorProviderError> {
    resources::checkpoint()?;
    let shared = IntegralFamily::new_with_limits(
        family.name(),
        family.loop_momenta().to_vec(),
        family.external_momenta().to_vec(),
        family.coefficient_context().clone(),
        family.dimension().clone(),
        family.denominators().to_vec(),
        family.external_gram().to_vec(),
        family.power_shifts().to_vec(),
        IntegralFamilyLimits::default(),
    )
    .map_err(error)?;
    if shared.fingerprint() != family.fingerprint() {
        return Err("rebuilt lower program family has a different authenticated identity".into());
    }
    Ok(Arc::new(shared))
}

fn kinematics(loops: usize) -> Result<(Vec<Atom>, Kinematics), FactorProviderError> {
    let momenta = (1..=loops).map(|axis| S.k.call(axis)).collect::<Vec<_>>();
    let kinematics = Kinematics::new()
        .with_momenta(momenta.iter().cloned())
        .map_err(error)?;
    Ok((momenta, kinematics))
}

/// Translate typed affine scalar-product rows, never coefficient text.
fn momentum_family(
    family: &IntegralFamily,
    slots: &[usize],
) -> Result<MomentumFamily, FactorProviderError> {
    resources::checkpoint()?;
    if family.external_count() != 0 {
        return Err("lower factor is not a vacuum family".into());
    }
    let (momenta, kinematics) = kinematics(family.loop_count())?;
    let products = family
        .coordinates()
        .iter()
        .map(|coordinate| match coordinate {
            ScalarProductCoordinate::LoopLoop { left, right } => kinematics
                .scalar_product(&momenta[*left], &momenta[*right])
                .map_err(error),
            _ => Err("lower vacuum factor has external scalar-product coordinates".into()),
        })
        .collect::<Result<Vec<_>, FactorProviderError>>()?;
    let denominators = slots
        .iter()
        .map(|slot| {
            let row = &family.denominators()[*slot];
            row.constant().to_expression()
                + row
                    .coefficients()
                    .iter()
                    .zip(&products)
                    .map(|(coefficient, product)| coefficient.to_expression() * product)
                    .sum::<Atom>()
        })
        .collect();
    MomentumFamily::new(momenta, Vec::new(), denominators, &kinematics).map_err(error)
}

fn physical_family(momenta: &[Atom], loops: usize) -> Result<MomentumFamily, FactorProviderError> {
    resources::checkpoint()?;
    let (basis, kinematics) = kinematics(loops)?;
    let denominators = momenta
        .iter()
        .map(|momentum| {
            kinematics
                .scalar_product(momentum, momentum)
                .map(|square| square - Atom::one())
                .map_err(error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    MomentumFamily::new(basis, Vec::new(), denominators, &kinematics).map_err(error)
}

fn mapped_key(
    source: &MomentumFamily,
    target: &MomentumFamily,
    source_slots: &[usize],
    integral: &IntegralKey,
    arity: usize,
) -> Result<Option<(IntegralKey, IntegralMapping)>, FactorProviderError> {
    // FeynKit's native search accepts only an exact denominator map
    // with real linear loop images and determinant +/-1. No U-only equality.
    // Its diagnostic count limit is absent; the surrounding request retains
    // its original RAM/time budget across this nested native operation.
    let Some(mapping) = source
        .find_mapping(target, None, ResourceLimits::default())
        .map_err(error)?
    else {
        return Ok(None);
    };
    let mut powers = vec![0; arity];
    for (source_position, target_position) in mapping.denominator_map().iter().enumerate() {
        powers[*target_position] = integral.powers()[source_slots[source_position]];
    }
    Ok(Some((
        IntegralKey::try_new(powers).map_err(error)?,
        mapping,
    )))
}

/// Replay the already verified physical witness over every source row, including
/// numerator auxiliaries. RustRed owns its exact affine expansion and guards.
fn signed_transport(
    source: &IntegralFamily,
    destination: Arc<IntegralFamily>,
    source_key: &IntegralKey,
    mapped_key: &IntegralKey,
    mapping: &IntegralMapping,
    limits: ReductionLimits,
) -> Result<symmetry::integral_transport::Prepared, FactorProviderError> {
    let loops = source.loop_count();
    let context = destination.coefficient_context();
    let momenta = (1..=loops).map(|axis| S.k.call(axis)).collect::<Vec<_>>();
    let variables = momenta.iter().collect::<Vec<_>>();
    let mut matrix = Vec::with_capacity(loops * loops);
    for (axis, (source_label, image)) in mapping.momentum_rules().iter().enumerate() {
        if source_label != &momenta[axis] {
            return Err("physical momentum witness has an unexpected source order".into());
        }
        let mut row = vec![context.zero(); loops];
        for (monomial, coefficient) in image.coefficient_list::<u16>(&variables) {
            let target_axis = momenta
                .iter()
                .position(|momentum| momentum == &monomial)
                .ok_or("physical momentum witness is not homogeneous linear syntax")?;
            let coefficient = coefficient
                .try_to_rational_polynomial(&Q, &Z, Some(context.one().get_variables().clone()))
                .map_err(error)?;
            context
                .validate_with_limits(&coefficient, limits.exact_algebra)
                .map_err(error)?;
            row[target_axis] = coefficient;
        }
        matrix.extend(row);
    }
    let momentum_map = symmetry::MomentumMap::new(
        symmetry::CoefficientMatrix::try_new(loops, loops, matrix).map_err(error)?,
        symmetry::CoefficientMatrix::try_new(loops, 0, []).map_err(error)?,
        symmetry::CoefficientMatrix::try_new(0, 0, []).map_err(error)?,
    );
    let verified =
        symmetry::verify(source, &destination, momentum_map, Default::default()).map_err(error)?;
    symmetry::integral_transport::compile(
        source,
        destination,
        Arc::new(verified),
        Mask::try_new(source_key.powers().iter().map(|power| *power > 0)).map_err(error)?,
        Mask::try_new(mapped_key.powers().iter().map(|power| *power > 0)).map_err(error)?,
        Default::default(),
    )
    .map_err(error)
}

fn result(
    family: Arc<IntegralFamily>,
    terms: &BTreeMap<IntegralKey, rustred::algebra::Coefficient>,
    nonzero_conditions: Vec<rustred::algebra::Coefficient>,
    statistics: ReductionStatistics,
) -> VacuumFactorReduction {
    VacuumFactorReduction {
        terms: terms
            .iter()
            .map(|(integral, coefficient)| VacuumFactorTerm {
                family: Arc::clone(&family),
                integral: integral.clone(),
                coefficient: coefficient.clone(),
            })
            .collect(),
        nonzero_conditions,
        work: VacuumFactorWork {
            rule_applications: statistics.rule_applications(),
            coalescing_additions: statistics.coalescing_additions(),
            cache_hits: statistics.cache_hits(),
        },
    }
}

#[derive(Debug)]
pub(super) enum FactorProviderError {
    Invalid(String),
    Native(VacuumProductReductionError),
    Resource(ResourceLimitError),
}
impl FactorProviderError {
    fn find_resource(error: &(dyn std::error::Error + 'static)) -> Option<ResourceLimitError> {
        super::RustRedEvaluationError::resource_cause(error).cloned()
    }
}
impl std::fmt::Display for FactorProviderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(error) => f.write_str(error),
            Self::Native(error) => error.fmt(f),
            Self::Resource(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for FactorProviderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Invalid(_) => None,
            Self::Native(error) => Some(error),
            Self::Resource(error) => Some(error),
        }
    }
}
impl From<ResourceLimitError> for FactorProviderError {
    fn from(error: ResourceLimitError) -> Self {
        Self::Resource(error)
    }
}
impl From<String> for FactorProviderError {
    fn from(error: String) -> Self {
        Self::Invalid(error)
    }
}
impl From<&str> for FactorProviderError {
    fn from(error: &str) -> Self {
        Self::Invalid(error.into())
    }
}
impl From<VacuumProductReductionError> for FactorProviderError {
    fn from(error: VacuumProductReductionError) -> Self {
        Self::Native(error)
    }
}

impl VacuumFactorProvider for ShippedFactorProvider {
    type Error = FactorProviderError;
    fn resource_limit_error(error: &Self::Error) -> Option<ResourceLimitError> {
        FactorProviderError::find_resource(error)
    }
    fn error_is_missing_coverage(error: &Self::Error) -> bool {
        matches!(
            error,
            FactorProviderError::Native(VacuumProductReductionError::UnresolvedFactor { .. })
        )
    }

    fn reduce(
        &mut self,
        family: Arc<IntegralFamily>,
        target: &IntegralKey,
        limits: ReductionLimits,
    ) -> Result<VacuumFactorReduction, FactorProviderError> {
        if target.powers().len() != family.denominator_count() {
            return Err("the lower-factor provider received a foreign-arity integral".into());
        }
        let has_numerator = target.powers().iter().any(|power| *power < 0);
        let slots = target
            .powers()
            .iter()
            .enumerate()
            .filter_map(|(slot, power)| (*power > 0).then_some(slot))
            .collect::<Vec<_>>();
        let source = momentum_family(&family, &slots)?;
        let data = self
            .data
            .clone()
            .map(Ok)
            .unwrap_or_else(RustRedData::vendored)
            .map_err(error)?;
        if let Some(shipped) = shipped_family_for_loop_count(family.loop_count()) {
            let snapshot = data.artifact(shipped).map_err(error)?;
            let artifact = snapshot.artifact();
            let shared = match self.closed_families.get(&family.loop_count()) {
                Some(shared) => Arc::clone(shared),
                None => {
                    let shared = share_family(artifact.family())?;
                    self.closed_families
                        .insert(family.loop_count(), Arc::clone(&shared));
                    shared
                }
            };
            let destination = momentum_family(
                &shared,
                &(0..shared.denominator_count()).collect::<Vec<_>>(),
            )?;
            let (mapped, mapping) = mapped_key(
                &source,
                &destination,
                &slots,
                target,
                shared.denominator_count(),
            )?
            .ok_or_else(|| "no exact momentum map to the shipped lower-loop program".to_owned())?;
            match snapshot.program().validate_entry(&mapped) {
                Ok(()) => {}
                Err(
                    rustred::reduction::ReductionError::OutsideCertifiedRootDomain { .. }
                    | rustred::reduction::ReductionError::OutsideCertifiedTotalExcessDomain {
                        ..
                    },
                ) => {
                    return Err(FactorProviderError::Native(
                        VacuumProductReductionError::UnresolvedFactor {
                            family: family.fingerprint().into(),
                            target: target.clone(),
                        },
                    ));
                }
                Err(refusal) => return Err(error(refusal)),
            }
            if has_numerator {
                let transport = signed_transport(
                    &family,
                    Arc::clone(&shared),
                    target,
                    &mapped,
                    &mapping,
                    limits,
                )?;
                let mut outside_endpoint = None;
                let mut callback_resource = None;
                let reduction = reduce_transported_vacuum_factor(
                    &transport,
                    target,
                    Arc::clone(&shared),
                    Default::default(),
                    limits,
                    |key, remaining| {
                        // Numerator transport can leave the certified entry
                        // domain even when its positive denominator map fits.
                        // The generic native callback has a Display boundary.
                        // Retain its typed operational cause before it becomes
                        // a diagnostic string, independently of entry coverage.
                        if let Err(refusal) = snapshot.program().validate_entry(key) {
                            callback_resource = FactorProviderError::find_resource(&refusal);
                            if matches!(
                                refusal,
                                rustred::reduction::ReductionError::OutsideCertifiedRootDomain { .. }
                                    | rustred::reduction::ReductionError::OutsideCertifiedTotalExcessDomain { .. }
                            ) {
                                outside_endpoint = Some(key.clone());
                            }
                            return Err(error(refusal));
                        }
                        let (reduced, conditions, work) =
                            match Reducer::reduce_artifact_unit_mass_guarded(
                                artifact, key, remaining,
                            ) {
                                Ok(reduced) => reduced,
                                Err(refusal) => {
                                    callback_resource =
                                        FactorProviderError::find_resource(&refusal);
                                    return Err(error(refusal));
                                }
                            };
                        Ok::<_, FactorProviderError>(result(
                            Arc::clone(&shared),
                            reduced.terms(),
                            conditions,
                            work,
                        ))
                    },
                );
                if let Some(resource) = callback_resource {
                    return Err(FactorProviderError::Resource(resource));
                }
                return match (reduction, outside_endpoint) {
                    (Err(VacuumProductReductionError::Provider { .. }), Some(target)) => {
                        Err(FactorProviderError::Native(
                            VacuumProductReductionError::UnresolvedFactor {
                                family: shared.fingerprint().into(),
                                target,
                            },
                        ))
                    }
                    (result, _) => result.map_err(FactorProviderError::from),
                };
            }
            let (reduced, conditions, work) =
                Reducer::reduce_artifact_unit_mass_guarded(artifact, &mapped, limits)
                    .map_err(error)?;
            return Ok(result(shared, reduced.terms(), conditions, work));
        }
        let Some(programs) = data.four_loop.as_ref() else {
            return Err(FactorProviderError::Native(
                VacuumProductReductionError::UnresolvedFactor {
                    family: family.fingerprint().into(),
                    target: target.clone(),
                },
            ));
        };
        let mut index = 0;
        while let Some(momenta) = programs.physical_momenta(index) {
            if family.loop_count() != 4 {
                break;
            }
            let destination = physical_family(momenta, family.loop_count())?;
            if let Some((mapped, mapping)) = mapped_key(&source, &destination, &slots, target, 10)?
            {
                let (native, _) = programs
                    .native_program(index)
                    .ok_or("missing shipped candidate program")?;
                if has_numerator {
                    let shared = native.family_shared();
                    let transport = signed_transport(
                        &family,
                        Arc::clone(&shared),
                        target,
                        &mapped,
                        &mapping,
                        limits,
                    )?;
                    return native
                        .reduce_transported_unit_mass_guarded(&transport, target, limits)
                        .map_err(FactorProviderError::from);
                }
                let (reduced, conditions, work) =
                    native.reduce_unit_mass_guarded(&mapped, limits)?;
                return Ok(result(
                    native.family_shared(),
                    reduced.terms(),
                    conditions,
                    work,
                ));
            }
            index += 1;
        }
        Err(FactorProviderError::Native(
            VacuumProductReductionError::UnresolvedFactor {
                family: family.fingerprint().into(),
                target: target.clone(),
            },
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lower_provider_resource_errors_are_not_coverage_and_allow_healthy_reuse() {
        use symbolica::utils::resources::ResourceBudget;

        crate::Vakint::initialize_vakint_symbols();
        let artifact = shipped_family_for_loop_count(1)
            .unwrap()
            .artifact()
            .unwrap();
        let family = share_family(artifact.family()).unwrap();
        let target = IntegralKey::try_new([3]).unwrap();
        let mut provider = ShippedFactorProvider::default();
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1),
                ..Default::default()
            },
            || Ok(2),
        );
        let refusal = budget
            .run(|| provider.reduce(Arc::clone(&family), &target, Default::default()))
            .unwrap_err();
        assert!(matches!(
            ShippedFactorProvider::resource_limit_error(&refusal),
            Some(ResourceLimitError::Memory {
                rss_bytes: 2,
                limit_bytes: 1,
                ..
            })
        ));
        assert!(!ShippedFactorProvider::error_is_missing_coverage(&refusal));
        // Inspect the original error after leaving the constrained scope; a
        // fresh global poll would succeed here and lose the actual cause.
        assert!(resources::checkpoint().is_ok());
        assert!(FactorProviderError::find_resource(&refusal).is_some());
        let recovered = provider
            .reduce(family, &target, Default::default())
            .unwrap();
        let direct =
            Reducer::reduce_artifact_unit_mass_guarded(artifact, &target, Default::default())
                .unwrap();
        assert_eq!(
            recovered
                .terms
                .iter()
                .map(|term| (term.integral.clone(), term.coefficient.clone()))
                .collect::<BTreeMap<_, _>>(),
            *direct.0.terms(),
        );
        assert_eq!(recovered.nonzero_conditions, direct.1);
    }

    #[test]
    fn lower_native_provider_preserves_exact_terms_and_guards() {
        crate::Vakint::initialize_vakint_symbols();
        let mut provider = ShippedFactorProvider::default();
        for (loops, powers) in [(1, vec![3]), (2, vec![2, 1, 1])] {
            let artifact = shipped_family_for_loop_count(loops)
                .unwrap()
                .artifact()
                .unwrap();
            let family = share_family(artifact.family()).unwrap();
            let target = IntegralKey::try_new(powers).unwrap();
            let direct = Reducer::reduce_artifact_unit_mass_guarded(
                artifact,
                &target,
                ReductionLimits::default(),
            )
            .unwrap();
            let reduced = provider
                .reduce(Arc::clone(&family), &target, ReductionLimits::default())
                .unwrap();
            let terms = reduced
                .terms
                .iter()
                .map(|term| {
                    assert_eq!(term.family.fingerprint(), family.fingerprint());
                    (term.integral.clone(), term.coefficient.clone())
                })
                .collect::<BTreeMap<_, _>>();
            assert_eq!(&terms, direct.0.terms());
            assert_eq!(reduced.nonzero_conditions, direct.1);
            let warm = provider
                .reduce(family, &target, ReductionLimits::default())
                .unwrap();
            assert_eq!(
                warm.work.cache_hits, reduced.work.cache_hits,
                "the provider itself must not retain an unaccounted native point cache"
            );
        }
    }

    #[test]
    fn lower_native_provider_honors_scoped_limits() {
        crate::Vakint::initialize_vakint_symbols();
        let artifact = shipped_family_for_loop_count(1)
            .unwrap()
            .artifact()
            .unwrap();
        let family = share_family(artifact.family()).unwrap();
        let target = IntegralKey::try_new([3]).unwrap();
        let mut provider = ShippedFactorProvider::default();
        assert!(
            provider
                .reduce(
                    family,
                    &target,
                    ReductionLimits {
                        max_rule_applications: 0,
                        ..Default::default()
                    }
                    .with_diagnostics()
                )
                .is_err()
        );
    }
}

#[cfg(test)]
#[path = "factor_provider_tests.rs"]
mod independent_tests;
