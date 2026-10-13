//! Select already shipped lower-loop programs by a verified momentum map.
//! RustRed owns guarded application, sum-of-products algebra and recursion.

use std::{collections::BTreeMap, sync::Arc};

use feynkit_graph::{IntegralFamily as MomentumFamily, IntegralMapping};
use feynkit_kinematics::Kinematics;
use rustred::family::{IntegralFamily, IntegralFamilyLimits, IntegralKey, ScalarProductCoordinate};
use rustred::reduction::terminal_normalization::{
    VacuumFactorProvider, VacuumFactorReduction, VacuumFactorTerm, VacuumFactorWork,
    reduce_transported_vacuum_factor,
};
use rustred::reduction::{Reducer, ReductionLimits, ReductionStatistics};
use rustred::sector::{Mask, symmetry};
use symbolica::atom::{Atom, AtomCore};
use symbolica::prelude::{Q, Z};

use super::{artifact::shipped_family_for_loop_count, four_loop};
use crate::symbols::S;

/// Only immutable family preparation is retained here. Scoped native caches are
/// cleared at both call boundaries; the enclosing RustRed session owns reuse.
#[derive(Default)]
pub(super) struct ShippedFactorProvider {
    closed_families: BTreeMap<usize, Arc<IntegralFamily>>,
}

fn error(value: impl std::fmt::Display) -> String {
    value.to_string()
}

fn share_family(family: &IntegralFamily) -> Result<Arc<IntegralFamily>, String> {
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

fn kinematics(loops: usize) -> Result<(Vec<Atom>, Kinematics), String> {
    let momenta = (1..=loops).map(|axis| S.k.call(axis)).collect::<Vec<_>>();
    let kinematics = Kinematics::new()
        .with_momenta(momenta.iter().cloned())
        .map_err(error)?;
    Ok((momenta, kinematics))
}

/// Translate typed affine scalar-product rows, never coefficient text.
fn momentum_family(family: &IntegralFamily, slots: &[usize]) -> Result<MomentumFamily, String> {
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
        .collect::<Result<Vec<_>, String>>()?;
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

fn physical_family(momenta: &[Atom], loops: usize) -> Result<MomentumFamily, String> {
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
) -> Result<Option<(IntegralKey, IntegralMapping)>, String> {
    // FeynKit's existing bounded search accepts only an exact denominator map
    // with real linear loop images and determinant +/-1. No U-only equality.
    let Some(mapping) = source.find_mapping(target, 1_000_000).map_err(error)? else {
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
) -> Result<symmetry::integral_transport::Prepared, String> {
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

impl VacuumFactorProvider for ShippedFactorProvider {
    type Error = String;

    fn reduce(
        &mut self,
        family: Arc<IntegralFamily>,
        target: &IntegralKey,
        limits: ReductionLimits,
    ) -> Result<VacuumFactorReduction, String> {
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
        if let Some(shipped) = shipped_family_for_loop_count(family.loop_count()) {
            let artifact = shipped.artifact().map_err(error)?;
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
            if has_numerator {
                let transport = signed_transport(
                    &family,
                    Arc::clone(&shared),
                    target,
                    &mapped,
                    &mapping,
                    limits,
                )?;
                return reduce_transported_vacuum_factor(
                    &transport,
                    target,
                    Arc::clone(&shared),
                    Default::default(),
                    limits,
                    |key, remaining| {
                        let (reduced, conditions, work) =
                            Reducer::reduce_artifact_unit_mass_guarded(artifact, key, remaining)
                                .map_err(error)?;
                        Ok::<_, String>(result(
                            Arc::clone(&shared),
                            reduced.terms(),
                            conditions,
                            work,
                        ))
                    },
                )
                .map_err(error);
            }
            let (reduced, conditions, work) =
                Reducer::reduce_artifact_unit_mass_guarded(artifact, &mapped, limits)
                    .map_err(error)?;
            return Ok(result(shared, reduced.terms(), conditions, work));
        }
        let mut index = 0;
        while let Some(momenta) = four_loop::physical_momenta(index)? {
            if family.loop_count() != 4 {
                break;
            }
            let destination = physical_family(momenta, family.loop_count())?;
            if let Some((mapped, mapping)) = mapped_key(&source, &destination, &slots, target, 10)?
            {
                let (native, _) =
                    four_loop::native_program(index)?.ok_or("missing shipped candidate program")?;
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
                    return native.reduce_transported_unit_mass_guarded(&transport, target, limits);
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
        Err(format!(
            "no shipped program matches the exact {}-loop factor",
            family.loop_count()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
                )
                .is_err()
        );
    }
}

#[cfg(test)]
#[path = "factor_provider_tests.rs"]
mod independent_tests;
