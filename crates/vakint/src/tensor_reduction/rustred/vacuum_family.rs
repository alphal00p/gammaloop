use std::collections::{HashMap, hash_map::Entry};

use ::rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use ::rustred::family::isp::IspCompletion;
use ::rustred::family::presentation::{
    AuxiliaryDenominator, CommonMassScale, DenominatorRole, FamilyConventions, FamilyPresentation,
    MetricConvention, MomentumCombination, MomentumRouting, PhysicalPropagator,
    PropagatorConvention,
};
use ::rustred::family::{AffineDenominator, IntegralFamilyError, IntegralKey};
use symbolica::atom::Atom;

use super::{TensorReductionError, is_exact_real_scalar};

/// One physical propagator after Vakint has matched and canonicalized a
/// single-scale vacuum topology.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct PhysicalVacuumPropagator {
    canonical_id: usize,
    loop_momentum: Vec<i64>,
    mass_squared: Atom,
    power: i64,
}

impl PhysicalVacuumPropagator {
    pub(super) fn new(
        canonical_id: usize,
        loop_momentum: Vec<i64>,
        mass_squared: Atom,
        power: i64,
    ) -> Self {
        Self {
            canonical_id,
            loop_momentum,
            mass_squared,
            power,
        }
    }
}

/// Topology-neutral input owned by Vakint's matcher boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct VacuumFamilySpec {
    loop_count: usize,
    propagators: Vec<PhysicalVacuumPropagator>,
}

impl VacuumFamilySpec {
    pub(super) fn new(loop_count: usize, propagators: Vec<PhysicalVacuumPropagator>) -> Self {
        Self {
            loop_count,
            propagators,
        }
    }
}

/// One invocation's authenticated structural families.
///
/// Concrete masses and integral powers deliberately do not participate in
/// this cache: they are validated and bound afresh for every term. The owner
/// lives beside Vakint's transactional expression clone and is discarded on
/// either success or failure.
#[derive(Debug, Default)]
pub(super) struct VacuumFamilyCache {
    families: HashMap<VacuumFamilyKey, VacuumFamilyStructure>,
}

impl VacuumFamilyCache {
    pub(super) fn bind<'cache>(
        &'cache mut self,
        term_index: usize,
        spec: VacuumFamilySpec,
    ) -> Result<VacuumFamilyBridge<'cache>, TensorReductionError> {
        validate_term_spec(term_index, &spec)?;
        let key = VacuumFamilyKey::from_spec(&spec);
        let structure = match self.families.entry(key) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let structure = VacuumFamilyStructure::try_new(term_index, entry.key())?;
                entry.insert(structure)
            }
        };
        VacuumFamilyBridge::bind(term_index, structure, spec)
    }

    #[cfg(test)]
    fn family_count(&self) -> usize {
        self.families.len()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct VacuumFamilyKey {
    loop_count: usize,
    physical_denominators: Box<[PhysicalDenominatorKey]>,
}

impl VacuumFamilyKey {
    fn from_spec(spec: &VacuumFamilySpec) -> Self {
        Self {
            loop_count: spec.loop_count,
            physical_denominators: spec
                .propagators
                .iter()
                .map(|propagator| PhysicalDenominatorKey {
                    canonical_id: propagator.canonical_id,
                    loop_momentum: propagator.loop_momentum.clone().into_boxed_slice(),
                })
                .collect(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct PhysicalDenominatorKey {
    canonical_id: usize,
    loop_momentum: Box<[i64]>,
}

/// Expensive authenticated data shared by structurally identical terms.
#[derive(Debug)]
struct VacuumFamilyStructure {
    presentation: FamilyPresentation,
    physical_denominator_ids: Box<[usize]>,
    appended_coordinate_ordinals: Box<[usize]>,
}

impl VacuumFamilyStructure {
    fn try_new(term_index: usize, key: &VacuumFamilyKey) -> Result<Self, TensorReductionError> {
        let propagator_count = key.physical_denominators.len();
        let coefficients = CoefficientContext::try_new(["d", "m2"])
            .map_err(|source| TensorReductionError::RustRedCoefficientContext { source })?;
        let dimension = coefficients
            .parameter("d")
            .expect("the adapter registers the dimension parameter");
        let abstract_mass = coefficients
            .parameter("m2")
            .expect("the adapter registers the mass parameter");
        let negative_mass = coefficients
            .try_neg(&abstract_mass, ExactAlgebraLimits::default())
            .map_err(|source| TensorReductionError::RustRedFamily {
                term: term_index,
                source: IntegralFamilyError::ExactAlgebra(source),
            })?;

        let scalar_product_count = key
            .loop_count
            .checked_add(1)
            .and_then(|successor| key.loop_count.checked_mul(successor))
            .and_then(|twice| twice.checked_div(2))
            .ok_or_else(|| TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: "vacuum scalar-product count overflowed".to_owned(),
            })?;
        let mut affine_denominators = Vec::new();
        let mut physical_roles = Vec::new();
        let mut physical_denominator_ids = Vec::new();
        affine_denominators.reserve(propagator_count);
        physical_roles.reserve(propagator_count);
        physical_denominator_ids.reserve(propagator_count);
        for propagator in &key.physical_denominators {
            let abstract_momentum = propagator
                .loop_momentum
                .iter()
                .map(|coefficient| coefficients.integer(*coefficient))
                .collect::<Vec<_>>();
            let row = momentum_squared_row(
                term_index,
                &coefficients,
                &abstract_momentum,
                scalar_product_count,
            )?;
            affine_denominators.push(AffineDenominator::new(negative_mass.clone(), row));
            physical_roles.push(DenominatorRole::Physical(PhysicalPropagator::new(
                format!("propagator-{}", propagator.canonical_id),
                MomentumCombination::new(abstract_momentum, Vec::new()),
                abstract_mass.clone(),
            )));
            physical_denominator_ids.push(propagator.canonical_id);
        }

        let loop_names = (1..=key.loop_count)
            .map(|loop_index| format!("k{loop_index}"))
            .collect::<Vec<_>>();
        let power_shifts = vec![coefficients.zero(); propagator_count];
        let completion = IspCompletion::try_new(
            format!("vakint-native-{}-loop-single-scale-vacuum", key.loop_count),
            loop_names.clone(),
            Vec::new(),
            coefficients,
            dimension,
            affine_denominators,
            Vec::new(),
            power_shifts,
        )
        .map_err(|source| TensorReductionError::RustRedIspCompletion {
            term: term_index,
            source,
        })?;
        let appended_coordinate_ordinals = completion
            .appended_coordinate_ordinals()
            .to_vec()
            .into_boxed_slice();
        let family = completion.into_family();
        let coefficients = family.coefficient_context();
        let mut denominator_roles = physical_roles;
        denominator_roles.extend(appended_coordinate_ordinals.iter().map(|ordinal| {
            DenominatorRole::Auxiliary(AuxiliaryDenominator::new(format!(
                "isp-coordinate-{ordinal}"
            )))
        }));
        let loop_linear = (0..key.loop_count)
            .map(|row| {
                (0..key.loop_count)
                    .map(|column| coefficients.integer(i64::from(row == column)))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let routing = MomentumRouting::new(
            loop_names,
            Vec::new(),
            loop_linear,
            vec![Vec::new(); key.loop_count],
            Vec::new(),
        );
        let common_mass_scale = CommonMassScale::new(
            family
                .coefficient_context()
                .parameter("m2")
                .expect("the adapter registers the mass parameter"),
        );
        let presentation = FamilyPresentation::try_new(
            family,
            denominator_roles,
            routing,
            FamilyConventions::new(
                MetricConvention::MinkowskiMostlyMinus,
                PropagatorConvention::MOMENTUM_SQUARED_MINUS_MASS_SQUARED,
            ),
            Some(common_mass_scale),
        )
        .map_err(|source| TensorReductionError::RustRedPresentation {
            term: term_index,
            source,
        })?;

        Ok(Self {
            presentation,
            physical_denominator_ids: physical_denominator_ids.into_boxed_slice(),
            appended_coordinate_ordinals,
        })
    }
}

/// Per-term binding to one cached structural family.
#[derive(Debug)]
pub(super) struct VacuumFamilyBridge<'cache> {
    structure: &'cache VacuumFamilyStructure,
    common_mass_squared: Atom,
    base_integral: IntegralKey,
}

impl<'cache> VacuumFamilyBridge<'cache> {
    fn bind(
        term_index: usize,
        structure: &'cache VacuumFamilyStructure,
        spec: VacuumFamilySpec,
    ) -> Result<Self, TensorReductionError> {
        let common_mass_squared = spec.propagators[0].mass_squared.clone();
        let mut base_powers = spec
            .propagators
            .into_iter()
            .map(|propagator| propagator.power)
            .collect::<Vec<_>>();
        base_powers.extend(std::iter::repeat_n(
            0,
            structure.appended_coordinate_ordinals.len(),
        ));
        let base_integral = IntegralKey::try_new(base_powers).map_err(|source| {
            TensorReductionError::RustRedTensor {
                term: term_index,
                source: source.into(),
            }
        })?;
        Ok(Self {
            structure,
            common_mass_squared,
            base_integral,
        })
    }

    pub(super) const fn presentation(&self) -> &FamilyPresentation {
        &self.structure.presentation
    }

    pub(super) const fn common_mass_squared(&self) -> &Atom {
        &self.common_mass_squared
    }

    pub(super) const fn base_integral(&self) -> &IntegralKey {
        &self.base_integral
    }

    pub(super) fn physical_denominator_ids(&self) -> &[usize] {
        &self.structure.physical_denominator_ids
    }

    pub(super) fn appended_coordinate_ordinals(&self) -> &[usize] {
        &self.structure.appended_coordinate_ordinals
    }
}

fn validate_term_spec(
    term_index: usize,
    spec: &VacuumFamilySpec,
) -> Result<(), TensorReductionError> {
    let propagator_count = spec.propagators.len();
    if spec.loop_count == 0 || propagator_count == 0 {
        return Err(TensorReductionError::RustRedUnsupportedFamily {
            term: term_index,
            loop_count: spec.loop_count,
            propagator_count,
        });
    }

    let mut previous_id = None;
    let common_mass_squared = &spec.propagators[0].mass_squared;
    for propagator in &spec.propagators {
        if propagator.loop_momentum.len() != spec.loop_count
            || propagator
                .loop_momentum
                .iter()
                .all(|coefficient| *coefficient == 0)
        {
            return Err(TensorReductionError::RustRedUnsupportedMomentum {
                term: term_index,
                propagator: propagator.canonical_id,
                momentum: format!("{:?}", propagator.loop_momentum),
            });
        }
        if previous_id.is_some_and(|previous| previous >= propagator.canonical_id) {
            return Err(TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: "physical propagator IDs are not in strict canonical order".to_owned(),
            });
        }
        previous_id = Some(propagator.canonical_id);
        if propagator.mass_squared.is_zero()
            || !is_exact_real_scalar(propagator.mass_squared.as_view())
            || propagator.mass_squared != *common_mass_squared
        {
            return Err(TensorReductionError::RustRedUnsupportedMass {
                term: term_index,
                propagator: propagator.canonical_id,
                mass: propagator.mass_squared.to_string(),
            });
        }
    }
    Ok(())
}

fn momentum_squared_row(
    term_index: usize,
    coefficients: &CoefficientContext,
    momentum: &[Coefficient],
    scalar_product_count: usize,
) -> Result<Vec<Coefficient>, TensorReductionError> {
    let mut row = Vec::with_capacity(scalar_product_count);
    for left in 0..momentum.len() {
        for right in left..momentum.len() {
            let mut coefficient = coefficients
                .try_mul(
                    &momentum[left],
                    &momentum[right],
                    ExactAlgebraLimits::default(),
                )
                .map_err(|source| TensorReductionError::RustRedFamily {
                    term: term_index,
                    source: IntegralFamilyError::ExactAlgebra(source),
                })?;
            if left != right {
                coefficient = coefficients
                    .try_mul(
                        &coefficient,
                        &coefficients.integer(2),
                        ExactAlgebraLimits::default(),
                    )
                    .map_err(|source| TensorReductionError::RustRedFamily {
                        term: term_index,
                        source: IntegralFamilyError::ExactAlgebra(source),
                    })?;
            }
            row.push(coefficient);
        }
    }
    debug_assert_eq!(row.len(), scalar_product_count);
    Ok(row)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn propagator(id: usize, momentum: &[i64], mass: i64, power: i64) -> PhysicalVacuumPropagator {
        PhysicalVacuumPropagator::new(id, momentum.to_vec(), Atom::num(mass), power)
    }

    #[test]
    fn complete_two_loop_basis_retains_physical_order_and_powers() {
        let mut cache = VacuumFamilyCache::default();
        let bridge = cache
            .bind(
                0,
                VacuumFamilySpec::new(
                    2,
                    vec![
                        propagator(1, &[1, 0], 1, 1),
                        propagator(2, &[0, 1], 1, 2),
                        propagator(3, &[1, 1], 1, -1),
                    ],
                ),
            )
            .unwrap();

        assert_eq!(bridge.physical_denominator_ids(), &[1, 2, 3]);
        assert!(bridge.appended_coordinate_ordinals().is_empty());
        assert_eq!(bridge.base_integral().powers(), &[1, 2, -1]);
        assert_eq!(bridge.common_mass_squared(), &Atom::num(1));
        assert!(
            bridge
                .presentation()
                .denominator_roles()
                .iter()
                .all(|role| matches!(role, DenominatorRole::Physical(_)))
        );
        let family = bridge.presentation().family();
        let coefficients = family.coefficient_context();
        assert_eq!(
            family.denominators()[2].coefficients(),
            &[
                coefficients.one(),
                coefficients.integer(2),
                coefficients.one()
            ]
        );
    }

    #[test]
    fn pinched_two_loop_family_appends_cross_isp_after_physical_rows() {
        let mut cache = VacuumFamilyCache::default();
        let bridge = cache
            .bind(
                0,
                VacuumFamilySpec::new(
                    2,
                    vec![propagator(1, &[1, 0], 1, 2), propagator(2, &[0, 1], 1, -1)],
                ),
            )
            .unwrap();

        assert_eq!(bridge.physical_denominator_ids(), &[1, 2]);
        assert_eq!(bridge.appended_coordinate_ordinals(), &[1]);
        assert_eq!(bridge.base_integral().powers(), &[2, -1, 0]);
        assert!(matches!(
            bridge.presentation().denominator_roles()[2],
            DenominatorRole::Auxiliary(_)
        ));
        let family = bridge.presentation().family();
        let coefficients = family.coefficient_context();
        let expected = [
            vec![coefficients.one(), coefficients.zero(), coefficients.zero()],
            vec![coefficients.zero(), coefficients.zero(), coefficients.one()],
            vec![coefficients.zero(), coefficients.one(), coefficients.zero()],
        ];
        for (denominator, expected_row) in family.denominators().iter().zip(expected) {
            assert_eq!(denominator.coefficients(), expected_row);
        }
    }

    #[test]
    fn cache_reuses_structure_but_rebinds_concrete_mass_and_powers() {
        let mut cache = VacuumFamilyCache::default();
        {
            let first = cache
                .bind(
                    0,
                    VacuumFamilySpec::new(
                        2,
                        vec![propagator(1, &[1, 0], 1, 1), propagator(2, &[0, 1], 1, 2)],
                    ),
                )
                .unwrap();
            assert_eq!(first.common_mass_squared(), &Atom::num(1));
            assert_eq!(first.base_integral().powers(), &[1, 2, 0]);
        }
        {
            let second = cache
                .bind(
                    1,
                    VacuumFamilySpec::new(
                        2,
                        vec![propagator(1, &[1, 0], 3, -2), propagator(2, &[0, 1], 3, 4)],
                    ),
                )
                .unwrap();
            assert_eq!(second.common_mass_squared(), &Atom::num(3));
            assert_eq!(second.base_integral().powers(), &[-2, 4, 0]);
        }
        let error = cache
            .bind(
                2,
                VacuumFamilySpec::new(
                    2,
                    vec![propagator(1, &[1, 0], 3, 1), propagator(2, &[0, 1], 4, 2)],
                ),
            )
            .unwrap_err();
        assert!(matches!(
            error,
            TensorReductionError::RustRedUnsupportedMass {
                term: 2,
                propagator: 2,
                ..
            }
        ));
        assert_eq!(cache.family_count(), 1);
    }
}
