use std::collections::{BTreeMap, BTreeSet};

use itertools::Itertools;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore};

use super::{IntegralFamily, IntegralFamilyError};

/// A verified, measure-preserving loop-momentum map into another family.
///
/// Source denominators embed injectively into the target; unused target entries
/// receive power zero. Scalar numerator substitutions are simultaneous.
#[derive(Clone, Debug)]
pub struct IntegralMapping {
    momentum_rules: Vec<(Atom, Atom)>,
    scalar_product_rules: BTreeMap<Atom, Atom>,
    propagators: PropagatorMapping,
}

/// A source-to-target propagator embedding without a momentum substitution.
///
/// Parametric mappings establish equality of Symanzik polynomials. They do not
/// establish integration prescriptions or supply transformations for arbitrary
/// numerators; express scalar numerator factors as signed propagator powers.
#[derive(Clone, Debug)]
pub struct PropagatorMapping {
    pub(super) denominator_map: Vec<usize>,
    pub(super) target_denominator_count: usize,
}

impl PropagatorMapping {
    /// Target denominator index for each source denominator, using zero-based indices.
    pub fn denominator_map(&self) -> &[usize] {
        &self.denominator_map
    }

    /// Move signed propagator powers into target-family order.
    pub fn map_powers(&self, powers: &[i32]) -> Result<Vec<i32>, IntegralFamilyError> {
        if powers.len() != self.denominator_map.len() {
            return Err(IntegralFamilyError::InvalidPowers);
        }
        let mut target = vec![0; self.target_denominator_count];
        for (source, destination) in self.denominator_map.iter().enumerate() {
            target[*destination] = powers[source];
        }
        Ok(target)
    }
}

impl IntegralMapping {
    /// Source loop momentum names paired with expressions in target coordinates.
    pub fn momentum_rules(&self) -> &[(Atom, Atom)] {
        &self.momentum_rules
    }

    /// Target denominator index for each source denominator, using zero-based indices.
    pub fn denominator_map(&self) -> &[usize] {
        self.propagators.denominator_map()
    }

    /// Move signed propagator powers into target-family order.
    pub fn map_powers(&self, powers: &[i32]) -> Result<Vec<i32>, IntegralFamilyError> {
        self.propagators.map_powers(powers)
    }

    /// Transform a scalar numerator in compact Spenso dot notation.
    ///
    /// Tensor indices must be contracted first. This uses the very same scalar
    /// substitutions that verified the denominator mapping.
    pub fn apply(&self, expression: &Atom) -> Atom {
        expression
            .replace_map(|view, _, output| {
                if let Some(value) = self.scalar_product_rules.get(&view.to_owned()) {
                    **output = value.clone();
                }
            })
            .expand()
    }
}

/// Exact decomposition of an inverse denominator into `scale * momentum² + remainder`.
///
/// The remainder is independent of every loop momentum. For an ordinary
/// massive propagator its squared mass is `-remainder / scale`.
#[derive(Clone, Debug)]
pub struct QuadraticMomentum {
    momentum: Atom,
    scale: Atom,
    remainder: Atom,
}

impl QuadraticMomentum {
    pub fn momentum(&self) -> &Atom {
        &self.momentum
    }

    pub fn scale(&self) -> &Atom {
        &self.scale
    }

    pub fn remainder(&self) -> &Atom {
        &self.remainder
    }
}

impl IntegralFamily {
    pub(super) fn compatible_kinematics(&self, target: &Self) -> Result<(), IntegralFamilyError> {
        if self.loop_momenta.len() != target.loop_momenta.len()
            || self.kinematics.dimension() != target.kinematics.dimension()
            || self.external_momenta.iter().collect::<BTreeSet<_>>()
                != target.external_momenta.iter().collect::<BTreeSet<_>>()
        {
            return Err(IntegralFamilyError::InvalidMapping(
                "loop counts, dimensions and external momentum names must agree".into(),
            ));
        }
        for (i, p) in self.external_momenta.iter().enumerate() {
            for q in &self.external_momenta[i..] {
                if !(self.kinematics.scalar_product(p, q)?
                    - target.kinematics.scalar_product(p, q)?)
                .together()
                .is_zero()
                {
                    return Err(IntegralFamilyError::InvalidMapping(
                        "external scalar-product assumptions differ".into(),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Verify explicit source-loop images in the target momentum coordinates.
    ///
    /// Returns `None` if the Jacobian determinant is not exactly +1 or -1, or
    /// coefficients cannot be proven real, or a source denominator has no
    /// distinct equal target denominator. External
    /// names and kinematics must agree. Target families may contain additional
    /// propagators, enabling subtopology embeddings.
    pub fn mapping_to(
        &self,
        target: &Self,
        loop_images: &[Atom],
    ) -> Result<Option<IntegralMapping>, IntegralFamilyError> {
        self.compatible_kinematics(target)?;
        if loop_images.len() != self.loop_momenta.len() {
            return Err(IntegralFamilyError::InvalidMapping(
                "supply one image per source loop momentum".into(),
            ));
        }
        let rep = Minkowski {}.new_rep(self.kinematics.dimension());
        let target_momenta = target
            .loop_momenta
            .iter()
            .chain(&target.external_momenta)
            .collect::<Vec<_>>();
        for image in loop_images {
            // Use the shared kinematics parser to validate linear combinations.
            let _ = target.kinematics.scalar_product(image, image)?;
            let terms = image.coefficient_list::<i32>(&target_momenta);
            if terms.iter().any(|(momentum, _)| {
                !target_momenta
                    .iter()
                    .any(|p| p.as_view() == momentum.as_view())
            }) {
                return Err(IntegralFamilyError::InvalidMapping(
                    "images must be linear combinations of target loop and external momenta".into(),
                ));
            }
            if terms
                .iter()
                .any(|(_, coefficient)| !coefficient.is_real().is_true())
            {
                return Ok(None);
            }
            if target
                .loop_momenta
                .iter()
                .any(|p| image.contains(rep.vector(p.as_view(), []).as_view()))
            {
                return Err(IntegralFamilyError::InvalidMapping(
                    "loop-dependent coefficients are not affine momentum shifts".into(),
                ));
            }
        }
        let (jacobian, _) = Self::affine_system(loop_images, &target.loop_momenta)
            .map_err(|e| IntegralFamilyError::InvalidMapping(e.to_string()))?;
        let determinant = jacobian
            .det()
            .map_err(|e| IntegralFamilyError::InvalidMapping(e.to_string()))?
            .to_expression()
            .together();
        if determinant != Atom::one() && determinant != Atom::num(-1) {
            return Ok(None);
        }
        let mut rules = BTreeMap::new();
        for (i, momentum) in self.loop_momenta.iter().enumerate() {
            for (j, other) in self.loop_momenta.iter().enumerate().skip(i) {
                rules.insert(
                    rep.inner_product(momentum, other),
                    target
                        .kinematics
                        .scalar_product(&loop_images[i], &loop_images[j])?,
                );
            }
            for other in &self.external_momenta {
                rules.insert(
                    rep.inner_product(momentum, other),
                    target.kinematics.scalar_product(&loop_images[i], other)?,
                );
            }
        }
        for (i, p) in self.external_momenta.iter().enumerate() {
            for q in &self.external_momenta[i..] {
                rules.insert(
                    rep.inner_product(p, q),
                    target.kinematics.scalar_product(p, q)?,
                );
            }
        }
        let mut mapping = IntegralMapping {
            momentum_rules: self
                .loop_momenta
                .iter()
                .cloned()
                .zip(loop_images.iter().cloned())
                .collect(),
            scalar_product_rules: rules,
            propagators: PropagatorMapping {
                denominator_map: Vec::new(),
                target_denominator_count: target.denominators.len(),
            },
        };
        for denominator in &self.denominators {
            let transformed = mapping.apply(denominator);
            let Some(index) = target.denominators.iter().enumerate().find_map(|(i, d)| {
                (!mapping.denominator_map().contains(&i) && (&transformed - d).together().is_zero())
                    .then_some(i)
            }) else {
                return Ok(None);
            };
            mapping.propagators.denominator_map.push(index);
        }
        Ok(Some(mapping))
    }

    /// Search exact affine loop shifts using independent quadratic propagators.
    ///
    /// This includes loop permutations, reversals, mixtures and external shifts.
    /// Candidates are derived from propagator quadratic forms, then verified by
    /// [`Self::mapping_to`]. Eikonal and auxiliary propagators participate in
    /// verification but cannot supply the quadratic search basis. This is not a
    /// parametric-polynomial equivalence test for identities without a loop shift.
    pub fn find_mapping(
        &self,
        target: &Self,
        max_candidates: usize,
    ) -> Result<Option<IntegralMapping>, IntegralFamilyError> {
        self.compatible_kinematics(target)?;
        if self.denominators.len() > target.denominators.len() {
            return Ok(None);
        }
        let source_quadratics = self.quadratic_momenta()?;
        let target_quadratics = target.quadratic_momenta()?;
        let mut basis = Vec::new();
        for quadratic in &source_quadratics {
            basis.push(quadratic);
            let momenta = basis.iter().map(|q| &q.momentum).collect::<Vec<_>>();
            let (matrix, _) = Atom::system_to_matrix::<u16, _, _>(&momenta, &self.loop_momenta)
                .map_err(|e| IntegralFamilyError::InvalidMapping(e.to_string()))?;
            if matrix.rank() < basis.len() {
                basis.pop();
            }
            if basis.len() == self.loop_momenta.len() {
                break;
            }
        }
        if basis.len() != self.loop_momenta.len() {
            return Err(IntegralFamilyError::NoQuadraticBasis);
        }
        let source_momenta = basis.iter().map(|q| &q.momentum).collect::<Vec<_>>();
        let (matrix, rhs) = Self::affine_system(&source_momenta, &self.loop_momenta)
            .map_err(|e| IntegralFamilyError::InvalidMapping(e.to_string()))?;
        let inverse = matrix
            .inv()
            .map_err(|e| IntegralFamilyError::InvalidMapping(e.to_string()))?
            .into_vec()
            .into_iter()
            .map(|v| v.to_expression())
            .collect::<Vec<_>>();
        let offsets = rhs
            .into_vec()
            .into_iter()
            .map(|v| v.to_expression())
            .collect::<Vec<_>>();
        let options = basis
            .iter()
            .map(|source| {
                target_quadratics
                    .iter()
                    .enumerate()
                    .filter(|(_, q)| (&source.remainder - &q.remainder).together().is_zero())
                    .flat_map(|(i, q)| {
                        let ratio = (&q.scale / &source.scale).pow(Atom::num(1) / 2);
                        let image = ratio * &q.momentum;
                        [(i, image.clone()), (i, -image)]
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        for (attempt, choice) in options.iter().multi_cartesian_product().enumerate() {
            if attempt >= max_candidates {
                return Err(IntegralFamilyError::MappingSearchLimit(max_candidates));
            }
            if choice
                .iter()
                .map(|(i, _)| *i)
                .collect::<BTreeSet<_>>()
                .len()
                != choice.len()
            {
                continue;
            }
            let images = (0..basis.len())
                .map(|i| {
                    choice
                        .iter()
                        .enumerate()
                        .map(|(j, (_, momentum))| {
                            &inverse[i * basis.len() + j] * (momentum + &offsets[j])
                        })
                        .sum::<Atom>()
                        .expand()
                })
                .collect::<Vec<_>>();
            if let Some(mapping) = self.mapping_to(target, &images)? {
                return Ok(Some(mapping));
            }
        }
        Ok(None)
    }

    fn quadratic_momenta(&self) -> Result<Vec<QuadraticMomentum>, IntegralFamilyError> {
        (0..self.denominators.len())
            .filter_map(|index| self.quadratic_denominator(index).transpose())
            .collect()
    }

    /// Decompose one ordered inverse denominator as a squared momentum plus a remainder.
    ///
    /// The normalization and external momentum shifts are retained exactly.
    /// Returns `None` for linear denominators or quadratic forms that are not
    /// the square of a single momentum. An out-of-range index is an error.
    /// This is the same decomposition used by the affine momentum-map search.
    pub fn quadratic_denominator(
        &self,
        index: usize,
    ) -> Result<Option<QuadraticMomentum>, IntegralFamilyError> {
        let rep = Minkowski {}.new_rep(self.kinematics.dimension());
        let denominator = self.denominators.get(index).ok_or_else(|| {
            IntegralFamilyError::InvalidBasis(format!("denominator index {index} is out of range"))
        })?;
        let coefficients = denominator
            .coefficient_list::<i32>(&self.scalar_products)
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let Some((pivot, scale)) = self.loop_momenta.iter().find_map(|p| {
            coefficients
                .get(&rep.inner_product(p, p))
                .filter(|c| !c.is_zero())
                .map(|c| (p, c))
        }) else {
            return Ok(None);
        };
        let mut momentum = pivot.clone();
        for other in self.loop_momenta.iter().chain(&self.external_momenta) {
            if other == pivot {
                continue;
            }
            if let Some(coefficient) = coefficients.get(&rep.inner_product(pivot, other)) {
                momentum += coefficient / (Atom::num(2) * scale) * other;
            }
        }
        let remainder = (denominator
            - scale * self.kinematics.scalar_product(&momentum, &momentum)?)
        .together();
        if self.rank_of(std::slice::from_ref(&remainder))? == 0 {
            Ok(Some(QuadraticMomentum {
                momentum,
                scale: scale.clone(),
                remainder,
            }))
        } else {
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_kinematics::Kinematics;
    use symbolica::parse;

    #[test]
    fn quadratic_denominators_retain_order_shifts_and_normalization() {
        let [k, p, mass, s] = ["quad::k", "quad::p", "quad::m2", "quad::s"]
            .map(|name| symbolica::symbol!(name).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, s)
            .unwrap();
        let momentum = &k + Atom::num(3) / 2 * &p;
        let denominator =
            Atom::num(2) * (kin.scalar_product(&momentum, &momentum).unwrap() - &mass);
        let family = IntegralFamily::new(
            vec![k.clone()],
            vec![p.clone()],
            vec![
                kin.scalar_product(&k, &p).unwrap(),
                denominator.clone(),
                mass.clone(),
            ],
            &kin,
        )
        .unwrap();
        assert!(family.quadratic_denominator(0).unwrap().is_none());
        let quadratic = family.quadratic_denominator(1).unwrap().unwrap();
        assert_eq!(quadratic.momentum(), &momentum);
        assert_eq!(quadratic.scale(), &Atom::num(2));
        assert!(
            (quadratic.remainder() + Atom::num(2) * &mass)
                .expand()
                .is_zero()
        );
        assert!(
            (quadratic.scale()
                * kin
                    .scalar_product(quadratic.momentum(), quadratic.momentum())
                    .unwrap()
                + quadratic.remainder()
                - denominator)
                .expand()
                .is_zero()
        );
        assert!(family.quadratic_denominator(2).unwrap().is_none());
        assert!(family.quadratic_denominator(3).is_err());
    }

    #[test]
    fn light_cone_shift_basis_maps_to_itself() {
        let [k, q, n, b] = [
            "map_lightcone::k",
            "map_lightcone::q",
            "map_lightcone::n",
            "map_lightcone::b",
        ]
        .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), q.clone(), n.clone(), b.clone()])
            .unwrap()
            .with_mass_squared(&n, Atom::Zero)
            .unwrap()
            .with_mass_squared(&b, Atom::Zero)
            .unwrap()
            .with_scalar_product(&n, &b, Atom::num(2))
            .unwrap();
        let mixed = &k - &q + &n - &b / 2;
        let momenta = [k.clone(), mixed.clone(), q.clone()];
        let family = IntegralFamily::new(
            vec![k.clone(), q.clone()],
            vec![n.clone(), b.clone()],
            momenta
                .iter()
                .map(|v| kin.scalar_product(v, v).unwrap())
                .collect(),
            &kin,
        )
        .unwrap();
        assert!(family.mapping_to(&family, &[k, q]).unwrap().is_some());
        assert!(family.find_mapping(&family, 100).unwrap().is_some());
    }

    #[test]
    fn reciprocal_loop_rescalings_preserve_the_measure() {
        let [k, q, l, r] = [
            "map_scale::k",
            "map_scale::q",
            "map_scale::l",
            "map_scale::r",
        ]
        .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), q.clone(), l.clone(), r.clone()])
            .unwrap();
        let source = IntegralFamily::new(
            vec![k.clone(), q.clone()],
            vec![],
            vec![
                kin.scalar_product(&k, &k).unwrap() - 1,
                kin.scalar_product(&q, &q).unwrap() - 2,
            ],
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![l.clone(), r.clone()],
            vec![],
            vec![
                4 * kin.scalar_product(&l, &l).unwrap() - 1,
                kin.scalar_product(&r, &r).unwrap() / 4 - 2,
            ],
            &kin,
        )
        .unwrap();
        assert!(
            source
                .mapping_to(&target, &[2 * &l, &r / 2])
                .unwrap()
                .is_some()
        );
        let mapping = source.find_mapping(&target, 100).unwrap().unwrap();
        for (i, j) in mapping.denominator_map().iter().enumerate() {
            assert!(
                (mapping.apply(&source.denominators[i]) - &target.denominators[*j])
                    .together()
                    .is_zero()
            );
        }
    }

    #[test]
    fn feyncalc_three_loop_topology_mapping() {
        // topo4 -> topo1 from the FCLoopFindTopologyMappings manual example.
        let [a, b, c, p] = ["map_fc3::a", "map_fc3::b", "map_fc3::c", "map_fc3::p"]
            .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([a.clone(), b.clone(), c.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let common = [c.clone(), b.clone(), a.clone(), &b + &c];
        let source_momenta = common.iter().cloned().chain([
            &a + &c,
            &b - &p,
            &a - &p,
            &a + &c - &p,
            &a + &b + &c - &p,
        ]);
        let target_momenta = common.iter().cloned().chain([
            &b - &p,
            &a - &p,
            &b + &c - &p,
            &a + &c - &p,
            &a + &b + &c - &p,
        ]);
        let source = IntegralFamily::new(
            vec![a.clone(), b.clone(), c.clone()],
            vec![p.clone()],
            source_momenta
                .map(|v| kin.scalar_product(&v, &v).unwrap())
                .collect(),
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![a.clone(), b.clone(), c.clone()],
            vec![p.clone()],
            target_momenta
                .map(|v| kin.scalar_product(&v, &v).unwrap())
                .collect(),
            &kin,
        )
        .unwrap();
        let explicit = source
            .mapping_to(&target, &[&p - &b, &p - &a, -&c])
            .unwrap()
            .unwrap();
        assert_eq!(explicit.denominator_map(), &[0, 5, 4, 7, 6, 2, 1, 3, 8]);
        let found = source.find_mapping(&target, 10_000).unwrap().unwrap();
        for (i, j) in found.denominator_map().iter().enumerate() {
            assert!(
                (found.apply(&source.denominators[i]) - &target.denominators[*j])
                    .together()
                    .is_zero()
            );
        }
    }

    #[test]
    fn bubble_shift_maps_denominators_powers_and_numerators() {
        let [k, l, p] = [
            parse!("map_bubble::k"),
            parse!("map_bubble::l"),
            parse!("map_bubble::p"),
        ];
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let square = |v: &Atom| kin.scalar_product(v, v).unwrap();
        let source = IntegralFamily::new(
            vec![k.clone()],
            vec![p.clone()],
            vec![square(&k) - parse!("m1"), square(&(&k + &p)) - parse!("m2")],
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![l.clone()],
            vec![p.clone()],
            vec![square(&l) - parse!("m2"), square(&(&l - &p)) - parse!("m1")],
            &kin,
        )
        .unwrap();
        let mapping = source.find_mapping(&target, 100).unwrap().unwrap();
        assert_eq!(mapping.denominator_map(), &[1, 0]);
        assert_eq!(mapping.map_powers(&[2, -1]).unwrap(), vec![-1, 2]);
        assert!(mapping.map_powers(&[1]).is_err());
        let kp = kin.scalar_product(&k, &p).unwrap();
        let expected = kin
            .scalar_product(&mapping.momentum_rules()[0].1, &p)
            .unwrap();
        assert!((mapping.apply(&kp) - expected).expand().is_zero());
        for (i, j) in mapping.denominator_map().iter().enumerate() {
            assert!(
                (mapping.apply(&source.denominators[i]) - &target.denominators[*j])
                    .together()
                    .is_zero()
            );
        }
        assert!(
            source
                .mapping_to(&target, &[Atom::num(2) * &l])
                .unwrap()
                .is_none()
        );
        assert!(matches!(
            source.find_mapping(&target, 0),
            Err(IntegralFamilyError::MappingSearchLimit(0))
        ));
    }

    #[test]
    fn mixed_loop_basis_and_subtopology_embedding() {
        let [k, q, l, r, p] = [
            "map_mixed::k",
            "map_mixed::q",
            "map_mixed::l",
            "map_mixed::r",
            "map_mixed::p",
        ]
        .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), q.clone(), l.clone(), r.clone(), p.clone()])
            .unwrap();
        let square = |v: &Atom| kin.scalar_product(v, v).unwrap();
        let source = IntegralFamily::new(
            vec![k.clone(), q.clone()],
            vec![p.clone()],
            vec![square(&k) - 1, square(&q) - 2, square(&(&k - &q + &p)) - 3],
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![l.clone(), r.clone()],
            vec![p.clone()],
            vec![
                square(&(&l + &p)) - 3,
                square(&(&l + &r + &p)) - 1,
                square(&(&r + &p)) - 2,
                square(&l),
            ],
            &kin,
        )
        .unwrap();
        let mapping = source.find_mapping(&target, 1000).unwrap().unwrap();
        assert_eq!(mapping.denominator_map(), &[1, 2, 0]);
        assert_eq!(mapping.map_powers(&[1, 2, 3]).unwrap(), vec![3, 1, 2, 0]);
        let images = [&l + &r + &p, &r + &p];
        let explicit = source.mapping_to(&target, &images).unwrap().unwrap();
        let kq = kin.scalar_product(&k, &q).unwrap();
        assert!(
            (explicit.apply(&kq) - kin.scalar_product(&images[0], &images[1]).unwrap())
                .expand()
                .is_zero()
        );
        assert!(target.find_mapping(&source, 100).unwrap().is_none());
    }

    #[test]
    fn eikonal_maps_are_verified_and_complex_contour_changes_are_rejected() {
        let [k, q, l, r, p] = [
            "map_eikonal::k",
            "map_eikonal::q",
            "map_eikonal::l",
            "map_eikonal::r",
            "map_eikonal::p",
        ]
        .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), q.clone(), l.clone(), r.clone(), p.clone()])
            .unwrap();
        let source = IntegralFamily::new(
            vec![k.clone()],
            vec![p.clone()],
            vec![kin.scalar_product(&k, &p).unwrap()],
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![l.clone()],
            vec![p.clone()],
            vec![kin.scalar_product(&(&l + &p), &p).unwrap()],
            &kin,
        )
        .unwrap();
        assert!(source.mapping_to(&target, &[&l + &p]).unwrap().is_some());
        assert!(matches!(
            source.find_mapping(&target, 100),
            Err(IntegralFamilyError::NoQuadraticBasis)
        ));
        let square = |v: &Atom| kin.scalar_product(v, v).unwrap();
        let source = IntegralFamily::new(
            vec![k.clone(), q.clone()],
            vec![],
            vec![square(&k), square(&q)],
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![l.clone(), r.clone()],
            vec![],
            vec![-square(&l), -square(&r)],
            &kin,
        )
        .unwrap();
        assert!(
            source
                .mapping_to(&target, &[parse!("𝑖") * &l, -parse!("𝑖") * &r])
                .unwrap()
                .is_none()
        );
    }
}
