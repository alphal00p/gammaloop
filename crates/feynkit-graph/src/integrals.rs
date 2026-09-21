//! Scalar-product bases for loop-integral families.
//!
//! Propagators are affine forms in the independent loop scalar products.
//! Symbolica owns the polynomial conversion, rank calculation and exact solve.

use std::collections::{BTreeMap, BTreeSet};

use feynkit_kinematics::{Kinematics, SymbolicKinematicsError};
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore, AtomView};
use thiserror::Error;

mod diagram;
mod mapping;
mod parametric;
mod partial_fraction;
pub use mapping::{IntegralMapping, PropagatorMapping};

#[derive(Clone, Debug, Error)]
pub enum IntegralFamilyError {
    #[error("an integral family needs at least one loop momentum")]
    NoLoops,
    #[error("loop and external momentum names must be distinct: {0}")]
    RepeatedMomentum(Atom),
    #[error(transparent)]
    Kinematics(#[from] SymbolicKinematicsError),
    #[error("loop scalar products cannot have on-shell assumptions: {0}")]
    LoopAssumption(Atom),
    #[error("invalid affine propagator basis: {0}")]
    InvalidBasis(String),
    #[error(
        "the {denominators} denominators have rank {rank}; partial-fraction them before completing the family"
    )]
    Dependent { denominators: usize, rank: usize },
    #[error(
        "the denominator basis has rank {rank}, but {required} loop scalar products are required"
    )]
    Incomplete { rank: usize, required: usize },
    #[error("supply one distinct symbol or labeled call per denominator")]
    InvalidLabels,
    #[error("supply one signed integer power per denominator")]
    InvalidPowers,
    #[error("partial fractioning exceeded its {0}-state budget")]
    PartialFractionLimit(usize),
    #[error("a propagator power overflowed during partial fractioning")]
    PowerOverflow,
    #[error("incompatible integral-family mapping: {0}")]
    InvalidMapping(String),
    #[error("mapping search exceeded its {0}-candidate budget")]
    MappingSearchLimit(usize),
    #[error(
        "automatic shift search requires quadratic propagators spanning all loop directions; supply explicit loop images instead"
    )]
    NoQuadraticBasis,
}

/// An ordered propagator family in ordinary Spenso scalar-product notation.
///
/// Denominators are inverse propagators, for example `(k+p)^2-m^2`, expanded
/// using [`Kinematics::scalar_product`]. Linear (eikonal) denominators are also
/// accepted. The external momentum list must be an independent basis; eliminate
/// dependent external momenta before constructing a family.
///
/// The family tracks algebraic completeness and verifies affine loop-momentum
/// maps. Parametric scaling certificates detect scaleless sectors in dimensional
/// regularization. Integration prescriptions and IBP solutions require additional
/// information beyond this algebraic representation.
#[derive(Clone, Debug)]
pub struct IntegralFamily {
    kinematics: Kinematics,
    loop_momenta: Vec<Atom>,
    external_momenta: Vec<Atom>,
    scalar_products: Vec<Atom>,
    denominators: Vec<Atom>,
    rank: usize,
}

impl IntegralFamily {
    /// Construct a family and compute its rank over the external invariants.
    pub fn new(
        loop_momenta: Vec<Atom>,
        external_momenta: Vec<Atom>,
        denominators: Vec<Atom>,
        kinematics: &Kinematics,
    ) -> Result<Self, IntegralFamilyError> {
        if loop_momenta.is_empty() {
            return Err(IntegralFamilyError::NoLoops);
        }
        let mut seen = BTreeSet::new();
        for momentum in loop_momenta.iter().chain(&external_momenta) {
            if !seen.insert(momentum) {
                return Err(IntegralFamilyError::RepeatedMomentum(momentum.clone()));
            }
        }
        let kinematics = kinematics.clone().with_momenta(seen.into_iter().cloned())?;
        let rep = Minkowski {}.new_rep(kinematics.dimension());
        let mut scalar_products = Vec::new();
        for (i, momentum) in loop_momenta.iter().enumerate() {
            for other in loop_momenta[i..].iter().chain(&external_momenta) {
                let product = rep.inner_product(momentum, other);
                if kinematics.apply(&product) != product {
                    return Err(IntegralFamilyError::LoopAssumption(product));
                }
                scalar_products.push(product);
            }
        }
        let denominators = denominators.iter().map(|d| kinematics.apply(d)).collect();
        let mut family = Self {
            kinematics,
            loop_momenta,
            external_momenta,
            scalar_products,
            denominators,
            rank: 0,
        };
        family.rank = family.rank_of(&family.denominators)?;
        Ok(family)
    }

    /// Integrated momentum names, in their original order.
    pub fn loop_momenta(&self) -> &[Atom] {
        &self.loop_momenta
    }

    /// Scoped kinematics, including the declared family momentum names.
    pub fn kinematics(&self) -> &Kinematics {
        &self.kinematics
    }

    /// Independent external momentum names, in their original order.
    pub fn external_momenta(&self) -> &[Atom] {
        &self.external_momenta
    }

    /// Loop-loop and loop-external products spanning the numerator space.
    pub fn scalar_products(&self) -> &[Atom] {
        &self.scalar_products
    }

    /// Ordered inverse propagators, including any completion terms.
    pub fn denominators(&self) -> &[Atom] {
        &self.denominators
    }

    /// Number of independent affine forms in the loop scalar products.
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Whether the inverse propagators span every loop scalar product.
    pub fn is_complete(&self) -> bool {
        self.rank == self.scalar_products.len()
    }

    /// Whether no denominator can be eliminated by an affine relation.
    pub fn is_independent(&self) -> bool {
        self.rank == self.denominators.len()
    }

    /// Retain the positive-power propagators of an integral's sector.
    ///
    /// Zero and negative powers are omitted; loop variables and external
    /// kinematics are retained. This describes sector support, not an algebraic
    /// removal of numerator factors from the original integrand.
    pub fn sector(&self, powers: &[i32]) -> Result<Self, IntegralFamilyError> {
        if powers.len() != self.denominators.len() {
            return Err(IntegralFamilyError::InvalidPowers);
        }
        Self::new(
            self.loop_momenta.clone(),
            self.external_momenta.clone(),
            self.denominators
                .iter()
                .zip(powers)
                .filter(|(_, power)| **power > 0)
                .map(|(denominator, _)| denominator.clone())
                .collect(),
            &self.kinematics,
        )
    }

    fn rank_of(&self, denominators: &[Atom]) -> Result<usize, IntegralFamilyError> {
        if denominators.is_empty() {
            return Ok(0);
        }
        let rep = Minkowski {}.new_rep(self.kinematics.dimension());
        let loop_vectors = self
            .loop_momenta
            .iter()
            .map(|p| rep.vector(p.as_view(), []))
            .collect::<Vec<_>>();
        // Reject hidden nonlinear dependence, including functions or inverse
        // powers of a loop scalar product treated as polynomial coefficients.
        for denominator in denominators {
            for (key, coefficient) in denominator.coefficient_list::<i32>(&self.scalar_products) {
                if !(key.is_one() || self.scalar_products.contains(&key))
                    || self
                        .scalar_products
                        .iter()
                        .any(|p| coefficient.contains(p.as_view()))
                    || loop_vectors
                        .iter()
                        .chain(&self.loop_momenta)
                        .any(|p| coefficient.contains(p.as_view()))
                {
                    return Err(IntegralFamilyError::InvalidBasis(format!(
                        "denominator is not affine in the declared loop scalar products: {denominator}"
                    )));
                }
            }
        }
        Atom::system_to_matrix::<u16, _, _>(denominators, &self.scalar_products)
            .map(|(matrix, _)| matrix.rank())
            .map_err(|error| IntegralFamilyError::InvalidBasis(error.to_string()))
    }

    /// Append irreducible scalar products until the propagators form a basis.
    ///
    /// Original propagators retain their positions. Added entries are auxiliary
    /// inverse propagators, whose nonpositive powers represent numerator factors.
    /// A dependent family must be partial-fractioned first.
    pub fn complete(&self) -> Result<Self, IntegralFamilyError> {
        self.require_independent()?;
        let mut result = self.clone();
        for product in &self.scalar_products {
            if result.is_complete() {
                break;
            }
            result.denominators.push(product.clone());
            let rank = result.rank_of(&result.denominators)?;
            if rank > result.rank {
                result.rank = rank;
            } else {
                result.denominators.pop();
            }
        }
        Ok(result)
    }

    fn require_independent(&self) -> Result<(), IntegralFamilyError> {
        if !self.is_independent() {
            return Err(IntegralFamilyError::Dependent {
                denominators: self.denominators.len(),
                rank: self.rank,
            });
        }
        Ok(())
    }

    /// Solve every loop scalar product in terms of labeled inverse propagators.
    ///
    /// Labels are ordinary Symbolica symbols or calls, such as `d1` or `D(1)`.
    /// The family must be independent and complete. The returned substitutions
    /// can also be applied with ordinary Symbolica operations.
    pub fn scalar_product_rules(
        &self,
        labels: &[Atom],
    ) -> Result<BTreeMap<Atom, Atom>, IntegralFamilyError> {
        self.require_independent()?;
        if !self.is_complete() {
            return Err(IntegralFamilyError::Incomplete {
                rank: self.rank,
                required: self.scalar_products.len(),
            });
        }
        self.validate_labels(labels)?;
        let equations = self
            .denominators
            .iter()
            .zip(labels)
            .map(|(d, label)| d - label)
            .collect::<Vec<_>>();
        let (matrix, rhs) = Atom::system_to_matrix::<u16, _, _>(&equations, &self.scalar_products)
            .map_err(|error| IntegralFamilyError::InvalidBasis(error.to_string()))?;
        let solution = matrix
            .solve(&rhs)
            .map_err(|error| IntegralFamilyError::InvalidBasis(error.to_string()))?;
        Ok(self
            .scalar_products
            .iter()
            .cloned()
            .zip(
                solution
                    .into_vec()
                    .into_iter()
                    .map(|value| value.to_expression()),
            )
            .collect())
    }

    fn validate_labels(&self, labels: &[Atom]) -> Result<(), IntegralFamilyError> {
        if labels.len() != self.denominators.len()
            || labels.iter().collect::<BTreeSet<_>>().len() != labels.len()
            || labels.iter().any(|label| {
                !matches!(label.as_view(), AtomView::Var(_) | AtomView::Fun(_))
                    || self
                        .scalar_products
                        .iter()
                        .any(|p| label.contains(p.as_view()))
                    || self
                        .denominators
                        .iter()
                        .any(|d| d.contains(label.as_view()))
            })
        {
            return Err(IntegralFamilyError::InvalidLabels);
        }
        Ok(())
    }

    /// Rewrite a numerator into inverse-propagator variables without integration.
    pub fn rewrite_numerator(
        &self,
        numerator: &Atom,
        labels: &[Atom],
    ) -> Result<Atom, IntegralFamilyError> {
        let rules = self.scalar_product_rules(labels)?;
        Ok(numerator
            .replace_map(|view, _, output| {
                if let Some(value) = rules.get(&view.to_owned()) {
                    **output = value.clone();
                }
            })
            .expand())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse;

    #[test]
    fn bubble_numerator_is_solved_in_propagator_variables() {
        let k = parse!("family_bubble::k");
        let p = parse!("family_bubble::p");
        let s = parse!("s");
        let m1 = parse!("m1sq");
        let m2 = parse!("m2sq");
        let kin = Kinematics::new()
            .with_momenta([k.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, s.clone())
            .unwrap();
        let kk = kin.scalar_product(&k, &k).unwrap();
        let kp = kin.scalar_product(&k, &p).unwrap();
        let shifted = &k + &p;
        let denominators = vec![
            &kk - &m1,
            kin.scalar_product(&shifted, &shifted).unwrap() - &m2,
        ];
        let family = IntegralFamily::new(vec![k], vec![p], denominators.clone(), &kin).unwrap();
        assert!(family.is_complete());
        assert!(family.is_independent());
        let labels = [parse!("d(1)"), parse!("d(2)")];
        let rules = family.scalar_product_rules(&labels).unwrap();
        assert!((&rules[&kk] - &labels[0] - &m1).expand().is_zero());
        let expected_kp = (&labels[1] - &labels[0] - s + m2 - m1) / 2;
        assert!((&rules[&kp] - &expected_kp).expand().is_zero());
        let reduced = family.rewrite_numerator(&kp.pow(2), &labels).unwrap();
        assert!((reduced - expected_kp.pow(2)).expand().is_zero());
        for product in family.scalar_products() {
            let reconstructed = rules[product].replace_map(|view, _, out| {
                if let Some(i) = labels.iter().position(|label| label.as_view() == view) {
                    **out = denominators[i].clone();
                }
            });
            assert!((reconstructed - product).expand().is_zero());
        }
    }

    #[test]
    fn incomplete_two_loop_family_gets_one_irreducible_product() {
        let [k, q, p] = [
            parse!("family_completion::k"),
            parse!("family_completion::q"),
            parse!("family_completion::p"),
        ];
        let kin = Kinematics::new()
            .with_momenta([k.clone(), q.clone(), p.clone()])
            .unwrap();
        let denominators = [&k, &q, &(&k - &p), &(&q - &p)]
            .map(|v| kin.scalar_product(v, v).unwrap())
            .to_vec();
        let family = IntegralFamily::new(
            vec![k.clone(), q.clone()],
            vec![p],
            denominators.clone(),
            &kin,
        )
        .unwrap();
        assert_eq!(family.rank(), 4);
        assert!(!family.is_complete());
        let completed = family.complete().unwrap();
        assert!(completed.is_complete());
        assert_eq!(&completed.denominators()[..4], denominators);
        assert_eq!(
            completed.denominators()[4],
            kin.scalar_product(&k, &q).unwrap()
        );
        assert_eq!(
            completed.complete().unwrap().denominators(),
            completed.denominators()
        );
    }

    #[test]
    fn dependence_eikonal_denominators_and_invalid_inputs() {
        let k = parse!("family_dependence::k");
        let p = parse!("family_dependence::p");
        let kin = Kinematics::new();
        let kk = kin.scalar_product(&k, &k).unwrap();
        let kp = kin.scalar_product(&k, &p).unwrap();
        let dependent = IntegralFamily::new(
            vec![k.clone()],
            vec![],
            vec![kk.clone(), &kk - parse!("m2")],
            &kin,
        )
        .unwrap();
        assert!(dependent.is_complete());
        assert!(!dependent.is_independent());
        assert!(matches!(
            dependent.complete(),
            Err(IntegralFamilyError::Dependent { .. })
        ));
        let eikonal = IntegralFamily::new(
            vec![k.clone()],
            vec![p.clone()],
            vec![kk.clone(), &kp + parse!("delta")],
            &kin,
        )
        .unwrap();
        assert!(eikonal.is_complete());
        assert!(
            eikonal
                .scalar_product_rules(&[parse!("d"), parse!("d")])
                .is_err()
        );
        for nonlinear in [kp.clone().pow(2), kp.clone().pow(-1)] {
            assert!(
                IntegralFamily::new(vec![k.clone()], vec![p.clone()], vec![nonlinear], &kin)
                    .is_err()
            );
        }
        let on_shell_loop = kin.clone().with_mass_squared(&k, Atom::Zero).unwrap();
        assert!(matches!(
            IntegralFamily::new(vec![k.clone()], vec![], vec![kk], &on_shell_loop),
            Err(IntegralFamilyError::LoopAssumption(_))
        ));
        assert!(IntegralFamily::new(vec![k.clone()], vec![k], vec![], &kin).is_err());
    }

    #[test]
    fn undeclared_external_directions_and_colliding_labels_are_rejected() {
        let [k, p, r] = [
            parse!("family_invalid::k"),
            parse!("family_invalid::p"),
            parse!("family_invalid::r"),
        ];
        let kin = Kinematics::new();
        let kk = kin.scalar_product(&k, &k).unwrap();
        let kr = kin.scalar_product(&k, &r).unwrap();
        assert!(IntegralFamily::new(vec![k.clone()], vec![p], vec![kk.clone(), kr], &kin).is_err());
        let m2 = parse!("family_invalid::m2");
        let family = IntegralFamily::new(vec![k], vec![], vec![kk - &m2], &kin).unwrap();
        assert!(matches!(
            family.scalar_product_rules(&[m2]),
            Err(IntegralFamilyError::InvalidLabels)
        ));
    }
}
