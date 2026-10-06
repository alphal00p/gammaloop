//! Scalar-product bases for loop-integral families.
//!
//! Propagators are affine forms in the independent loop scalar products.
//! Symbolica owns the polynomial conversion, rank calculation and exact solve.

use std::collections::{BTreeMap, BTreeSet};

use feynkit_kinematics::{Kinematics, SymbolicKinematicsError};
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::{
    algebraic::AlgebraicExtension,
    rational::{Q, RationalField},
    rational_polynomial::{RationalPolynomial, RationalPolynomialField},
};
use symbolica::poly::{PolyVariable, polynomial::MultivariatePolynomial};
use symbolica::tensors::matrix::Matrix;
use thiserror::Error;

mod diagram;
mod mapping;
mod parametric;
mod partial_fraction;
pub use mapping::{IntegralMapping, PropagatorMapping, QuadraticMomentum};

type GaussianField = AlgebraicExtension<RationalField>;
type FamilyCoefficient = RationalPolynomial<GaussianField, u16>;
type FamilyField = RationalPolynomialField<GaussianField, u16>;
type CoefficientMatrix = Matrix<FamilyField>;
type ExactPolynomial = MultivariatePolynomial<FamilyField, u16>;

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
/// information beyond this algebraic representation. Coefficients are exact
/// rational functions over the native Gaussian-rational field, with `i^2 = -1`.
/// Algebraic coefficient support does not authorize complex contour changes.
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
    /// Group over the exact native Gaussian-rational field, retaining the original row
    /// denominator. Default expression-field conversion uses statistical zero
    /// tests, which can discard tiny symbolic coefficients after expansion.
    fn exact_polynomial(
        expression: &impl AtomCore,
        variables: &[Atom],
    ) -> Result<ExactPolynomial, IntegralFamilyError> {
        let variables = variables
            .iter()
            .cloned()
            .map(PolyVariable::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(IntegralFamilyError::InvalidBasis)?;
        let field = GaussianField::complex(Q);
        expression
            .try_to_rational_polynomial(&field, &field, None)
            .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?
            .to_polynomial(&variables, false)
            .map_err(|e| IntegralFamilyError::InvalidBasis(e.into()))
    }

    fn affine_coefficients(
        expression: &impl AtomCore,
        variables: &[Atom],
    ) -> Result<Vec<FamilyCoefficient>, IntegralFamilyError> {
        let polynomial = Self::exact_polynomial(expression, variables)?;
        let field = GaussianField::complex(Q);
        let mut row =
            vec![RationalPolynomial::new(&field, Default::default()); variables.len() + 1];
        for term in &polynomial {
            let mut powers = term.exponents.iter().enumerate().filter(|(_, p)| **p != 0);
            let column = match (powers.next(), powers.next()) {
                (None, None) => variables.len(),
                (Some((index, 1)), None) => index,
                _ => {
                    return Err(IntegralFamilyError::InvalidBasis(
                        "expected affine expressions".into(),
                    ));
                }
            };
            row[column] = term.coefficient.clone();
        }
        Ok(row)
    }

    fn scalar_product_coefficients(
        &self,
        expression: &Atom,
    ) -> Result<BTreeMap<Atom, Atom>, IntegralFamilyError> {
        self.scalar_products
            .iter()
            .cloned()
            .chain(std::iter::once(Atom::one()))
            .zip(Self::affine_coefficients(
                expression,
                &self.scalar_products,
            )?)
            .map(|(key, coefficient)| {
                let coefficient = coefficient.to_expression();
                self.require_external_coefficient(&coefficient)?;
                Ok((key, coefficient))
            })
            .collect()
    }

    /// Preserve affine coefficients, including rational row normalization.
    /// Symbolica's equation-to-matrix conversion may clear denominators, which
    /// preserves solutions and ranks but changes determinants and relations
    /// between the original expressions. Extract those coefficients explicitly;
    /// Symbolica still owns all polynomial and matrix operations.
    fn affine_system(
        expressions: &[impl AtomCore],
        variables: &[Atom],
    ) -> Result<(CoefficientMatrix, CoefficientMatrix), IntegralFamilyError> {
        let mut entries: Vec<FamilyCoefficient> = Vec::new();
        for expression in expressions {
            entries.extend(Self::affine_coefficients(expression, variables)?);
        }
        // Every entry must use the same variable ordering for exact arithmetic.
        if let Some((first, rest)) = entries.split_first_mut() {
            for _ in 0..2 {
                for entry in &mut *rest {
                    first.unify_variables(entry);
                }
            }
        }
        let field = RationalPolynomialField::new(GaussianField::complex(Q));
        let mut matrix = Vec::new();
        let mut rhs = Vec::new();
        for row in entries.chunks(variables.len() + 1) {
            matrix.extend_from_slice(&row[..variables.len()]);
            rhs.push(-row[variables.len()].clone());
        }
        Ok((
            Matrix::from_linear(
                matrix,
                expressions.len() as u32,
                variables.len() as u32,
                field.clone(),
            )
            .map_err(IntegralFamilyError::InvalidBasis)?,
            Matrix::new_vec(rhs, field),
        ))
    }

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

    fn require_external_coefficient(&self, coefficient: &Atom) -> Result<(), IntegralFamilyError> {
        let rep = Minkowski {}.new_rep(self.kinematics.dimension());
        if self
            .scalar_products
            .iter()
            .chain(&self.loop_momenta)
            .any(|p| coefficient.contains(p.as_view()))
            || self
                .loop_momenta
                .iter()
                .any(|p| coefficient.contains(rep.vector(p.as_view(), []).as_view()))
        {
            return Err(IntegralFamilyError::InvalidBasis(format!(
                "coefficient depends on loop momenta: {coefficient}"
            )));
        }
        Ok(())
    }

    fn rank_of(&self, denominators: &[Atom]) -> Result<usize, IntegralFamilyError> {
        if denominators.is_empty() {
            return Ok(0);
        }
        let (matrix, constants) = Self::affine_system(denominators, &self.scalar_products)?;
        // Native rational-polynomial conversion treats functions and powers as
        // external atoms. Reject any hidden loop dependence in those atoms too.
        for coefficient in matrix
            .row_iter()
            .flatten()
            .chain(constants.row_iter().flatten())
        {
            self.require_external_coefficient(&coefficient.to_expression())?;
        }
        Ok(matrix.rank())
    }

    /// Append independent candidates, then scalar products, to form a basis.
    ///
    /// Original propagators retain their positions. Added entries are auxiliary
    /// inverse propagators, whose nonpositive powers represent numerator factors.
    /// A dependent family must be partial-fractioned first.
    /// Candidates are affine inverse propagators in this family's momenta,
    /// tried in order after applying its kinematics. Redundant candidates are
    /// skipped; bare scalar products fill any directions still missing.
    pub fn complete(&self, candidates: &[Atom]) -> Result<Self, IntegralFamilyError> {
        self.require_independent()?;
        let candidates = candidates
            .iter()
            .map(|candidate| self.kinematics.apply(candidate).expand())
            .collect::<Vec<_>>();
        // Validate the entire supplied pool, including entries beyond the first
        // complete basis, so invalid candidates are never silently accepted.
        if !candidates.is_empty() {
            self.rank_of(&candidates)?;
        }
        let mut result = self.clone();
        for product in candidates.iter().chain(&self.scalar_products) {
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
        let (matrix, rhs) = Self::affine_system(&equations, &self.scalar_products)?;
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
        let completed = family.complete(&[]).unwrap();
        assert!(completed.is_complete());
        assert_eq!(&completed.denominators()[..4], denominators);
        assert_eq!(
            completed.denominators()[4],
            kin.scalar_product(&k, &q).unwrap()
        );
        assert_eq!(
            completed.complete(&[]).unwrap().denominators(),
            completed.denominators()
        );
        let preferred = kin.scalar_product(&(&k - &q), &(&k - &q)).unwrap();
        let pool = [denominators[0].clone(), preferred.clone()];
        let selected = family.complete(&pool).unwrap();
        assert_eq!(&selected.denominators()[..4], denominators);
        assert_eq!(selected.denominators()[4], preferred);
        assert_eq!(
            selected.complete(&pool).unwrap().denominators(),
            selected.denominators()
        );
        // Complete with the physical mixed propagator rather than a bare k.q.
        let labels = [
            parse!("c1"),
            parse!("c2"),
            parse!("c3"),
            parse!("c4"),
            parse!("c5"),
        ];
        let reduced = selected
            .rewrite_numerator(&kin.scalar_product(&k, &q).unwrap(), &labels)
            .unwrap();
        assert!(
            (reduced - (&labels[0] + &labels[1] - &labels[4]) / 2)
                .expand()
                .is_zero()
        );
        assert_eq!(
            family
                .complete(&[denominators[0].clone()])
                .unwrap()
                .denominators(),
            completed.denominators()
        );
        let nonlinear = kin.scalar_product(&k, &q).unwrap().pow(2);
        assert!(family.complete(&[preferred, nonlinear.clone()]).is_err());
        assert!(selected.complete(&[nonlinear]).is_err());
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
            dependent.complete(&[]),
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

    #[test]
    fn gaussian_affine_rank_row_scales_and_rules_are_exact() {
        let k = parse!("gaussian_affine::k");
        let p = parse!("gaussian_affine::p");
        let kin = Kinematics::new();
        let x = kin.scalar_product(&k, &k).unwrap();
        let y = kin.scalar_product(&k, &p).unwrap();
        let imaginary = parse!("𝑖");
        let dependent = IntegralFamily::new(
            vec![k.clone()],
            vec![p.clone()],
            vec![&x + &imaginary * &y, &imaginary * &x - &y],
            &kin,
        )
        .unwrap();
        assert_eq!(
            dependent.rank(),
            1,
            "the coefficient field must enforce i squared = -1"
        );
        let scale = parse!("(gaussian_affine::a+𝑖*gaussian_affine::b)/10^1000");
        let denominators: [Atom; 2] = [
            &scale * (&x + (Atom::one() + &imaginary) * &y / 3 - (Atom::num(2) - &imaginary) / 5),
            (&imaginary * &x + &y + 7) / (Atom::num(3) - &imaginary),
        ];
        for expanded in [false, true] {
            let input = denominators
                .iter()
                .map(|d| if expanded { d.expand() } else { d.clone() })
                .collect::<Vec<_>>();
            let family =
                IntegralFamily::new(vec![k.clone()], vec![p.clone()], input.clone(), &kin).unwrap();
            let variables = [x.clone(), y.clone()];
            let (matrix, rhs) = IntegralFamily::affine_system(&input, &variables).unwrap();
            for (i, row) in matrix.row_iter().enumerate() {
                let reconstructed = row
                    .iter()
                    .zip(&variables)
                    .map(|(c, v)| c.to_expression() * v)
                    .sum::<Atom>()
                    - rhs[(i as u32, 0)].to_expression();
                assert!((reconstructed - &input[i]).together().is_zero());
            }
            assert!(
                (matrix[(0, 0)].to_expression() - &scale)
                    .together()
                    .is_zero()
            );
            let labels = [parse!("gaussian_affine::d1"), parse!("gaussian_affine::d2")];
            let rules = family.scalar_product_rules(&labels).unwrap();
            for variable in &variables {
                let reconstructed = rules[variable].replace_map(|view, _, output| {
                    if let Some(i) = labels.iter().position(|label| label.as_view() == view) {
                        **output = input[i].clone();
                    }
                });
                assert!((reconstructed - variable).together().is_zero());
            }
        }
    }

    #[test]
    fn exact_affine_extraction_preserves_original_rational_row_scales() {
        let k = parse!("family_exact_scale::k");
        let p = parse!("family_exact_scale::p");
        let kin = Kinematics::new();
        let x = kin.scalar_product(&k, &k).unwrap();
        let y = kin.scalar_product(&k, &p).unwrap();
        let q = parse!("(family_exact_scale::a+family_exact_scale::b)/10^1000");
        let expressions = [
            (&x + Atom::num(2) * &y - 7) * &q / 3,
            (Atom::num(3) * &x - &y + 5) / (Atom::num(11) * &q),
        ];
        let expected = [
            vec![&q / 3, Atom::num(2) * &q / 3, Atom::num(7) * &q / 3],
            vec![
                Atom::num(3) / (Atom::num(11) * &q),
                -Atom::one() / (Atom::num(11) * &q),
                -Atom::num(5) / (Atom::num(11) * &q),
            ],
        ];
        for expanded in [false, true] {
            let input = expressions
                .iter()
                .map(|d| if expanded { d.expand() } else { d.clone() })
                .collect::<Vec<_>>();
            let (matrix, rhs) =
                IntegralFamily::affine_system(&input, &[x.clone(), y.clone()]).unwrap();
            for (i, row) in matrix.row_iter().enumerate() {
                for (j, coefficient) in row.iter().enumerate() {
                    assert!(
                        (coefficient.to_expression() - &expected[i][j])
                            .together()
                            .is_zero()
                    );
                }
                assert!(
                    (rhs[(i as u32, 0)].to_expression() - &expected[i][2])
                        .together()
                        .is_zero()
                );
            }
            let family =
                IntegralFamily::new(vec![k.clone()], vec![p.clone()], input, &kin).unwrap();
            assert_eq!(family.rank(), 2);
            let labels = [
                parse!("family_exact_scale::d1"),
                parse!("family_exact_scale::d2"),
            ];
            for (expression, label) in expressions.iter().zip(&labels) {
                let rewritten = family.rewrite_numerator(expression, &labels).unwrap();
                assert!((rewritten - label).together().is_zero());
            }
        }
    }

    #[test]
    fn tiny_hidden_loop_dependence_and_nonlinear_terms_are_rejected_exactly() {
        let k = parse!("family_exact_guard::k");
        let p = parse!("family_exact_guard::p");
        let r = parse!("family_exact_guard::r");
        let kin = Kinematics::new();
        let x = kin.scalar_product(&k, &k).unwrap();
        let y = kin.scalar_product(&k, &p).unwrap();
        let hidden = symbolica::symbol!("family_exact_guard::hidden");
        let family = IntegralFamily::new(
            vec![k.clone()],
            vec![p.clone()],
            vec![x.clone(), y.clone()],
            &kin,
        )
        .unwrap();
        for scale in [
            parse!("(family_exact_guard::a+family_exact_guard::b)/10^1000"),
            parse!("(family_exact_guard::a+𝑖*family_exact_guard::b)/10^1000"),
            parse!("(family_exact_guard::a+family_exact_guard::b)*10^1000"),
        ] {
            for bad in [
                x.clone().pow(2),
                x.clone().pow(-1),
                hidden.call(&x),
                hidden.call(&k),
                kin.scalar_product(&k, &r).unwrap(),
            ] {
                for expanded in [false, true] {
                    let bad = &scale * &bad;
                    let bad = if expanded { bad.expand() } else { bad };
                    assert!(matches!(
                        IntegralFamily::new(
                            vec![k.clone()],
                            vec![p.clone()],
                            vec![&y + &bad],
                            &kin
                        ),
                        Err(IntegralFamilyError::InvalidBasis(_))
                    ));
                    assert!(matches!(
                        family.require_external_coefficient(&bad),
                        Err(IntegralFamilyError::InvalidBasis(_))
                    ));
                }
            }
            family.require_external_coefficient(&scale).unwrap();
        }
        let gaussian =
            IntegralFamily::new(vec![k], vec![p], vec![parse!("𝑖") * x - 1], &kin).unwrap();
        assert_eq!(gaussian.rank(), 1);
    }
}
