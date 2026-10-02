//! Scoped scalar-product assumptions; no global Symbolica state is mutated.

use spenso::{
    network::library::symbolic::ETS,
    structure::{
        dimension::{Dimension, DimensionError},
        representation::{Minkowski, RepName},
    },
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
};
use thiserror::Error;

#[derive(Clone, Debug, Error)]
pub enum SymbolicKinematicsError {
    #[error("expected an unindexed momentum symbol or labeled call, got {0}")]
    InvalidMomentum(Atom),
    #[error("{0} is not linear in the declared momenta; declare momenta before combining them")]
    NonlinearMomentum(Atom),
    #[error("phase-space densities require four-dimensional kinematics")]
    PhaseSpaceDimension,
}

/// Symbolic Lorentz kinematics in Spenso's scalar-product notation.
///
/// Momentum names are unindexed symbols or labeled calls such as `Q(1)`.
/// Assumptions belong to this value and can be reused across expressions.
/// Tensor and Dirac contractions remain the responsibility of Idenso; apply
/// these substitutions after converting contractions to compact scalar products.
#[derive(Clone, Debug)]
pub struct Kinematics {
    dimension: Dimension,
    products: BTreeMap<Atom, Atom>,
    momenta: BTreeSet<Atom>,
}

impl Default for Kinematics {
    fn default() -> Self {
        Self {
            dimension: 4.into(),
            products: BTreeMap::new(),
            momenta: BTreeSet::new(),
        }
    }
}

impl Kinematics {
    /// Construct an empty set of four-dimensional scalar-product assumptions.
    pub fn new() -> Self {
        Self::default()
    }

    /// Lorentz dimension carried by this context's scalar products.
    pub fn dimension(&self) -> Dimension {
        self.dimension
    }

    /// Start with no assumptions in a symbolic or nonnegative integer dimension.
    ///
    /// Use a symbol such as `D` during dimensional regularization, then
    /// substitute `D = 4 - 2 eps` in the contracted result with Symbolica.
    pub fn in_dimension(dimension: &Atom) -> Result<Self, DimensionError> {
        Ok(Self {
            dimension: Dimension::try_from(dimension.as_view())?,
            products: BTreeMap::new(),
            momenta: BTreeSet::new(),
        })
    }

    /// Declare momentum names for linear combinations with symbolic coefficients.
    ///
    /// Setting a scalar product also declares both momentum names. Symbols not
    /// declared as momenta are scalar coefficients inside a linear combination.
    pub fn with_momenta(
        mut self,
        momenta: impl IntoIterator<Item = Atom>,
    ) -> Result<Self, SymbolicKinematicsError> {
        for momentum in momenta {
            if !matches!(momentum.as_view(), AtomView::Var(_) | AtomView::Fun(_)) {
                return Err(SymbolicKinematicsError::InvalidMomentum(momentum));
            }
            self.momenta.insert(momentum);
        }
        Ok(self)
    }

    /// Set one scalar product. Reversed arguments use the same symmetric key.
    pub fn with_scalar_product(
        self,
        left: &Atom,
        right: &Atom,
        value: Atom,
    ) -> Result<Self, SymbolicKinematicsError> {
        let mut result = self.with_momenta([left.clone(), right.clone()])?;
        let rep = Minkowski {}.new_rep(result.dimension);
        result
            .products
            .insert(rep.inner_product(left, right), value.clone());
        // Dirac simplification may retain metric shorthand instead of dots.
        result.products.insert(
            function!(
                ETS.metric,
                rep.vector(left.as_view(), []),
                rep.vector(right.as_view(), [])
            ),
            value,
        );
        Ok(result)
    }

    /// Set the on-shell relation `p.p = mass_squared`.
    pub fn with_mass_squared(
        self,
        momentum: &Atom,
        mass_squared: Atom,
    ) -> Result<Self, SymbolicKinematicsError> {
        self.with_scalar_product(momentum, momentum, mass_squared)
    }

    /// Initial-state denominator for scattering, or a decay in its rest frame.
    ///
    /// Two momenta give `4 sqrt((p1.p2)^2 - p1^2 p2^2)`. One momentum gives
    /// `2 sqrt(p^2)` for a rest-frame decay. Divide the squared matrix element
    /// times the phase-space measure by this value. Physical momenta are future
    /// directed and on shell; square roots retain Symbolica's branch semantics.
    pub fn flux(
        &self,
        first: &Atom,
        second: Option<&Atom>,
    ) -> Result<Atom, SymbolicKinematicsError> {
        let flux = match second {
            Some(second) => crate::InitialStateFlux::Scattering {
                momentum_dot: self.scalar_product(first, second)?,
                mass_squared_product: self.scalar_product(first, first)?
                    * self.scalar_product(second, second)?,
            },
            None => crate::InitialStateFlux::Decay {
                energy: self.scalar_product(first, first)?.sqrt(),
            },
        };
        Ok(flux.denominator(|value| value.sqrt()))
    }

    /// Four-dimensional `dPhi_2 / dOmega` for two outgoing on-shell momenta.
    ///
    /// The measure includes `(2 pi)^4 delta^4(P-p1-p2)` and each
    /// `d^3p / ((2 pi)^3 2E)`. The solid angle is measured in the pair rest
    /// frame. This is the physical above-threshold branch, without an
    /// identical-particle factor or initial-state flux. Symbolica retains
    /// roots unless positivity assumptions justify their simplification.
    pub fn two_body_phase_space(
        &self,
        first: &Atom,
        second: &Atom,
    ) -> Result<Atom, SymbolicKinematicsError> {
        if self.dimension != Dimension::from(4) {
            return Err(SymbolicKinematicsError::PhaseSpaceDimension);
        }
        let invariant = self.scalar_product(first, first)?
            + Atom::num(2) * self.scalar_product(first, second)?
            + self.scalar_product(second, second)?;
        Ok(self.flux(first, Some(second))?
            / (Atom::num(64) * Atom::var(symbolica::atom::Symbol::PI).pow(2) * invariant))
    }

    /// Four-dimensional `dPhi_3 / (ds12 ds23)`, with overall orientation integrated.
    ///
    /// Here `sij=(pi+pj)^2`. The normalization includes `(2 pi)^4 delta^4`
    /// and `d^3p / ((2 pi)^3 2E)` for each final particle. The density is
    /// `1/(128 pi^3 P^2)`, where `P=p1+p2+p3`, inside the physical Dalitz
    /// region. Masses determine that region, not the density. This method does
    /// not impose its boundaries or include flux or identical-particle factors.
    /// Use an orientation-independent or orientation-averaged squared amplitude.
    pub fn three_body_phase_space(
        &self,
        first: &Atom,
        second: &Atom,
        third: &Atom,
    ) -> Result<Atom, SymbolicKinematicsError> {
        if self.dimension != Dimension::from(4) {
            return Err(SymbolicKinematicsError::PhaseSpaceDimension);
        }
        let momenta = [first, second, third];
        let mut invariant = Atom::Zero;
        for (i, p) in momenta.iter().enumerate() {
            invariant += self.scalar_product(p, p)?;
            for q in &momenta[i + 1..] {
                invariant += Atom::num(2) * self.scalar_product(p, q)?;
            }
        }
        Ok(Atom::one()
            / (Atom::num(128) * Atom::var(symbolica::atom::Symbol::PI).pow(3) * invariant))
    }

    /// Scalar products for `p1 + p2 -> p3 + p4`.
    ///
    /// Mass arguments are squared masses. The invariants obey
    /// `s=(p1+p2)^2`, `t=(p1-p3)^2`, `u=(p1-p4)^2` and
    /// `s+t+u=sum(mass_squared)`. This constructor does not eliminate an
    /// invariant; use ordinary Symbolica substitution to choose which to remove.
    pub fn mandelstam(
        momenta: [&Atom; 4],
        mass_squared: [Atom; 4],
        invariants: [Atom; 3],
    ) -> Result<Self, SymbolicKinematicsError> {
        let [s, t, u] = invariants;
        let mut result = Self::new();
        for i in 0..4 {
            result = result.with_mass_squared(momenta[i], mass_squared[i].clone())?;
        }
        for (i, j, invariant, sign) in [
            (0, 1, &s, 1),
            (2, 3, &s, 1),
            (0, 2, &t, -1),
            (1, 3, &t, -1),
            (0, 3, &u, -1),
            (1, 2, &u, -1),
        ] {
            let value = Atom::num(sign) * (invariant - &mass_squared[i] - &mass_squared[j]) / 2;
            result = result.with_scalar_product(momenta[i], momenta[j], value.expand())?;
        }
        Ok(result)
    }

    /// Expand a bilinear scalar product and apply the known assumptions.
    ///
    /// Single momentum names need no declaration. For sums and products, declare
    /// the momentum names first; all remaining symbols are scalar coefficients.
    pub fn scalar_product(
        &self,
        left: &Atom,
        right: &Atom,
    ) -> Result<Atom, SymbolicKinematicsError> {
        let left = self.linear_terms(left)?;
        let right = self.linear_terms(right)?;
        let rep = Minkowski {}.new_rep(self.dimension);
        let mut result = Atom::Zero;
        for (p, a) in &left {
            for (q, b) in &right {
                let dot = rep.inner_product(p, q);
                result += a * b * self.products.get(&dot).cloned().unwrap_or(dot);
            }
        }
        Ok(result.expand())
    }

    fn linear_terms(&self, momentum: &Atom) -> Result<Vec<(Atom, Atom)>, SymbolicKinematicsError> {
        if momentum.is_zero() {
            return Ok(Vec::new());
        }
        if matches!(momentum.as_view(), AtomView::Var(_) | AtomView::Fun(_)) {
            return Ok(vec![(momentum.clone(), Atom::one())]);
        }
        let variables = self.momenta.iter().collect::<Vec<_>>();
        let terms = momentum.coefficient_list::<i32>(&variables);
        if terms.iter().any(|(key, coefficient)| {
            !self.momenta.contains(key)
                || variables.iter().any(|p| coefficient.contains(p.as_view()))
        }) {
            return Err(SymbolicKinematicsError::NonlinearMomentum(momentum.clone()));
        }
        Ok(terms)
    }

    /// Substitute known dot products and compact metric contractions.
    ///
    /// Substitutions are simultaneous: a value containing another assumed
    /// scalar product is retained verbatim. This avoids recursive rewrite cycles.
    pub fn apply(&self, expression: &Atom) -> Atom {
        // Exact atoms, rather than wildcard patterns, keep assumptions local to
        // their momenta and preserve unknown scalar products.
        expression.replace_map(|view, _, output| {
            if let Some(value) = self.products.get(&view.to_owned()) {
                **output = value.clone();
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse;

    #[test]
    fn two_body_measure_has_physical_normalization_and_mass_symmetry() {
        let p = parse!("p");
        let q = parse!("q");
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let kin = Kinematics::new()
            .with_mass_squared(&p, Atom::num(9))
            .unwrap()
            .with_mass_squared(&q, Atom::num(25))
            .unwrap()
            .with_scalar_product(&p, &q, Atom::num(113))
            .unwrap();
        assert_eq!(kin.flux(&p, Some(&q)).unwrap(), Atom::num(448));
        assert_eq!(kin.flux(&p, None).unwrap(), Atom::num(6));
        let measure = kin.two_body_phase_space(&p, &q).unwrap();
        assert!(
            (measure.clone() * Atom::num(260) * pi.clone().pow(2) - Atom::num(7))
                .expand()
                .is_zero()
        );
        assert_eq!(measure, kin.two_body_phase_space(&q, &p).unwrap());
        let threshold = kin.with_scalar_product(&p, &q, Atom::num(15)).unwrap();
        assert!(threshold.two_body_phase_space(&p, &q).unwrap().is_zero());
        assert!(
            Kinematics::in_dimension(&parse!("D"))
                .unwrap()
                .two_body_phase_space(&p, &q)
                .is_err()
        );
        // Momentum names also work without assumptions or prior declarations.
        assert!(Kinematics::new().two_body_phase_space(&p, &q).is_ok());
    }

    #[test]
    fn dalitz_density_has_massless_volume_and_permutation_symmetry() {
        let [p, q, r, s, s12, s23] = [
            parse!("p"),
            parse!("q"),
            parse!("r"),
            parse!("s"),
            parse!("s12"),
            parse!("s23"),
        ];
        let mut kin = Kinematics::new();
        for momentum in [&p, &q, &r] {
            kin = kin.with_mass_squared(momentum, Atom::Zero).unwrap();
        }
        kin = kin
            .with_scalar_product(&p, &q, &s12 / 2)
            .unwrap()
            .with_scalar_product(&q, &r, &s23 / 2)
            .unwrap()
            .with_scalar_product(&p, &r, (&s - &s12 - &s23) / 2)
            .unwrap();
        let density = kin.three_body_phase_space(&p, &q, &r).unwrap();
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        // The massless Dalitz triangle has area s^2/2. Its integrated volume
        // also follows from two sequential two-body phase spaces.
        assert!(
            (density.clone() * s.clone().pow(2) / 2 - &s / (Atom::num(256) * pi.clone().pow(3)))
                .together()
                .is_zero()
        );
        assert_eq!(density, kin.three_body_phase_space(&r, &p, &q).unwrap());
        let massive = kin.clone().with_mass_squared(&p, Atom::num(9)).unwrap();
        assert!(
            (massive.three_body_phase_space(&p, &q, &r).unwrap()
                * (s + Atom::num(9))
                * Atom::num(128)
                * pi.pow(3)
                - Atom::one())
            .together()
            .is_zero()
        );
        assert!(matches!(
            Kinematics::in_dimension(&parse!("D"))
                .unwrap()
                .three_body_phase_space(&p, &q, &r),
            Err(SymbolicKinematicsError::PhaseSpaceDimension)
        ));
        assert!(kin.three_body_phase_space(&(&p * &q), &q, &r).is_err());
    }

    #[test]
    fn bilinear_momenta_use_symbolica_coefficients() {
        let p = parse!("Q(1)");
        let q = parse!("Q(2)");
        let a = parse!("a");
        let b = parse!("b");
        let kin = Kinematics::new()
            .with_scalar_product(&p, &p, parse!("m2"))
            .unwrap()
            .with_scalar_product(&q, &q, Atom::Zero)
            .unwrap()
            .with_scalar_product(&p, &q, parse!("s/2"))
            .unwrap();
        let result = kin
            .scalar_product(&(&a * &p + &b * &q), &(&p - &q))
            .unwrap();
        let expected = parse!("a*m2 + (b-a)*s/2");
        assert!((result - expected).expand().is_zero());
        assert!(kin.scalar_product(&Atom::Zero, &p).unwrap().is_zero());
        assert!(kin.scalar_product(&(&p * &q), &p).is_err());
        assert!(kin.scalar_product(&parse!("Q(1)/(1+Q(2))"), &p).is_err());
        assert!(kin.scalar_product(&(&p + Atom::one()), &p).is_err());
        assert!(Kinematics::new().with_momenta([&p + &q]).is_err());
    }

    #[test]
    fn declarations_expand_unknown_products_without_global_assumptions() {
        let p = parse!("p");
        let q = parse!("q");
        let kin = Kinematics::new()
            .with_momenta([p.clone(), q.clone()])
            .unwrap();
        let difference = kin.scalar_product(&(&p + &q), &(&p - &q)).unwrap();
        let expected = kin.scalar_product(&p, &p).unwrap() - kin.scalar_product(&q, &q).unwrap();
        assert_eq!(difference, expected);
        assert!(Kinematics::new().scalar_product(&(&p + &q), &p).is_err());
    }

    #[test]
    fn unequal_mass_mandelstam_respects_momentum_conservation() {
        let p = [parse!("p1"), parse!("p2"), parse!("p3"), parse!("p4")];
        let masses = [
            parse!("m1sq"),
            parse!("m2sq"),
            parse!("m3sq"),
            parse!("m4sq"),
        ];
        let s = parse!("s");
        let t = parse!("t");
        let u = masses.iter().cloned().sum::<Atom>() - &s - &t;
        let kin = Kinematics::mandelstam(
            p.each_ref(),
            masses.clone(),
            [s.clone(), t.clone(), u.clone()],
        )
        .unwrap();
        for i in 0..4 {
            assert_eq!(kin.scalar_product(&p[i], &p[i]).unwrap(), masses[i]);
            let conserved = kin.scalar_product(&p[i], &p[0]).unwrap()
                + kin.scalar_product(&p[i], &p[1]).unwrap()
                - kin.scalar_product(&p[i], &p[2]).unwrap()
                - kin.scalar_product(&p[i], &p[3]).unwrap();
            assert!(conserved.expand().is_zero(), "{conserved}");
        }
        assert!(
            (kin.scalar_product(&p[0], &p[1]).unwrap() * 2 + &masses[0] + &masses[1] - s)
                .expand()
                .is_zero()
        );
        assert!(
            (&masses[0] + &masses[2] - kin.scalar_product(&p[0], &p[2]).unwrap() * 2 - t)
                .expand()
                .is_zero()
        );
        assert!(
            (&masses[0] + &masses[3] - kin.scalar_product(&p[0], &p[3]).unwrap() * 2 - u)
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn substitutions_are_scoped_symmetric_and_simultaneous() {
        let rep = Minkowski {}.new_rep(4);
        let p = parse!("p");
        let q = parse!("q");
        let r = parse!("r");
        let pq = rep.inner_product(&p, &q);
        let qr = rep.inner_product(&q, &r);
        let kin = Kinematics::new()
            .with_scalar_product(&p, &q, qr.clone())
            .unwrap()
            .with_scalar_product(&q, &r, Atom::num(7))
            .unwrap();
        assert_eq!(kin.scalar_product(&q, &p).unwrap(), qr);
        assert_eq!(kin.apply(&pq), qr);
        assert_eq!(Kinematics::new().apply(&pq), pq);
        let shorthand = function!(
            ETS.metric,
            rep.vector(p.as_view(), []),
            rep.vector(q.as_view(), [])
        );
        assert_eq!(kin.apply(&shorthand), qr);
        let unknown = rep.inner_product(&p, &r);
        assert_eq!(kin.apply(&unknown), unknown);
    }

    #[test]
    fn symbolic_dimensions_do_not_mix_with_four_dimensions() {
        let d = parse!("D");
        let p = parse!("p");
        let kin = Kinematics::in_dimension(&d)
            .unwrap()
            .with_mass_squared(&p, Atom::num(7))
            .unwrap();
        assert_eq!(kin.scalar_product(&p, &p).unwrap(), Atom::num(7));
        let four_dimensional = Kinematics::new().scalar_product(&p, &p).unwrap();
        assert_eq!(kin.apply(&four_dimensional), four_dimensional);
        assert!(Kinematics::in_dimension(&parse!("4-2*eps")).is_err());
    }
}
