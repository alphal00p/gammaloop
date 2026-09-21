//! External-state completeness relations shared by FeynKit and GammaLoop.

use feynkit_model::{Model, Particle};
use idenso::{dirac::AGS, representations::Bispinor};
use spenso::network::library::symbolic::ETS;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use thiserror::Error;

/// Invalid or unsupported external-state polarization sum.
#[derive(Debug, Error)]
pub enum SpinSumError {
    #[error("polarization sums for UFO spin {0} are not implemented")]
    UnsupportedSpin(i64),
    #[error("expected UFO spin {expected}, got {actual}")]
    WrongSpin { expected: i64, actual: i64 },
    #[error("the axial reference must have nonzero momentum contraction")]
    OrthogonalReference,
    #[error("a momentum must be an unindexed symbol or function call")]
    InvalidMomentum,
}

/// An axial-gauge reference contracted with the two open Lorentz slots.
///
/// The reference need not be null. All contractions use the same Lorentz
/// metric as the momentum and open slots.
pub struct AxialReference {
    pub left: Atom,
    pub right: Atom,
    pub momentum_dot_reference: Atom,
    pub norm_squared: Atom,
}

/// Completeness relations for scalars, Dirac fermions, and vector particles.
///
/// Inputs and outputs are ordinary Symbolica expressions in Spenso notation.
/// Kinematics, on-shell assumptions, and the choice of axial reference stay
/// with the caller. A covariant massless vector sum assumes a gauge-invariant
/// amplitude; a physical sum requires an explicit reference.
pub struct SpinSum {
    spin: i64,
    mass: Atom,
    massless: bool,
    antiparticle: bool,
    average: bool,
    covariant: bool,
}

impl SpinSum {
    pub fn new(particle: &Particle, model: &Model) -> Result<Self, SpinSumError> {
        if !(1..=3).contains(&particle.spin) {
            return Err(SpinSumError::UnsupportedSpin(particle.spin));
        }
        let massless = particle.is_massless(model);
        Ok(Self {
            spin: particle.spin,
            mass: particle.symbolic_mass(model)
                .replace(symbolica::symbol!("UFO::ZERO")).with(Atom::Zero),
            massless,
            antiparticle: particle.is_antiparticle(),
            average: false,
            covariant: false,
        })
    }

    /// Divide by the number of physical spin states in four dimensions.
    pub fn averaged(mut self, average: bool) -> Self {
        self.average = average;
        self
    }

    /// Use `-g` for vector sums, including massive vectors in Feynman gauge.
    pub fn covariant(mut self, covariant: bool) -> Self {
        self.covariant = covariant;
        self
    }

    /// Build a four-dimensional completeness tensor from unindexed momenta.
    ///
    /// The two indices are bare Symbolica indices, not Spenso slots. Fermions
    /// use `bis(4, index)` and vectors use `mink(4, index)`. A momentum may be
    /// a symbol (`p`) or a labeled function (`Q(1)`). For massless vectors,
    /// omitting the reference selects the covariant sum.
    pub fn expression(
        &self,
        momentum: &Atom,
        indices: [Atom; 2],
        reference: Option<&Atom>,
    ) -> Result<Atom, SpinSumError> {
        use symbolica::atom::AtomView;
        if !matches!(momentum.as_view(), AtomView::Var(_) | AtomView::Fun(_))
            || reference
                .is_some_and(|r| !matches!(r.as_view(), AtomView::Var(_) | AtomView::Fun(_)))
        {
            return Err(SpinSumError::InvalidMomentum);
        }
        let lorentz = Minkowski {}.new_rep(4);
        match self.spin {
            1 => Ok(Atom::one()),
            2 => {
                let spinor = Bispinor {}.new_rep(4);
                let [left, right] = indices.map(|index| spinor.pattern(index));
                let slash = function!(
                    AGS.gamma,
                    &left,
                    &right,
                    lorentz.vector(momentum.as_view(), [])
                );
                self.fermion(slash, &left, &right)
            }
            3 => {
                let components = indices
                    .clone()
                    .map(|index| lorentz.vector(momentum.as_view(), [index]));
                let reference = reference.map(|reference| AxialReference {
                    left: lorentz.vector(reference.as_view(), [indices[0].clone()]),
                    right: lorentz.vector(reference.as_view(), [indices[1].clone()]),
                    momentum_dot_reference: lorentz.inner_product(momentum, reference),
                    norm_squared: lorentz.inner_product(reference, reference),
                });
                let [left, right] = indices.map(|index| lorentz.pattern(index));
                self.vector(
                    [&components[0], &components[1]],
                    &left,
                    &right,
                    reference.as_ref(),
                )
            }
            spin => Err(SpinSumError::UnsupportedSpin(spin)),
        }
    }

    pub fn averaging_factor(&self) -> Atom {
        let states = match (self.average, self.spin, self.massless) {
            (false, _, _) | (_, 1, _) => 1,
            (_, 2, _) | (_, 3, true) => 2,
            _ => 3,
        };
        Atom::one() / Atom::num(states)
    }

    /// Return `(slash(p) +/- m identity)` with optional spin averaging.
    ///
    /// `slashed_momentum` must have the supplied open spinor slots. Passing
    /// an already contracted slash lets callers use their own dummy indices.
    pub fn fermion(
        &self,
        slashed_momentum: Atom,
        left: &Atom,
        right: &Atom,
    ) -> Result<Atom, SpinSumError> {
        if self.spin != 2 {
            return Err(SpinSumError::WrongSpin {
                expected: 2,
                actual: self.spin,
            });
        }
        let sign = if self.antiparticle { -1 } else { 1 };
        Ok(
            (slashed_momentum + Atom::num(sign) * &self.mass * function!(ETS.metric, left, right))
                * self.averaging_factor(),
        )
    }

    /// Return the covariant massless, massive Proca, or axial projector.
    ///
    /// `momentum` contains the momentum components in `left` and `right`.
    /// A supplied axial reference is used only for a massless particle.
    pub fn vector(
        &self,
        momentum: [&Atom; 2],
        left: &Atom,
        right: &Atom,
        reference: Option<&AxialReference>,
    ) -> Result<Atom, SpinSumError> {
        if self.spin != 3 {
            return Err(SpinSumError::WrongSpin {
                expected: 3,
                actual: self.spin,
            });
        }
        let mut projector = -function!(ETS.metric, left, right);
        if !self.covariant {
            if !self.massless {
                projector += momentum[0] * momentum[1] / self.mass.clone().pow(2);
            } else if let Some(reference) = reference {
                if reference.momentum_dot_reference.is_zero() {
                    return Err(SpinSumError::OrthogonalReference);
                }
                projector += (momentum[0] * &reference.right + &reference.left * momentum[1])
                    / &reference.momentum_dot_reference
                    - &reference.norm_squared * momentum[0] * momentum[1]
                        / reference.momentum_dot_reference.clone().pow(2);
            }
        }
        Ok(projector * self.averaging_factor())
    }
}
