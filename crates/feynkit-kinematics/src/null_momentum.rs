//! Constructive future-null states, distinct from arbitrary off-shell vectors.

use numerica::domains::float::Complex;
use thiserror::Error;

use crate::{
    Energy, FourMomentum, Helicity, KinematicScalar, Sign, ThreeMomentum, Wavefunction,
    WavefunctionError, WavefunctionKind,
};

#[cfg(test)]
mod tests;

/// Failure to construct a future-null momentum from represented spatial data.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum FutureNullMomentumError {
    #[error("a future-null state requires nonzero spatial momentum")]
    ZeroSpatialMomentum,
    #[error(transparent)]
    Wavefunction(#[from] WavefunctionError),
}

/// An immutable momentum whose defining energy is the positive spatial norm.
///
/// This owns the construction `E = hypot(hypot(px, py), pz)`. It does not certify
/// that rounded components obey `p.p == 0` bit for bit, or that an event conserves
/// momentum, and it does not bind a model mass or particle species. Use this same
/// owned vector for coefficients and other consumers;
/// do not attach its wavefunctions to a different numerical four-vector.
///
/// There is no raw four-vector import, mutable accessor, deserialization,
/// scalar-type conversion or tag-preserving boost/crossing operation. Rebuild
/// from the original source when changing precision. A copied or consumed raw
/// vector loses this construction contract. Physical incoming and outgoing
/// particles both use positive energy; species/role crossing belongs to the
/// caller's amplitude convention.
///
/// ```compile_fail
/// use feynkit_kinematics::{FutureNullMomentum, ThreeMomentum};
/// let state = FutureNullMomentum::from_spatial(ThreeMomentum::new(1.0, 2.0, 3.0)).unwrap();
/// state.as_momentum().temporal.value = 7.0;
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct FutureNullMomentum<T: KinematicScalar> {
    momentum: FourMomentum<T>,
}

impl<T: KinematicScalar> FutureNullMomentum<T> {
    /// Construct from finite, real, nonzero represented spatial components.
    ///
    /// A supplied zero vector is an invalid source. If exact nonzero input was
    /// already lost before this call, the source reconstruction layer must
    /// detect that loss; these represented values cannot recover its provenance.
    pub fn from_spatial(spatial: ThreeMomentum<T>) -> Result<Self, FutureNullMomentumError> {
        if spatial.components().iter().any(|value| !value.is_finite()) {
            return Err(WavefunctionError::NonFinite.into());
        }
        if spatial
            .components()
            .iter()
            .all(|value| **value == value.zero())
        {
            return Err(FutureNullMomentumError::ZeroSpatialMomentum);
        }
        let energy = spatial.px.hypot(&spatial.py).hypot(&spatial.pz);
        if !energy.is_finite() {
            return Err(WavefunctionError::NonFinite.into());
        }
        if energy <= energy.zero() {
            return Err(WavefunctionError::InsufficientPrecision {
                bits: energy.get_precision(),
                reason: "nonzero future-null source lost its positive spatial norm",
            }
            .into());
        }
        Ok(Self {
            momentum: FourMomentum::new(Energy::new(energy), spatial),
        })
    }

    /// Borrow the same immutable four-vector used by the owned external states.
    pub fn as_momentum(&self) -> &FourMomentum<T> {
        &self.momentum
    }

    /// Drop the construction contract and return a generic, mutable four-vector.
    pub fn into_momentum(self) -> FourMomentum<T> {
        self.momentum
    }

    /// Construct a physical massless state using the existing MadGraph phases.
    ///
    /// Fermions use the checked helicity spinors `xi±` and `A=sqrt(E)*sqrt(2)`:
    /// `u+=(0,A xi+)`, `u-=(A xi-,0)`, `v+=(-A xi-,0)`, `v-=(0,-A xi+)`.
    /// Each entry denotes a Weyl block. The opposite block is identically zero
    /// by construction, never by clipping a computed mass-shell residual.
    /// Barred states reuse the existing Dirac/vector adjoint.
    ///
    /// Vectors admit only transverse helicities; scalars admit helicity zero.
    /// Massive states, past-directed states and raw off-shell imports are outside
    /// this API. The generic `FourMomentum::wavefunction` remains unchanged.
    /// A near-south nonaxis chart can still require more working precision;
    /// reconstruct the source after its typed precision failure.
    pub fn wavefunction(
        &self,
        kind: WavefunctionKind,
        helicity: Helicity,
    ) -> Result<Wavefunction<T>, WavefunctionError> {
        if matches!(
            kind,
            WavefunctionKind::EpsilonBar | WavefunctionKind::UBar | WavefunctionKind::VBar
        ) {
            return Ok(self.wavefunction(kind.bar(), helicity)?.bar());
        }
        if kind == WavefunctionKind::Epsilon && helicity.is_zero() {
            return Err(WavefunctionError::UndefinedLongitudinal);
        }
        if !kind.is_spinor() {
            return self.momentum.wavefunction(kind, helicity);
        }
        let sign = if helicity.is_plus() {
            Sign::Positive
        } else if helicity.is_minus() {
            Sign::Negative
        } else {
            return Err(WavefunctionError::Helicity(kind));
        };
        let xi = self
            .momentum
            .helicity_spinor(if kind == WavefunctionKind::U {
                sign
            } else {
                -sign
            })?;
        let energy = &self.momentum.temporal.value;
        // Forming 2*E first can overflow although its square root fits.
        let amplitude = energy.sqrt() * energy.from_usize(2).sqrt();
        if !amplitude.is_finite() {
            return Err(WavefunctionError::NonFinite);
        }
        if amplitude <= amplitude.zero() {
            return Err(WavefunctionError::InsufficientPrecision {
                bits: energy.get_precision(),
                reason: "positive future-null spinor normalization was lost",
            });
        }
        let phase = if kind == WavefunctionKind::U {
            amplitude
        } else {
            -amplitude
        };
        let active = [
            Self::scale_component(&xi[0], &phase)?,
            Self::scale_component(&xi[1], &phase)?,
        ];
        let zero = Complex::new(energy.zero(), energy.zero());
        let components = match (kind, sign) {
            (WavefunctionKind::U, Sign::Positive) | (WavefunctionKind::V, Sign::Negative) => {
                vec![zero.clone(), zero, active[0].clone(), active[1].clone()]
            }
            _ => vec![active[0].clone(), active[1].clone(), zero.clone(), zero],
        };
        Wavefunction::from_components(kind, components)
    }

    // Keep scalar multiplication explicit: each nonzero Cartesian component
    // must survive, including after the checked dimensionless chart succeeds.
    fn scale_component(value: &Complex<T>, factor: &T) -> Result<Complex<T>, WavefunctionError> {
        let result = Complex::new(
            value.re.clone() * factor.clone(),
            value.im.clone() * factor.clone(),
        );
        if !result.re.is_finite() || !result.im.is_finite() {
            return Err(WavefunctionError::NonFinite);
        }
        if factor != &factor.zero()
            && ((value.re != value.re.zero() && result.re == result.re.zero())
                || (value.im != value.im.zero() && result.im == result.im.zero()))
        {
            return Err(WavefunctionError::InsufficientPrecision {
                bits: factor.get_precision(),
                reason: "nonzero future-null wavefunction component underflowed",
            });
        }
        Ok(result)
    }
}
