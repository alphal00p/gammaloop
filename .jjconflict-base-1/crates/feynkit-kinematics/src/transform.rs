use serde::{Deserialize, Deserializer, Serialize, de};
use thiserror::Error;

use crate::{FourMomentum, KinematicScalar, ThreeMomentum};

/// Cartesian axis for an exact quarter-turn rotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    /// X axis.
    X,
    /// Y axis.
    Y,
    /// Z axis.
    Z,
}

/// An active spatial rotation with Euler matrix `Rz(gamma) * Ry(beta) * Rx(alpha)`.
#[derive(Debug, Default, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum Rotation<T> {
    /// No rotation.
    #[default]
    Identity,
    /// Euler angles `(alpha, beta, gamma)` without an intermediate precision conversion.
    Euler {
        /// First Euler angle.
        alpha: T,
        /// Second Euler angle.
        beta: T,
        /// Third Euler angle.
        gamma: T,
    },
    /// An exact positive `pi/2` rotation around a Cartesian axis.
    QuarterTurn {
        /// Rotation axis.
        axis: Axis,
    },
}

impl<T> Rotation<T> {
    /// Construct an Euler-angle rotation.
    pub const fn euler(alpha: T, beta: T, gamma: T) -> Self {
        Self::Euler { alpha, beta, gamma }
    }

    /// Construct an exact positive quarter turn.
    pub const fn quarter_turn(axis: Axis) -> Self {
        Self::QuarterTurn { axis }
    }

    /// Return whether this is the identity transformation.
    pub const fn is_identity(&self) -> bool {
        matches!(self, Self::Identity)
    }

    /// Transform the angle scalar type.
    pub fn map<U>(self, mut map: impl FnMut(T) -> U) -> Rotation<U> {
        match self {
            Self::Identity => Rotation::Identity,
            Self::Euler { alpha, beta, gamma } => {
                Rotation::euler(map(alpha), map(beta), map(gamma))
            }
            Self::QuarterTurn { axis } => Rotation::quarter_turn(axis),
        }
    }
}

impl<T: KinematicScalar> Rotation<T> {
    /// Rotate a spatial momentum.
    pub fn rotate_three(&self, momentum: &ThreeMomentum<T>) -> ThreeMomentum<T> {
        match self {
            Self::Identity => momentum.clone(),
            Self::QuarterTurn { axis: Axis::X } => ThreeMomentum::new(
                momentum.px.clone(),
                -momentum.pz.clone(),
                momentum.py.clone(),
            ),
            Self::QuarterTurn { axis: Axis::Y } => ThreeMomentum::new(
                momentum.pz.clone(),
                momentum.py.clone(),
                -momentum.px.clone(),
            ),
            Self::QuarterTurn { axis: Axis::Z } => ThreeMomentum::new(
                -momentum.py.clone(),
                momentum.px.clone(),
                momentum.pz.clone(),
            ),
            Self::Euler { alpha, beta, gamma } => {
                let sin_alpha = alpha.sin();
                let cos_alpha = alpha.cos();
                let sin_beta = beta.sin();
                let cos_beta = beta.cos();
                let sin_gamma = gamma.sin();
                let cos_gamma = gamma.cos();

                let px = momentum.px.clone();
                let py = momentum.py.clone();
                let pz = momentum.pz.clone();

                ThreeMomentum::new(
                    cos_gamma.clone() * cos_beta.clone() * px.clone()
                        + (-cos_alpha.clone() * sin_gamma.clone()
                            + sin_alpha.clone() * sin_beta.clone() * cos_gamma.clone())
                            * py.clone()
                        + (sin_alpha.clone() * sin_gamma.clone()
                            + cos_alpha.clone() * sin_beta.clone() * cos_gamma)
                            * pz.clone(),
                    sin_gamma.clone() * cos_beta.clone() * px.clone()
                        + (cos_alpha.clone() * gamma.cos()
                            + sin_alpha.clone() * sin_beta.clone() * sin_gamma.clone())
                            * py.clone()
                        + (-sin_alpha.clone() * gamma.cos()
                            + cos_alpha.clone() * sin_beta * sin_gamma)
                            * pz.clone(),
                    -beta.sin() * px
                        + cos_beta.clone() * sin_alpha * py
                        + cos_alpha * cos_beta * pz,
                )
            }
        }
    }

    /// Apply the inverse rotation to a spatial momentum.
    pub fn inverse_rotate_three(&self, momentum: &ThreeMomentum<T>) -> ThreeMomentum<T> {
        match self {
            Self::Identity => momentum.clone(),
            Self::QuarterTurn { axis: Axis::X } => ThreeMomentum::new(
                momentum.px.clone(),
                momentum.pz.clone(),
                -momentum.py.clone(),
            ),
            Self::QuarterTurn { axis: Axis::Y } => ThreeMomentum::new(
                -momentum.pz.clone(),
                momentum.py.clone(),
                momentum.px.clone(),
            ),
            Self::QuarterTurn { axis: Axis::Z } => ThreeMomentum::new(
                momentum.py.clone(),
                -momentum.px.clone(),
                momentum.pz.clone(),
            ),
            Self::Euler { alpha, beta, gamma } => {
                let sin_alpha = alpha.sin();
                let cos_alpha = alpha.cos();
                let sin_beta = beta.sin();
                let cos_beta = beta.cos();
                let sin_gamma = gamma.sin();
                let cos_gamma = gamma.cos();

                let px = momentum.px.clone();
                let py = momentum.py.clone();
                let pz = momentum.pz.clone();

                ThreeMomentum::new(
                    cos_gamma.clone() * cos_beta.clone() * px.clone()
                        + sin_gamma.clone() * cos_beta.clone() * py.clone()
                        - sin_beta.clone() * pz.clone(),
                    (-cos_alpha.clone() * sin_gamma.clone()
                        + sin_alpha.clone() * sin_beta.clone() * cos_gamma.clone())
                        * px.clone()
                        + (cos_alpha.clone() * cos_gamma.clone()
                            + sin_alpha.clone() * sin_beta.clone() * sin_gamma.clone())
                            * py.clone()
                        + cos_beta.clone() * sin_alpha * pz.clone(),
                    (alpha.sin() * sin_gamma
                        + cos_alpha.clone() * sin_beta.clone() * cos_gamma.clone())
                        * px
                        + (-alpha.sin() * cos_gamma + cos_alpha.clone() * sin_beta * gamma.sin())
                            * py
                        + cos_alpha * cos_beta * pz,
                )
            }
        }
    }

    /// Rotate only the spatial component of a four-momentum.
    pub fn rotate_four(&self, momentum: &FourMomentum<T>) -> FourMomentum<T> {
        FourMomentum::new(
            momentum.temporal.clone(),
            self.rotate_three(&momentum.spatial),
        )
    }

    /// Apply the inverse rotation to only the spatial component of a four-momentum.
    pub fn inverse_rotate_four(&self, momentum: &FourMomentum<T>) -> FourMomentum<T> {
        FourMomentum::new(
            momentum.temporal.clone(),
            self.inverse_rotate_three(&momentum.spatial),
        )
    }
}

/// A physical Lorentz boost parameterized by the three-velocity `beta`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Boost<T> {
    beta: ThreeMomentum<T>,
}

#[derive(Deserialize)]
struct BoostData<T> {
    beta: ThreeMomentum<T>,
}

/// Invalid Lorentz-boost parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum BoostError {
    /// A beta component is non-finite.
    #[error("boost velocity contains a non-finite component")]
    NonFinite,
    /// Physical boosts require `|beta| < 1`.
    #[error("boost velocity must satisfy |beta| < 1")]
    Superluminal,
}

impl<'de, T> Deserialize<'de> for Boost<T>
where
    T: KinematicScalar + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let BoostData { beta } = BoostData::deserialize(deserializer)?;
        Self::new(beta).map_err(de::Error::custom)
    }
}

impl<T: KinematicScalar> Boost<T> {
    /// Validate and construct a boost.
    pub fn new(beta: ThreeMomentum<T>) -> Result<Self, BoostError> {
        if beta
            .components()
            .into_iter()
            .any(|component| !component.is_finite())
        {
            return Err(BoostError::NonFinite);
        }

        let beta_squared = beta.norm_squared();
        if beta_squared >= beta_squared.one() {
            return Err(BoostError::Superluminal);
        }
        Ok(Self { beta })
    }

    /// Borrow the boost velocity.
    pub fn beta(&self) -> &ThreeMomentum<T> {
        &self.beta
    }

    /// Apply the boost using the mostly-minus metric convention.
    pub fn apply(&self, momentum: &FourMomentum<T>) -> FourMomentum<T> {
        let beta_squared = self.beta.norm_squared();
        let one = beta_squared.one();
        let zero = beta_squared.zero();
        let gamma = (one.clone() - beta_squared.clone()).sqrt().inv();
        let beta_dot_momentum = momentum.spatial.dot(&self.beta);
        let gamma_two = if beta_squared > zero {
            (gamma.clone() - one) / beta_squared
        } else {
            zero
        };
        let spatial_factor =
            gamma_two * beta_dot_momentum.clone() + gamma.clone() * momentum.temporal.value.clone();

        FourMomentum::from_args(
            gamma * (momentum.temporal.value.clone() + beta_dot_momentum),
            self.beta.px.clone() * spatial_factor.clone() + momentum.spatial.px.clone(),
            self.beta.py.clone() * spatial_factor.clone() + momentum.spatial.py.clone(),
            self.beta.pz.clone() * spatial_factor + momentum.spatial.pz.clone(),
        )
    }

    /// Apply the inverse boost without reconstructing configuration constants.
    pub fn apply_inverse(&self, momentum: &FourMomentum<T>) -> FourMomentum<T> {
        let inverse = Self { beta: -&self.beta };
        inverse.apply(momentum)
    }
}

#[cfg(test)]
mod tests {
    use numerica::domains::float::{FloatLike, Real};

    use super::*;

    fn assert_close(lhs: f64, rhs: f64) {
        assert!((lhs - rhs).abs() < 1.0e-12, "{lhs} != {rhs}");
    }

    fn assert_momentum_close(lhs: &FourMomentum<f64>, rhs: &FourMomentum<f64>) {
        for (left, right) in lhs.into_iter().zip(rhs) {
            assert_close(*left, *right);
        }
    }

    #[test]
    fn exact_quarter_turn_matches_euler_rotation() {
        let zero = 0.0;
        let half_pi = zero.pi() / zero.from_usize(2);
        let momentum = FourMomentum::from_args(2.0, 3.0, 1.0, 2.0);

        for (axis, euler) in [
            (Axis::X, Rotation::euler(half_pi, zero, zero)),
            (Axis::Y, Rotation::euler(zero, half_pi, zero)),
            (Axis::Z, Rotation::euler(zero, zero, half_pi)),
        ] {
            let exact = Rotation::quarter_turn(axis).rotate_four(&momentum);
            let rotated = euler.rotate_four(&momentum);
            assert_momentum_close(&exact, &rotated);
            assert_momentum_close(&euler.inverse_rotate_four(&rotated), &momentum);
        }
    }

    #[test]
    fn inverse_restores_a_general_euler_rotation() {
        let momentum = FourMomentum::from_args(7.0, 3.0, -2.0, 5.0);
        let rotation = Rotation::euler(0.37, -0.28, 0.91);
        let rotated = rotation.rotate_four(&momentum);

        assert_close(rotated.mass_squared(), momentum.mass_squared());
        assert_momentum_close(&rotation.inverse_rotate_four(&rotated), &momentum);
    }

    #[test]
    fn boost_and_inverse_preserve_invariant_mass() {
        let momentum = FourMomentum::from_args(5.0, 0.0, 0.0, 0.0);
        let boost = Boost::new(ThreeMomentum::new(0.6, 0.0, 0.0)).unwrap();
        let boosted = boost.apply(&momentum);

        assert_close(boosted.temporal.value, 6.25);
        assert_close(boosted.spatial.px, 3.75);
        assert_close(boosted.mass_squared(), momentum.mass_squared());
        assert_momentum_close(&boost.apply_inverse(&boosted), &momentum);
    }

    #[test]
    fn rejects_lightlike_or_non_finite_boost_velocity() {
        assert_eq!(
            Boost::new(ThreeMomentum::new(1.0, 0.0, 0.0)).unwrap_err(),
            BoostError::Superluminal
        );
        assert_eq!(
            Boost::new(ThreeMomentum::new(f64::NAN, 0.0, 0.0)).unwrap_err(),
            BoostError::NonFinite
        );
    }

    #[test]
    fn deserialization_rejects_superluminal_boost_velocity() {
        let error = serde_json::from_str::<Boost<f64>>(r#"{"beta":{"px":1.0,"py":0.0,"pz":0.0}}"#)
            .unwrap_err();

        assert!(error.to_string().contains("|beta| < 1"));
    }
}
