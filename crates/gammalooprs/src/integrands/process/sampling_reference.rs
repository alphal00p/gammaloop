//! Reference functions used by process-level sampling acceptance tests.
//!
//! These functions deliberately live below the process parameterization layer: an
//! acceptance run can replace the physical graph value while retaining the real
//! graph routing, parameterization and Jacobian.

use crate::momentum::sample::LoopMomenta;
use crate::utils::{F, FloatLike};
use color_eyre::Result;
use eyre::eyre;

/// A normalized isotropic Gaussian on the raw loop-momentum coordinates.
///
/// `center` is flattened as `(p_0x, p_0y, p_0z, p_1x, ...)`.  The function is
/// normalized with respect to the unconstrained spatial loop-momentum volume,
/// making it suitable for testing a real process parameterization independently
/// of its physical numerator and denominator.
#[derive(Clone, Debug, PartialEq)]
pub struct GaussianReferenceFunction {
    width: f64,
    center: Vec<f64>,
}

impl GaussianReferenceFunction {
    /// Construct a Gaussian with a positive width and a center of the given dimension.
    pub fn new(width: f64, center: Vec<f64>) -> Result<Self> {
        if !width.is_finite() || width <= 0.0 {
            return Err(eyre!(
                "Gaussian reference width must be finite and positive"
            ));
        }
        if center.iter().any(|component| !component.is_finite()) {
            return Err(eyre!("Gaussian reference center must be finite"));
        }
        if center.len() % 3 != 0 {
            return Err(eyre!(
                "Gaussian reference center has dimension {}, expected a multiple of three",
                center.len()
            ));
        }
        Ok(Self { width, center })
    }

    /// Construct a centered Gaussian for `n_loop_momenta` spatial loop momenta.
    pub fn centered(width: f64, n_loop_momenta: usize) -> Result<Self> {
        Self::new(width, vec![0.0; 3 * n_loop_momenta])
    }

    pub fn width(&self) -> f64 {
        self.width
    }

    pub fn center(&self) -> &[f64] {
        &self.center
    }

    pub(crate) fn evaluate<T: FloatLike>(&self, loop_momenta: &LoopMomenta<F<T>>) -> Result<F<T>> {
        let dimension = 3 * loop_momenta.0.len();
        if self.center.len() != dimension {
            return Err(eyre!(
                "Gaussian reference has dimension {}, but the process has {} loop-momentum coordinates",
                self.center.len(),
                dimension
            ));
        }

        let zero = loop_momenta
            .0
            .first()
            .map(|momentum| momentum.px.zero())
            .ok_or_else(|| eyre!("Gaussian reference cannot evaluate a zero-loop process"))?;
        let width = F::<T>::from_f64(self.width);
        let width_squared = width.square();
        let two = zero.from_i64(2);
        let normalization = (two.clone() * zero.PI()).sqrt().inv() / width;
        let normalization = (0..dimension).fold(zero.one(), |acc, _| acc * &normalization);
        let mut distance_squared = zero.zero();
        for (momentum, center) in loop_momenta.0.iter().zip(self.center.chunks_exact(3)) {
            let dx = momentum.px.clone() - F::from_f64(center[0]);
            let dy = momentum.py.clone() - F::from_f64(center[1]);
            let dz = momentum.pz.clone() - F::from_f64(center[2]);
            distance_squared += dx.square() + dy.square() + dz.square();
        }
        let exponent = -distance_squared / (two * width_squared);
        Ok(normalization * F(exponent.0.exp()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::momentum::ThreeMomentum;

    #[test]
    fn centered_gaussian_is_normalized_at_the_origin() {
        let reference = GaussianReferenceFunction::centered(2.0, 1).unwrap();
        let momenta = LoopMomenta(vec![ThreeMomentum::new(F(0.0), F(0.0), F(0.0))]);
        let value = reference.evaluate(&momenta).unwrap().0;
        let expected = (2.0 * std::f64::consts::PI * 4.0).powf(-1.5);
        assert!((value - expected).abs() < 1.0e-14);
    }

    #[test]
    fn gaussian_rejects_wrong_loop_dimension() {
        let reference = GaussianReferenceFunction::centered(1.0, 1).unwrap();
        let momenta = LoopMomenta(vec![
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
        ]);
        let error = reference.evaluate(&momenta).unwrap_err().to_string();
        assert!(error.contains("dimension"));
    }
}
