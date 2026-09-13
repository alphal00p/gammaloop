//! Reference functions used by process-level sampling acceptance tests.
//!
//! These functions deliberately live below the process parameterization layer: an
//! acceptance run can replace the physical graph value while retaining the real
//! graph routing, parameterization and Jacobian.

use crate::integrands::evaluation::EvaluationResult;
use crate::momentum::sample::LoopMomenta;
use crate::utils::{F, FloatLike};
use color_eyre::Result;
use eyre::eyre;

/// The result of evaluating a reference function through a process map,
/// together with the raw loop momenta seen by that map.  Keeping the latter
/// available lets acceptance tests check moments without bypassing the
/// parameterisation and its Jacobian.
#[derive(Clone, Debug)]
pub struct ReferenceSampleEvaluation {
    pub evaluation: EvaluationResult,
    pub loop_momenta: LoopMomenta<F<f64>>,
}

/// Summary of a reference-function acceptance run.
///
/// The estimates include both the sampling-grid weight and the process
/// parameterisation Jacobian, exactly as the production integrator does.  For
/// a normalized Gaussian, `normalization` should approach one and
/// `second_moment` should approach `|center|^2 + D * width^2`.
#[derive(Clone, Debug, Default)]
pub struct ReferenceSamplingReport {
    pub sample_count: usize,
    pub finite_sample_count: usize,
    pub normalization: f64,
    pub normalization_squared: f64,
    pub normalization_stderr: f64,
    pub second_moment: f64,
    pub expected_second_moment: f64,
    pub second_moment_stderr: f64,
    pub jacobian_min: f64,
    pub jacobian_max: f64,
    pub jacobian_abs_max: f64,
}

impl ReferenceSamplingReport {
    pub fn from_evaluations(
        evaluations: impl IntoIterator<Item = ReferenceSampleEvaluation>,
        reference: &GaussianReferenceFunction,
    ) -> Self {
        let mut report = Self {
            jacobian_min: f64::INFINITY,
            jacobian_max: f64::NEG_INFINITY,
            ..Self::default()
        };
        let mut normalization_sum = 0.0;
        let mut normalization_square_sum = 0.0;
        let mut moment_sum = 0.0;
        let mut moment_square_sum = 0.0;
        let expected = reference
            .center()
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            + reference.center().len() as f64 * reference.width().powi(2);
        report.expected_second_moment = expected;

        for sample in evaluations {
            report.sample_count += 1;
            let Some(jacobian) = sample.evaluation.parameterization_jacobian else {
                continue;
            };
            let jacobian = jacobian.0;
            report.jacobian_min = report.jacobian_min.min(jacobian);
            report.jacobian_max = report.jacobian_max.max(jacobian);
            report.jacobian_abs_max = report.jacobian_abs_max.max(jacobian.abs());
            let weight = sample.evaluation.integrator_weight.0;
            let value = sample.evaluation.integrand_result.re.0 * jacobian * weight;
            if !value.is_finite() {
                continue;
            }
            report.finite_sample_count += 1;
            normalization_sum += value;
            normalization_square_sum += value * value;
            let radius_squared = sample
                .loop_momenta
                .0
                .iter()
                .flat_map(|momentum| [&momentum.px.0, &momentum.py.0, &momentum.pz.0])
                .map(|component| component * component)
                .sum::<f64>();
            let moment = value * radius_squared;
            moment_sum += moment;
            moment_square_sum += moment * moment;
        }

        if report.finite_sample_count > 0 {
            let count = report.finite_sample_count as f64;
            report.normalization = normalization_sum;
            report.normalization_squared = normalization_square_sum;
            report.normalization_stderr =
                ((normalization_square_sum / count - (normalization_sum / count).powi(2)).max(0.0)
                    / count)
                    .sqrt();
            report.second_moment = moment_sum;
            report.second_moment_stderr =
                ((moment_square_sum / count - (moment_sum / count).powi(2)).max(0.0) / count)
                    .sqrt();
        }
        report
    }
}

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
    use crate::integrands::evaluation::EvaluationResult;
    use crate::momentum::ThreeMomentum;
    use spenso::algebra::complex::Complex;

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

    #[test]
    fn reference_report_records_gaussian_target_moment() {
        let reference = GaussianReferenceFunction::new(2.0, vec![1.0, 0.0, 0.0]).unwrap();
        let mut evaluation = EvaluationResult::zero();
        evaluation.integrand_result = Complex::new_re(F(1.0));
        evaluation.parameterization_jacobian = Some(F(1.0));
        evaluation.integrator_weight = F(1.0);
        let sample = ReferenceSampleEvaluation {
            evaluation,
            loop_momenta: LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0))]),
        };
        let report = ReferenceSamplingReport::from_evaluations([sample], &reference);
        assert_eq!(report.sample_count, 1);
        assert_eq!(report.normalization, 1.0);
        assert_eq!(report.expected_second_moment, 13.0);
    }
}
