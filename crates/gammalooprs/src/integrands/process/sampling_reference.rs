//! Reference functions used by process-level sampling acceptance tests.
//!
//! These functions deliberately live below the process parameterization layer: an
//! acceptance run can replace the physical graph value while retaining the real
//! graph routing, parameterization and Jacobian.

use super::ProcessIntegrand;
use crate::graph::FeynmanGraph;
use crate::integrands::evaluation::{EvaluationResult, PreciseEvaluationResult};
use crate::momentum::sample::LoopMomenta;
use crate::utils::{F, FloatLike};
use color_eyre::Result;
use eyre::eyre;
use schemars::JsonSchema;
use serde::{Deserialize, Deserializer, Serialize};
use symbolica::numerical_integration::StatisticsAccumulator;

/// The result of evaluating a reference function through a process map,
/// together with moments evaluated at every actual mapped raw point. These
/// retain the same channel factors as the value; one outer sample can contain
/// several points, whose contributions must be summed before statistical squaring.
#[derive(Clone, Debug)]
pub struct ReferenceSampleEvaluation {
    pub evaluation: EvaluationResult,
    pub moments: ReferenceMoments<f64>,
}

/// Reference contributions carry the same factors as the scalar; any remaining
/// parameterization Jacobian and the outer grid weight are still outside.
/// Jacobian extrema describe the effective sampling factors at the individual
/// mapped points, since an explicit channel sum has no single map Jacobian.
#[derive(Clone, Debug)]
pub struct ReferenceMoments<T: FloatLike> {
    pub second_moment: F<T>,
    pub jacobian_min: f64,
    pub jacobian_max: f64,
}

impl<T: FloatLike> ReferenceMoments<T> {
    pub(crate) fn rescale(&mut self, factor: &F<T>) {
        self.second_moment *= factor;
        // Extrema are reporting diagnostics, not inputs to a map or precision rescue.
        let factor = factor.clone().into_ff64().0;
        let a = self.jacobian_min * factor;
        let b = self.jacobian_max * factor;
        self.jacobian_min = a.min(b);
        self.jacobian_max = a.max(b);
    }

    pub(crate) fn merge_in_place(&mut self, other: Self) {
        self.second_moment += other.second_moment;
        self.jacobian_min = self.jacobian_min.min(other.jacobian_min);
        self.jacobian_max = self.jacobian_max.max(other.jacobian_max);
    }

    pub(crate) fn try_into_f64(self) -> Result<ReferenceMoments<f64>> {
        let narrowed = self.second_moment.into_ff64();
        if !self.second_moment.0.is_finite()
            || !narrowed.0.is_finite()
            || (narrowed == narrowed.zero() && self.second_moment != self.second_moment.zero())
            || !self.jacobian_min.is_finite()
            || !self.jacobian_max.is_finite()
        {
            return Err(eyre!(
                "native reference moment {} or its Jacobian diagnostics cannot be represented by the f64 acceptance reporting boundary",
                self.second_moment
            ));
        }
        Ok(ReferenceMoments {
            second_moment: narrowed,
            jacobian_min: self.jacobian_min,
            jacobian_max: self.jacobian_max,
        })
    }
}

/// Summary of a reference-function acceptance run.
///
/// The estimates include both the sampling-grid weight and the process
/// parameterisation Jacobian, exactly as the production integrator does, with
/// one average over the requested draws. For
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
        let mut normalization = StatisticsAccumulator::<F<f64>>::new();
        let mut second_moment = StatisticsAccumulator::<F<f64>>::new();
        let mut last_value = 0.0;
        let mut last_moment = 0.0;
        let mut normalization_square_sum = 0.0;
        let expected = reference
            .center()
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            + reference.center().len() as f64 * reference.width().powi(2);
        report.expected_second_moment = expected;

        for sample in evaluations {
            report.sample_count += 1;
            let jacobian = sample
                .evaluation
                .parameterization_jacobian
                .unwrap_or(F(f64::NAN));
            let mut moments = sample.moments;
            moments.rescale(&jacobian);
            report.jacobian_min = report.jacobian_min.min(moments.jacobian_min);
            report.jacobian_max = report.jacobian_max.max(moments.jacobian_max);
            report.jacobian_abs_max = report
                .jacobian_abs_max
                .max(moments.jacobian_min.abs())
                .max(moments.jacobian_max.abs());
            let weight = sample.evaluation.integrator_weight.0;
            let value = sample.evaluation.integrand_result.re.0 * jacobian.0 * weight;
            let moment = moments.second_moment.0 * weight;
            last_value = value;
            last_moment = moment;
            if !value.is_finite()
                || !moment.is_finite()
                || !(value * value).is_finite()
                || !(moment * moment).is_finite()
                || !jacobian.0.is_finite()
                || jacobian.0 <= 0.0
                || !moments.jacobian_min.is_finite()
                || !moments.jacobian_max.is_finite()
            {
                // A failed draw remains part of the requested sample count.
                // Its loss is structural, never a reason to renormalize the finite subset.
                normalization.add_sample(F(0.0), None);
                second_moment.add_sample(F(0.0), None);
                continue;
            }
            report.finite_sample_count += 1;
            normalization.add_sample(F(value), None);
            second_moment.add_sample(F(moment), None);
            normalization_square_sum += value * value;
        }

        if report.sample_count == 1 && report.finite_sample_count == 1 {
            // The production statistics owner requires two draws for a variance.
            report.normalization = last_value;
            report.second_moment = last_moment;
            report.normalization_squared = last_value * last_value;
        } else if report.sample_count > 1 {
            normalization.update_iter(false);
            second_moment.update_iter(false);
            report.normalization = normalization.avg.0;
            report.normalization_squared = normalization_square_sum / report.sample_count as f64;
            report.normalization_stderr = normalization.err.0;
            report.second_moment = second_moment.avg.0;
            report.second_moment_stderr = second_moment.err.0;
        }
        if report.finite_sample_count != report.sample_count {
            // Counts/extrema remain diagnostic, but no finite central estimate
            // may make a structurally invalid acceptance run appear successful.
            report.normalization = f64::NAN;
            report.normalization_squared = f64::NAN;
            report.normalization_stderr = f64::NAN;
            report.second_moment = f64::NAN;
            report.second_moment_stderr = f64::NAN;
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
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct GaussianReferenceFunction {
    width: f64,
    #[schemars(default)]
    center: Vec<f64>,
}

impl<'de> Deserialize<'de> for GaussianReferenceFunction {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> std::result::Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Input {
            width: f64,
            #[serde(default)]
            center: Vec<f64>,
        }
        let input = Input::deserialize(deserializer)?;
        Self::new(input.width, input.center).map_err(serde::de::Error::custom)
    }
}

impl std::str::FromStr for GaussianReferenceFunction {
    type Err = serde_json::Error;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        serde_json::from_str(value)
    }
}

impl GaussianReferenceFunction {
    /// Workspace reporting uses two dimensionless observables, both of mean one.
    pub const INTEGRATION_CONVENTION: &str = "gaussian_and_normalized_raw_second_moment_v1";

    /// Resolve centered shorthand against the actual generated process. A single
    /// raw-frame Gaussian cannot describe graph groups of different dimensions.
    pub fn for_integrand(&self, integrand: &ProcessIntegrand) -> Result<Self> {
        let loops = match integrand {
            ProcessIntegrand::Amplitude(integrand) => integrand
                .data
                .graph_terms
                .iter()
                .map(|term| term.graph.get_loop_number())
                .collect::<std::collections::BTreeSet<_>>(),
            ProcessIntegrand::CrossSection(integrand) => integrand
                .data
                .graph_terms
                .iter()
                .map(|term| term.graph.get_loop_number())
                .collect::<std::collections::BTreeSet<_>>(),
        };
        if loops.len() != 1 || loops.contains(&0) {
            return Err(eyre!(
                "reference integration requires one nonzero loop dimension per process slot; found {loops:?}"
            ));
        }
        let n_loops = *loops.first().unwrap();
        let reference = if self.center.is_empty() {
            Self::centered(self.width, n_loops)?
        } else {
            self.clone()
        };
        if reference.center.len() != 3 * n_loops {
            return Err(eyre!(
                "reference center has {} coordinates, but the process has {}",
                reference.center.len(),
                3 * n_loops
            ));
        }
        // This is the f64 CLI/statistics metadata boundary; native reference APIs
        // retain their wider-range input and moment contracts.
        let scale = reference.expected_second_moment::<f64>().0;
        if !scale.is_finite() || scale <= 0.0 {
            return Err(eyre!(
                "reference integration requires a finite positive raw second-moment scale for f64 reporting"
            ));
        }
        let settings = integrand.get_settings();
        if settings.selectors.values().any(|selector| selector.active)
            || !settings.observables.is_empty()
        {
            return Err(eyre!(
                "reference integration substitutes the physical event value; disable physical selectors and observables for this acceptance run"
            ));
        }
        Ok(reference)
    }

    /// Construct the known raw moment in the active precision from input settings.
    /// An empty center is only a CLI descriptor shorthand, resolved before evaluation.
    pub fn expected_second_moment<T: FloatLike>(&self) -> F<T> {
        let width = F::<T>::from_f64(self.width);
        self.center.iter().fold(width.zero(), |sum, value| {
            sum + F::<T>::from_f64(*value).square()
        }) + width.square() * width.from_usize(self.center.len())
    }

    /// Report the reference through the ordinary integrator's existing complex
    /// accumulators: Re is normalization and Im is the normalized raw moment.
    /// Both have already received native map factors; the grid weight stays outside.
    pub(crate) fn integration_report(
        &self,
        result: PreciseEvaluationResult,
    ) -> Result<EvaluationResult> {
        macro_rules! report {
            ($result:expr) => {{
                let mut result = $result;
                if result.evaluation_metadata.is_nan {
                    return Err(eyre!(
                        "reference integration cannot report a nonfinite finalized value"
                    ));
                }
                let moment = result
                    .reference_moments
                    .take()
                    .ok_or_else(|| eyre!("reference integration produced no mapped moment"))?;
                let scale = self.expected_second_moment();
                let normalized = &moment.second_moment / &scale;
                if scale.is_nan()
                    || scale.is_infinite()
                    || scale <= scale.zero()
                    || normalized.is_nan()
                    || normalized.is_infinite()
                    || (normalized == normalized.zero()
                        && moment.second_moment != moment.second_moment.zero())
                {
                    return Err(eyre!(
                        "native normalized reference moment is not representable"
                    ));
                }
                if let Some(absolute) = &mut result.absolute_integrand_result {
                    absolute.im = normalized.abs();
                }
                result.integrand_result.im = normalized;
                result.try_into_f64()
            }};
        }
        match result {
            PreciseEvaluationResult::Double(result) => report!(result),
            PreciseEvaluationResult::Quad(result) => report!(result),
            PreciseEvaluationResult::Arb(result) => report!(result),
        }
    }

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
        if !center.len().is_multiple_of(3) {
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
        for (momentum, center) in loop_momenta
            .0
            .iter()
            .zip(self.center.as_chunks::<3>().0.iter())
        {
            let dx = momentum.px.clone() - F::from_f64(center[0]);
            let dy = momentum.py.clone() - F::from_f64(center[1]);
            let dz = momentum.pz.clone() - F::from_f64(center[2]);
            distance_squared += dx.square() + dy.square() + dz.square();
        }
        let exponent = -distance_squared / (two * width_squared);
        // Native map rescue does not certify arbitrarily small Gaussian tails:
        // this existing body may underflow before weighting. Log-domain reference
        // weighting and heavier-tailed normalized probes remain separate work.
        Ok(normalization * F(exponent.0.exp()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrands::evaluation::{EvaluationResult, PreciseEvaluationResult};
    use crate::momentum::ThreeMomentum;
    use spenso::algebra::complex::Complex;

    #[test]
    fn native_reference_moments_reject_loss_at_reporting_boundary() {
        use crate::utils::ArbPrec;
        let one = F::<ArbPrec>::default().one();
        for exponent in [-400, 400] {
            let moments = ReferenceMoments {
                second_moment: one.from_i64(10).powi(exponent),
                jacobian_min: 1.0,
                jacobian_max: 1.0,
            };
            assert!(!moments.second_moment.is_nan() && !moments.second_moment.is_infinite());
            assert!(
                moments
                    .try_into_f64()
                    .unwrap_err()
                    .to_string()
                    .contains("f64 acceptance reporting boundary")
            );
        }
        for second_moment in [one.zero(), one.from_i64(10).powi(-200), one.clone()] {
            let moments = ReferenceMoments {
                second_moment: second_moment.clone(),
                jacobian_min: 0.5,
                jacobian_max: 2.0,
            }
            .try_into_f64()
            .unwrap();
            assert_eq!(moments.second_moment, second_moment.into_ff64());
        }
    }

    #[test]
    fn reference_descriptor_validates_serde_and_normalizes_native_moment_before_reporting() {
        use crate::integrands::evaluation::{EvaluationMetaData, GenericEvaluationResult};
        use crate::utils::ArbPrec;
        let shorthand: GaussianReferenceFunction = serde_json::from_str(r#"{"width":2}"#).unwrap();
        assert!(shorthand.center().is_empty());
        for invalid in [
            r#"{"width":0}"#,
            r#"{"width":-2}"#,
            r#"{"width":2,"center":[1]}"#,
            r#"{"width":2,"centre":[0,0,0]}"#,
        ] {
            assert!(serde_json::from_str::<GaussianReferenceFunction>(invalid).is_err());
        }
        assert!(toml::from_str::<GaussianReferenceFunction>("width=2\ncenter=[inf,0,0]").is_err());
        let schema = schemars::schema_for!(GaussianReferenceFunction);
        assert_eq!(schema.as_value()["required"], serde_json::json!(["width"]));

        // The raw moment is outside binary64 range, but the normalized
        // observable is finite. Narrow only after native normalization.
        let reference = GaussianReferenceFunction::centered(5.0e153, 1).unwrap();
        let one = F::<ArbPrec>::default().one();
        let moment = reference.expected_second_moment::<ArbPrec>() * one.from_usize(4);
        assert!(moment.into_ff64().0.is_infinite());
        let native = GenericEvaluationResult {
            reference_moments: Some(ReferenceMoments {
                second_moment: moment,
                jacobian_min: 1.0,
                jacobian_max: 1.0,
            }),
            integrand_result: Complex::new_re(one.clone()),
            absolute_integrand_result: Some(Complex::new_re(one.clone())),
            parameterization_jacobian: Some(one.clone()),
            integrator_weight: one.from_usize(7),
            event_groups: Default::default(),
            evaluation_metadata: EvaluationMetaData::new_empty(),
        };
        let output = reference
            .integration_report(PreciseEvaluationResult::Arb(native.clone()))
            .unwrap();
        assert_eq!(output.integrand_result, Complex::new(F(1.0), F(4.0)));
        assert_eq!(output.integrator_weight, F(7.0));
        let mut failed = native;
        failed.evaluation_metadata.is_nan = true;
        assert!(
            reference
                .integration_report(PreciseEvaluationResult::Arb(failed))
                .is_err()
        );
    }

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
            moments: ReferenceMoments {
                second_moment: F(1.0),
                jacobian_min: 1.0,
                jacobian_max: 1.0,
            },
        };
        let report = ReferenceSamplingReport::from_evaluations([sample], &reference);
        assert_eq!(report.sample_count, 1);
        assert_eq!(report.normalization, 1.0);
        assert_eq!(report.expected_second_moment, 13.0);
    }

    #[test]
    fn reference_report_averages_weighted_draws_once_and_keeps_failed_draws() {
        let reference = GaussianReferenceFunction::centered(1.0, 1).unwrap();
        let samples = [1.0, 2.0, 3.0, 4.0].map(|weight| {
            let mut evaluation = EvaluationResult::zero();
            evaluation.integrand_result = Complex::new_re(F(1.0));
            evaluation.parameterization_jacobian = Some(F(2.0));
            evaluation.integrator_weight = F(weight);
            ReferenceSampleEvaluation {
                evaluation,
                moments: ReferenceMoments {
                    second_moment: F(3.0),
                    jacobian_min: 1.0,
                    jacobian_max: 1.0,
                },
            }
        });
        let report = ReferenceSamplingReport::from_evaluations(samples.clone(), &reference);
        assert_eq!(report.normalization, 5.0);
        assert_eq!(report.normalization_squared, 30.0);
        assert!((report.normalization_stderr - (5.0f64 / 3.0).sqrt()).abs() < 1.0e-14);
        assert_eq!(report.second_moment, 15.0);
        assert!((report.second_moment_stderr - 15.0f64.sqrt()).abs() < 1.0e-14);
        let mut failed = samples.clone();
        failed[3].moments.second_moment = F(f64::NAN);
        let report = ReferenceSamplingReport::from_evaluations(failed, &reference);
        assert_eq!(report.sample_count, 4);
        assert_eq!(report.finite_sample_count, 3);
        assert!(report.normalization.is_nan());
        assert!(report.second_moment.is_nan());
        let mut overflow = samples;
        overflow[0].evaluation.integrator_weight = F(1.0e200);
        let report = ReferenceSamplingReport::from_evaluations(overflow, &reference);
        assert_eq!(report.finite_sample_count, 3);
        assert!(report.normalization.is_nan());
    }

    #[test]
    fn reference_channel_moments_are_scaled_and_summed_before_statistical_squaring() {
        use crate::integrands::evaluation::GraphEvaluationResult;
        let reference = GaussianReferenceFunction::centered(1.0, 1).unwrap();
        let samples = [(1.0, 3.0), (3.0, 1.0)].map(|(a, b)| {
            let mut total = GraphEvaluationResult::zero(F(0.0));
            for (value, factor, radius_squared) in [(a, 2.0, 1.0), (b, 2.0, 4.0)] {
                let mut contribution = GraphEvaluationResult::zero(F(0.0));
                contribution.integrand_result = Complex::new_re(F(value));
                contribution.reference_moments = Some(ReferenceMoments {
                    second_moment: F(value * radius_squared),
                    jacobian_min: 1.0,
                    jacobian_max: 1.0,
                });
                contribution.apply_sampling_factor(F(factor));
                total.merge_in_place(contribution);
            }
            let mut evaluation = EvaluationResult::zero();
            evaluation.integrand_result = total.integrand_result;
            evaluation.parameterization_jacobian = Some(F(1.0));
            evaluation.integrator_weight = F(1.0);
            ReferenceSampleEvaluation {
                evaluation,
                moments: total.reference_moments.unwrap(),
            }
        });
        let report = ReferenceSamplingReport::from_evaluations(samples, &reference);
        assert_eq!(report.normalization, 8.0);
        assert_eq!(report.normalization_squared, 64.0);
        assert_eq!(report.normalization_stderr, 0.0);
        assert_eq!(report.second_moment, 20.0);
        assert!((report.second_moment_stderr - 6.0).abs() < 1.0e-14);
    }
}
