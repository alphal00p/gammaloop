use std::time::Duration;

use crate::evaluation_result::EvaluationMetaData;
use crate::evaluation_result::EvaluationResult;
use crate::integrands::*;
use crate::utils;
use crate::utils::FloatLike;
use crate::ParameterizationMapping;
use crate::Precision;
use crate::Settings;
use num::traits::ToPrimitive;
use num::Complex;
use serde::{Deserialize, Serialize};
use symbolica::numerical_integration::{ContinuousGrid, Grid, Sample};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct HFunctionTestSettings {
    pub h_function: crate::HFunctionSettings,
}

#[derive(Clone)]
pub struct HFunctionTestIntegrand {
    pub settings: Settings,
    pub n_dim: usize,
    pub integrand_settings: HFunctionTestSettings,
}

#[allow(unused)]
impl HFunctionTestIntegrand {
    pub fn new(
        settings: Settings,
        integrand_settings: HFunctionTestSettings,
    ) -> HFunctionTestIntegrand {
        let n_dim = 1;
        HFunctionTestIntegrand {
            settings,
            n_dim,
            integrand_settings,
        }
    }

    fn evaluate_sample_generic<T: FloatLike>(&self, xs: &[T]) -> (Complex<T>, Duration, Duration) {
        let e_cm = Into::<T>::into(self.settings.kinematics.e_cm);
        let mut jac = T::one();

        let parameterization_start = std::time::Instant::now();

        let t = match self.settings.parameterization.mapping {
            ParameterizationMapping::Log => {
                // r = e_cm * ln(1 + b*x/(1-x))
                let x = xs[0];
                let b: T = Into::<T>::into(self.settings.parameterization.b);
                let r = e_cm * (T::one() + b * x / (T::one() - x)).ln();
                jac *= e_cm * b / (T::one() - x) / (T::one() + x * (b - T::one()));

                r
            }
            ParameterizationMapping::Linear => {
                // r = e_cm * b * x/(1-x)
                let b = Into::<T>::into(self.settings.parameterization.b);
                let radius = e_cm * b * xs[0] / (T::one() - xs[0]);
                jac *= <T as num::traits::Float>::powi(e_cm * b + radius, 2) / e_cm / b;
                radius
            }
        };

        let parameterization_time = parameterization_start.elapsed();

        let evaluation_time = std::time::Instant::now();
        let h = utils::h(t, None, None, &self.integrand_settings.h_function);
        let evaluation_time = evaluation_time.elapsed();

        (
            Complex::new(h * jac, T::zero()),
            parameterization_time,
            evaluation_time,
        )
    }
}

#[allow(unused)]
impl HasIntegrand for HFunctionTestIntegrand {
    fn create_grid(&self) -> Grid<f64> {
        Grid::Continuous(ContinuousGrid::new(
            self.n_dim,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
            self.settings.integrator.bin_number_evolution.clone(),
            self.settings.integrator.train_on_avg,
        ))
    }

    fn get_n_dim(&self) -> usize {
        self.n_dim
    }

    fn evaluate_sample(
        &self,
        sample: &Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
        _max_eval: f64,
    ) -> EvaluationResult {
        let start_evaluate_sample = std::time::Instant::now();

        let xs = match sample {
            Sample::Continuous(_w, v) => v,
            _ => panic!("Wrong sample type"),
        };

        let mut sample_xs = vec![];
        sample_xs.extend(xs);
        if self.settings.general.debug > 1 {
            println!(
                "Sampled x-space : ( {} )",
                sample_xs
                    .iter()
                    .map(|&x| format!("{:.16}", x))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            println!("Integrator weight : {:+.16e}", wgt);
        }

        // TODO implement stability check
        let (integration_result, parameterization_timing, evaluation_timing, precision) =
            if use_f128 {
                let sample_xs_f128 = sample_xs
                    .iter()
                    .map(|x| Into::<f128::f128>::into(*x))
                    .collect::<Vec<_>>();
                if self.settings.general.debug > 1 {
                    println!(
                        "f128 Upcasted x-space sample : ( {} )",
                        sample_xs_f128
                            .iter()
                            .map(|&x| format!("{:+e}", x))
                            .collect::<Vec<_>>()
                            .join(", ")
                    );
                }
                let (res, parameterization_timing, evaluation_timing) =
                    self.evaluate_sample_generic(sample_xs_f128.as_slice());
                (
                    Complex::new(
                        f128::f128::to_f64(&res.re).unwrap(),
                        f128::f128::to_f64(&res.im).unwrap(),
                    ),
                    parameterization_timing,
                    evaluation_timing,
                    Precision::Quad,
                )
            } else {
                let (res, parameterization_timing, evaluation_timing) =
                    self.evaluate_sample_generic(sample_xs.as_slice());

                (
                    res,
                    parameterization_timing,
                    evaluation_timing,
                    Precision::Double,
                )
            };

        let is_nan = integration_result.re.is_nan() || integration_result.im.is_nan();

        let evaluation_metadata = EvaluationMetaData {
            total_timing: start_evaluate_sample.elapsed(),
            rep3d_evaluation_time: evaluation_timing,
            parameterization_time: parameterization_timing,
            relative_instability_error: Complex::new(0.0, 0.0),
            highest_precision: precision,
            is_nan,
        };

        EvaluationResult {
            integrand_result: integration_result,
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata,
        }
    }
}
