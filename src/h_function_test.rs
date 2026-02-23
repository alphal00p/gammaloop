use std::time::Duration;

use crate::evaluation_result::EvaluationResult;
use crate::evaluation_result::{EvaluationMetaData, StabilityEvaluation};
use crate::integrands::*;
use crate::model::Model;
use crate::settings::RuntimeSettings;
use crate::settings::runtime::HFunctionSettings;
use crate::settings::runtime::ParameterizationMapping;
use crate::settings::runtime::Precision;
use crate::settings::runtime::SamplingSettings;
use crate::utils;
use crate::utils::F;
use crate::utils::FloatLike;
use crate::utils::f128;
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::domains::float::FloatLike as SymFloatLike;
use symbolica::numerical_integration::{ContinuousGrid, Grid, Sample};
use tracing::info;

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Default, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
// #[trait_decode(trait= GammaLoopContext)]
pub struct HFunctionTestSettings {
    pub h_function: HFunctionSettings,
}

#[derive(Clone)]
pub struct HFunctionTestIntegrand {
    pub settings: RuntimeSettings,
    pub n_dim: usize,
    pub integrand_settings: HFunctionTestSettings,
}

#[allow(unused)]
impl HFunctionTestIntegrand {
    pub(crate) fn new(
        settings: RuntimeSettings,
        integrand_settings: HFunctionTestSettings,
    ) -> HFunctionTestIntegrand {
        let n_dim = 1;
        HFunctionTestIntegrand {
            settings,
            n_dim,
            integrand_settings,
        }
    }

    fn evaluate_sample_generic<T: FloatLike>(
        &self,
        xs: &[F<T>],
    ) -> (Complex<F<T>>, Duration, Duration) {
        let e_cm: F<T> = F::<T>::from_f64(self.settings.kinematics.e_cm);
        let one = e_cm.one();
        let mut jac = e_cm.one();

        let parameterization_start = std::time::Instant::now();

        let t = match &self.settings.sampling {
            SamplingSettings::Default(parameterization_settings) => {
                match parameterization_settings.mapping {
                    ParameterizationMapping::Log => {
                        // r = e_cm * ln(1 + b*x/(1-x))
                        let x = xs[0].clone();
                        let b: F<T> = F::<T>::from_f64(parameterization_settings.b);
                        let r = &e_cm * (&one + &b * &x / (&one - &x)).ln();
                        jac *= &e_cm * &b / (&one - &x) / (&one + &x * (&b - &one));

                        r
                    }
                    ParameterizationMapping::Linear => {
                        // r = e_cm * b * x/(1-x)
                        let b: F<T> = F::<T>::from_f64(parameterization_settings.b);
                        let radius = &e_cm * &b * &xs[0] / (&one - &xs[0]);
                        jac *= (&e_cm * &b + &radius).powi(2) / &e_cm / &b;
                        radius
                    }
                }
            }
            _ => {
                panic!("Unsupported sampling type");
            }
        };

        let parameterization_time = parameterization_start.elapsed();

        let evaluation_time = std::time::Instant::now();
        let h = utils::h(&t, None, None, &self.integrand_settings.h_function);
        let evaluation_time = evaluation_time.elapsed();

        ((h * jac).into(), parameterization_time, evaluation_time)
    }
}

#[allow(unused)]
impl HasIntegrand for HFunctionTestIntegrand {
    fn name(&self) -> String {
        "HFunctionTestIntegrand".to_string()
    }
    fn create_grid(&self) -> Grid<F<f64>> {
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
        &mut self,
        sample: &Sample<F<f64>>,
        _model: &Model,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        _max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        let start_evaluate_sample = std::time::Instant::now();

        let xs = match sample {
            Sample::Continuous(_w, v) => v,
            _ => panic!("Wrong sample type"),
        };

        let mut sample_xs: Vec<F<f64>> = vec![];
        sample_xs.extend(xs);

        info!(
            "Sampled x-space : ( {} )",
            sample_xs
                .iter()
                .map(|&x| format!("{:.16}", x))
                .collect::<Vec<_>>()
                .join(", ")
        );
        info!("Integrator weight : {:+.16e}", wgt);

        // TODO implement stability check
        let (integration_result, parameterization_timing, evaluation_timing, precision) =
            if use_f128 {
                let sample_xs_f128 = sample_xs
                    .iter()
                    .map(|x| x.higher())
                    .collect::<Vec<F<f128>>>();
                info!(
                    "f128 Upcasted x-space sample : ( {} )",
                    sample_xs_f128
                        .iter()
                        .map(|x| format!("{:+e}", x))
                        .collect::<Vec<_>>()
                        .join(", ")
                );

                let (res, parameterization_timing, evaluation_timing) =
                    self.evaluate_sample_generic(sample_xs_f128.as_slice());
                (
                    Complex::new(res.re.lower(), res.im.lower()),
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
            relative_instability_error: Complex::new_zero(),
            highest_precision: precision,
            is_nan,
            final_is_stable: !is_nan,
            loop_momenta_escalation: None,
            stability_evaluations: vec![StabilityEvaluation {
                precision,
                result: integration_result,
                parameterization_time: parameterization_timing,
                ltd_evaluation_time: evaluation_timing,
                is_stable: !is_nan,
                instability_reason: None,
                rotated_results: Vec::new(),
            }],
        };

        Ok(EvaluationResult {
            integrand_result: integration_result,
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata,
        })
    }
}
