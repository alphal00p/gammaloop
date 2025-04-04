use std::time::Duration;

use bincode::{Decode, Encode};
use colored::Colorize;
use log::info;
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use symbolica::domains::float::ConstructibleFloat;

use crate::observables::Event;
use crate::{
    utils::{format_evaluation_time, F},
    Precision,
};

/// The result of an evaluation of the integrand
#[derive(Clone)]
pub struct EvaluationResult {
    pub integrand_result: Complex<F<f64>>,
    pub integrator_weight: F<f64>,
    pub event_buffer: Vec<Event>,
    pub evaluation_metadata: EvaluationMetaData,
}

impl EvaluationResult {
    pub fn zero() -> Self {
        Self {
            integrand_result: Complex::new_zero(),
            integrator_weight: F(0.0),
            event_buffer: Vec::new(),
            evaluation_metadata: EvaluationMetaData::new_empty(),
        }
    }
}

/// Useful metadata generated during the evaluation, this may be expanded in the future to include more information
#[derive(Copy, Clone)]
pub struct EvaluationMetaData {
    pub total_timing: Duration,
    pub rep3d_evaluation_time: Duration,
    pub parameterization_time: Duration,
    pub relative_instability_error: Complex<F<f64>>,
    pub highest_precision: Precision,
    pub is_nan: bool,
}

impl EvaluationMetaData {
    pub fn new_empty() -> Self {
        Self {
            total_timing: Duration::ZERO,
            rep3d_evaluation_time: Duration::ZERO,
            parameterization_time: Duration::ZERO,
            relative_instability_error: Complex::new_zero(),
            highest_precision: Precision::Double,
            is_nan: false,
        }
    }
}

/// This struct merges the evaluation metadata of many evaluations into a single struct
#[derive(Copy, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct StatisticsCounter {
    pub num_evals: usize,
    sum_rep3d_evaluation_time: Duration,
    sum_parameterization_time: Duration,
    sum_total_evaluation_time: Duration,
    sum_relative_instability_error: (F<f64>, F<f64>),
    num_double_precision_evals: usize,
    num_quadruple_precision_evals: usize,
    num_nan_evals: usize,
}

impl StatisticsCounter {
    /// Turn a slice of evaluation results into a statistics counter
    pub fn from_evaluation_results(data: &[EvaluationResult]) -> Self {
        let statistics_counter = data.iter().fold(
            StatisticsCounter::new_empty(),
            |mut accumulator, data_entry| {
                accumulator.sum_rep3d_evaluation_time +=
                    data_entry.evaluation_metadata.rep3d_evaluation_time;
                accumulator.sum_parameterization_time +=
                    data_entry.evaluation_metadata.parameterization_time;
                accumulator.sum_relative_instability_error.0 +=
                    data_entry.evaluation_metadata.relative_instability_error.re;
                accumulator.sum_relative_instability_error.1 +=
                    data_entry.evaluation_metadata.relative_instability_error.im;
                accumulator.sum_total_evaluation_time +=
                    data_entry.evaluation_metadata.total_timing;

                accumulator.num_evals += 1;
                match data_entry.evaluation_metadata.highest_precision {
                    Precision::Double => accumulator.num_double_precision_evals += 1,
                    Precision::Quad => accumulator.num_quadruple_precision_evals += 1,
                    _ => (),
                }

                if data_entry.evaluation_metadata.is_nan {
                    accumulator.num_nan_evals += 1;
                }

                accumulator
            },
        );

        statistics_counter
    }

    /// Merge two statistics counters into a single one, but keeping the original ones unchanged
    pub fn merged(&self, other: &Self) -> Self {
        Self {
            sum_rep3d_evaluation_time: self.sum_rep3d_evaluation_time
                + other.sum_rep3d_evaluation_time,
            sum_parameterization_time: self.sum_parameterization_time
                + other.sum_parameterization_time,
            sum_relative_instability_error: (
                self.sum_relative_instability_error.0 + other.sum_relative_instability_error.0,
                self.sum_relative_instability_error.1 + other.sum_relative_instability_error.1,
            ),
            num_evals: self.num_evals + other.num_evals,
            num_double_precision_evals: (self.num_double_precision_evals
                + other.num_double_precision_evals),
            num_quadruple_precision_evals: (self.num_quadruple_precision_evals
                + other.num_quadruple_precision_evals),
            sum_total_evaluation_time: self.sum_total_evaluation_time
                + other.sum_total_evaluation_time,
            num_nan_evals: self.num_nan_evals + other.num_nan_evals,
        }
    }

    pub fn new_empty() -> Self {
        Self {
            sum_rep3d_evaluation_time: Duration::ZERO,
            sum_parameterization_time: Duration::ZERO,
            sum_relative_instability_error: (F(0.0), F(0.0)),
            sum_total_evaluation_time: Duration::ZERO,
            num_evals: 0,
            num_double_precision_evals: 0,
            num_quadruple_precision_evals: 0,
            num_nan_evals: 0,
        }
    }

    /// Merge a list of statistics counters into a single one.
    pub fn merge_list(list: Vec<Self>) -> Self {
        if let Some(merged) = list.into_iter().reduce(|acc, x| acc.merged(&x)) {
            merged
        } else {
            Self::new_empty()
        }
    }

    /// Compute the average time spent in the evaluate_sample function.
    pub fn get_avg_total_timing(&self) -> Duration {
        let avg_total_timing = self.sum_total_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_total_timing)
    }

    /// Compute the average time spent in a single evaluation of the three-dimensional representation.
    /// Note that a single evaluation contains at least two evaluations of the three-dimensional representation.
    pub fn get_avg_rep3d_timing(&self) -> Duration {
        let avg_rep3d_timing = self.sum_rep3d_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_rep3d_timing)
    }

    /// Compute the average time spent in the parameterization of the integrand. Especaially useful for monitoring the performance of tropical sampling.
    pub fn get_avg_param_timing(&self) -> Duration {
        let avg_param_timing = self.sum_parameterization_time.as_secs_f64() / self.num_evals as f64;

        Duration::from_secs_f64(avg_param_timing)
    }

    /// Get the average relative error computed during instability checks
    pub fn get_avg_instabillity_error(&self) -> (F<f64>, F<f64>) {
        (
            self.sum_relative_instability_error.0 / F::<f64>::new_from_usize(self.num_evals),
            self.sum_relative_instability_error.1 / F::<f64>::new_from_usize(self.num_evals),
        )
    }

    /// Get the percentage of evaluations that were done in double precision and were stable.
    pub fn get_percentage_f64(&self) -> f64 {
        self.num_double_precision_evals as f64 / self.num_evals as f64 * 100.0
    }

    /// Get the percentage of evaluations that went to quadruple precision and were stable.
    pub fn get_percentage_f128(&self) -> f64 {
        self.num_quadruple_precision_evals as f64 / self.num_evals as f64 * 100.0
    }

    pub fn get_percentage_nan(&self) -> f64 {
        self.num_nan_evals as f64 / self.num_evals as f64 * 100.0
    }

    #[allow(clippy::format_in_format_args)]
    pub fn display_status(&self) {
        let time_ltd_formatted = format_evaluation_time(self.get_avg_rep3d_timing());
        let param_time_formatted = format_evaluation_time(self.get_avg_param_timing());
        let total_time = format_evaluation_time(self.get_avg_total_timing());

        info!(
            "|  {}  | {} {} | {} {} | {} {}",
            format!("{:-7}", "timing").blue().bold(),
            format!("{:-7}", "total:"),
            format!("{:-9}", total_time).green(),
            format!("{:-7}", "param:"),
            format!("{:-9}", param_time_formatted).green(),
            format!("{:-7}", "ltd:"),
            format!("{:-9}", time_ltd_formatted).green(),
        );

        info!(
            "|  {}  | {} {} | {} {} | {} {}",
            format!("{:-7}", "evals").blue().bold(),
            format!("{:-7}", "f64:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_f64())).green(),
            format!("{:-7}", "f128:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_f128())).green(),
            format!("{:-7}", "nan:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_nan())).green(),
        );
    }
}
