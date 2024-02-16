use std::time::Duration;

use num::Complex;
use serde::{Deserialize, Serialize};

use crate::{observables::Event, Precision};

pub struct EvaluationResult {
    pub integrand_result: Complex<f64>,
    pub integrator_weight: f64,
    pub event_buffer: Vec<Event>,
    pub evaluation_metadata: EvaluationMetaData,
}

pub struct EvaluationMetaData {
    pub total_timing: Duration,
    pub rep3d_evaluation_time: Duration,
    pub parameterization_time: Duration,
    pub relative_instability_error: Complex<f64>,
    pub highest_precision: Precision,
}

#[derive(Copy, Clone)]
pub struct StatisticsCounter {
    pub num_evals: usize,
    sum_rep3d_evaluation_time: Duration,
    sum_parameterization_time: Duration,
    sum_total_evaluation_time: Duration,
    sum_relative_instability_error: Complex<f64>,
    num_double_precision_evals: usize,
    num_quadruple_precision_evals: usize,
}

#[derive(Serialize, Deserialize)]
pub struct SerializableMeteDataStatistics {
    num_evals: usize,
    sum_rep3d_evaluation_time: u128,
    sum_parameterization_time: u128,
    sum_total_evaluation_time: u128,
    sum_relative_instability_error: (f64, f64),
    num_double_precision_evals: usize,
    num_quadruple_precision_evals: usize,
}

impl SerializableMeteDataStatistics {
    pub fn from_metadata_statistics(metadata: StatisticsCounter) -> Self {
        Self {
            sum_rep3d_evaluation_time: metadata.sum_rep3d_evaluation_time.as_nanos(),
            sum_parameterization_time: metadata.sum_parameterization_time.as_nanos(),
            sum_total_evaluation_time: metadata.sum_total_evaluation_time.as_nanos(),
            sum_relative_instability_error: (
                metadata.sum_relative_instability_error.re,
                metadata.sum_relative_instability_error.im,
            ),
            num_evals: metadata.num_evals,
            num_double_precision_evals: metadata.num_double_precision_evals,
            num_quadruple_precision_evals: metadata.num_quadruple_precision_evals,
        }
    }

    pub fn into_metadata_statistics(self) -> StatisticsCounter {
        StatisticsCounter {
            sum_rep3d_evaluation_time: Duration::from_nanos(self.sum_rep3d_evaluation_time as u64), // this cast might ruin the statistics on large runs
            sum_parameterization_time: Duration::from_nanos(self.sum_parameterization_time as u64),
            sum_total_evaluation_time: Duration::from_nanos(self.sum_total_evaluation_time as u64),
            sum_relative_instability_error: Complex::new(
                self.sum_relative_instability_error.0,
                self.sum_relative_instability_error.1,
            ),
            num_evals: self.num_evals,
            num_double_precision_evals: self.num_double_precision_evals,
            num_quadruple_precision_evals: self.num_quadruple_precision_evals,
        }
    }
}

impl StatisticsCounter {
    pub fn from_evaluation_results(data: &[EvaluationResult]) -> Self {
        let statistics_counter = data.iter().fold(
            StatisticsCounter::new_empty(),
            |mut accumulator, data_entry| {
                accumulator.sum_rep3d_evaluation_time +=
                    data_entry.evaluation_metadata.rep3d_evaluation_time;
                accumulator.sum_parameterization_time +=
                    data_entry.evaluation_metadata.parameterization_time;
                accumulator.sum_relative_instability_error +=
                    data_entry.evaluation_metadata.relative_instability_error;
                accumulator.sum_total_evaluation_time +=
                    data_entry.evaluation_metadata.total_timing;

                accumulator.num_evals += 1;
                match data_entry.evaluation_metadata.highest_precision {
                    Precision::Double => accumulator.num_double_precision_evals += 1,
                    Precision::Quad => accumulator.num_quadruple_precision_evals += 1,
                    _ => (),
                }

                accumulator
            },
        );

        statistics_counter
    }

    // this isn't completly correct if one of them is emtpy
    pub fn merged(&self, other: &Self) -> Self {
        Self {
            sum_rep3d_evaluation_time: self.sum_rep3d_evaluation_time
                + other.sum_rep3d_evaluation_time,
            sum_parameterization_time: self.sum_parameterization_time
                + other.sum_parameterization_time,
            sum_relative_instability_error: (self.sum_relative_instability_error
                + other.sum_relative_instability_error),
            num_evals: self.num_evals + other.num_evals,
            num_double_precision_evals: (self.num_double_precision_evals
                + other.num_double_precision_evals),
            num_quadruple_precision_evals: (self.num_quadruple_precision_evals
                + other.num_quadruple_precision_evals),
            sum_total_evaluation_time: self.sum_total_evaluation_time
                + other.sum_total_evaluation_time,
        }
    }

    pub fn new_empty() -> Self {
        Self {
            sum_rep3d_evaluation_time: Duration::ZERO,
            sum_parameterization_time: Duration::ZERO,
            sum_relative_instability_error: Complex::new(0.0, 0.0),
            sum_total_evaluation_time: Duration::ZERO,
            num_evals: 0,
            num_double_precision_evals: 0,
            num_quadruple_precision_evals: 0,
        }
    }

    pub fn merge_list(list: Vec<Self>) -> Self {
        if let Some(merged) = list.into_iter().reduce(|acc, x| acc.merged(&x)) {
            merged
        } else {
            Self::new_empty()
        }
    }

    pub fn get_avg_total_timing(&self) -> Duration {
        let avg_total_timing = self.sum_total_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_total_timing)
    }

    pub fn get_avg_rep3d_timing(&self) -> Duration {
        let avg_rep3d_timing = self.sum_rep3d_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_rep3d_timing)
    }

    pub fn get_avg_param_timing(&self) -> Duration {
        let avg_param_timing = self.sum_parameterization_time.as_secs_f64() / self.num_evals as f64;

        Duration::from_secs_f64(avg_param_timing)
    }

    pub fn get_avg_instabillity_error(&self) -> Complex<f64> {
        self.sum_relative_instability_error / self.num_evals as f64
    }

    pub fn get_percentage_f64(&self) -> f64 {
        self.num_double_precision_evals as f64 / self.num_evals as f64 * 100.0
    }

    pub fn get_percentage_f128(&self) -> f64 {
        self.num_quadruple_precision_evals as f64 / self.num_evals as f64 * 100.0
    }
}
