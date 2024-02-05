use std::time::Duration;

use colored::Colorize;
use log::info;
use num::Complex;

use crate::{observables::Event, utils::format_evaluation_time, Precision};

const MAX_QUAD_PREC_FOR_GREEN: f64 = 10.0;
const MAX_REL_PREC_FOR_GREEN: f64 = 1e-6;

pub struct EvaluationResult {
    pub integrand_result: Complex<f64>,
    pub integrator_weight: f64,
    pub event_buffer: Vec<Event>,
    pub evaluation_metadata: EvaluationMetaData,
}

pub struct EvaluationMetaData {
    pub rep3d_evaluation_time: Duration,
    pub parameterization_time: Duration,
    pub relative_instability_error: Complex<f64>,
    pub highest_precision: Precision,
}

pub struct MetaDataStatistics {
    pub avg_rep3d_evaluation_time: Duration,
    pub avg_parameterization_time: Duration,
    pub avg_relative_instability_error: Complex<f64>,
    pub num_evals: usize,
    pub num_double_precision_evals: usize,
    pub num_quadruple_precision_evals: usize,
}

struct TempStatisticsCounter {
    total_rep3d_evaluation_time: Duration,
    total_parameterization_time: Duration,
    total_relative_instability_error: Complex<f64>,
    num_double_precision_evals: usize,
    num_quadruple_precision_evals: usize,
}

impl TempStatisticsCounter {
    fn new() -> Self {
        Self {
            total_rep3d_evaluation_time: Duration::ZERO,
            total_parameterization_time: Duration::ZERO,
            total_relative_instability_error: Complex::new(0.0, 0.0),
            num_double_precision_evals: 0,
            num_quadruple_precision_evals: 0,
        }
    }

    fn to_metadata_statistics(&self, len: usize) -> MetaDataStatistics {
        MetaDataStatistics {
            avg_rep3d_evaluation_time: self.total_rep3d_evaluation_time / len as u32,
            avg_parameterization_time: self.total_parameterization_time / len as u32,
            avg_relative_instability_error: self.total_relative_instability_error / len as f64,
            num_evals: len,
            num_double_precision_evals: self.num_double_precision_evals,
            num_quadruple_precision_evals: self.num_quadruple_precision_evals,
        }
    }
}

impl MetaDataStatistics {
    pub fn from_evaluation_results(data: &[EvaluationResult]) -> Self {
        let len = data.len();

        let temp_statistics_counter = data.iter().fold(
            TempStatisticsCounter::new(),
            |mut accumulator, data_entry| {
                accumulator.total_rep3d_evaluation_time +=
                    data_entry.evaluation_metadata.rep3d_evaluation_time;
                accumulator.total_parameterization_time +=
                    data_entry.evaluation_metadata.parameterization_time;
                accumulator.total_relative_instability_error +=
                    data_entry.evaluation_metadata.relative_instability_error;
                match data_entry.evaluation_metadata.highest_precision {
                    Precision::Double => accumulator.num_double_precision_evals += 1,
                    Precision::Quad => accumulator.num_quadruple_precision_evals += 1,
                    _ => (),
                }

                accumulator
            },
        );

        temp_statistics_counter.to_metadata_statistics(len)
    }

    pub fn merge(self, other: Self) -> Self {
        Self {
            avg_rep3d_evaluation_time: (self.avg_rep3d_evaluation_time
                + other.avg_rep3d_evaluation_time)
                / 2,
            avg_parameterization_time: (self.avg_parameterization_time
                + other.avg_parameterization_time)
                / 2,
            avg_relative_instability_error: (self.avg_relative_instability_error
                + other.avg_relative_instability_error)
                / 2.0,
            num_evals: self.num_evals + other.num_evals,
            num_double_precision_evals: (self.num_double_precision_evals
                + other.num_double_precision_evals),
            num_quadruple_precision_evals: (self.num_quadruple_precision_evals
                + other.num_quadruple_precision_evals),
        }
    }

    fn new_empty() -> Self {
        Self {
            avg_rep3d_evaluation_time: Duration::ZERO,
            avg_parameterization_time: Duration::ZERO,
            avg_relative_instability_error: Complex::new(0.0, 0.0),
            num_evals: 0,
            num_double_precision_evals: 0,
            num_quadruple_precision_evals: 0,
        }
    }

    pub fn merge_list(list: Vec<Self>) -> Self {
        if let Some(merged) = list.into_iter().reduce(|acc, x| acc.merge(x)) {
            merged
        } else {
            Self::new_empty()
        }
    }

    pub fn print_stats(&self) {
        let precentage_double =
            self.num_double_precision_evals as f64 / self.num_evals as f64 * 100.0;
        let precentage_quad =
            self.num_quadruple_precision_evals as f64 / self.num_evals as f64 * 100.0;

        let percentage_double_str = format!("{:.2}%", precentage_double);
        let percentage_quad_str = format!("{:.2}%", precentage_quad);

        let (percentage_double_str_formatted, percentage_quad_str_formatted) =
            if precentage_quad > MAX_QUAD_PREC_FOR_GREEN {
                (percentage_double_str.red(), percentage_quad_str.red())
            } else {
                (percentage_double_str.green(), percentage_quad_str.green())
            };

        let (avg_rel_error_str_re, avg_rel_error_str_im) = (
            format!("{:+.2e}", self.avg_relative_instability_error.re),
            format!("{:+.2e}", self.avg_relative_instability_error.im),
        );

        let avg_rel_error_str_re_formatted =
            if self.avg_relative_instability_error.re > MAX_REL_PREC_FOR_GREEN {
                avg_rel_error_str_re.red()
            } else {
                avg_rel_error_str_re.green()
            };

        let avg_rel_error_str_im_formatted =
            if self.avg_relative_instability_error.im > MAX_REL_PREC_FOR_GREEN {
                avg_rel_error_str_im.red()
            } else {
                avg_rel_error_str_im.green()
            };

        let param_time_str = format_evaluation_time(self.avg_parameterization_time);
        let rep3d_time_str = format_evaluation_time(self.avg_rep3d_evaluation_time);

        info!("Average parameterization time: {}", param_time_str,);
        info!("Average rep3d evaluation time: {}", rep3d_time_str,);
        info!(
            "Average relative instability error: ( {}, {} )",
            avg_rel_error_str_re_formatted, avg_rel_error_str_im_formatted
        );
        info!(
            "Percentage of double precision evaluations: {}",
            percentage_double_str_formatted
        );
        info!(
            "Percentage of quadruple precision evaluations: {}",
            percentage_quad_str_formatted
        );
    }
}
