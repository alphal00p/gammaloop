use std::fmt::Display;
use std::time::Duration;

use bincode::{Decode, Encode};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::domains::float::Constructible;
use tracing::info;

use crate::observables::{Event, GenericEvent};
use crate::{
    settings::runtime::Precision,
    utils::{F, FloatLike, format_evaluation_time},
};

#[derive(Clone, Debug)]
pub struct GraphEvaluationResult<T: FloatLike> {
    pub integrand_result: Complex<F<T>>,
    pub generated_events: Vec<GenericEvent<T>>,
    pub event_processing_time: Duration,
    pub generated_event_count: usize,
    pub accepted_event_count: usize,
}

impl<T: FloatLike> GraphEvaluationResult<T> {
    pub fn zero(zero: F<T>) -> Self {
        Self {
            integrand_result: Complex::new_re(zero),
            generated_events: Vec::new(),
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
        }
    }

    pub fn merge_in_place(&mut self, mut other: Self) {
        self.integrand_result += other.integrand_result;
        self.generated_events.append(&mut other.generated_events);
        self.event_processing_time += other.event_processing_time;
        self.generated_event_count += other.generated_event_count;
        self.accepted_event_count += other.accepted_event_count;
    }

    pub fn into_f64(self) -> GraphEvaluationResult<f64> {
        GraphEvaluationResult {
            integrand_result: Complex::new(
                self.integrand_result.re.into_ff64(),
                self.integrand_result.im.into_ff64(),
            ),
            generated_events: self
                .generated_events
                .into_iter()
                .map(|event| event.to_f64())
                .collect(),
            event_processing_time: self.event_processing_time,
            generated_event_count: self.generated_event_count,
            accepted_event_count: self.accepted_event_count,
        }
    }
}

/// The result of an evaluation of the integrand
#[derive(Clone, Serialize, Debug)]
pub struct EvaluationResult {
    pub integrand_result: Complex<F<f64>>,
    pub integrator_weight: F<f64>,
    pub event_buffer: Vec<Event>,
    pub event_processing_time: Duration,
    pub generated_event_count: usize,
    pub accepted_event_count: usize,
    pub evaluation_metadata: EvaluationMetaData,
}

impl EvaluationResult {
    pub fn zero() -> Self {
        Self {
            integrand_result: Complex::new_zero(),
            integrator_weight: F(0.0),
            event_buffer: Vec::new(),
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
            evaluation_metadata: EvaluationMetaData::new_empty(),
        }
    }
}

/// Per-precision evaluation details produced during stability checks.
#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub enum StabilityFailureReason {
    ErrorThreshold,
    ZeroError,
    WeightThreshold,
}

#[derive(Clone, Debug, Serialize)]
pub struct RotatedEvaluation {
    pub rotation: String,
    pub result: Complex<F<f64>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct StabilityEvaluation {
    pub precision: Precision,
    pub result: Complex<F<f64>>,
    pub parameterization_time: Duration,
    pub integrand_evaluation_time: Duration,
    pub evaluator_evaluation_time: Duration,
    pub is_stable: bool,
    pub instability_reason: Option<StabilityFailureReason>,
    pub rotated_results: Vec<RotatedEvaluation>,
}

#[derive(Clone, Debug, Serialize)]
pub struct LoopMomentaEscalationMetrics {
    pub sum_norm: f64,
    pub threshold: f64,
}

/// Useful metadata generated during the evaluation, this may be expanded in the future to include more information
#[derive(Clone, Serialize, Debug)]
pub struct EvaluationMetaData {
    pub total_timing: Duration,
    pub integrand_evaluation_time: Duration,
    pub evaluator_evaluation_time: Duration,
    pub parameterization_time: Duration,
    pub event_time: Duration,
    pub relative_instability_error: Complex<F<f64>>,
    pub highest_precision: Precision,
    pub is_nan: bool,
    pub final_is_stable: bool,
    pub loop_momenta_escalation: Option<LoopMomentaEscalationMetrics>,
    pub stability_evaluations: Vec<StabilityEvaluation>,
}

impl Display for EvaluationMetaData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let stability_summary = if self.stability_evaluations.is_empty() {
            "none".to_string()
        } else {
            self.stability_evaluations
                .iter()
                .map(|evaluation| {
                    format!(
                        "{}:{}",
                        evaluation.precision,
                        if evaluation.is_stable {
                            "stable"
                        } else {
                            "unstable"
                        }
                    )
                })
                .collect::<Vec<_>>()
                .join(", ")
        };
        write!(
            f,
            "EvaluationMetaData {{ total_timing: {:?}, integrand_evaluation_time: {:?}, evaluator_evaluation_time: {:?}, parameterization_time: {:?}, event_time: {:?}, relative_instability_error: {:?}, highest_precision: {:?}, is_nan: {}, final_is_stable: {}, loop_momenta_escalation: {:?}, stability_evaluations: {} }}",
            self.total_timing,
            self.integrand_evaluation_time,
            self.evaluator_evaluation_time,
            self.parameterization_time,
            self.event_time,
            self.relative_instability_error,
            self.highest_precision,
            self.is_nan,
            self.final_is_stable,
            self.loop_momenta_escalation,
            stability_summary
        )
    }
}

impl EvaluationMetaData {
    pub(crate) fn new_empty() -> Self {
        Self {
            total_timing: Duration::ZERO,
            integrand_evaluation_time: Duration::ZERO,
            evaluator_evaluation_time: Duration::ZERO,
            parameterization_time: Duration::ZERO,
            event_time: Duration::ZERO,
            relative_instability_error: Complex::new_zero(),
            highest_precision: Precision::Double,
            is_nan: false,
            final_is_stable: true,
            loop_momenta_escalation: None,
            stability_evaluations: Vec::new(),
        }
    }
}

/// This struct merges the evaluation metadata of many evaluations into a single struct
#[derive(Debug, Copy, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct StatisticsCounter {
    pub num_evals: usize,
    sum_integrand_evaluation_time: Duration,
    sum_evaluator_evaluation_time: Duration,
    sum_parameterization_time: Duration,
    sum_event_time: Duration,
    sum_total_evaluation_time: Duration,
    sum_relative_instability_error: (F<f64>, F<f64>),
    num_double_precision_evals: usize,
    num_quadruple_precision_evals: usize,
    num_arb_precision_evals: usize,
    num_nan_evals: usize,
    sum_generated_event_count: usize,
    sum_accepted_event_count: usize,
}

impl StatisticsCounter {
    /// Turn a slice of evaluation results into a statistics counter
    pub(crate) fn from_evaluation_results(data: &[EvaluationResult]) -> Self {
        data.iter().fold(
            StatisticsCounter::new_empty(),
            |mut accumulator, data_entry| {
                accumulator.sum_integrand_evaluation_time +=
                    data_entry.evaluation_metadata.integrand_evaluation_time;
                accumulator.sum_evaluator_evaluation_time +=
                    data_entry.evaluation_metadata.evaluator_evaluation_time;
                accumulator.sum_parameterization_time +=
                    data_entry.evaluation_metadata.parameterization_time;
                accumulator.sum_event_time += data_entry.evaluation_metadata.event_time;
                accumulator.sum_relative_instability_error.0 +=
                    data_entry.evaluation_metadata.relative_instability_error.re;
                accumulator.sum_relative_instability_error.1 +=
                    data_entry.evaluation_metadata.relative_instability_error.im;
                accumulator.sum_total_evaluation_time +=
                    data_entry.evaluation_metadata.total_timing;
                accumulator.sum_generated_event_count += data_entry.generated_event_count;
                accumulator.sum_accepted_event_count += data_entry.accepted_event_count;

                accumulator.num_evals += 1;
                match data_entry.evaluation_metadata.highest_precision {
                    Precision::Double => accumulator.num_double_precision_evals += 1,
                    Precision::Quad => accumulator.num_quadruple_precision_evals += 1,
                    Precision::Arb => accumulator.num_arb_precision_evals += 1,
                    // _ => (),
                }

                if data_entry.evaluation_metadata.is_nan {
                    accumulator.num_nan_evals += 1;
                }

                accumulator
            },
        )
    }

    /// Merge two statistics counters into a single one, but keeping the original ones unchanged
    pub(crate) fn merged(&self, other: &Self) -> Self {
        Self {
            sum_integrand_evaluation_time: self.sum_integrand_evaluation_time
                + other.sum_integrand_evaluation_time,
            sum_evaluator_evaluation_time: self.sum_evaluator_evaluation_time
                + other.sum_evaluator_evaluation_time,
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
            num_arb_precision_evals: self.num_arb_precision_evals + other.num_arb_precision_evals,
            sum_total_evaluation_time: self.sum_total_evaluation_time
                + other.sum_total_evaluation_time,
            num_nan_evals: self.num_nan_evals + other.num_nan_evals,
            sum_event_time: self.sum_event_time + other.sum_event_time,
            sum_generated_event_count: self.sum_generated_event_count
                + other.sum_generated_event_count,
            sum_accepted_event_count: self.sum_accepted_event_count
                + other.sum_accepted_event_count,
        }
    }

    pub(crate) fn new_empty() -> Self {
        Self {
            sum_integrand_evaluation_time: Duration::ZERO,
            sum_evaluator_evaluation_time: Duration::ZERO,
            sum_parameterization_time: Duration::ZERO,
            sum_event_time: Duration::ZERO,
            sum_relative_instability_error: (F(0.0), F(0.0)),
            sum_total_evaluation_time: Duration::ZERO,
            num_evals: 0,
            num_double_precision_evals: 0,
            num_quadruple_precision_evals: 0,
            num_arb_precision_evals: 0,
            num_nan_evals: 0,
            sum_generated_event_count: 0,
            sum_accepted_event_count: 0,
        }
    }

    /// Merge a list of statistics counters into a single one.
    #[allow(dead_code)]
    pub(crate) fn merge_list(list: Vec<Self>) -> Self {
        if let Some(merged) = list.into_iter().reduce(|acc, x| acc.merged(&x)) {
            merged
        } else {
            Self::new_empty()
        }
    }

    /// Compute the average time spent in the evaluate_sample function.
    pub(crate) fn get_avg_total_timing(&self) -> Duration {
        let avg_total_timing = self.sum_total_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_total_timing)
    }

    /// Compute the average time spent in the original integrand evaluation call.
    pub(crate) fn get_avg_integrand_timing(&self) -> Duration {
        let avg_integrand_timing =
            self.sum_integrand_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_integrand_timing)
    }

    /// Compute the average time spent inside Symbolica evaluator function calls.
    pub(crate) fn get_avg_evaluator_timing(&self) -> Duration {
        let avg_evaluator_timing =
            self.sum_evaluator_evaluation_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_evaluator_timing)
    }

    /// Compute the average time spent in the parameterization of the integrand. Especaially useful for monitoring the performance of tropical sampling.
    pub(crate) fn get_avg_param_timing(&self) -> Duration {
        let avg_param_timing = self.sum_parameterization_time.as_secs_f64() / self.num_evals as f64;

        Duration::from_secs_f64(avg_param_timing)
    }

    /// Get the average relative error computed during instability checks
    #[allow(dead_code)]
    pub(crate) fn get_avg_instabillity_error(&self) -> (F<f64>, F<f64>) {
        (
            self.sum_relative_instability_error.0 / F::<f64>::new_from_usize(self.num_evals),
            self.sum_relative_instability_error.1 / F::<f64>::new_from_usize(self.num_evals),
        )
    }

    /// Get the percentage of evaluations that were done in double precision and were stable.
    pub(crate) fn get_percentage_f64(&self) -> f64 {
        self.num_double_precision_evals as f64 / self.num_evals as f64 * 100.0
    }

    /// Get the percentage of evaluations that went to quadruple precision and were stable.
    pub(crate) fn get_percentage_f128(&self) -> f64 {
        self.num_quadruple_precision_evals as f64 / self.num_evals as f64 * 100.0
    }

    pub(crate) fn get_percentage_arb(&self) -> f64 {
        self.num_arb_precision_evals as f64 / self.num_evals as f64 * 100.0
    }

    pub(crate) fn get_percentage_nan(&self) -> f64 {
        self.num_nan_evals as f64 / self.num_evals as f64 * 100.0
    }

    #[allow(clippy::format_in_format_args)]
    pub(crate) fn display_status(&self) {
        let time_integrand_formatted = format_evaluation_time(self.get_avg_integrand_timing());
        let time_evaluators_formatted = format_evaluation_time(self.get_avg_evaluator_timing());
        let param_time_formatted = format_evaluation_time(self.get_avg_param_timing());
        let total_time = format_evaluation_time(self.get_avg_total_timing());

        info!(
            "|  {}  | {} {} | {} {} | {} {} | {} {}",
            format!("{:-7}", "timing").blue().bold(),
            format!("{:-7}", "total:"),
            format!("{:-9}", total_time).green(),
            format!("{:-7}", "param:"),
            format!("{:-9}", param_time_formatted).green(),
            format!("{:-7}", "itg:"),
            format!("{:-9}", time_integrand_formatted).green(),
            format!("{:-11}", "evaluators:"),
            format!("{:-9}", time_evaluators_formatted).green(),
        );

        info!(
            "|  {}  | {} {} | {} {} | {} {} | {} {}",
            format!("{:-7}", "evals").blue().bold(),
            format!("{:-7}", "f64:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_f64())).green(),
            format!("{:-7}", "f128:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_f128())).green(),
            format!("{:-7}", "arb:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_arb())).green(),
            format!("{:-7}", "nan:"),
            format!("{:-9}", format!("{:.2}%", self.get_percentage_nan())).green(),
        );
    }
}
