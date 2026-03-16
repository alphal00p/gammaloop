use std::fmt::Display;
use std::time::Duration;

use bincode::{Decode, Encode};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use tabled::{Table, Tabled, settings::Style};
use tracing::info;

use crate::observables::{
    EventGroupList, GenericEventGroupList, ObservableSnapshotBundle,
    events::{format_complex_generic, format_optional_real_generic, format_real_generic},
};
use crate::{
    settings::runtime::Precision,
    utils::{ArbPrec, F, FloatLike, f128, format_evaluation_time},
};

#[derive(Clone, Debug)]
pub struct GraphEvaluationResult<T: FloatLike> {
    pub integrand_result: Complex<F<T>>,
    pub event_groups: GenericEventGroupList<T>,
    pub event_processing_time: Duration,
    pub generated_event_count: usize,
    pub accepted_event_count: usize,
}

impl<T: FloatLike> GraphEvaluationResult<T> {
    pub fn zero(zero: F<T>) -> Self {
        Self {
            integrand_result: Complex::new_re(zero),
            event_groups: GenericEventGroupList::default(),
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
        }
    }

    pub fn merge_in_place(&mut self, mut other: Self) {
        self.integrand_result += other.integrand_result;
        self.event_groups.append(&mut other.event_groups);
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
            event_groups: self.event_groups.to_f64(),
            event_processing_time: self.event_processing_time,
            generated_event_count: self.generated_event_count,
            accepted_event_count: self.accepted_event_count,
        }
    }
}

/// The result of an evaluation of the integrand
#[derive(Clone, Serialize, Debug)]
pub struct EvaluationResult {
    /// Integrand value before any parameterization Jacobian is applied.
    pub integrand_result: Complex<F<f64>>,
    pub parameterization_jacobian: Option<F<f64>>,
    /// Monte Carlo sample weight supplied by the integrator/grid, excluding the parameterization Jacobian.
    pub integrator_weight: F<f64>,
    pub event_groups: EventGroupList,
    pub evaluation_metadata: EvaluationMetaData,
}

#[derive(Clone, Serialize, Debug)]
pub struct EvaluationResultOutput {
    pub integrand_result: Complex<F<f64>>,
    pub parameterization_jacobian: Option<F<f64>>,
    pub integrator_weight: F<f64>,
    pub event_groups: EventGroupList,
    pub evaluation_metadata: Option<EvaluationMetaData>,
}

#[derive(Clone, Debug)]
pub struct GenericEvaluationResult<T: FloatLike> {
    pub integrand_result: Complex<F<T>>,
    pub parameterization_jacobian: Option<F<T>>,
    pub integrator_weight: F<T>,
    pub event_groups: GenericEventGroupList<T>,
    pub evaluation_metadata: EvaluationMetaData,
}

#[derive(Clone, Debug)]
pub struct GenericEvaluationResultOutput<T: FloatLike> {
    pub integrand_result: Complex<F<T>>,
    pub parameterization_jacobian: Option<F<T>>,
    pub integrator_weight: F<T>,
    pub event_groups: GenericEventGroupList<T>,
    pub evaluation_metadata: Option<EvaluationMetaData>,
}

impl<T: FloatLike> GenericEvaluationResult<T> {
    pub fn into_output(self, minimal_output: bool) -> GenericEvaluationResultOutput<T> {
        GenericEvaluationResultOutput {
            integrand_result: self.integrand_result,
            parameterization_jacobian: self.parameterization_jacobian,
            integrator_weight: self.integrator_weight,
            event_groups: self.event_groups,
            evaluation_metadata: (!minimal_output).then_some(self.evaluation_metadata),
        }
    }
}

#[derive(Clone, Debug)]
pub enum PreciseEvaluationResultOutput {
    Double(GenericEvaluationResultOutput<f64>),
    Quad(GenericEvaluationResultOutput<f128>),
    Arb(GenericEvaluationResultOutput<ArbPrec>),
}

#[derive(Clone, Debug)]
pub enum PreciseEvaluationResult {
    Double(GenericEvaluationResult<f64>),
    Quad(GenericEvaluationResult<f128>),
    Arb(GenericEvaluationResult<ArbPrec>),
}

impl PreciseEvaluationResult {
    pub fn into_output(self, minimal_output: bool) -> PreciseEvaluationResultOutput {
        match self {
            PreciseEvaluationResult::Double(result) => {
                PreciseEvaluationResultOutput::Double(result.into_output(minimal_output))
            }
            PreciseEvaluationResult::Quad(result) => {
                PreciseEvaluationResultOutput::Quad(result.into_output(minimal_output))
            }
            PreciseEvaluationResult::Arb(result) => {
                PreciseEvaluationResultOutput::Arb(result.into_output(minimal_output))
            }
        }
    }
}

impl PreciseEvaluationResultOutput {
    pub fn precision(&self) -> Precision {
        match self {
            PreciseEvaluationResultOutput::Double(_) => Precision::Double,
            PreciseEvaluationResultOutput::Quad(_) => Precision::Quad,
            PreciseEvaluationResultOutput::Arb(_) => Precision::Arb,
        }
    }

    pub fn evaluation_metadata(&self) -> Option<&EvaluationMetaData> {
        match self {
            PreciseEvaluationResultOutput::Double(result) => result.evaluation_metadata.as_ref(),
            PreciseEvaluationResultOutput::Quad(result) => result.evaluation_metadata.as_ref(),
            PreciseEvaluationResultOutput::Arb(result) => result.evaluation_metadata.as_ref(),
        }
    }
}

impl EvaluationResult {
    pub fn zero() -> Self {
        Self {
            integrand_result: Complex::new_zero(),
            parameterization_jacobian: None,
            integrator_weight: F(0.0),
            event_groups: EventGroupList::default(),
            evaluation_metadata: EvaluationMetaData::new_empty(),
        }
    }

    pub fn into_output(self, minimal_output: bool) -> EvaluationResultOutput {
        EvaluationResultOutput {
            integrand_result: self.integrand_result,
            parameterization_jacobian: self.parameterization_jacobian,
            integrator_weight: self.integrator_weight,
            event_groups: self.event_groups,
            evaluation_metadata: (!minimal_output).then_some(self.evaluation_metadata),
        }
    }
}

fn fmt_evaluation_result_output<T: FloatLike>(
    f: &mut std::fmt::Formatter<'_>,
    precision: Option<Precision>,
    integrand_result: &Complex<F<T>>,
    parameterization_jacobian: Option<&F<T>>,
    integrator_weight: &F<T>,
    event_groups: &GenericEventGroupList<T>,
    evaluation_metadata: Option<&EvaluationMetaData>,
) -> std::fmt::Result {
    let mut summary_rows = Vec::with_capacity(5);
    if let Some(precision) = precision {
        summary_rows.push(EvaluationSummaryRow {
            field: "precision".to_string(),
            value: precision.to_string(),
        });
    }
    summary_rows.extend([
        EvaluationSummaryRow {
            field: "integrand result".to_string(),
            value: format_complex_generic(integrand_result),
        },
        EvaluationSummaryRow {
            field: "parameterization jacobian".to_string(),
            value: format_optional_real_generic(parameterization_jacobian),
        },
        EvaluationSummaryRow {
            field: "integrator weight".to_string(),
            value: format_real_generic(integrator_weight),
        },
        EvaluationSummaryRow {
            field: "event groups".to_string(),
            value: format_count(event_groups.len()),
        },
    ]);

    writeln!(f, "{}", "Evaluation result".bold().bright_green())?;
    writeln!(f, "{}", Table::new(summary_rows).with(Style::rounded()))?;
    if let Some(metadata) = evaluation_metadata {
        writeln!(f)?;
        write!(f, "{metadata}")?;
    }

    if !event_groups.is_empty() {
        writeln!(f)?;
        writeln!(
            f,
            "{}",
            format!(
                "Generated {} event group(s)",
                format_count(event_groups.len())
            )
            .bold()
            .bright_cyan()
        )?;
        write!(f, "{event_groups}")?;
    }

    Ok(())
}

fn fmt_sample_result<D: Display>(
    f: &mut std::fmt::Formatter<'_>,
    evaluation: &D,
    observables: &ObservableSnapshotBundle,
) -> std::fmt::Result {
    write!(f, "{evaluation}")?;

    if let Some(observables_table) = summarize_observables(observables) {
        writeln!(f)?;
        writeln!(f)?;
        writeln!(f, "{}", "Observable snapshots".bold().bright_magenta())?;
        write!(f, "{observables_table}")?;
    }

    Ok(())
}

fn fmt_batch_result<D: Display>(
    f: &mut std::fmt::Formatter<'_>,
    samples: &[D],
    observables: &ObservableSnapshotBundle,
) -> std::fmt::Result {
    let summary_rows = [EvaluationSummaryRow {
        field: "samples".to_string(),
        value: format_count(samples.len()),
    }];

    writeln!(f, "{}", "Batch evaluation result".bold().bright_green())?;
    writeln!(f, "{}", Table::new(summary_rows).with(Style::rounded()))?;

    if let Some(observables_table) = summarize_observables(observables) {
        writeln!(f)?;
        writeln!(f, "{}", "Observable snapshots".bold().bright_magenta())?;
        writeln!(f, "{observables_table}")?;
    }

    for (sample_index, sample) in samples.iter().enumerate() {
        writeln!(f)?;
        if sample_index > 0 {
            writeln!(f)?;
        }
        writeln!(
            f,
            "{}",
            format!("Sample {sample_index}").bold().bright_cyan()
        )?;
        write!(f, "{sample}")?;
    }

    Ok(())
}

impl Display for EvaluationResultOutput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_evaluation_result_output(
            f,
            None,
            &self.integrand_result,
            self.parameterization_jacobian.as_ref(),
            &self.integrator_weight,
            &self.event_groups,
            self.evaluation_metadata.as_ref(),
        )
    }
}

impl<T: FloatLike> Display for GenericEvaluationResultOutput<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_evaluation_result_output(
            f,
            None,
            &self.integrand_result,
            self.parameterization_jacobian.as_ref(),
            &self.integrator_weight,
            &self.event_groups,
            self.evaluation_metadata.as_ref(),
        )
    }
}

impl Display for PreciseEvaluationResultOutput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PreciseEvaluationResultOutput::Double(result) => fmt_evaluation_result_output(
                f,
                Some(Precision::Double),
                &result.integrand_result,
                result.parameterization_jacobian.as_ref(),
                &result.integrator_weight,
                &result.event_groups,
                result.evaluation_metadata.as_ref(),
            ),
            PreciseEvaluationResultOutput::Quad(result) => fmt_evaluation_result_output(
                f,
                Some(Precision::Quad),
                &result.integrand_result,
                result.parameterization_jacobian.as_ref(),
                &result.integrator_weight,
                &result.event_groups,
                result.evaluation_metadata.as_ref(),
            ),
            PreciseEvaluationResultOutput::Arb(result) => fmt_evaluation_result_output(
                f,
                Some(Precision::Arb),
                &result.integrand_result,
                result.parameterization_jacobian.as_ref(),
                &result.integrator_weight,
                &result.event_groups,
                result.evaluation_metadata.as_ref(),
            ),
        }
    }
}

impl Display for SampleEvaluationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.evaluation)
    }
}

#[derive(Clone, Serialize, Debug)]
pub struct SampleEvaluationResult {
    pub evaluation: EvaluationResultOutput,
}

#[derive(Clone, Debug)]
pub struct PreciseSampleEvaluationResult {
    pub evaluation: PreciseEvaluationResultOutput,
}

#[derive(Clone, Serialize, Debug)]
pub struct SingleSampleEvaluationResult {
    pub sample: SampleEvaluationResult,
    pub observables: ObservableSnapshotBundle,
}

#[derive(Clone, Serialize, Debug)]
pub struct BatchSampleEvaluationResult {
    pub samples: Vec<SampleEvaluationResult>,
    pub observables: ObservableSnapshotBundle,
}

#[derive(Clone, Debug)]
pub struct PreciseSingleSampleEvaluationResult {
    pub sample: PreciseSampleEvaluationResult,
    pub observables: ObservableSnapshotBundle,
}

#[derive(Clone, Debug)]
pub struct PreciseBatchSampleEvaluationResult {
    pub samples: Vec<PreciseSampleEvaluationResult>,
    pub observables: ObservableSnapshotBundle,
}

impl Display for PreciseSampleEvaluationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.evaluation)
    }
}

impl Display for SingleSampleEvaluationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_sample_result(f, &self.sample.evaluation, &self.observables)
    }
}

impl Display for BatchSampleEvaluationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_batch_result(f, &self.samples, &self.observables)
    }
}

impl Display for PreciseSingleSampleEvaluationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_sample_result(f, &self.sample.evaluation, &self.observables)
    }
}

impl Display for PreciseBatchSampleEvaluationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_batch_result(f, &self.samples, &self.observables)
    }
}

#[derive(Tabled)]
struct EvaluationSummaryRow {
    field: String,
    value: String,
}

#[derive(Tabled)]
struct HistogramSummaryRow {
    name: String,
    bins: usize,
    #[tabled(rename = "in-range")]
    in_range_entries: String,
    phase: String,
    range: String,
    underflow: String,
    overflow: String,
}

#[derive(Tabled)]
struct StabilitySummaryRow {
    level: String,
    relative_accuracy: String,
    time: String,
    status: String,
}

fn format_duration(duration: Duration) -> String {
    format_evaluation_time(duration)
}

fn format_count(value: usize) -> String {
    if value < 1_000 {
        return value.to_string();
    }

    let value = value as f64;
    for (scale, suffix) in [
        (1_000_000_000_f64, "B"),
        (1_000_000_f64, "M"),
        (1_000_f64, "K"),
    ] {
        if value >= scale {
            let scaled = value / scale;
            let precision = if scaled >= 100.0 {
                0
            } else if scaled >= 10.0 {
                1
            } else {
                2
            };
            return format!("{scaled:.precision$}{suffix}");
        }
    }

    value.round().to_string()
}

fn format_percentage(value: f64, significant_digits: usize) -> String {
    if !value.is_finite() {
        return "None".red().to_string();
    }

    if value == 0.0 {
        let decimals = significant_digits.saturating_sub(1);
        return format!("{:.*}%", decimals, 0.0);
    }

    let abs_value = value.abs();
    let exponent = abs_value.log10().floor() as i32;
    if exponent < -2 || exponent >= significant_digits as i32 {
        return format!("{:.*e}%", significant_digits.saturating_sub(1), value);
    }

    let decimals = (significant_digits as i32 - exponent - 1).max(0) as usize;
    format!("{value:.decimals$}%")
}

fn summarize_observables(observables: &ObservableSnapshotBundle) -> Option<String> {
    if observables.histograms.is_empty() {
        return None;
    }

    let rows = observables
        .histograms
        .iter()
        .map(|(name, histogram)| HistogramSummaryRow {
            name: name.clone(),
            bins: histogram.bins.len(),
            in_range_entries: format_count(histogram.statistics.in_range_entry_count),
            phase: format!("{:?}", histogram.phase).to_lowercase(),
            range: format!("[{:+.16e}, {:+.16e}]", histogram.x_min, histogram.x_max),
            underflow: format_count(histogram.underflow_bin.entry_count),
            overflow: format_count(histogram.overflow_bin.entry_count),
        })
        .collect::<Vec<_>>();

    Some(Table::new(rows).with(Style::rounded()).to_string())
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
pub struct StabilityResult {
    pub precision: Precision,
    pub estimated_relative_accuracy: Option<F<f64>>,
    pub accepted_as_stable: bool,
    pub total_time: Duration,
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
    pub event_processing_time: Duration,
    pub generated_event_count: usize,
    pub accepted_event_count: usize,
    pub relative_instability_error: Complex<F<f64>>,
    pub is_nan: bool,
    pub loop_momenta_escalation: Option<LoopMomentaEscalationMetrics>,
    pub stability_results: Vec<StabilityResult>,
}

impl Display for EvaluationMetaData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let rows = vec![
            EvaluationSummaryRow {
                field: "nan".to_string(),
                value: self.is_nan.to_string(),
            },
            EvaluationSummaryRow {
                field: "generated events".to_string(),
                value: format_count(self.generated_event_count),
            },
            EvaluationSummaryRow {
                field: "accepted events".to_string(),
                value: format_count(self.accepted_event_count),
            },
            EvaluationSummaryRow {
                field: "parameterization time".to_string(),
                value: format_duration(self.parameterization_time),
            },
            EvaluationSummaryRow {
                field: "integrand evaluation time".to_string(),
                value: format_duration(self.integrand_evaluation_time),
            },
            EvaluationSummaryRow {
                field: "evaluator evaluation time".to_string(),
                value: format_duration(self.evaluator_evaluation_time),
            },
            EvaluationSummaryRow {
                field: "event processing time".to_string(),
                value: format_duration(self.event_processing_time),
            },
            EvaluationSummaryRow {
                field: "total evaluation time".to_string(),
                value: format_duration(self.total_timing),
            },
        ];
        writeln!(f, "{}", "Evaluation metadata".bold().bright_yellow())?;
        writeln!(f, "{}", Table::new(rows).with(Style::rounded()))?;

        if !self.stability_results.is_empty() {
            let stability_rows = self
                .stability_results
                .iter()
                .map(|result| StabilitySummaryRow {
                    level: result.precision.to_string(),
                    relative_accuracy: format_optional_real_generic(
                        result.estimated_relative_accuracy.as_ref(),
                    ),
                    time: format_duration(result.total_time),
                    status: if result.accepted_as_stable {
                        "stable".green().to_string()
                    } else {
                        "unstable".red().to_string()
                    },
                })
                .collect::<Vec<_>>();
            writeln!(f)?;
            writeln!(f, "{}", "Stability results".bold().bright_magenta())?;
            write!(f, "{}", Table::new(stability_rows).with(Style::rounded()))?;
        }

        Ok(())
    }
}

impl EvaluationMetaData {
    pub(crate) fn new_empty() -> Self {
        Self {
            total_timing: Duration::ZERO,
            integrand_evaluation_time: Duration::ZERO,
            evaluator_evaluation_time: Duration::ZERO,
            parameterization_time: Duration::ZERO,
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
            relative_instability_error: Complex::new_zero(),
            is_nan: false,
            loop_momenta_escalation: None,
            stability_results: Vec::new(),
        }
    }

    pub(crate) fn final_precision(&self) -> Option<Precision> {
        self.stability_results.last().map(|result| result.precision)
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
    num_nan_or_unstable_evals: usize,
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
                accumulator.sum_event_time += data_entry.evaluation_metadata.event_processing_time;
                accumulator.sum_relative_instability_error.0 +=
                    data_entry.evaluation_metadata.relative_instability_error.re;
                accumulator.sum_relative_instability_error.1 +=
                    data_entry.evaluation_metadata.relative_instability_error.im;
                accumulator.sum_total_evaluation_time +=
                    data_entry.evaluation_metadata.total_timing;
                accumulator.sum_generated_event_count +=
                    data_entry.evaluation_metadata.generated_event_count;
                accumulator.sum_accepted_event_count +=
                    data_entry.evaluation_metadata.accepted_event_count;

                accumulator.num_evals += 1;
                match data_entry
                    .evaluation_metadata
                    .final_precision()
                    .unwrap_or(Precision::Double)
                {
                    Precision::Double => accumulator.num_double_precision_evals += 1,
                    Precision::Quad => accumulator.num_quadruple_precision_evals += 1,
                    Precision::Arb => accumulator.num_arb_precision_evals += 1,
                    // _ => (),
                }

                if data_entry.evaluation_metadata.is_nan {
                    accumulator.num_nan_evals += 1;
                }
                if data_entry.evaluation_metadata.is_nan
                    || data_entry
                        .evaluation_metadata
                        .stability_results
                        .last()
                        .map(|result| !result.accepted_as_stable)
                        .unwrap_or(false)
                {
                    accumulator.num_nan_or_unstable_evals += 1;
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
            num_nan_or_unstable_evals: self.num_nan_or_unstable_evals
                + other.num_nan_or_unstable_evals,
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
            num_nan_or_unstable_evals: 0,
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

    pub(crate) fn get_avg_event_timing(&self) -> Duration {
        let avg_event_timing = self.sum_event_time.as_secs_f64() / self.num_evals as f64;
        Duration::from_secs_f64(avg_event_timing)
    }

    /// Get the average relative error computed during instability checks
    #[allow(dead_code)]
    pub(crate) fn get_avg_instabillity_error(&self) -> (F<f64>, F<f64>) {
        (
            self.sum_relative_instability_error.0
                / self
                    .sum_relative_instability_error
                    .0
                    .from_usize(self.num_evals),
            self.sum_relative_instability_error.1
                / self
                    .sum_relative_instability_error
                    .1
                    .from_usize(self.num_evals),
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

    pub(crate) fn get_percentage_nan_or_unstable(&self) -> f64 {
        self.num_nan_or_unstable_evals as f64 / self.num_evals as f64 * 100.0
    }

    pub(crate) fn selection_efficiency_percentage(&self) -> Option<f64> {
        if self.sum_generated_event_count == 0 {
            None
        } else {
            Some(
                self.sum_accepted_event_count as f64 / self.sum_generated_event_count as f64
                    * 100.0,
            )
        }
    }

    #[allow(clippy::format_in_format_args)]
    pub(crate) fn display_status(&self) {
        let time_integrand_formatted = format_evaluation_time(self.get_avg_integrand_timing());
        let time_evaluators_formatted = format_evaluation_time(self.get_avg_evaluator_timing());
        let param_time_formatted = format_evaluation_time(self.get_avg_param_timing());
        let event_time_formatted = format_evaluation_time(self.get_avg_event_timing());
        let total_time = format_evaluation_time(self.get_avg_total_timing());
        let selection_efficiency = self.selection_efficiency_percentage();
        let selection_efficiency_display = selection_efficiency
            .map(|value| format_percentage(value, 3).green().to_string())
            .unwrap_or_else(|| "None".red().to_string());
        let nan_or_unstable = self.get_percentage_nan_or_unstable();
        let nan_or_unstable_display = if nan_or_unstable > 0.0 {
            format_percentage(nan_or_unstable, 2).red().to_string()
        } else {
            format_percentage(nan_or_unstable, 2).green().to_string()
        };

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

        info!(
            "|  {}  | {} {} | {} {} | {} {} | {} {}",
            format!("{:-7}", "events").blue().bold(),
            format!("{:-7}", "Evts:"),
            format!("{:-9}", format_count(self.sum_generated_event_count)).green(),
            format!("{:-7}", "Sel. %:"),
            format!("{:-9}", selection_efficiency_display),
            format!("{:-7}", "t. obs:"),
            format!("{:-9}", event_time_formatted).green(),
            format!("{:-16}", "NaNs/Unstable %:"),
            format!("{:-9}", nan_or_unstable_display),
        );
    }
}
