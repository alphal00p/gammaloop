use std::{
    fs,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use clap::Args;
use color_eyre::eyre::{eyre, Context};
use color_eyre::Result;
use colored::Colorize;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use tabled::{
    builder::Builder,
    settings::{
        object::{Columns, Object, Rows},
        style::HorizontalLine,
        Alignment, Modify, Style,
    },
};
use tracing::info;

use crate::{
    commands::evaluate_samples::{evaluate_sample, evaluate_samples, EvaluateSamples},
    completion::CompletionArgExt,
    state::{ProcessRef, State},
};
use gammalooprs::{
    integrands::evaluation::EvaluationResultOutput,
    settings::{
        runtime::{RotationSetting, StabilityLevelSetting},
        RuntimeSettings,
    },
    utils::F,
};

const INSPECT_BENCH_WARMUP_SAMPLES: usize = 10;
const INSPECT_BENCH_DEFAULT_BATCHES: usize = 10;
const INSPECT_BENCH_OTHER_RELATIVE_THRESHOLD: f64 = 5.0e-3;
const INSPECT_BENCH_OTHER_ABSOLUTE_THRESHOLD_SECONDS: f64 = 1.0e-9;

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Inspect {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,
    /// The integrand name to inspect
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,
    /// The point to inspect (x y) or (p0 px ...)
    #[arg(
        short = 'x',
        long = "point",
        num_args = 2..,
        value_name = "POINT",
        value_delimiter = ',',
        allow_negative_numbers = true,
    )]
    pub point: Vec<f64>,

    /// Evaluate in f128 precision
    #[arg(short = 'f', long = "use_arb_prec")]
    pub use_arb_prec: bool,

    /// Interpret point as momentum-space coordinates
    #[arg(short = 'm', long)]
    pub momentum_space: bool,

    /// The discrete dimensions of the sample
    #[arg(
        short = 'd',
        long = "discrete-dim",
        value_name = "DIMS",
        num_args = 1..,
        value_delimiter = ',',
        conflicts_with_all = ["graph_id", "orientation_id"],
    )]
    pub discrete_dim: Vec<usize>,

    /// Select a specific graph in momentum-space inspect
    #[arg(long = "graph-id", value_name = "GRAPH_ID")]
    pub graph_id: Option<usize>,

    /// Select a specific orientation of the selected graph in momentum-space inspect
    #[arg(
        long = "orientation-id",
        value_name = "ORIENTATION_ID",
        requires = "graph_id",
        conflicts_with = "discrete_dim"
    )]
    pub orientation_id: Option<usize>,

    /// Benchmark repeated inspect evaluations for approximately this duration in seconds.
    #[arg(long, value_name = "SECONDS")]
    pub bench: Option<String>,

    /// Number of batches used to estimate inspect benchmark uncertainty.
    #[arg(
        long = "n_batches",
        value_name = "N",
        requires = "bench",
        value_parser = parse_positive_usize
    )]
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub n_batches: Option<usize>,

    /// Temporarily benchmark only the raw integrand path by disabling stability rotations,
    /// caches, generated events, additional event weights, selectors, and observables.
    #[arg(long = "minimal-integrand", requires = "bench")]
    #[serde(default)]
    pub minimal_integrand: bool,

    /// Write the inspect or benchmark result as pretty JSON to this path.
    #[arg(long = "json-output", value_name = "PATH")]
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub json_output: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct InspectBenchBatchTiming {
    parameterization: f64,
    integrand: f64,
    event_processing: f64,
    evaluator: f64,
    other: f64,
    total: f64,
}

#[derive(Debug, Clone)]
struct InspectBenchRow {
    category: &'static str,
    timings: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct InspectBenchSummaryRow {
    category: String,
    mean_seconds_per_sample: f64,
    standard_error_seconds_per_sample: f64,
    percentage_of_total: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct InspectBenchOutput {
    process_id: usize,
    integrand_name: String,
    momentum_space: bool,
    point: Vec<f64>,
    discrete_dim: Vec<usize>,
    graph_name: Option<String>,
    orientation_id: Option<usize>,
    target_seconds: f64,
    warmup_samples: usize,
    warmup_seconds: f64,
    total_samples: usize,
    n_batches: usize,
    minimal_integrand: bool,
    parameterization_jacobian: Option<f64>,
    displayed_result: Complex<f64>,
    batches: Vec<InspectBenchBatchTiming>,
    summary: Vec<InspectBenchSummaryRow>,
}

impl InspectBenchRow {
    fn mean_and_standard_error(&self) -> (f64, f64) {
        mean_and_standard_error(&self.timings)
    }
}

impl Inspect {
    fn validate_selector_mode(&self) -> Result<()> {
        if !self.momentum_space && (self.graph_id.is_some() || self.orientation_id.is_some()) {
            return Err(eyre!(
                "Graph and orientation selectors are only supported in momentum-space inspect."
            ));
        }
        Ok(())
    }

    fn resolve_graph_name(
        &self,
        state: &mut State,
        process_id: usize,
        integrand_name: &str,
    ) -> Result<Option<String>> {
        let Some(graph_id) = self.graph_id else {
            return Ok(None);
        };
        let integrand = state
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        let graph_name = integrand.graph_name_by_id(graph_id).ok_or_else(|| {
            eyre!(
                "Graph id {} is out of range for integrand '{}'; it has {} graphs.",
                graph_id,
                integrand_name,
                integrand.graph_count()
            )
        })?;
        Ok(Some(graph_name.to_string()))
    }

    fn validate_x_space_point(
        &self,
        state: &mut State,
        process_id: usize,
        integrand_name: &str,
    ) -> Result<()> {
        let integrand = state
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        let expected_dimension = integrand.expected_x_space_dimension(&self.discrete_dim)?;
        if self.point.len() != expected_dimension {
            return Err(eyre!(
                "Expected {} x-space coordinates for this integrand selection, got {}.",
                expected_dimension,
                self.point.len()
            ));
        }
        Ok(())
    }

    fn repeated_point_array(&self, rows: usize) -> Result<Array2<f64>> {
        let mut values = Vec::with_capacity(rows * self.point.len());
        for _ in 0..rows {
            values.extend_from_slice(&self.point);
        }
        Array2::from_shape_vec((rows, self.point.len()), values).map_err(Into::into)
    }

    fn repeated_discrete_dim_array(&self, rows: usize) -> Result<Array2<usize>> {
        let mut values = Vec::with_capacity(rows * self.discrete_dim.len());
        for _ in 0..rows {
            values.extend_from_slice(&self.discrete_dim);
        }
        Array2::from_shape_vec((rows, self.discrete_dim.len()), values).map_err(Into::into)
    }

    fn evaluate_samples_request<'a>(
        &'a self,
        process_id: usize,
        integrand_name: String,
        points: &'a Array2<f64>,
        discrete_dims: &'a Array2<usize>,
        graph_name: &Option<String>,
        return_generated_events: bool,
    ) -> EvaluateSamples<'a> {
        let batch_len = points.nrows();
        EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name),
            use_arb_prec: self.use_arb_prec,
            minimal_output: false,
            return_generated_events: Some(return_generated_events),
            momentum_space: self.momentum_space,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete_dims.view()),
            graph_names: graph_name
                .as_ref()
                .map(|name| vec![Some(name.clone()); batch_len]),
            orientations: self
                .orientation_id
                .map(|orientation| vec![Some(orientation); batch_len]),
        }
    }

    fn displayed_result_from_evaluation(
        &self,
        evaluation: &EvaluationResultOutput,
    ) -> Result<(Option<f64>, Complex<f64>)> {
        let raw_result = evaluation.integrand_result;
        let jacobian = evaluation.parameterization_jacobian.map(|jac| jac.0);
        let displayed_result = if self.momentum_space {
            raw_result
        } else {
            let jacobian = jacobian.ok_or_else(|| {
                eyre!("x-space inspect requires a parameterization jacobian, but none was returned")
            })?;
            raw_result.map(|entry| entry * F(jacobian))
        };
        Ok((jacobian, displayed_result.map(|entry| entry.0)))
    }

    fn write_json_output<T: Serialize>(&self, value: &T) -> Result<()> {
        let Some(path) = self.json_output.as_deref() else {
            return Ok(());
        };
        let path = resolve_json_output_path(path)?;
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Could not create {}", parent.display()))?;
        }
        fs::write(&path, serde_json::to_vec_pretty(value)?)
            .with_context(|| format!("Could not write inspect JSON output {}", path.display()))?;
        info!(
            "Wrote inspect JSON output to {}",
            path.display().to_string().green()
        );
        Ok(())
    }

    fn run_bench(
        &self,
        state: &mut State,
        process_id: usize,
        integrand_name: String,
        graph_name: Option<String>,
    ) -> Result<(Option<f64>, Complex<f64>)> {
        let target = parse_bench_target_duration(
            self.bench
                .as_deref()
                .expect("bench target should be present on bench path"),
        )?;
        let n_batches = self.n_batches.unwrap_or(INSPECT_BENCH_DEFAULT_BATCHES);

        let warmup_points = self.repeated_point_array(INSPECT_BENCH_WARMUP_SAMPLES)?;
        let warmup_discrete_dims =
            self.repeated_discrete_dim_array(INSPECT_BENCH_WARMUP_SAMPLES)?;
        let warmup_start = Instant::now();
        let warmup = evaluate_samples(
            state,
            &self.evaluate_samples_request(
                process_id,
                integrand_name.clone(),
                &warmup_points,
                &warmup_discrete_dims,
                &graph_name,
                false,
            ),
        )?;
        let warmup_timing = warmup_start.elapsed();
        let warmup_samples = warmup.samples.len().max(1);
        let warmup_per_sample = warmup_timing.as_secs_f64() / warmup_samples as f64;
        let total_samples = if warmup_per_sample > 0.0 {
            ((target.as_secs_f64() / warmup_per_sample).ceil() as usize).max(n_batches)
        } else {
            n_batches
        };

        let batch_sizes = split_samples_into_batches(total_samples, n_batches);
        let progress_bar = inspect_bench_progress_bar(n_batches);
        progress_bar.set_message(
            format!(
                "{} {}",
                "running benchmark batches".cyan(),
                format!("0/{n_batches}").blue()
            )
            .to_string(),
        );

        let mut batch_timings = Vec::with_capacity(n_batches);
        for (batch_index, batch_size) in batch_sizes.into_iter().enumerate() {
            let points = self.repeated_point_array(batch_size)?;
            let discrete_dims = self.repeated_discrete_dim_array(batch_size)?;
            let batch_start = Instant::now();
            let batch = evaluate_samples(
                state,
                &self.evaluate_samples_request(
                    process_id,
                    integrand_name.clone(),
                    &points,
                    &discrete_dims,
                    &graph_name,
                    false,
                ),
            )?;
            let total_timing = batch_start.elapsed();
            batch_timings.push(inspect_bench_batch_timing(total_timing, &batch)?);
            update_inspect_bench_progress(
                &progress_bar,
                batch_index + 1,
                n_batches,
                &batch_timings,
                total_samples,
            );
        }
        progress_bar.finish_and_clear();
        drop(progress_bar);

        let warmup_evaluation = warmup
            .samples
            .first()
            .ok_or_else(|| eyre!("inspect --bench warmup returned no sample"))?;
        let (jacobian, displayed_result) =
            self.displayed_result_from_evaluation(&warmup_evaluation.evaluation)?;
        let bench_output = InspectBenchOutput {
            process_id,
            integrand_name: integrand_name.clone(),
            momentum_space: self.momentum_space,
            point: self.point.clone(),
            discrete_dim: self.discrete_dim.clone(),
            graph_name: graph_name.clone(),
            orientation_id: self.orientation_id,
            target_seconds: target.as_secs_f64(),
            warmup_samples,
            warmup_seconds: warmup_timing.as_secs_f64(),
            total_samples,
            n_batches,
            minimal_integrand: self.minimal_integrand,
            parameterization_jacobian: jacobian,
            displayed_result,
            batches: batch_timings.clone(),
            summary: inspect_bench_summary_rows(&batch_timings),
        };
        self.write_json_output(&bench_output)?;

        info!(
            "\n{}\n{}",
            format!(
                "Inspect benchmark for integrand '{}' at fixed input point: target {}, warmup {} sample(s), profiled {} sample(s) in {} batch(es).",
                integrand_name.green(),
                format_duration_smart(target).yellow(),
                warmup_samples.to_string().blue(),
                total_samples.to_string().blue(),
                n_batches.to_string().blue(),
            ),
            render_inspect_bench_table(&batch_timings)
        );

        Ok((jacobian, displayed_result))
    }

    pub fn run(&self, state: &mut State) -> Result<(Option<f64>, Complex<f64>)> {
        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;
        self.validate_selector_mode()?;
        if self.minimal_integrand && self.bench.is_none() {
            return Err(eyre!(
                "inspect --minimal-integrand is only supported together with --bench"
            ));
        }
        if !self.momentum_space {
            self.validate_x_space_point(state, process_id, &integrand_name)?;
        }
        let graph_name = self.resolve_graph_name(state, process_id, &integrand_name)?;
        if self.bench.is_some() {
            let integrand_name_for_restore = integrand_name.clone();
            return self.with_optional_minimal_integrand_settings(
                state,
                process_id,
                &integrand_name_for_restore,
                |state| self.run_bench(state, process_id, integrand_name, graph_name),
            );
        }
        let points = Array2::from_shape_vec((1, self.point.len()), self.point.clone())?;
        let discrete_dims =
            Array2::from_shape_vec((1, self.discrete_dim.len()), self.discrete_dim.clone())?;
        let result = evaluate_sample(
            state,
            &EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(integrand_name.clone()),
                use_arb_prec: self.use_arb_prec,
                minimal_output: false,
                return_generated_events: Some(true),
                momentum_space: self.momentum_space,
                points: points.view(),
                integrator_weights: None,
                discrete_dims: Some(discrete_dims.view()),
                graph_names: graph_name.map(|name| vec![Some(name)]),
                orientations: self
                    .orientation_id
                    .map(|orientation| vec![Some(orientation)]),
            },
        )?;
        self.write_json_output(&result.sample)?;
        let evaluation = &result.sample.evaluation;

        let raw_result = evaluation.integrand_result;
        let jacobian = evaluation.parameterization_jacobian.map(|jac| jac.0);
        let displayed_result = if self.momentum_space {
            raw_result
        } else {
            let jacobian = jacobian.ok_or_else(|| {
                eyre!("x-space inspect requires a parameterization jacobian, but none was returned")
            })?;
            raw_result.map(|entry| entry * F(jacobian))
        };

        let point_label = if self.momentum_space {
            "Input point in momentum space"
        } else {
            "Input point in unit hypercube xs"
        };
        info!(
            "\n{}:\n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
            point_label,
            format!(
                "( {} )",
                self.point
                    .iter()
                    .map(|x| format!("{x:.16}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
            .blue(),
            integrand_name.green(),
            format!(
                "( {:+.16e}, {:+.16e} i)",
                displayed_result.re, displayed_result.im
            )
            .blue(),
        );

        if let Some(jacobian) = jacobian {
            info!(
                "Parameterization jacobian for this point: {:+.16e}",
                jacobian
            );
        }
        info!("\n{}", result);

        Ok((jacobian, displayed_result.map(|entry| entry.0)))
    }

    fn with_optional_minimal_integrand_settings<R>(
        &self,
        state: &mut State,
        process_id: usize,
        integrand_name: &str,
        run: impl FnOnce(&mut State) -> Result<R>,
    ) -> Result<R> {
        if !self.minimal_integrand {
            return run(state);
        }

        let original_settings = {
            let integrand = state
                .process_list
                .get_integrand_mut(process_id, integrand_name)?;
            let original_settings = integrand.get_settings().clone();
            apply_minimal_integrand_settings(integrand.get_mut_settings());
            original_settings
        };

        let result = run(state);
        let integrand = state
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        *integrand.get_mut_settings() = original_settings;
        result
    }
}

fn apply_minimal_integrand_settings(settings: &mut RuntimeSettings) {
    settings.general.enable_cache = false;
    settings.general.debug_cache = false;
    settings.general.generate_events = false;
    settings.general.store_additional_weights_in_event = false;
    settings.observables = Default::default();
    settings.selectors = Default::default();
    settings.stability.rotation_axis = vec![RotationSetting::None {}];
    settings.stability.levels = vec![StabilityLevelSetting::default_double()];
    settings.stability.check_on_norm = false;
    settings.stability.escalate_if_exact_zero = false;
    settings.stability.loop_momenta_norm_escalation_factor = -1.0;
    settings.stability.recording = None;
}

fn parse_bench_target_duration(input: &str) -> Result<Duration> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(eyre!("inspect --bench requires a non-empty duration"));
    }

    let lower = trimmed.to_ascii_lowercase();
    let (number, multiplier) = if let Some(number) = lower.strip_suffix("microseconds") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("microsecond") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("usec") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("us") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("µs") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("milliseconds") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("millisecond") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("msec") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("ms") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("seconds") {
        (number, 1.0)
    } else if let Some(number) = lower.strip_suffix("second") {
        (number, 1.0)
    } else if let Some(number) = lower.strip_suffix("sec") {
        (number, 1.0)
    } else if let Some(number) = lower.strip_suffix('s') {
        (number, 1.0)
    } else {
        (trimmed, 1.0)
    };
    let seconds = number
        .trim()
        .parse::<f64>()
        .with_context(|| format!("Could not parse inspect --bench duration '{input}'"))?
        * multiplier;
    if !seconds.is_finite() || seconds <= 0.0 {
        return Err(eyre!(
            "inspect --bench duration must be positive and finite"
        ));
    }
    Ok(Duration::from_secs_f64(seconds))
}

fn parse_positive_usize(input: &str) -> std::result::Result<usize, String> {
    let value = input
        .parse::<usize>()
        .map_err(|err| format!("expected a positive integer: {err}"))?;
    if value == 0 {
        return Err("expected a positive integer".to_string());
    }
    Ok(value)
}

fn resolve_json_output_path(path: &Path) -> Result<PathBuf> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}

fn split_samples_into_batches(total_samples: usize, batches: usize) -> Vec<usize> {
    let base = total_samples / batches;
    let remainder = total_samples % batches;
    (0..batches)
        .map(|index| base + usize::from(index < remainder))
        .collect()
}

fn inspect_bench_batch_timing(
    total_timing: Duration,
    batch: &gammalooprs::integrands::evaluation::BatchSampleEvaluationResult,
) -> Result<InspectBenchBatchTiming> {
    let sample_count = batch.samples.len();
    if sample_count == 0 {
        return Err(eyre!("inspect --bench batch returned no samples"));
    }

    let mut parameterization = 0.0;
    let mut integrand_inclusive = 0.0;
    let mut event_processing = 0.0;
    let mut evaluator = 0.0;
    for sample in &batch.samples {
        let metadata = sample
            .evaluation
            .evaluation_metadata
            .as_ref()
            .ok_or_else(|| eyre!("inspect --bench requires evaluation metadata"))?;
        parameterization += metadata.parameterization_time.as_secs_f64();
        integrand_inclusive += metadata.integrand_evaluation_time.as_secs_f64();
        event_processing += metadata.event_processing_time.as_secs_f64();
        evaluator += metadata.evaluator_evaluation_time.as_secs_f64();
    }

    let inv_sample_count = 1.0 / sample_count as f64;
    let parameterization = parameterization * inv_sample_count;
    let evaluator = evaluator * inv_sample_count;
    // integrand_evaluation_time is inclusive of evaluator calls. Keep the
    // displayed rows disjoint by showing the non-evaluator integrand residual.
    let integrand = (integrand_inclusive * inv_sample_count - evaluator).max(0.0);
    let event_processing = event_processing * inv_sample_count;
    // Keep Total as the wall-clock time around evaluate_samples. The metadata
    // rows below are subtracted from that wrapper timing so other/overhead
    // captures unclassified evaluation work plus evaluate_samples overhead.
    let total = total_timing.as_secs_f64() * inv_sample_count;
    let known = parameterization + integrand + event_processing + evaluator;
    let other = (total - known).max(0.0);

    Ok(InspectBenchBatchTiming {
        parameterization,
        integrand,
        event_processing,
        evaluator,
        other,
        total,
    })
}

fn inspect_bench_summary_rows(
    batch_timings: &[InspectBenchBatchTiming],
) -> Vec<InspectBenchSummaryRow> {
    let rows = inspect_bench_rows(batch_timings);
    let total_values = batch_timings
        .iter()
        .map(|timing| timing.total)
        .collect::<Vec<_>>();
    let (total_mean, total_standard_error) = mean_and_standard_error(&total_values);

    let mut summary = rows
        .into_iter()
        .map(|row| {
            let (mean, standard_error) = row.mean_and_standard_error();
            InspectBenchSummaryRow {
                category: row.category.to_string(),
                mean_seconds_per_sample: mean,
                standard_error_seconds_per_sample: standard_error,
                percentage_of_total: (total_mean > 0.0).then_some(mean / total_mean * 100.0),
            }
        })
        .collect::<Vec<_>>();
    summary.push(InspectBenchSummaryRow {
        category: "Total".to_string(),
        mean_seconds_per_sample: total_mean,
        standard_error_seconds_per_sample: total_standard_error,
        percentage_of_total: (total_mean > 0.0).then_some(100.0),
    });
    summary
}

fn inspect_bench_progress_bar(n_batches: usize) -> ProgressBar {
    let bar = ProgressBar::new(n_batches as u64);
    let style = ProgressStyle::with_template(
        "{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len}",
    )
    .unwrap_or_else(|_| ProgressStyle::default_bar());
    bar.set_style(style);
    bar
}

fn update_inspect_bench_progress(
    progress_bar: &ProgressBar,
    completed_batches: usize,
    n_batches: usize,
    batch_timings: &[InspectBenchBatchTiming],
    total_samples: usize,
) {
    let total_values = batch_timings
        .iter()
        .map(|timing| timing.total)
        .collect::<Vec<_>>();
    let (total_mean, total_standard_error) = mean_and_standard_error(&total_values);
    let projected_total = Duration::from_secs_f64((total_mean * total_samples as f64).max(0.0));
    let message = format!(
        "{} {} | {} {} | {} {}",
        "batch".cyan(),
        format!("{completed_batches}/{n_batches}").blue(),
        "total/sample".cyan(),
        format_timing_with_uncertainty(total_mean, total_standard_error).yellow(),
        "projected total".cyan(),
        format_duration_smart(projected_total).yellow(),
    );
    progress_bar.set_position(completed_batches as u64);
    progress_bar.set_message(message);
}

fn inspect_bench_rows(batch_timings: &[InspectBenchBatchTiming]) -> Vec<InspectBenchRow> {
    let mut rows = vec![
        InspectBenchRow {
            category: "parameterization",
            timings: batch_timings
                .iter()
                .map(|timing| timing.parameterization)
                .collect(),
        },
        InspectBenchRow {
            category: "integrand",
            timings: batch_timings
                .iter()
                .map(|timing| timing.integrand)
                .collect(),
        },
        InspectBenchRow {
            category: "event processing",
            timings: batch_timings
                .iter()
                .map(|timing| timing.event_processing)
                .collect(),
        },
    ];

    let total_values = batch_timings
        .iter()
        .map(|timing| timing.total)
        .collect::<Vec<_>>();
    let (total_mean, _) = mean_and_standard_error(&total_values);
    let other_values = batch_timings
        .iter()
        .map(|timing| timing.other)
        .collect::<Vec<_>>();
    let (other_mean, _) = mean_and_standard_error(&other_values);
    if other_mean
        > (total_mean * INSPECT_BENCH_OTHER_RELATIVE_THRESHOLD)
            .max(INSPECT_BENCH_OTHER_ABSOLUTE_THRESHOLD_SECONDS)
    {
        rows.push(InspectBenchRow {
            category: "other/overhead",
            timings: other_values,
        });
    }

    rows.push(InspectBenchRow {
        category: "evaluators",
        timings: batch_timings
            .iter()
            .map(|timing| timing.evaluator)
            .collect(),
    });
    rows
}

fn render_inspect_bench_table(batch_timings: &[InspectBenchBatchTiming]) -> String {
    let rows = inspect_bench_rows(batch_timings);
    let total_values = batch_timings
        .iter()
        .map(|timing| timing.total)
        .collect::<Vec<_>>();
    let (total_mean, _) = mean_and_standard_error(&total_values);

    let total_row = InspectBenchRow {
        category: "Total",
        timings: total_values,
    };
    let total_record_index = rows.len() + 1;

    let mut builder = Builder::new();
    builder.push_record([
        "category".bold().cyan().to_string(),
        "time / sample".bold().cyan().to_string(),
        "% total".bold().cyan().to_string(),
    ]);
    for row in &rows {
        let (mean, standard_error) = row.mean_and_standard_error();
        builder.push_record([
            row.category.cyan().to_string(),
            format_timing_with_uncertainty(mean, standard_error)
                .yellow()
                .to_string(),
            format_percentage_cell(mean, total_mean).blue().to_string(),
        ]);
    }

    let (total_mean, total_standard_error) = total_row.mean_and_standard_error();
    builder.push_record([
        total_row.category.bold().green().to_string(),
        format_timing_with_uncertainty(total_mean, total_standard_error)
            .bold()
            .green()
            .to_string(),
        "100.0%".bold().green().to_string(),
    ]);

    let mut table = builder.build();
    table.with(
        Style::rounded().horizontals([
            (
                1,
                HorizontalLine::new('─')
                    .intersection('┼')
                    .left('├')
                    .right('┤'),
            ),
            (
                total_record_index,
                HorizontalLine::new('─')
                    .intersection('┼')
                    .left('├')
                    .right('┤'),
            ),
        ]),
    );
    table.with(Modify::new(Rows::first()).with(Alignment::center()));
    table.with(Modify::new(Rows::new(1..).intersect(Columns::one(0))).with(Alignment::left()));
    table.with(Modify::new(Rows::new(1..).intersect(Columns::new(1..))).with(Alignment::right()));
    table.to_string()
}

fn mean_and_standard_error(values: &[f64]) -> (f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    if values.len() < 2 {
        return (mean, 0.0);
    }
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    (mean, variance.sqrt() / (values.len() as f64).sqrt())
}

fn format_percentage_cell(value: f64, total: f64) -> String {
    if total <= 0.0 {
        return "-".to_string();
    }
    format!("{:.1}%", value / total * 100.0)
}

fn format_duration_smart(duration: Duration) -> String {
    let (scaled, unit) = scale_seconds_for_display(duration.as_secs_f64());
    format!("{} {unit}", format_significant(scaled, 3))
}

fn format_timing_with_uncertainty(mean_seconds: f64, standard_error_seconds: f64) -> String {
    let (scaled_mean, unit) = scale_seconds_for_display(mean_seconds);
    let scale = unit_scale(unit);
    let scaled_error = standard_error_seconds * scale;
    format!(
        "{}±{} {unit}",
        format_significant(scaled_mean, 3),
        format_significant(scaled_error, 2)
    )
}

fn scale_seconds_for_display(seconds: f64) -> (f64, &'static str) {
    let abs_seconds = seconds.abs();
    if abs_seconds < 1.0e-6 {
        (seconds * 1.0e9, "ns")
    } else if abs_seconds < 1.0e-3 {
        (seconds * 1.0e6, "µs")
    } else if abs_seconds < 1.0 {
        (seconds * 1.0e3, "ms")
    } else {
        (seconds, "s")
    }
}

fn unit_scale(unit: &str) -> f64 {
    match unit {
        "ns" => 1.0e9,
        "µs" => 1.0e6,
        "ms" => 1.0e3,
        "s" => 1.0,
        _ => 1.0,
    }
}

fn format_significant(value: f64, significant_digits: usize) -> String {
    if !value.is_finite() {
        return value.to_string();
    }
    if value == 0.0 {
        return format!("{:.*}", significant_digits.saturating_sub(1), 0.0);
    }
    let decimals = (significant_digits as i32 - 1 - value.abs().log10().floor() as i32).max(0);
    format!("{:.*}", decimals as usize, value)
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use gammalooprs::{
        integrands::evaluation::{
            BatchSampleEvaluationResult, EvaluationResult, SampleEvaluationResult,
        },
        observables::ObservableSnapshotBundle,
        settings::runtime::{RotationSetting, StabilityRecordingSettings},
    };

    use super::{
        apply_minimal_integrand_settings, format_timing_with_uncertainty,
        inspect_bench_batch_timing, parse_bench_target_duration, split_samples_into_batches,
        RuntimeSettings, StabilityLevelSetting,
    };

    #[test]
    fn inspect_bench_duration_parser_accepts_common_units() {
        assert_eq!(
            parse_bench_target_duration("5").unwrap(),
            Duration::from_secs(5)
        );
        assert_eq!(
            parse_bench_target_duration("250ms").unwrap(),
            Duration::from_millis(250)
        );
        assert_eq!(
            parse_bench_target_duration("10 us").unwrap(),
            Duration::from_micros(10)
        );
    }

    #[test]
    fn inspect_bench_batches_cover_all_samples() {
        let batches = split_samples_into_batches(23, 10);
        assert_eq!(batches.len(), 10);
        assert_eq!(batches.iter().sum::<usize>(), 23);
        assert_eq!(batches[0], 3);
        assert_eq!(batches[9], 2);
    }

    #[test]
    fn inspect_bench_timing_formatter_uses_shared_unit() {
        assert_eq!(
            format_timing_with_uncertainty(1.234e-6, 4.5e-8),
            "1.23±0.045 µs"
        );
    }

    #[test]
    fn inspect_bench_table_includes_total_row_after_separator() {
        let table = super::render_inspect_bench_table(&[
            super::InspectBenchBatchTiming {
                parameterization: 1.0e-7,
                integrand: 4.0e-4,
                event_processing: 0.0,
                evaluator: 1.7e-2,
                other: 1.0e-5,
                total: 1.74101e-2,
            },
            super::InspectBenchBatchTiming {
                parameterization: 1.2e-7,
                integrand: 4.2e-4,
                event_processing: 0.0,
                evaluator: 1.8e-2,
                other: 1.0e-5,
                total: 1.843012e-2,
            },
        ]);

        assert!(table.contains("│ Total"));
        assert!(table.contains("100.0%"));
        assert!(table.lines().any(|line| line.starts_with('├')));
        assert!(table
            .lines()
            .last()
            .is_some_and(|line| line.starts_with('╰')));
    }

    #[test]
    fn inspect_bench_timing_uses_wrapper_total_and_metadata_components() {
        let batch = BatchSampleEvaluationResult {
            samples: vec![
                sample_with_metadata(
                    Duration::from_micros(10),
                    Duration::from_micros(1),
                    Duration::from_micros(5),
                    Duration::from_micros(3),
                    Duration::from_micros(1),
                ),
                sample_with_metadata(
                    Duration::from_micros(14),
                    Duration::from_micros(2),
                    Duration::from_micros(6),
                    Duration::from_micros(4),
                    Duration::from_micros(2),
                ),
            ],
            observables: ObservableSnapshotBundle::default(),
            numerical_stability: None,
        };

        let timing = inspect_bench_batch_timing(Duration::from_micros(30), &batch).unwrap();

        assert_close(timing.total, 15.0e-6);
        assert_close(timing.parameterization, 1.5e-6);
        assert_close(timing.integrand, 2.0e-6);
        assert_close(timing.event_processing, 1.5e-6);
        assert_close(timing.evaluator, 3.5e-6);
        assert_close(timing.other, 6.5e-6);
    }

    fn sample_with_metadata(
        total: Duration,
        parameterization: Duration,
        integrand: Duration,
        evaluator: Duration,
        event_processing: Duration,
    ) -> SampleEvaluationResult {
        let mut evaluation = EvaluationResult::zero();
        evaluation.evaluation_metadata.total_timing = total;
        evaluation.evaluation_metadata.parameterization_time = parameterization;
        evaluation.evaluation_metadata.integrand_evaluation_time = integrand;
        evaluation.evaluation_metadata.evaluator_evaluation_time = evaluator;
        evaluation.evaluation_metadata.event_processing_time = event_processing;
        SampleEvaluationResult {
            evaluation: evaluation.into_output(false),
        }
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "expected {expected:e}, got {actual:e}"
        );
    }

    #[test]
    fn minimal_integrand_settings_disable_benchmark_overhead_without_touching_subtraction() {
        let mut settings = RuntimeSettings::default();
        settings.general.enable_cache = true;
        settings.general.debug_cache = true;
        settings.general.generate_events = true;
        settings.general.store_additional_weights_in_event = true;
        settings.stability.rotation_axis = vec![RotationSetting::Pi2X {}, RotationSetting::Pi2Y {}];
        settings.stability.levels = vec![
            StabilityLevelSetting::default_double(),
            StabilityLevelSetting::default_quad(),
            StabilityLevelSetting::default_arb(),
        ];
        settings.stability.check_on_norm = true;
        settings.stability.escalate_if_exact_zero = true;
        settings.stability.loop_momenta_norm_escalation_factor = 2.0;
        settings.stability.recording = Some(StabilityRecordingSettings {
            record_rotated_results: true,
            record_all_stability_levels: true,
            record_loop_momenta_escalation: true,
        });
        let original_subtraction = settings.subtraction.clone();

        apply_minimal_integrand_settings(&mut settings);

        assert!(!settings.general.enable_cache);
        assert!(!settings.general.debug_cache);
        assert!(!settings.general.generate_events);
        assert!(!settings.general.store_additional_weights_in_event);
        assert!(settings.observables.is_empty());
        assert!(settings.selectors.values().all(|selector| !selector.active));
        assert_eq!(
            settings.stability.rotation_axis,
            vec![RotationSetting::None {}]
        );
        assert_eq!(
            settings.stability.levels,
            vec![StabilityLevelSetting::default_double()]
        );
        assert!(!settings.stability.check_on_norm);
        assert!(!settings.stability.escalate_if_exact_zero);
        assert_eq!(settings.stability.loop_momenta_norm_escalation_factor, -1.0);
        assert_eq!(settings.stability.recording, None);
        assert_eq!(settings.subtraction, original_subtraction);
    }
}
