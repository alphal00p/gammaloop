use std::{
    collections::BTreeMap,
    env, fs,
    io::{self, IsTerminal},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::{
    integrands::{evaluation::EvaluationResult, process::ProcessIntegrand},
    model::Model,
    observables::events::AdditionalWeightKey,
    utils::F,
};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::{prelude::*, ThreadPoolBuilder};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use tracing::info;

use crate::{
    commands::evaluate_samples::{build_havana_sample, build_momentum_input},
    completion::CompletionArgExt,
    state::{ProcessRef, State},
    CLISettings,
};

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Approach {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// The integrand name to approach
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// The midpoint to inspect (x y) or (p0 px ...)
    #[arg(
        short = 'x',
        long = "point",
        num_args = 1..,
        value_name = "POINT",
        value_delimiter = ',',
        allow_negative_numbers = true,
    )]
    pub point: Vec<f64>,

    /// Direction vector for one approach axis, as comma-separated components
    #[arg(
        long = "approach-axis",
        value_name = "AXIS",
        num_args = 1..,
        allow_negative_numbers = true,
    )]
    pub approach_axes: Vec<String>,

    /// Number of points on each side of the midpoint
    #[arg(long = "n-points", value_name = "N")]
    pub n_points: usize,

    /// Use linear spacing in t
    #[arg(long, conflicts_with = "logarithmic")]
    pub linear: bool,

    /// Use logarithmic spacing in |t|
    #[arg(long, conflicts_with = "linear")]
    pub logarithmic: bool,

    /// Smallest non-zero |t| for logarithmic spacing
    #[arg(
        long = "min-abs-t",
        default_value_t = 1.0e-6,
        allow_negative_numbers = true
    )]
    pub min_abs_t: f64,

    /// Number of worker threads for approach evaluations
    #[arg(long = "n-cores", value_name = "N")]
    pub n_cores: Option<usize>,

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

    /// Select a specific graph in momentum-space approach
    #[arg(long = "graph-id", value_name = "GRAPH_ID")]
    pub graph_id: Option<usize>,

    /// Select a specific orientation of the selected graph in momentum-space approach
    #[arg(
        long = "orientation-id",
        value_name = "ORIENTATION_ID",
        requires = "graph_id",
        conflicts_with = "discrete_dim"
    )]
    pub orientation_id: Option<usize>,

    /// Path to the JSON output file
    #[arg(long = "output-results", value_hint = clap::ValueHint::FilePath)]
    pub output_results: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum ApproachSpacing {
    Linear,
    Logarithmic,
}

#[derive(Debug, Clone)]
struct ApproachJob {
    index: usize,
    axis_index: usize,
    axis_point_index: usize,
    t: f64,
    point: Vec<f64>,
    skip_reason: Option<String>,
}

#[derive(Debug, Serialize)]
struct ApproachOutput {
    schema_version: u32,
    command: ApproachCommandJson,
    process: ApproachProcessJson,
    integrand: ApproachIntegrandJson,
    space: String,
    base_point: Vec<f64>,
    axes: Vec<Vec<f64>>,
    spacing: ApproachSpacingJson,
    n_cores: usize,
    points_per_axis: usize,
    evaluated_points: usize,
    skipped_points: usize,
    points: Vec<ApproachPointRecord>,
}

#[derive(Debug, Serialize)]
struct ApproachCommandJson {
    name: &'static str,
    use_arb_prec: bool,
    graph_id: Option<usize>,
    orientation_id: Option<usize>,
    discrete_dim: Vec<usize>,
}

#[derive(Debug, Serialize)]
struct ApproachProcessJson {
    id: usize,
    name: String,
}

#[derive(Debug, Serialize)]
struct ApproachIntegrandJson {
    name: String,
    kind: String,
}

#[derive(Debug, Serialize)]
struct ApproachSpacingJson {
    kind: ApproachSpacing,
    n_points: usize,
    min_abs_t: Option<f64>,
    t_values: Vec<f64>,
}

#[derive(Debug, Serialize)]
struct ApproachPointRecord {
    index: usize,
    axis_index: usize,
    axis_point_index: usize,
    t: f64,
    point: Vec<f64>,
    status: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    skip_reason: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    evaluation: Option<ApproachEvaluationRecord>,
}

#[derive(Debug, Serialize)]
struct ApproachEvaluationRecord {
    integrand_result: ComplexJson,
    parameterization_jacobian: Option<f64>,
    integrator_weight: f64,
    total_weight: ComplexJson,
    event_weight_sum: ComplexJson,
    additional_weight_sums: BTreeMap<String, ComplexJson>,
    contributions: Vec<ContributionRecord>,
    events: Vec<EventRecord>,
    metadata: EvaluationMetadataJson,
}

#[derive(Debug, Serialize)]
struct EvaluationMetadataJson {
    generated_event_count: usize,
    accepted_event_count: usize,
    is_nan: bool,
    total_time_seconds: f64,
    parameterization_time_seconds: f64,
    integrand_evaluation_time_seconds: f64,
    evaluator_evaluation_time_seconds: f64,
    event_processing_time_seconds: f64,
}

#[derive(Debug, Serialize)]
struct EventRecord {
    event_group_index: usize,
    event_index: usize,
    graph_id: usize,
    graph_name: Option<String>,
    graph_group_id: Option<usize>,
    orientation_id: Option<usize>,
    cut_id: usize,
    cut_edges: Vec<usize>,
    lmb_channel_id: Option<usize>,
    lmb_sample_id: Option<usize>,
    weight: ComplexJson,
    additional_weights: BTreeMap<String, ComplexJson>,
}

#[derive(Debug, Serialize)]
struct ContributionRecord {
    label: String,
    contribution: String,
    graph_id: usize,
    graph_name: Option<String>,
    graph_group_id: Option<usize>,
    orientation_id: Option<usize>,
    cut_id: usize,
    cut_edges: Vec<usize>,
    lmb_sample_id: Option<usize>,
    weight: ComplexJson,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ContributionKey {
    contribution: String,
    graph_id: usize,
    graph_name: Option<String>,
    graph_group_id: Option<usize>,
    orientation_id: Option<usize>,
    cut_id: usize,
    cut_edges: Vec<usize>,
    lmb_sample_id: Option<usize>,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct ComplexJson {
    re: f64,
    im: f64,
}

impl Approach {
    pub fn run(&self, state: &mut State, cli_settings: &CLISettings) -> Result<PathBuf> {
        self.validate_cli()?;
        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;
        let process_name = state.process_list.processes[process_id]
            .definition
            .folder_name
            .clone();
        let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;
        let base_integrand = state
            .process_list
            .get_integrand(process_id, &integrand_name)?
            .require_generated()?
            .clone();

        if !self.momentum_space && (self.graph_id.is_some() || self.orientation_id.is_some()) {
            return Err(eyre!(
                "Graph and orientation selectors are only supported in momentum-space approach."
            ));
        }
        if !self.momentum_space {
            let expected_dimension =
                base_integrand.expected_x_space_dimension(&self.discrete_dim)?;
            if self.point.len() != expected_dimension {
                return Err(eyre!(
                    "Expected {} x-space coordinates for this integrand selection, got {}.",
                    expected_dimension,
                    self.point.len()
                ));
            }
        }

        let graph_name = self.resolve_graph_name(&base_integrand, &integrand_name)?;
        let axes = self.parse_axes()?;
        let spacing = self.spacing()?;
        let t_values = self.t_values(spacing)?;
        let jobs = self.build_jobs(&axes, &t_values);
        let n_cores = self
            .n_cores
            .unwrap_or(cli_settings.global.n_cores.integrate);
        if n_cores == 0 {
            return Err(eyre!("--n-cores must be at least 1."));
        }

        let total_jobs = jobs.len();
        let progress = approach_progress_bar(total_jobs as u64);
        let pool = ThreadPoolBuilder::new()
            .num_threads(n_cores)
            .build()
            .context("failed to build approach evaluation thread pool")?;
        let worker_count = n_cores.min(total_jobs.max(1));
        let worker_inputs = (0..worker_count)
            .map(|worker_index| {
                let worker_jobs = jobs
                    .iter()
                    .skip(worker_index)
                    .step_by(worker_count)
                    .cloned()
                    .collect::<Vec<_>>();
                let mut integrand = base_integrand.clone();
                force_event_output(&mut integrand);
                (integrand, model.clone(), worker_jobs)
            })
            .collect::<Vec<_>>();
        let completed_jobs = Arc::new(AtomicUsize::new(0));
        let completed_workers = Arc::new(AtomicUsize::new(0));
        let results = pool.install(|| {
            worker_inputs
                .into_par_iter()
                .map(
                    |(mut integrand, model, worker_jobs)| -> Result<Vec<ApproachPointRecord>> {
                        integrand
                            .warm_up(&model)
                            .context("failed to warm up approach worker integrand")?;
                        let mut records = Vec::with_capacity(worker_jobs.len());
                        for job in worker_jobs {
                            let record = self.evaluate_job(
                                &job,
                                &mut integrand,
                                &model,
                                graph_name.as_deref(),
                                axes.len(),
                                t_values.len(),
                            );
                            let completed_total =
                                completed_jobs.fetch_add(1, Ordering::Relaxed) + 1;
                            progress.inc(1);
                            progress.set_message(progress_message(
                                job.axis_index,
                                axes.len(),
                                job.axis_point_index,
                                t_values.len(),
                                completed_total,
                                total_jobs,
                            ));
                            records.push(record?);
                        }
                        let completed_worker_count =
                            completed_workers.fetch_add(1, Ordering::Relaxed) + 1;
                        progress.set_message(format!(
                            "{} {}",
                            "collecting worker results".bright_magenta().bold(),
                            format!("{completed_worker_count:>2}/{worker_count:<2}")
                                .bright_white()
                                .bold()
                        ));
                        if completed_worker_count == worker_count {
                            progress.println(
                                "Approach evaluations complete; collecting worker results"
                                    .bright_blue()
                                    .bold()
                                    .to_string(),
                            );
                        }
                        Ok(records)
                    },
                )
                .collect::<Result<Vec<_>>>()
        });
        progress.finish_with_message(
            "approach evaluations complete; finalizing output"
                .bright_green()
                .bold()
                .to_string(),
        );

        info!(
            "{} {}",
            "Finalizing approach results".bright_blue().bold(),
            "(sorting records and preparing JSON)".bright_black()
        );
        let mut point_records = results?
            .into_iter()
            .flatten()
            .collect::<Vec<ApproachPointRecord>>();
        point_records.sort_by_key(|record| record.index);

        let evaluated_points = point_records
            .iter()
            .filter(|record| record.status == "evaluated")
            .count();
        let skipped_points = point_records.len() - evaluated_points;
        let output = ApproachOutput {
            schema_version: 1,
            command: ApproachCommandJson {
                name: "approach",
                use_arb_prec: self.use_arb_prec,
                graph_id: self.graph_id,
                orientation_id: self.orientation_id,
                discrete_dim: self.discrete_dim.clone(),
            },
            process: ApproachProcessJson {
                id: process_id,
                name: process_name,
            },
            integrand: ApproachIntegrandJson {
                name: integrand_name,
                kind: base_integrand.kind_name().to_string(),
            },
            space: if self.momentum_space {
                "momentum".to_string()
            } else {
                "coordinate".to_string()
            },
            base_point: self.point.clone(),
            axes,
            spacing: ApproachSpacingJson {
                kind: spacing,
                n_points: self.n_points,
                min_abs_t: matches!(spacing, ApproachSpacing::Logarithmic)
                    .then_some(self.min_abs_t),
                t_values,
            },
            n_cores,
            points_per_axis: 2 * self.n_points + 1,
            evaluated_points,
            skipped_points,
            points: point_records,
        };

        let output_path = self
            .output_results
            .clone()
            .unwrap_or_else(|| PathBuf::from("approach_result.json"));
        let relative_output = relative_path_display(&output_path);
        info!(
            "{} {}",
            "Writing approach JSON to".bright_blue().bold(),
            relative_output.bright_cyan()
        );
        if let Some(parent) = output_path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent).with_context(|| {
                    format!(
                        "failed to create approach output directory '{}'",
                        parent.display()
                    )
                })?;
            }
        }
        let file = fs::File::create(&output_path).with_context(|| {
            format!(
                "failed to create approach output file '{}'",
                output_path.display()
            )
        })?;
        serde_json::to_writer_pretty(file, &output).with_context(|| {
            format!(
                "failed to serialize approach output to '{}'",
                output_path.display()
            )
        })?;

        let pdf_path = output_path.with_extension("pdf");
        let relative_pdf = relative_path_display(&pdf_path);
        let plot_command = format!(
            "python3 assets/plot_approach_result.py {} --output {}",
            shell_quote(&relative_output),
            shell_quote(&relative_pdf)
        );
        info!(
            "{} {}",
            "Approach results written to".bright_green().bold(),
            relative_output.bright_cyan()
        );
        info!(
            "{} {}",
            "Plot with:".bright_green().bold(),
            plot_command.bright_cyan()
        );

        Ok(output_path)
    }

    fn validate_cli(&self) -> Result<()> {
        if self.point.is_empty() {
            return Err(eyre!("approach requires a midpoint supplied with --point."));
        }
        if self.approach_axes.is_empty() {
            return Err(eyre!(
                "approach requires at least one --approach-axis value."
            ));
        }
        if self.n_points == 0 {
            return Err(eyre!("--n-points must be at least 1."));
        }
        if self.logarithmic && !(self.min_abs_t > 0.0 && self.min_abs_t <= 1.0) {
            return Err(eyre!(
                "--min-abs-t must be greater than 0 and at most 1 for logarithmic spacing."
            ));
        }
        if self.n_cores == Some(0) {
            return Err(eyre!("--n-cores must be at least 1."));
        }
        Ok(())
    }

    fn spacing(&self) -> Result<ApproachSpacing> {
        if self.linear && self.logarithmic {
            return Err(eyre!("--linear and --logarithmic are mutually exclusive."));
        }
        Ok(if self.logarithmic {
            ApproachSpacing::Logarithmic
        } else {
            ApproachSpacing::Linear
        })
    }

    fn parse_axes(&self) -> Result<Vec<Vec<f64>>> {
        self.approach_axes
            .iter()
            .map(|raw_axis| parse_axis(raw_axis, self.point.len()))
            .collect()
    }

    fn resolve_graph_name(
        &self,
        integrand: &ProcessIntegrand,
        integrand_name: &str,
    ) -> Result<Option<String>> {
        let Some(graph_id) = self.graph_id else {
            return Ok(None);
        };
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

    fn t_values(&self, spacing: ApproachSpacing) -> Result<Vec<f64>> {
        let magnitudes = match spacing {
            ApproachSpacing::Linear => (1..=self.n_points)
                .map(|index| index as f64 / self.n_points as f64)
                .collect::<Vec<_>>(),
            ApproachSpacing::Logarithmic => {
                if self.n_points == 1 {
                    vec![1.0]
                } else {
                    let log_min = self.min_abs_t.ln();
                    let denom = (self.n_points - 1) as f64;
                    (0..self.n_points)
                        .map(|index| (log_min * (1.0 - index as f64 / denom)).exp())
                        .collect()
                }
            }
        };
        let mut t_values = magnitudes
            .iter()
            .rev()
            .map(|value| -*value)
            .collect::<Vec<_>>();
        t_values.push(0.0);
        t_values.extend(magnitudes);
        Ok(t_values)
    }

    fn build_jobs(&self, axes: &[Vec<f64>], t_values: &[f64]) -> Vec<ApproachJob> {
        let mut jobs = Vec::with_capacity(axes.len() * t_values.len());
        for (axis_index, axis) in axes.iter().enumerate() {
            for (axis_point_index, t) in t_values.iter().copied().enumerate() {
                let point = self
                    .point
                    .iter()
                    .zip(axis)
                    .map(|(center, direction)| center + t * direction)
                    .collect::<Vec<_>>();
                let skip_reason = (!self.momentum_space)
                    .then(|| coordinate_skip_reason(&point))
                    .flatten();
                jobs.push(ApproachJob {
                    index: jobs.len(),
                    axis_index,
                    axis_point_index,
                    t,
                    point,
                    skip_reason,
                });
            }
        }
        jobs
    }

    fn evaluate_job(
        &self,
        job: &ApproachJob,
        integrand: &mut ProcessIntegrand,
        model: &Model,
        graph_name: Option<&str>,
        _axis_count: usize,
        _axis_point_count: usize,
    ) -> Result<ApproachPointRecord> {
        if let Some(skip_reason) = &job.skip_reason {
            return Ok(ApproachPointRecord {
                index: job.index,
                axis_index: job.axis_index,
                axis_point_index: job.axis_point_index,
                t: job.t,
                point: job.point.clone(),
                status: "skipped",
                skip_reason: Some(skip_reason.clone()),
                evaluation: None,
            });
        }

        let mut samples = if self.momentum_space {
            let input = build_momentum_input(
                integrand,
                &job.point,
                1.0,
                &self.discrete_dim,
                graph_name,
                self.orientation_id,
            )?;
            integrand
                .evaluate_momentum_configurations_raw(model, &[input], self.use_arb_prec)?
                .samples
        } else {
            let sample = build_havana_sample(integrand, &job.point, &self.discrete_dim, 1.0)?;
            integrand
                .evaluate_samples_raw(
                    model,
                    &[sample],
                    1,
                    self.use_arb_prec,
                    false,
                    Complex::new_zero(),
                )?
                .samples
        };
        let evaluation = samples
            .pop()
            .ok_or_else(|| eyre!("approach evaluation did not return a sample"))?;
        Ok(ApproachPointRecord {
            index: job.index,
            axis_index: job.axis_index,
            axis_point_index: job.axis_point_index,
            t: job.t,
            point: job.point.clone(),
            status: "evaluated",
            skip_reason: None,
            evaluation: Some(approach_evaluation_record(integrand, evaluation)?),
        })
    }
}

fn parse_axis(raw_axis: &str, expected_dimension: usize) -> Result<Vec<f64>> {
    let trimmed = raw_axis
        .trim()
        .trim_start_matches('[')
        .trim_end_matches(']');
    let axis = trimmed
        .split(',')
        .map(str::trim)
        .filter(|component| !component.is_empty())
        .map(|component| {
            component.parse::<f64>().map_err(|err| {
                eyre!(
                    "Could not parse approach-axis component '{}' as a float: {}",
                    component,
                    err
                )
            })
        })
        .collect::<Result<Vec<_>>>()?;
    if axis.len() != expected_dimension {
        return Err(eyre!(
            "Approach axis has {} components, but the point has {} components.",
            axis.len(),
            expected_dimension
        ));
    }
    if axis.iter().all(|component| *component == 0.0) {
        return Err(eyre!("Approach axis must not be the zero vector."));
    }
    Ok(axis)
}

fn coordinate_skip_reason(point: &[f64]) -> Option<String> {
    let invalid = point
        .iter()
        .enumerate()
        .filter(|(_, value)| !(0.0..=1.0).contains(*value))
        .map(|(index, value)| format!("{index}:{value:+.16e}"))
        .collect::<Vec<_>>();
    (!invalid.is_empty()).then(|| {
        format!(
            "coordinate-space point lies outside the unit hypercube at component(s) {}",
            invalid.join(", ")
        )
    })
}

fn force_event_output(integrand: &mut ProcessIntegrand) {
    let settings = integrand.get_mut_settings();
    settings.general.generate_events = true;
    settings.general.store_additional_weights_in_event = true;
}

fn approach_evaluation_record(
    integrand: &ProcessIntegrand,
    evaluation: EvaluationResult,
) -> Result<ApproachEvaluationRecord> {
    let parameterization_jacobian = evaluation.parameterization_jacobian.map(|value| value.0);
    let integrator_weight = evaluation.integrator_weight.0;
    let total_scale = parameterization_jacobian.unwrap_or(1.0) * integrator_weight;
    let total_weight = Complex::new(
        F(evaluation.integrand_result.re.0 * total_scale),
        F(evaluation.integrand_result.im.0 * total_scale),
    );

    let mut event_weight_sum = Complex::new(F(0.0), F(0.0));
    let mut additional_weight_sums = BTreeMap::<String, Complex<F<f64>>>::new();
    let mut contribution_sums = BTreeMap::<ContributionKey, Complex<F<f64>>>::new();
    let mut events = Vec::new();
    let parameterization_settings = integrand
        .get_settings()
        .sampling
        .get_parameterization_settings()
        .unwrap_or_default();

    for (event_group_index, event_group) in evaluation.event_groups.iter().enumerate() {
        for (event_index, event) in event_group.iter().enumerate() {
            event_weight_sum += event.weight;
            let graph_id = event.cut_info.graph_id;
            let graph_name = integrand
                .graph_name_by_id(graph_id)
                .map(ToString::to_string);
            let graph_group_id = integrand.graph_group_id_by_graph_id(graph_id);
            let cut_id = event.cut_info.cut_id;
            let cut_edges = integrand.cut_edge_ids(graph_id, cut_id).unwrap_or_default();
            let lmb_channel_id = event.cut_info.lmb_channel_id;
            let lmb_sample_id = lmb_channel_id
                .map(|channel_id| {
                    integrand.lmb_sample_id_for_channel(
                        graph_id,
                        channel_id,
                        &parameterization_settings,
                    )
                })
                .transpose()?
                .flatten();

            let event_key = ContributionKey {
                contribution: "event_weight".to_string(),
                graph_id,
                graph_name: graph_name.clone(),
                graph_group_id,
                orientation_id: event.cut_info.orientation_id,
                cut_id,
                cut_edges: cut_edges.clone(),
                lmb_sample_id,
            };
            add_contribution(&mut contribution_sums, event_key, event.weight);

            let mut additional_weights = BTreeMap::new();
            for (key, value) in &event.additional_weights.weights {
                let label = additional_weight_key_label(*key);
                *additional_weight_sums
                    .entry(label.clone())
                    .or_insert_with(zero_complex) += *value;
                additional_weights.insert(label.clone(), complex_json(*value));
                let contribution_key = ContributionKey {
                    contribution: label,
                    graph_id,
                    graph_name: graph_name.clone(),
                    graph_group_id,
                    orientation_id: event.cut_info.orientation_id,
                    cut_id,
                    cut_edges: cut_edges.clone(),
                    lmb_sample_id,
                };
                add_contribution(&mut contribution_sums, contribution_key, *value);
            }

            events.push(EventRecord {
                event_group_index,
                event_index,
                graph_id,
                graph_name,
                graph_group_id,
                orientation_id: event.cut_info.orientation_id,
                cut_id,
                cut_edges,
                lmb_channel_id,
                lmb_sample_id,
                weight: complex_json(event.weight),
                additional_weights,
            });
        }
    }

    Ok(ApproachEvaluationRecord {
        integrand_result: complex_json(evaluation.integrand_result),
        parameterization_jacobian,
        integrator_weight,
        total_weight: complex_json(total_weight),
        event_weight_sum: complex_json(event_weight_sum),
        additional_weight_sums: additional_weight_sums
            .into_iter()
            .map(|(key, value)| (key, complex_json(value)))
            .collect(),
        contributions: contribution_sums
            .into_iter()
            .map(|(key, weight)| contribution_record(key, weight))
            .collect(),
        events,
        metadata: EvaluationMetadataJson {
            generated_event_count: evaluation.evaluation_metadata.generated_event_count,
            accepted_event_count: evaluation.evaluation_metadata.accepted_event_count,
            is_nan: evaluation.evaluation_metadata.is_nan,
            total_time_seconds: evaluation.evaluation_metadata.total_timing.as_secs_f64(),
            parameterization_time_seconds: evaluation
                .evaluation_metadata
                .parameterization_time
                .as_secs_f64(),
            integrand_evaluation_time_seconds: evaluation
                .evaluation_metadata
                .integrand_evaluation_time
                .as_secs_f64(),
            evaluator_evaluation_time_seconds: evaluation
                .evaluation_metadata
                .evaluator_evaluation_time
                .as_secs_f64(),
            event_processing_time_seconds: evaluation
                .evaluation_metadata
                .event_processing_time
                .as_secs_f64(),
        },
    })
}

fn add_contribution(
    contribution_sums: &mut BTreeMap<ContributionKey, Complex<F<f64>>>,
    key: ContributionKey,
    value: Complex<F<f64>>,
) {
    *contribution_sums.entry(key).or_insert_with(zero_complex) += value;
}

fn contribution_record(key: ContributionKey, weight: Complex<F<f64>>) -> ContributionRecord {
    ContributionRecord {
        label: contribution_label(&key),
        contribution: key.contribution,
        graph_id: key.graph_id,
        graph_name: key.graph_name,
        graph_group_id: key.graph_group_id,
        orientation_id: key.orientation_id,
        cut_id: key.cut_id,
        cut_edges: key.cut_edges,
        lmb_sample_id: key.lmb_sample_id,
        weight: complex_json(weight),
    }
}

fn contribution_label(key: &ContributionKey) -> String {
    let graph = key
        .graph_name
        .as_deref()
        .map(ToString::to_string)
        .unwrap_or_else(|| format!("#{}", key.graph_id));
    let edge_label = if key.cut_edges.is_empty() {
        "[]".to_string()
    } else {
        format!(
            "[{}]",
            key.cut_edges
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )
    };
    let mut parts = vec![
        key.contribution.clone(),
        format!("graph={graph}"),
        format!("cut={} edges={edge_label}", key.cut_id),
    ];
    if let Some(orientation_id) = key.orientation_id {
        parts.push(format!("orientation={orientation_id}"));
    }
    if let Some(lmb_sample_id) = key.lmb_sample_id {
        parts.push(format!("lmb_sample={lmb_sample_id}"));
    }
    parts.join(" ")
}

fn additional_weight_key_label(key: AdditionalWeightKey) -> String {
    match key {
        AdditionalWeightKey::FullMultiplicativeFactor => "full_multiplicative_factor".to_string(),
        AdditionalWeightKey::Original => "original".to_string(),
        AdditionalWeightKey::ThresholdCounterterm { subset_index } => {
            format!("threshold_counterterm_{subset_index}")
        }
        AdditionalWeightKey::AmplitudeThresholdCounterterm {
            esurface_id,
            overlap_group,
        } => format!("ct_{esurface_id}_{overlap_group}"),
    }
}

fn zero_complex() -> Complex<F<f64>> {
    Complex::new(F(0.0), F(0.0))
}

fn complex_json(value: Complex<F<f64>>) -> ComplexJson {
    ComplexJson {
        re: value.re.0,
        im: value.im.0,
    }
}

fn approach_progress_bar(len: u64) -> Arc<ProgressBar> {
    let bar = if io::stderr().is_terminal() {
        ProgressBar::new(len)
    } else {
        ProgressBar::hidden()
    };
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise} | ETA: {eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
        )
        .expect("approach progress bar template should be valid"),
    );
    Arc::new(bar)
}

fn progress_message(
    axis_index: usize,
    axis_count: usize,
    axis_point_index: usize,
    axis_point_count: usize,
    completed_total: usize,
    total: usize,
) -> String {
    format!(
        "{} {} {} {} {} {}",
        "axis".bright_cyan().bold(),
        format!("{:>2}/{:<2}", axis_index + 1, axis_count)
            .bright_white()
            .bold(),
        "point".bright_cyan().bold(),
        format!("{:>5}/{:<5}", axis_point_index + 1, axis_point_count)
            .bright_white()
            .bold(),
        "total".bright_cyan().bold(),
        format!("{:>5}/{:<5}", completed_total, total)
            .bright_green()
            .bold()
    )
}

fn relative_path_display(path: &Path) -> String {
    if path.is_absolute() {
        if let Ok(current_dir) = env::current_dir() {
            if let Ok(relative) = path.strip_prefix(current_dir) {
                return relative.display().to_string();
            }
        }
    }
    path.display().to_string()
}

fn shell_quote(value: &str) -> String {
    if value
        .chars()
        .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '/' | '.' | '_' | '-' | '='))
    {
        value.to_string()
    } else {
        format!("'{}'", value.replace('\'', "'\\''"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn command_with_spacing(n_points: usize, logarithmic: bool, min_abs_t: f64) -> Approach {
        Approach {
            point: vec![0.5, 0.5],
            approach_axes: vec!["0.1,0.0".to_string()],
            n_points,
            logarithmic,
            min_abs_t,
            ..Default::default()
        }
    }

    #[test]
    fn linear_t_values_are_ordered_around_zero() {
        let command = command_with_spacing(2, false, 1.0e-6);
        assert_eq!(
            command.t_values(ApproachSpacing::Linear).unwrap(),
            vec![-1.0, -0.5, 0.0, 0.5, 1.0]
        );
    }

    #[test]
    fn logarithmic_t_values_use_min_abs_t() {
        let command = command_with_spacing(3, true, 1.0e-4);
        let t_values = command.t_values(ApproachSpacing::Logarithmic).unwrap();
        assert_eq!(t_values.len(), 7);
        assert!((t_values[0] + 1.0).abs() < 1.0e-14);
        assert!((t_values[2] + 1.0e-4).abs() < 1.0e-14);
        assert_eq!(t_values[3], 0.0);
        assert!((t_values[4] - 1.0e-4).abs() < 1.0e-14);
        assert!((t_values[6] - 1.0).abs() < 1.0e-14);
    }

    #[test]
    fn parse_axis_rejects_dimension_mismatch() {
        let err = parse_axis("1.0,0.0,2.0", 2).unwrap_err();
        assert!(format!("{err:#}").contains("Approach axis has 3 components"));
    }

    #[test]
    fn coordinate_jobs_skip_points_outside_unit_hypercube() {
        let command = Approach {
            point: vec![0.95, 0.5],
            approach_axes: vec!["0.1,0.0".to_string()],
            n_points: 1,
            ..Default::default()
        };
        let axes = command.parse_axes().unwrap();
        let t_values = command.t_values(ApproachSpacing::Linear).unwrap();
        let jobs = command.build_jobs(&axes, &t_values);
        assert!(jobs[2].skip_reason.is_some());
        assert!(jobs[1].skip_reason.is_none());
    }

    #[test]
    fn zero_n_cores_is_rejected() {
        let command = Approach {
            point: vec![0.5, 0.5],
            approach_axes: vec!["0.1,0.0".to_string()],
            n_points: 1,
            n_cores: Some(0),
            ..Default::default()
        };
        let err = command.validate_cli().unwrap_err();
        assert!(format!("{err:#}").contains("--n-cores must be at least 1"));
    }

    #[test]
    fn invalid_log_min_abs_t_is_rejected() {
        let command = Approach {
            point: vec![0.5, 0.5],
            approach_axes: vec!["0.1,0.0".to_string()],
            n_points: 1,
            logarithmic: true,
            min_abs_t: 0.0,
            ..Default::default()
        };
        let err = command.validate_cli().unwrap_err();
        assert!(format!("{err:#}").contains("--min-abs-t must be greater than 0"));
    }

    #[test]
    fn simulated_parallel_completion_is_sorted_back_to_job_order() {
        let command = Approach {
            point: vec![0.5, 0.5],
            approach_axes: vec!["0.1,0.0".to_string(), "0.0,0.2".to_string()],
            n_points: 2,
            ..Default::default()
        };
        let axes = command.parse_axes().unwrap();
        let t_values = command.t_values(ApproachSpacing::Linear).unwrap();
        let jobs = command.build_jobs(&axes, &t_values);
        let mut records = (0..2)
            .flat_map(|worker_index| jobs.iter().skip(worker_index).step_by(2))
            .map(|job| ApproachPointRecord {
                index: job.index,
                axis_index: job.axis_index,
                axis_point_index: job.axis_point_index,
                t: job.t,
                point: job.point.clone(),
                status: "skipped",
                skip_reason: None,
                evaluation: None,
            })
            .collect::<Vec<_>>();

        assert_ne!(
            records
                .iter()
                .map(|record| record.index)
                .collect::<Vec<_>>(),
            (0..records.len()).collect::<Vec<_>>(),
            "strided worker completion should be out of original order in this simulation",
        );
        records.sort_by_key(|record| record.index);
        assert_eq!(
            records
                .iter()
                .map(|record| (record.index, record.axis_index, record.axis_point_index))
                .collect::<Vec<_>>(),
            jobs.iter()
                .map(|job| (job.index, job.axis_index, job.axis_point_index))
                .collect::<Vec<_>>(),
        );
    }
}
