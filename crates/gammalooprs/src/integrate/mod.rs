#![allow(dead_code)]

//! This module steers the integration process.
//! It contains two main ways of integrating the integrand.
//! The havana_integrate function is mostly used for local runs.
//! The master node in combination with batch_integrate is for distributed runs.

use bincode::Decode;
use bincode::Encode;
use color_eyre::{Report, Result};
use colored::Colorize;
use itertools::Itertools;
use itertools::izip;
use rayon::ThreadPoolBuilder;
use rayon::iter::repeat_n;
use serde::Deserialize;
use serde::Serialize;
use spenso::algebra::algebraic_traits::IsZero;
use symbolica::domains::float::Constructible;
use symbolica::numerical_integration::{Grid, MonteCarloRng, Sample, StatisticsAccumulator};

use crate::INTERRUPTED;
use crate::Integrand;
use crate::integrands::HasIntegrand;
use crate::integrands::evaluation::EvaluationResult;
use crate::integrands::evaluation::StatisticsCounter;
use crate::model::Model;
use crate::observables::{EventGroupList, ObservableAccumulatorBundle, ObservableFileFormat};
use crate::settings::IntegratorSettings;
use crate::settings::RuntimeSettings;
use crate::settings::runtime::{IntegratedPhase, IntegrationResult};
use crate::utils;
use crate::utils::F;
use crate::utils::normalize_tabled_separator_rows;
use crate::{is_interrupted, set_interrupted};
use rayon::prelude::*;
use spenso::algebra::complex::Complex;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::time::Duration;
use std::time::Instant;
use tabled::{
    Table, Tabled,
    builder::Builder,
    settings::{
        Alignment, Modify, Panel, Span,
        object::{Cell, Rows},
        style::{HorizontalLine, On, Style, VerticalLine},
        themes::BorderCorrection,
        width::Width,
    },
};
#[allow(unused_imports)]
use tracing::{debug, error, info, trace, warn};

#[derive(Clone, Copy)]
enum ObservableFlushReason {
    Iteration(usize),
    Final,
    Interrupted,
}

// const N_INTEGRAND_ACCUMULATORS: usize = 2;

/// Intended to be a copy of the integrand for each core.
pub struct UserData {
    pub integrand: Vec<Integrand>,
}

#[derive(Tabled)]
pub struct IntegralResult {
    id: String,
    n_samples: String,
    #[tabled(rename = "n_samples[%]")]
    n_samples_perc: String,
    #[tabled(rename = "<I>")]
    integral: String,
    #[tabled(rename = "sqrt(σ)")]
    variance: String,
    err: String,
    #[tabled(rename = "err[%]")]
    err_perc: String,
    #[tabled(rename = "PDF")]
    pdf: String,
}

#[derive(Serialize, Deserialize, Encode, Decode)]
pub struct ComplexAccumulator {
    #[bincode(with_serde)]
    pub re: StatisticsAccumulator<F<f64>>,
    #[bincode(with_serde)]
    pub im: StatisticsAccumulator<F<f64>>,
}

impl ComplexAccumulator {
    pub(crate) fn new() -> Self {
        Self {
            re: StatisticsAccumulator::new(),
            im: StatisticsAccumulator::new(),
        }
    }

    /// used to escalate the stability test on high evaluation values
    pub(crate) fn get_worst_case(&self) -> Complex<F<f64>> {
        Complex::new(
            self.re
                .max_eval_positive
                .abs()
                .max(self.re.max_eval_negative.abs()),
            self.im
                .max_eval_positive
                .abs()
                .max(self.im.max_eval_negative.abs()),
        )
    }

    /// add a sample to the accumulator
    pub(crate) fn add_sample(
        &mut self,
        result: Complex<F<f64>>,
        sample_weight: F<f64>,
        sample: Option<&Sample<F<f64>>>,
    ) {
        self.re.add_sample(result.re * sample_weight, sample);
        self.im.add_sample(result.im * sample_weight, sample);
    }

    pub(crate) fn merge(&mut self, other: &Self) {
        self.re.merge_samples_no_reset(&other.re);
        self.im.merge_samples_no_reset(&other.im);
    }

    pub(crate) fn update_iter(&mut self, use_weighted_average: bool) {
        self.re.update_iter(use_weighted_average);
        self.im.update_iter(use_weighted_average);
    }

    pub(crate) fn max_weight_rows(&self, i_itg: usize) -> Vec<[String; 3]> {
        let max_evals = [
            &self.re.max_eval_positive,
            &self.re.max_eval_negative,
            &self.im.max_eval_positive,
            &self.im.max_eval_negative,
        ];

        let max_eval_samples = [
            &self.re.max_eval_positive_xs,
            &self.re.max_eval_negative_xs,
            &self.im.max_eval_positive_xs,
            &self.im.max_eval_negative_xs,
        ];

        let sign_strs = ["+", "-", "+", "-"];
        let phase_strs = ["re", "re", "im", "im"];

        izip!(
            max_evals.iter(),
            max_eval_samples.iter(),
            sign_strs.iter(),
            phase_strs.iter()
        )
        .filter_map(|(max_eval, max_eval_sample, sign_str, phase_str)| {
            if max_eval.is_non_zero() {
                Some([
                    format!(
                        "itg #{:-3} {} [{}] ",
                        format!("{:<3}", i_itg + 1),
                        format!("{:<2}", phase_str).blue(),
                        format!("{:<1}", sign_str).blue()
                    ),
                    format!("{:+.16e}", max_eval),
                    if let Some(sample) = max_eval_sample {
                        format_max_eval_sample(sample)
                    } else {
                        "N/A".to_string()
                    },
                ])
            } else {
                None
            }
        })
        .collect()
    }
}

struct IntegralResultCells {
    integrand: String,
    value: String,
    relative_error: String,
    chi_sq: String,
    delta_sigma: Option<String>,
    delta_percent: Option<String>,
    mwi: String,
}

const MAX_SHARED_TABLE_WIDTH: usize = 150;

fn result_group_separator() -> VerticalLine<On, On, ()> {
    VerticalLine::new('│').top('┬').bottom('┴')
}

fn format_max_eval_coordinate(value: F<f64>) -> String {
    let formatted = format!("{:.16e}", value.0);
    let Some((mantissa, exponent)) = formatted.rsplit_once('e') else {
        return formatted;
    };
    let exponent = exponent.parse::<i32>().unwrap_or_default();
    format!("{mantissa}e{exponent:+03}")
}

fn format_max_eval_coordinates(xs: &[F<f64>]) -> String {
    format!(
        "[ {} ]",
        xs.iter()
            .map(|value| format_max_eval_coordinate(*value))
            .join(" ")
    )
}

fn format_max_eval_sample(sample: &Sample<F<f64>>) -> String {
    match sample {
        Sample::Continuous(_, xs) => {
            format!("xs: {}", format_max_eval_coordinates(xs))
        }
        Sample::Discrete(_, graph_index, Some(nested_sample)) => match nested_sample.as_ref() {
            Sample::Continuous(_, xs) => {
                format!(
                    "graph: {graph_index}, xs: {}",
                    format_max_eval_coordinates(xs)
                )
            }
            Sample::Discrete(_, channel_index, Some(nested_cont_sample)) => {
                match nested_cont_sample.as_ref() {
                    Sample::Continuous(_, xs) => format!(
                        "graph: {graph_index}, channel: {channel_index}, xs: {}",
                        format_max_eval_coordinates(xs)
                    ),
                    _ => String::from("N/A"),
                }
            }
            _ => String::from("N/A"),
        },
        _ => String::from("N/A"),
    }
}

fn format_iteration_points(points: usize) -> String {
    format!("{:.0}K", points as f64 / 1000.0)
}

fn format_total_points(points: usize) -> String {
    if points >= 10_000_000 {
        format!("{:.0}M", points as f64 / 1_000_000.0)
    } else {
        format!("{:.0}K", points as f64 / 1000.0)
    }
}

fn build_integral_result_cells(
    itg: &StatisticsAccumulator<F<f64>>,
    i_itg: usize,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) -> IntegralResultCells {
    let relative_error = if itg.avg.is_non_zero() {
        if (itg.err / itg.avg).abs().0 > 0.01 {
            format!(
                "{:-8}",
                format!("{:.3}%", (itg.err / itg.avg).abs().0 * 100.).red()
            )
        } else {
            format!(
                "{:-8}",
                format!("{:.3}%", (itg.err / itg.avg).abs().0 * 100.).green()
            )
        }
    } else {
        format!("{:-8}", "")
    };

    let chi_sq = if itg.chi_sq / F::<f64>::new_from_usize(i_iter) > F(5.) {
        format!("{:-6.3} χ²/dof", itg.chi_sq.0 / (i_iter as f64)).red()
    } else {
        format!("{:-6.3} χ²/dof", itg.chi_sq.0 / (i_iter as f64)).normal()
    };

    let (delta_sigma, delta_percent) = if i_itg == 1 {
        if let Some(t) = trgt {
            let delta_in_sigmas = (t - itg.avg).abs().0 / itg.err.0;
            let delta_in_percent = if t.abs() > F(0.) {
                (t - itg.avg).abs().0 / t.abs().0 * 100.
            } else {
                0.
            };
            let is_outside_target = delta_in_sigmas > 5.
                || (t.abs().is_non_zero() && ((t - itg.avg).abs() / t.abs()).0 > 0.01);
            let sigma_display = if is_outside_target {
                delta_in_sigmas
            } else if t.is_non_zero() && itg.avg.is_non_zero() {
                delta_in_sigmas
            } else {
                0.
            };
            let sigma_text = format!("Δ = {:.3}σ", sigma_display);
            let percent_text = format!("Δ = {:.3}%", delta_in_percent);
            if is_outside_target {
                (
                    Some(sigma_text.red().to_string()),
                    Some(percent_text.red().to_string()),
                )
            } else {
                (
                    Some(sigma_text.green().to_string()),
                    Some(percent_text.green().to_string()),
                )
            }
        } else {
            (None, None)
        }
    } else {
        (None, None)
    };

    let mwi = if itg.avg.abs().0 != 0. {
        let mwi = itg.max_eval_negative.abs().max(itg.max_eval_positive.abs())
            / (itg.avg.abs() * F::<f64>::new_from_usize(itg.processed_samples));
        if mwi > F(1.) {
            format!("  mwi: {:<10.4e}", mwi.0).red().to_string()
        } else {
            format!("  mwi: {:<10.4e}", mwi.0).normal().to_string()
        }
    } else {
        format!("  mwi: {:<10.4e}", 0.).normal().to_string()
    };

    IntegralResultCells {
        integrand: format!(
            "itg #{:-3} {}:",
            format!("{:<3}", i_itg),
            format!("{:-2}", tag).blue().bold(),
        ),
        value: utils::format_uncertainty(itg.avg, itg.err)
            .blue()
            .bold()
            .to_string(),
        relative_error: relative_error.to_string(),
        chi_sq: chi_sq.to_string(),
        delta_sigma,
        delta_percent,
        mwi,
    }
}

fn build_iteration_status_header_left(elapsed_time: Duration, iter: usize) -> String {
    format!(
        "[ {} ] {}",
        format!(
            "{:^7}",
            utils::format_wdhms(elapsed_time.as_secs() as usize)
        )
        .bold(),
        format!("Iteration #{:-4}", iter).bold().green(),
    )
}

fn build_iteration_status_header_middle(cur_points: usize, total_points: usize) -> String {
    format!(
        "n_pts = {} {}",
        format_iteration_points(cur_points),
        format!("n_tot = {}", format_total_points(total_points))
            .bold()
            .green(),
    )
}

fn build_iteration_status_header_tail(
    cores: usize,
    elapsed_time: Duration,
    n_samples_evaluated: usize,
) -> String {
    let average_sample_time = if n_samples_evaluated == 0 {
        "N/A".red().to_string()
    } else {
        utils::format_evaluation_time_from_f64(
            elapsed_time.as_secs_f64() / (n_samples_evaluated as f64) * (cores as f64),
        )
        .bold()
        .blue()
        .to_string()
    };

    format!("{average_sample_time} /sample/core")
}

fn build_iteration_results_table(
    integration_state: &IntegrationState,
    cores: usize,
    elapsed_time: Duration,
    cur_points: usize,
    n_samples_evaluated: usize,
    target: &Option<Complex<F<f64>>>,
) -> Table {
    let has_target_columns = target.is_some();
    let n_columns = if has_target_columns { 7 } else { 5 };
    let mut builder = Builder::new();
    let mut first_row = vec![
        build_iteration_status_header_left(elapsed_time, integration_state.iter),
        String::new(),
        build_iteration_status_header_middle(cur_points, integration_state.num_points),
    ];
    first_row.resize(n_columns - 1, String::new());
    first_row.push(build_iteration_status_header_tail(
        cores,
        elapsed_time,
        n_samples_evaluated,
    ));
    builder.push_record(first_row);

    let mut rows_without_target = Vec::new();
    let mut row_index = 1usize;

    for (i_integrand, integrand) in integration_state.all_integrals.iter().enumerate() {
        for (tag, accumulator, target_component) in [
            (
                "re",
                &integrand.re,
                if i_integrand == 0 {
                    target.map(|value| value.re)
                } else {
                    None
                },
            ),
            (
                "im",
                &integrand.im,
                if i_integrand == 0 {
                    target.map(|value| value.im)
                } else {
                    None
                },
            ),
        ] {
            let cells = build_integral_result_cells(
                accumulator,
                i_integrand + 1,
                integration_state.iter,
                tag,
                target_component,
            );
            if has_target_columns {
                if let (Some(delta_sigma), Some(delta_percent)) =
                    (cells.delta_sigma, cells.delta_percent)
                {
                    builder.push_record(vec![
                        cells.integrand,
                        cells.value,
                        cells.relative_error,
                        cells.chi_sq,
                        delta_sigma,
                        delta_percent,
                        cells.mwi,
                    ]);
                } else {
                    builder.push_record(vec![
                        cells.integrand,
                        cells.value,
                        cells.relative_error,
                        cells.chi_sq,
                        cells.mwi,
                        String::new(),
                        String::new(),
                    ]);
                    rows_without_target.push(row_index);
                }
            } else {
                builder.push_record(vec![
                    cells.integrand,
                    cells.value,
                    cells.relative_error,
                    cells.chi_sq,
                    cells.mwi,
                ]);
            }
            row_index += 1;
        }
    }

    let mut table = builder.build();
    table.modify((0, 0), Span::column(2));
    table.modify((0, 2), Span::column((n_columns - 3) as isize));
    for row in rows_without_target {
        table.modify((row, 4), Span::column(3));
    }

    if has_target_columns {
        table.with(
            Style::rounded()
                .remove_horizontals()
                .verticals([
                    (3, result_group_separator()),
                    (4, result_group_separator()),
                    (6, result_group_separator()),
                ])
                .horizontals([(
                    1,
                    HorizontalLine::new('─')
                        .intersection('┬')
                        .left('├')
                        .right('┤'),
                )])
                .remove_vertical(),
        );
    } else {
        table.with(
            Style::rounded()
                .remove_horizontals()
                .verticals([(3, result_group_separator()), (4, result_group_separator())])
                .horizontals([(
                    1,
                    HorizontalLine::new('─')
                        .intersection('┬')
                        .left('├')
                        .right('┤'),
                )])
                .remove_vertical(),
        );
    }

    table.with(BorderCorrection::span());
    table.with(Modify::new(Rows::new(0..)).with(Alignment::left()));
    table.with(Modify::new(Cell::new(0, 2)).with(Alignment::center()));
    table.with(Modify::new(Cell::new(0, n_columns - 1)).with(Alignment::center()));
    table
}

fn build_max_weight_details_table(all_integrals: &[ComplexAccumulator]) -> Table {
    let mut builder = Builder::new();
    builder.push_record([
        "Integrand".bold().blue().to_string(),
        "Max eval".bold().blue().to_string(),
        "Max eval coordinates".bold().blue().to_string(),
    ]);

    for (i_itg, integral) in all_integrals.iter().enumerate() {
        for row in integral.max_weight_rows(i_itg) {
            builder.push_record(row);
        }
    }

    let mut table = builder.build();
    table.with(Panel::header(
        "Maximum weight details".bold().green().to_string(),
    ));
    table.with(Style::rounded().remove_horizontals().horizontals([
        (1, HorizontalLine::full('─', '┬', '├', '┤')),
        (2, HorizontalLine::full('─', '┼', '├', '┤')),
    ]));
    table.with(BorderCorrection::span());
    table.with(Modify::new(Rows::new(0..2)).with(Alignment::center()));
    table.with(Modify::new(Rows::new(2..)).with(Alignment::left()));
    table
}

fn render_max_weight_details_table(all_integrals: &[ComplexAccumulator]) -> String {
    normalize_tabled_separator_rows(&build_max_weight_details_table(all_integrals).to_string())
}

fn suppress_iteration_header_tail_separator(rendered: &str) -> String {
    rendered
        .lines()
        .enumerate()
        .map(|(line_index, line)| {
            if line_index != 1 {
                return line.to_string();
            }

            let vertical_positions = line.match_indices('│').map(|(idx, _)| idx).collect_vec();
            if vertical_positions.len() < 3 {
                return line.to_string();
            }

            let suppressed_index = vertical_positions[vertical_positions.len() - 2];
            let mut updated = line.to_string();
            updated.replace_range(suppressed_index..suppressed_index + '│'.len_utf8(), " ");
            updated
        })
        .join("\n")
}

fn render_tables_with_shared_width(mut tables: Vec<Table>) -> String {
    let max_width = tables
        .iter()
        .map(Table::total_width)
        .max()
        .unwrap_or(0)
        .min(MAX_SHARED_TABLE_WIDTH);

    for table in &mut tables {
        if table.total_width() < max_width {
            table.with(Width::increase(max_width));
        }
    }

    tables
        .into_iter()
        .enumerate()
        .map(|(table_index, table)| {
            let rendered = table.to_string();
            let rendered = if table_index == 0 {
                suppress_iteration_header_tail_separator(&rendered)
            } else {
                rendered
            };
            normalize_tabled_separator_rows(&rendered)
        })
        .collect_vec()
        .join("\n")
}

/// struct to keep track of state, used in the havana_integrate function
/// the idea is to save this to disk after each iteration, so that the integration can be resumed
#[derive(Serialize, Deserialize, Encode, Decode)]
pub struct IntegrationState {
    pub num_points: usize,
    #[bincode(with_serde)]
    pub integral: ComplexAccumulator,
    #[bincode(with_serde)]
    pub all_integrals: Vec<ComplexAccumulator>,
    pub stats: StatisticsCounter,
    pub max_stream_id: usize,
    #[bincode(with_serde)]
    pub grid: Grid<F<f64>>,
    pub iter: usize,
}

impl IntegrationState {
    fn new_from_settings<GridGenerator>(create_grid: GridGenerator) -> Self
    where
        GridGenerator: Fn() -> Grid<F<f64>>,
    {
        let num_points = 0;
        let iter = 0;
        //let rng = MonteCarloRng::new(settings.integrator.seed, 0); // in havana_integrate, samples are generated on a single thread
        let grid = create_grid();
        let integral = ComplexAccumulator::new();
        let all_integrals = vec![ComplexAccumulator::new()];
        let stats = StatisticsCounter::new_empty();
        let max_stream_id = 0;

        Self {
            num_points,
            integral,
            all_integrals,
            stats,
            max_stream_id,
            grid,
            iter,
        }
    }

    fn update_iter(&mut self, use_weighted_average: bool) {
        self.integral.update_iter(use_weighted_average);
        self.all_integrals
            .iter_mut()
            .for_each(|acc| acc.update_iter(use_weighted_average));
    }

    // this thing is an unholy mess
    fn display_orientation_results(&self, settings: &RuntimeSettings) {
        if settings.sampling.sample_orientations()
            && let Grid::Discrete(graph_grids) = &self.grid
        {
            for (i_graph, graph_grid) in graph_grids.bins.iter().enumerate() {
                info!("results for graph #{}", i_graph);
                info!("-------------------------");
                if let Grid::Discrete(orientation_grids) = graph_grid.sub_grid.as_ref().unwrap() {
                    let mut sum_all_positive = F(0.0);
                    let mut sum_all_negative = F(0.0);
                    for (i_orientation, orientation_grid) in
                        orientation_grids.bins.iter().enumerate()
                    {
                        print_integral_result(
                            &orientation_grid.accumulator,
                            1,
                            self.iter,
                            &format!("orientation_{i_orientation}"),
                            None,
                        );

                        if orientation_grid.accumulator.avg.positive() {
                            sum_all_positive += orientation_grid.accumulator.avg;
                        } else {
                            sum_all_negative += orientation_grid.accumulator.avg;
                        }
                    }

                    info!("sum_all_positive: {}", sum_all_positive);
                    info!("sum_all_negative: {}", sum_all_negative);
                }
            }
        }
    }
}

pub struct CoreResult {
    pub stats: StatisticsCounter,
    pub integral: ComplexAccumulator,
    pub grid: Grid<F<f64>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IntegrationStatusKind {
    Iteration,
    Final,
}

/// Integrate function used for local runs
pub fn havana_integrate<T, S>(
    settings: &RuntimeSettings,
    model: &Model,
    user_data_generator: T,
    target: Option<Complex<F<f64>>>,
    state: Option<IntegrationState>,
    workspace: Option<PathBuf>,
    show_max_wgt_info_per_iteration: bool,
    show_integration_statistics_per_iteration: bool,
    mut status_emitter: S,
) -> Result<IntegrationResult>
where
    T: Fn(&RuntimeSettings) -> UserData,
    S: FnMut(IntegrationStatusKind, String) -> Result<()>,
{
    let mut user_data = user_data_generator(settings);

    let mut integration_state = if let Some(integration_state) = state {
        integration_state
    } else {
        IntegrationState::new_from_settings(|| user_data.integrand[0].create_grid())
    };

    let sampling_str = settings.sampling.describe_settings();
    let dimension = user_data.integrand[0].get_n_dim();
    let discrete_depth = settings.sampling.discrete_depth();
    let is_tropical_sampling = settings.sampling.get_parameterization_settings().is_none();

    let cont_dim_str = if is_tropical_sampling {
        format!("a median continious dimension of {}", dimension)
    } else {
        format!("{} continuous dimensions", dimension)
    };

    let graph_string = if discrete_depth > 0 {
        let num_graphs = match &integration_state.grid {
            Grid::Discrete(g) => g.bins.len(),
            _ => unreachable!(),
        };
        format!(
            "{discrete_depth} nested discrete grids with {} {} and ",
            num_graphs,
            if num_graphs > 1 { "graphs" } else { "graph" }
        )
    } else {
        String::new()
    };

    let grid_str = format!("{graph_string}{cont_dim_str} using {sampling_str}");

    let cores = user_data.integrand.len();

    let t_start = Instant::now();

    info!(
        "Integrating using {} ltd with {} {} over {} ...",
        if settings.general.use_ltd {
            "naive"
        } else {
            "cff"
        },
        cores,
        if cores > 1 { "cores" } else { "core" },
        grid_str
    );
    info!("");

    let pool = ThreadPoolBuilder::new().num_threads(cores).build().unwrap();

    let mut n_samples_evaluated = 0;
    let mut interrupted = false;
    'integrateLoop: while integration_state.num_points < settings.integrator.n_max {
        // ensure we do not overshoot
        let cur_points = {
            let cur_points_not_final_iter = settings.integrator.n_start
                + settings.integrator.n_increase * integration_state.iter;
            if cur_points_not_final_iter + integration_state.num_points > settings.integrator.n_max
            {
                settings.integrator.n_max - integration_state.num_points
            } else {
                cur_points_not_final_iter
            }
        };

        // the number of points per core is the same for all cores, except for the last one
        let target_points_per_core = (cur_points - 1) / cores + 1;
        let n_points_per_core = repeat_n(target_points_per_core, cores - 1).chain(
            rayon::iter::once(cur_points - target_points_per_core * (cores - 1)),
        );

        let current_max_evals = integration_state.integral.get_worst_case();

        let grids = repeat_n(integration_state.grid.clone_without_samples(), cores);

        let core_results: Vec<Result<CoreResult>> = pool.install(|| {
            user_data
                .integrand
                .par_iter_mut()
                .enumerate()
                .zip(grids)
                .zip(n_points_per_core)
                .map(
                    |(((core_id, integrand), mut grid), n_points)| -> Result<CoreResult> {
                        let mut rng = MonteCarloRng::new(
                            settings.integrator.seed + integration_state.iter as u64,
                            0,
                        );

                        for _ in 0..(target_points_per_core * core_id) {
                            let mut sample = Sample::new();
                            grid.sample(&mut rng, &mut sample);
                        }

                        let samples = (0..n_points)
                            .map(|_| {
                                let mut sample = Sample::new();
                                grid.sample(&mut rng, &mut sample);
                                sample
                            })
                            .collect_vec();

                        let mut core_accumulator = ComplexAccumulator::new();

                        let mut results = Vec::new();
                        for s in samples.iter() {
                            if INTERRUPTED.load(std::sync::atomic::Ordering::Relaxed) {
                                break;
                            }

                            let mut result = integrand.evaluate_sample(
                                s,
                                model,
                                s.get_weight(),
                                integration_state.iter,
                                false,
                                current_max_evals,
                            )?;

                            integrand.process_evaluation_result(&result);
                            maybe_discard_generated_events(&mut result, settings);
                            let jacobian = result.parameterization_jacobian.unwrap_or(F(1.0));
                            let effective_integrand_result =
                                result.integrand_result * Complex::new_re(jacobian);

                            core_accumulator.add_sample(
                                effective_integrand_result,
                                s.get_weight(),
                                Some(s),
                            );

                            let training_eval = match settings.integrator.integrated_phase {
                                IntegratedPhase::Real => effective_integrand_result.re,
                                IntegratedPhase::Imag => effective_integrand_result.im,
                                IntegratedPhase::Both => unimplemented!(),
                            };

                            if let Err(err) = grid.add_training_sample(s, training_eval) {
                                println!("Error adding training sample to grid: {}", err);
                            };

                            results.push(result);
                        }

                        let evaluation_statistics =
                            StatisticsCounter::from_evaluation_results(&results);

                        Ok(CoreResult {
                            stats: evaluation_statistics,
                            integral: core_accumulator,
                            grid,
                        })
                    },
                )
                .collect()
        });
        let core_results: Vec<CoreResult> = core_results.into_iter().collect::<Result<_>>()?;

        if is_interrupted() {
            warn!("{}", "Integration iterrupted by user".yellow());
            interrupted = true;
            break 'integrateLoop;
        }

        for core_result in core_results.iter() {
            integration_state.stats = integration_state.stats.merged(&core_result.stats);
            integration_state.integral.merge(&core_result.integral);
            integration_state.all_integrals[0].merge(&core_result.integral);
        }

        for core_result in core_results.iter() {
            integration_state
                .grid
                .merge(&core_result.grid)
                .expect("could not merge grids");
        }

        integration_state.grid.update(
            F(settings.integrator.discrete_dim_learning_rate),
            F(settings.integrator.continuous_dim_learning_rate),
        );

        integration_state.update_iter(false);

        integration_state.iter += 1;

        integration_state.num_points += cur_points;
        n_samples_evaluated += cur_points;

        // now merge all statistics and observables into the first
        let (first, others) = user_data.integrand.split_at_mut(1);
        for other_itg in others {
            first[0].merge_runtime_results(other_itg)?;
        }

        // now write the observables to disk
        if let Some(itg) = user_data.integrand.first_mut() {
            itg.update_runtime_results(integration_state.iter);
            if settings.integrator.observables_output.per_iteration {
                write_observables_output(
                    itg,
                    settings,
                    workspace.as_deref(),
                    ObservableFlushReason::Iteration(integration_state.iter),
                );
            }
        }

        status_emitter(
            IntegrationStatusKind::Iteration,
            render_iteration_status_block(
                &integration_state,
                cores,
                t_start.elapsed(),
                cur_points,
                n_samples_evaluated,
                &target,
                show_max_wgt_info_per_iteration && settings.integrator.show_max_wgt_info,
                show_integration_statistics_per_iteration,
            ),
        )?;

        // now write the integration state to disk if a workspace has been provided.
        if let Some(ref workspace_path) = workspace {
            let integration_state_path = workspace_path.join("integration_state");

            match fs::write(
                integration_state_path,
                bincode::encode_to_vec(&integration_state, bincode::config::standard())
                    .unwrap_or_else(|_| panic!("failed to serialize the integration state")),
            ) {
                Ok(_) => {}
                Err(_) => warn!("Warning: failed to write integration state to disk"),
            }

            // write the settings to the workspace as well
            let settings_path = workspace_path.join("settings.yaml");
            let settings_string = serde_yaml::to_string(settings)
                .unwrap_or_else(|_| panic!("failed to serialize the settings to a yaml string"));

            match fs::write(settings_path, settings_string) {
                Ok(_) => {}
                Err(_) => warn!("Warning: failed to write settings to disk"),
            }
        }
    }
    // Reset the interrupted flag
    set_interrupted(false);

    if let Some(itg) = user_data.integrand.first() {
        write_observables_output(
            itg,
            settings,
            workspace.as_deref(),
            if interrupted {
                ObservableFlushReason::Interrupted
            } else {
                ObservableFlushReason::Final
            },
        );
    }

    if integration_state.num_points > 0 {
        status_emitter(
            IntegrationStatusKind::Final,
            render_iteration_status_block(
                &integration_state,
                cores,
                t_start.elapsed(),
                0,
                n_samples_evaluated,
                &target,
                true,
                true,
            ),
        )?;
    } else {
        info!("");
        warn!(
            "{}",
            "No final integration results to display since no iteration completed.".yellow()
        );
        info!("");
    }

    integration_state.display_orientation_results(settings);

    Ok(IntegrationResult {
        neval: integration_state.integral.re.processed_samples,
        real_zero: integration_state.integral.re.num_zero_evaluations,
        im_zero: integration_state.integral.im.num_zero_evaluations,
        result: Complex::new(
            integration_state.integral.re.avg,
            integration_state.integral.im.avg,
        ),
        error: Complex::new(
            integration_state.integral.re.err,
            integration_state.integral.im.err,
        ),
        real_chisq: integration_state.integral.re.chi_sq,
        im_chisq: integration_state.integral.im.chi_sq,
    })
}

/// Batch integrate function used for distributed runs, used by the worker nodes.
/// Evaluates a batch of points and returns the results in a manner specified by the user.
pub(crate) fn batch_integrate(
    integrand: &mut Integrand,
    model: &Model,
    input: BatchIntegrateInput,
) -> Result<BatchResult> {
    let samples = match input.samples {
        SampleInput::SampleList { samples } => samples,
        SampleInput::Grid {
            mut grid,
            num_points,
            seed,
            thread_id,
        } => {
            let mut rng = MonteCarloRng::new(seed, thread_id);

            (0..num_points)
                .map(|_| {
                    let mut sample = Sample::new();
                    grid.sample(&mut rng, &mut sample);
                    sample
                })
                .collect_vec()
        }
    };

    let (evaluation_results, metadata_statistics) = evaluate_sample_list(
        integrand,
        &samples,
        model,
        input.num_cores,
        input.iter,
        input.max_eval,
    )?;

    let integrand_output = generate_integrand_output(
        integrand,
        &evaluation_results,
        &samples,
        input.integrand_output_settings,
        input.settings.integrator.integrated_phase,
    );

    let event_output = generate_event_output(
        integrand,
        evaluation_results,
        input.event_output_settings,
        input.settings,
    );

    Ok(BatchResult {
        statistics: metadata_statistics,
        integrand_data: integrand_output,
        event_data: event_output,
    })
}

/// Map the evaluation result on to the right output specified by the user.
fn generate_integrand_output(
    integrand: &Integrand,
    evaluation_results: &[EvaluationResult],
    samples: &[Sample<F<f64>>],
    integrand_output_settings: IntegralOutputSettings,
    integrated_phase: IntegratedPhase,
) -> BatchIntegrateOutput {
    fn effective_integrand_result(result: &EvaluationResult) -> Complex<F<f64>> {
        let jacobian = result.parameterization_jacobian.unwrap_or(F(1.0));
        result.integrand_result * Complex::new_re(jacobian)
    }

    match integrand_output_settings {
        IntegralOutputSettings::Default => {
            let integrand_values = evaluation_results
                .iter()
                .map(effective_integrand_result)
                .collect();

            BatchIntegrateOutput::Default(integrand_values, samples.to_vec())
        }
        IntegralOutputSettings::Accumulator => {
            let mut real_accumulator = StatisticsAccumulator::new();
            let mut imag_accumulator = StatisticsAccumulator::new();
            let mut grid = integrand.create_grid();

            for (result, sample) in evaluation_results.iter().zip(samples.iter()) {
                let effective_result = effective_integrand_result(result);
                real_accumulator
                    .add_sample(effective_result.re * sample.get_weight(), Some(sample));
                imag_accumulator
                    .add_sample(effective_result.im * sample.get_weight(), Some(sample));

                match integrated_phase {
                    IntegratedPhase::Real => {
                        grid.add_training_sample(sample, effective_result.re)
                            .unwrap();
                    }
                    IntegratedPhase::Imag => {
                        grid.add_training_sample(sample, effective_result.im)
                            .unwrap();
                    }
                    IntegratedPhase::Both => {
                        unimplemented!()
                    }
                }
            }

            BatchIntegrateOutput::Accumulator(
                Box::new((real_accumulator, imag_accumulator)),
                Box::new(grid),
            )
        }
    }
}

/// Process events into histograms or event lists, as specified by the user.
fn generate_event_output(
    integrand: &Integrand,
    evaluation_results: Vec<EvaluationResult>,
    event_output_settings: EventOutputSettings,
    _settings: &RuntimeSettings,
) -> EventOutput {
    match event_output_settings {
        EventOutputSettings::None => EventOutput::None,

        EventOutputSettings::EventList => {
            let mut event_groups = EventGroupList::default();
            for mut result in evaluation_results {
                event_groups.append(&mut result.event_groups);
            }
            EventOutput::EventList { event_groups }
        }
        EventOutputSettings::Histogram => integrand
            .observable_accumulator_bundle()
            .map(|histograms| EventOutput::Histogram { histograms })
            .unwrap_or(EventOutput::None),
    }
}

fn maybe_discard_generated_events(result: &mut EvaluationResult, settings: &RuntimeSettings) {
    if !settings.should_return_generated_events() {
        result.event_groups.clear();
    }
}

fn maybe_discard_generated_events_for_integrand(
    result: &mut EvaluationResult,
    integrand: &Integrand,
) {
    if let Integrand::ProcessIntegrand(process_integrand) = integrand {
        maybe_discard_generated_events(result, process_integrand.get_settings());
    }
}

fn observable_output_path(
    workspace: &Path,
    format: ObservableFileFormat,
    reason: ObservableFlushReason,
) -> Option<PathBuf> {
    let extension = match format {
        ObservableFileFormat::None => return None,
        ObservableFileFormat::Hwu => "hwu",
        ObservableFileFormat::Json => "json",
    };

    let file_name = match reason {
        ObservableFlushReason::Iteration(iter) => format!("observables_iter_{iter:04}.{extension}"),
        ObservableFlushReason::Final => format!("observables_final.{extension}"),
        ObservableFlushReason::Interrupted => format!("observables_interrupted.{extension}"),
    };

    Some(workspace.join(file_name))
}

fn write_observables_output(
    integrand: &Integrand,
    settings: &RuntimeSettings,
    workspace: Option<&Path>,
    reason: ObservableFlushReason,
) {
    let format = settings.integrator.observables_output.format;
    let Some(workspace) = workspace else {
        return;
    };
    let Some(path) = observable_output_path(workspace, format, reason) else {
        return;
    };

    if let Err(err) = integrand.write_observable_snapshots(&path, format) {
        warn!(
            "failed to write observables output to {}: {}",
            path.display(),
            err
        );
    }
}

/// This function actually evaluates the list of samples in parallel.
fn evaluate_sample_list(
    integrand: &mut Integrand,
    samples: &[Sample<F<f64>>],
    model: &Model,
    num_cores: usize,
    iter: usize,
    max_eval: Complex<F<f64>>,
) -> Result<(Vec<EvaluationResult>, StatisticsCounter)> {
    // todo!()
    let list_size = samples.len();
    let nvec_per_core = (list_size - 1) / num_cores + 1;

    let sample_chunks = samples.par_chunks(nvec_per_core);
    let integrands = (0..nvec_per_core)
        .map(|_| integrand.clone())
        .collect_vec()
        .into_par_iter();

    let evaluation_results_per_core: Vec<Result<(Vec<EvaluationResult>, Integrand)>> =
        sample_chunks
            .zip(integrands)
            .map(|(chunk, mut integrand)| {
                let evaluation_results = chunk
                    .iter()
                    .map(|sample| {
                        integrand
                            .evaluate_sample(
                                sample,
                                model,
                                sample.get_weight(),
                                iter,
                                false,
                                max_eval,
                            )
                            .map(|mut result| {
                                integrand.process_evaluation_result(&result);
                                maybe_discard_generated_events_for_integrand(
                                    &mut result,
                                    &integrand,
                                );
                                result
                            })
                    })
                    .collect::<Result<Vec<_>>>()?;

                Ok((evaluation_results, integrand))
            })
            .collect();
    let mut evaluation_results_per_core: Vec<(Vec<EvaluationResult>, Integrand)> =
        evaluation_results_per_core
            .into_iter()
            .collect::<Result<_>>()?;

    for (_, worker_integrand) in evaluation_results_per_core.iter_mut() {
        integrand.merge_runtime_results(worker_integrand)?;
    }

    let evaluation_results = evaluation_results_per_core
        .into_iter()
        .flat_map(|(results, _)| results)
        .collect_vec();

    let meta_data_statistics = StatisticsCounter::from_evaluation_results(&evaluation_results);

    Ok((evaluation_results, meta_data_statistics))
}

/// Different ways of passing samples to the batch_integrate function
/// One is simply a list of samples, the other is a grid and the desired number of samples
/// The worker then generates the samples itself.
#[derive(Serialize, Deserialize)]
pub enum SampleInput {
    SampleList {
        samples: Vec<Sample<F<f64>>>,
    },
    Grid {
        grid: Box<Grid<F<f64>>>,
        num_points: usize,
        seed: u64,
        thread_id: usize,
    },
}

type AccumulatorPair = (StatisticsAccumulator<F<f64>>, StatisticsAccumulator<F<f64>>);

#[derive(Serialize, Deserialize)]
pub enum BatchIntegrateOutput {
    Default(Vec<Complex<F<f64>>>, Vec<Sample<F<f64>>>),
    Accumulator(Box<AccumulatorPair>, Box<Grid<F<f64>>>),
}

/// Different ways of processing events, EventList is a list of events, Histogram does accumulation of events on the worker nodes, so the
/// master node only has to merge the histograms.
#[derive(Serialize, Deserialize)]
pub enum EventOutput {
    None,
    EventList {
        event_groups: EventGroupList,
    },
    Histogram {
        histograms: ObservableAccumulatorBundle,
    },
}

/// The result of evaluating a batch of points
#[derive(Serialize, Deserialize, Encode, Decode)]
pub struct BatchResult {
    pub statistics: StatisticsCounter,
    #[bincode(with_serde)]
    pub integrand_data: BatchIntegrateOutput,
    #[bincode(with_serde)]
    pub event_data: EventOutput,
}

/// Input for the batch_integrate function, created by the master node
pub struct BatchIntegrateInput<'a> {
    // global run info:
    pub max_eval: Complex<F<f64>>,
    pub iter: usize,
    pub settings: &'a RuntimeSettings,
    // input data:
    pub samples: SampleInput,
    pub integrand_output_settings: IntegralOutputSettings,
    pub event_output_settings: EventOutputSettings,
    pub num_cores: usize,
}

/// Choose whether to output the integrand results, or a statistics accumulator
#[derive(Serialize, Deserialize)]
pub enum IntegralOutputSettings {
    Default,
    Accumulator,
}

/// Choose whether to output an eventlist, a histogram, or nothing
#[derive(Serialize, Deserialize)]
pub enum EventOutputSettings {
    None,
    EventList,
    Histogram,
}

#[derive(Serialize, Deserialize, Encode, Decode)]
pub struct SerializableBatchIntegrateInput {
    #[bincode(with_serde)]
    pub max_eval: Complex<F<f64>>,
    pub iter: usize,
    #[bincode(with_serde)]
    pub samples: SampleInput,
    #[bincode(with_serde)]
    pub integrand_output_settings: IntegralOutputSettings,
    #[bincode(with_serde)]
    pub event_output_settings: EventOutputSettings,
    pub num_cores: usize,
}

impl SerializableBatchIntegrateInput {
    pub(crate) fn into_batch_integrate_input(
        self,
        settings: &RuntimeSettings,
    ) -> BatchIntegrateInput<'_> {
        BatchIntegrateInput {
            max_eval: self.max_eval,
            iter: self.iter,
            samples: self.samples,
            integrand_output_settings: self.integrand_output_settings,
            event_output_settings: self.event_output_settings,
            num_cores: self.num_cores,
            settings,
        }
    }
}

/// Master node which accumulates data from jobs, and creates the input for jobs
#[derive(Clone)]
pub struct MasterNode {
    grid: Grid<F<f64>>,
    integrator_settings: IntegratorSettings,
    master_accumulator_re: StatisticsAccumulator<F<f64>>,
    master_accumulator_im: StatisticsAccumulator<F<f64>>,
    statistics: StatisticsCounter,
    observable_accumulators: Option<ObservableAccumulatorBundle>,
    event_groups: EventGroupList,
    current_iter: usize,
}

impl MasterNode {
    pub(crate) fn new(grid: Grid<F<f64>>, integrator_settings: IntegratorSettings) -> Self {
        MasterNode {
            grid,
            integrator_settings,
            master_accumulator_im: StatisticsAccumulator::new(),
            master_accumulator_re: StatisticsAccumulator::new(),
            statistics: StatisticsCounter::new_empty(),
            observable_accumulators: None,
            event_groups: EventGroupList::default(),
            current_iter: 0,
        }
    }

    /// Update the grid with the data from another grid.
    fn update_grid_with_grid(&mut self, other_grid: &Grid<F<f64>>) -> Result<(), String> {
        self.grid.merge(other_grid)
    }

    /// Update the grid with the data from a set of samples.
    fn update_grid_with_samples(
        &mut self,
        samples_points: &[Sample<F<f64>>],
        results: &[Complex<F<f64>>],
    ) -> Result<(), String> {
        let integrated_phase = self.integrator_settings.integrated_phase;

        for (sample_point, result) in samples_points.iter().zip(results.iter()) {
            match integrated_phase {
                IntegratedPhase::Real => self.grid.add_training_sample(sample_point, result.re)?,
                IntegratedPhase::Imag => self.grid.add_training_sample(sample_point, result.im)?,
                IntegratedPhase::Both => {
                    unimplemented!("integrated phase both not yet implemented")
                }
            }
        }

        Ok(())
    }

    /// Update the accumulators with the data from another set of accumulators.
    pub(crate) fn update_accumulators_with_accumulators(
        &mut self,
        mut real_accumulator: StatisticsAccumulator<F<f64>>,
        mut imaginary_accumulator: StatisticsAccumulator<F<f64>>,
    ) {
        self.master_accumulator_re
            .merge_samples(&mut real_accumulator);

        self.master_accumulator_im
            .merge_samples(&mut imaginary_accumulator);
    }

    /// Update the accumulators with the data from a set of samples.
    fn update_accumuators_with_samples(
        &mut self,
        sample_points: &[Sample<F<f64>>],
        results: &[Complex<F<f64>>],
    ) {
        for (sample_point, result) in sample_points.iter().zip(results.iter()) {
            self.master_accumulator_re
                .add_sample(result.re * sample_point.get_weight(), Some(sample_point));

            self.master_accumulator_im
                .add_sample(result.im * sample_point.get_weight(), Some(sample_point));
        }
    }

    /// Update the metadata statistics with the data from another set of metadata statistics.
    fn update_metadata_statistics(&mut self, statistics: StatisticsCounter) {
        self.statistics = self.statistics.merged(&statistics);
    }

    /// Finish the current iteration. This should be called after all jobs have been processed.
    pub(crate) fn update_iter(&mut self) {
        self.grid.update(
            F(self.integrator_settings.discrete_dim_learning_rate),
            F(self.integrator_settings.continuous_dim_learning_rate),
        );
        self.master_accumulator_re.update_iter(false);
        self.master_accumulator_im.update_iter(false);
        if let Some(observable_accumulators) = self.observable_accumulators.as_mut() {
            observable_accumulators.update_results();
        }

        self.current_iter += 1;
    }

    /// Write the input for a batch job to a file.
    pub(crate) fn write_batch_input(
        &mut self,
        num_cores: usize,
        num_samples: usize,
        export_grid: bool,
        output_accumulator: bool,
        workspace_path: &str,
        job_id: usize,
    ) -> Result<(), Report> {
        let max_eval = Complex::new(
            self.master_accumulator_re
                .max_eval_positive
                .max(self.master_accumulator_re.max_eval_negative),
            self.master_accumulator_im
                .max_eval_positive
                .max(self.master_accumulator_im.max_eval_negative),
        );

        let samples = if export_grid {
            SampleInput::Grid {
                grid: Box::new(self.grid.clone()),
                num_points: num_samples,
                seed: self.integrator_settings.seed,
                thread_id: job_id,
            }
        } else {
            let mut rng = rand::rng();
            let mut samples_temp = vec![Sample::new(); num_samples];
            for sample in samples_temp.iter_mut() {
                self.grid.sample(&mut rng, sample);
            }
            SampleInput::SampleList {
                samples: samples_temp,
            }
        };

        let integrand_output_settings = if output_accumulator {
            IntegralOutputSettings::Accumulator
        } else {
            IntegralOutputSettings::Default
        };

        let input = SerializableBatchIntegrateInput {
            num_cores,
            max_eval,
            iter: self.current_iter,
            samples,
            integrand_output_settings,
            event_output_settings: EventOutputSettings::None,
        };

        let input_bytes = bincode::encode_to_vec(&input, bincode::config::standard())?;
        let job_name = format!("job_{}", job_id);
        let job_path = std::path::Path::new(workspace_path).join(job_name);

        std::fs::write(job_path, input_bytes)?;

        Ok(())
    }

    /// Process the output of a batch job.
    pub(crate) fn process_batch_output(&mut self, output: BatchResult) -> Result<(), String> {
        self.update_metadata_statistics(output.statistics);

        match output.integrand_data {
            BatchIntegrateOutput::Default(results, samples) => {
                self.update_accumuators_with_samples(&samples, &results);
                self.update_grid_with_samples(&samples, &results)?;
            }
            BatchIntegrateOutput::Accumulator(accumulators, grid) => {
                let (real_accumulator, imag_accumulator) = *accumulators;
                self.update_accumulators_with_accumulators(real_accumulator, imag_accumulator);
                self.update_grid_with_grid(&grid)?;
            }
        }

        match output.event_data {
            EventOutput::None => {}
            EventOutput::EventList { mut event_groups } => {
                self.event_groups.append(&mut event_groups);
            }
            EventOutput::Histogram { mut histograms } => {
                if let Some(existing) = self.observable_accumulators.as_mut() {
                    existing
                        .merge_samples(&mut histograms)
                        .map_err(|err| err.to_string())?;
                } else {
                    self.observable_accumulators = Some(histograms);
                }
            }
        }

        Ok(())
    }

    /// Display the current status of the integration. Usually called after each iteration.
    pub(crate) fn display_status(&self) {
        let status_block = [
            render_integral_result(
                &self.master_accumulator_re,
                1,
                self.current_iter,
                "re",
                None,
            ),
            render_integral_result(
                &self.master_accumulator_im,
                1,
                self.current_iter,
                "im",
                None,
            ),
            self.statistics.render_status_table(),
        ]
        .join("\n");

        info!("\n{status_block}");
    }
}

pub fn emit_integration_status_via_tracing(
    kind: IntegrationStatusKind,
    status_block: impl AsRef<str>,
) -> Result<()> {
    let status_block = status_block.as_ref();
    match kind {
        IntegrationStatusKind::Iteration => {
            info!("\n{status_block}");
            info!("");
        }
        IntegrationStatusKind::Final => {
            info!("");
            info!("{}", "Final integration results:".bold().green());
            info!("");
            info!("\n{status_block}");
            info!("");
        }
    }

    Ok(())
}

fn render_iteration_status_block(
    integration_state: &IntegrationState,
    cores: usize,
    elapsed_time: Duration,
    cur_points: usize,
    n_samples_evaluated: usize,
    target: &Option<Complex<F<f64>>>,
    show_max_wgt_info: bool,
    show_statistics: bool,
) -> String {
    let mut tables = vec![build_iteration_results_table(
        integration_state,
        cores,
        elapsed_time,
        cur_points,
        n_samples_evaluated,
        target,
    )];

    if show_max_wgt_info {
        tables.push(build_max_weight_details_table(
            &integration_state.all_integrals,
        ));
    }

    if show_statistics {
        tables.push(integration_state.stats.build_status_table());
    }

    render_tables_with_shared_width(tables)
}

pub fn print_integral_result(
    itg: &StatisticsAccumulator<F<f64>>,
    i_itg: usize,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) {
    info!("{}", render_integral_result(itg, i_itg, i_iter, tag, trgt));
}

#[allow(clippy::format_in_format_args)]
fn render_integral_result(
    itg: &StatisticsAccumulator<F<f64>>,
    i_itg: usize,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) -> String {
    let cells = build_integral_result_cells(itg, i_itg, i_iter, tag, trgt);
    let delta = match (cells.delta_sigma.as_ref(), cells.delta_percent.as_ref()) {
        (Some(delta_sigma), Some(delta_percent)) => format!("{delta_sigma}, {delta_percent}"),
        _ => String::new(),
    };

    format!(
        "|  {} {} {} {} {} {}",
        cells.integrand, cells.value, cells.relative_error, cells.chi_sq, delta, cells.mwi
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::numerical_integration::ContinuousGrid;

    fn make_accumulator(
        re_avg: f64,
        re_err: f64,
        re_chi_sq: f64,
        im_avg: f64,
        im_err: f64,
        im_chi_sq: f64,
    ) -> ComplexAccumulator {
        let mut accumulator = ComplexAccumulator::new();
        accumulator.re.avg = F(re_avg);
        accumulator.re.err = F(re_err);
        accumulator.re.chi_sq = F(re_chi_sq);
        accumulator.re.processed_samples = 100_000;
        accumulator.re.max_eval_positive = F(1.0);
        accumulator.im.avg = F(im_avg);
        accumulator.im.err = F(im_err);
        accumulator.im.chi_sq = F(im_chi_sq);
        accumulator.im.processed_samples = 100_000;
        accumulator.im.max_eval_positive = F(1.0);
        accumulator
    }

    fn make_integration_state() -> IntegrationState {
        let mut state = IntegrationState::new_from_settings(|| {
            Grid::Continuous(ContinuousGrid::new(1, 64, 100, None, false))
        });
        state.iter = 1;
        state.num_points = 100_000;
        let mut evaluation = EvaluationResult::zero();
        evaluation.evaluation_metadata.integrand_evaluation_time = Duration::from_micros(414);
        evaluation.evaluation_metadata.evaluator_evaluation_time = Duration::from_micros(279);
        evaluation.evaluation_metadata.parameterization_time = Duration::from_nanos(5_800);
        evaluation.evaluation_metadata.total_timing = Duration::from_micros(462);
        state.stats = StatisticsCounter::from_evaluation_results(&[evaluation]);
        state.all_integrals = vec![
            make_accumulator(7.5e-5, 9.8e-5, 0.394, 3.2e-5, 1.5e-5, 0.378),
            make_accumulator(2.1e-5, 2.0e-6, 0.221, -1.7e-5, 3.0e-6, 0.187),
        ];
        state
    }

    #[test]
    fn max_weight_details_table_renders_titled_table_without_wrapping_parentheses() {
        let mut accumulator = ComplexAccumulator::new();
        accumulator.re.max_eval_positive = F(2.3840672847728);
        accumulator.re.max_eval_positive_xs = Some(Sample::Discrete(
            F(1.0),
            0,
            Some(Box::new(Sample::Discrete(
                F(1.0),
                0,
                Some(Box::new(Sample::Continuous(
                    F(1.0),
                    vec![
                        F(0.9695746085826609),
                        F(0.5327714835649985),
                        F(0.003410284786539597),
                    ],
                ))),
            ))),
        ));

        let rendered = render_max_weight_details_table(&[accumulator]);

        assert!(rendered.contains("Maximum weight details"), "{rendered}");
        assert!(rendered.contains("Integrand"), "{rendered}");
        assert!(rendered.contains("Max eval"), "{rendered}");
        assert!(rendered.contains("Max eval coordinates"), "{rendered}");
        assert!(
            rendered.contains("graph: 0, channel: 0, xs: [ "),
            "{rendered}"
        );
        assert!(rendered.contains("e-01"), "{rendered}");
        assert!(rendered.contains("e-03 ]"), "{rendered}");
        assert!(!rendered.contains(", 0.532"), "{rendered}");
        assert!(!rendered.contains("( graph:"), "{rendered}");
    }

    #[test]
    fn iteration_status_block_uses_compact_header_and_optional_statistics() {
        let state = make_integration_state();
        let rendered = render_iteration_status_block(
            &state,
            4,
            Duration::from_secs(1),
            100_000,
            100_000,
            &Some(Complex::new(F(1.0e-4), F(2.0e-5))),
            false,
            false,
        );

        assert!(rendered.contains("Iteration #   1"), "{rendered}");
        assert!(rendered.contains("n_pts = 100K n_tot = 100K"), "{rendered}");
        assert!(rendered.contains("/sample/core"), "{rendered}");
        assert!(rendered.contains("Δ ="), "{rendered}");
        assert!(!rendered.contains("Integration statistics"), "{rendered}");
        let lines = rendered.lines().collect::<Vec<_>>();
        assert!(
            lines.first().is_some_and(|line| !line.contains('┬')),
            "{rendered}"
        );
        assert_eq!(lines[1].matches('│').count(), 2, "{rendered}");
        assert!(lines[2].contains('┬'), "{rendered}");
        assert!(
            lines.last().is_some_and(|line| line.contains('┴')),
            "{rendered}"
        );
    }

    #[test]
    fn iteration_status_block_omits_delta_columns_when_no_target_is_provided() {
        let state = make_integration_state();
        let rendered = render_iteration_status_block(
            &state,
            4,
            Duration::from_secs(1),
            100_000,
            100_000,
            &None,
            false,
            true,
        );

        assert!(!rendered.contains("Δ="), "{rendered}");
        assert!(rendered.contains("Integration statistics"), "{rendered}");
        assert!(rendered.contains("mwi:"), "{rendered}");
    }
}

#[test]
fn test_threading() {
    use symbolica::numerical_integration::ContinuousGrid;

    fn test_fn(x: f64) -> f64 {
        x * x * (x * 6.).sin()
    }

    let mut acc_1 = StatisticsAccumulator::<f64>::new();
    let mut acc_2 = StatisticsAccumulator::<f64>::new();

    let samples_per_sample = 4;
    let samples_per_iter = 100000;
    let n_iter = 10;

    let mut rng = MonteCarloRng::new(42, 0);

    let mut grid = Grid::<f64>::Continuous(ContinuousGrid::new(1, 64, 100, None, false));

    for _i_iter in 0..n_iter {
        let mut multiplice_accs = vec![StatisticsAccumulator::<f64>::new(); samples_per_sample];
        for _i_sample in 0..samples_per_iter {
            let n_samples = (0..samples_per_sample)
                .map(|_| {
                    let mut sample = Sample::new();
                    grid.sample(&mut rng, &mut sample);
                    sample
                })
                .collect_vec();

            let n_evals = n_samples
                .iter()
                .map(|s| match s {
                    Sample::Continuous(_, xs) => test_fn(xs[0]),
                    _ => unreachable!(),
                })
                .collect_vec();

            for (i_eval, (sample, eval)) in n_samples.iter().zip(&n_evals).enumerate() {
                acc_1.add_sample(eval * sample.get_weight(), Some(sample));
                multiplice_accs[i_eval].add_sample(eval * sample.get_weight(), Some(sample));
            }
        }

        for acc in multiplice_accs {
            acc_2.merge_samples(&mut acc.clone());
        }
        acc_1.update_iter(false);
        acc_2.update_iter(false);
    }

    println!("acc1: {:?}", acc_1);
    println!("acc2: {:?}", acc_2);
}
