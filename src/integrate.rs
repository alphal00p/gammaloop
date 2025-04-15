//! This module steers the integration process.
//! It contains two main ways of integrating the integrand.
//! The havana_integrate function is mostly used for local runs.
//! The master node in combination with batch_integrate is for distributed runs.

use bincode::Decode;
use bincode::Encode;
use color_eyre::Report;
use colored::Colorize;
use itertools::izip;
use itertools::Itertools;
use rayon::iter::repeatn;
use serde::Deserialize;
use serde::Serialize;
use spenso::contraction::IsZero;
use symbolica::domains::float::ConstructibleFloat;
use symbolica::domains::float::Real;
use symbolica::domains::float::RealNumberLike;
use symbolica::domains::float::SingleFloat;
use symbolica::numerical_integration::{Grid, MonteCarloRng, Sample, StatisticsAccumulator};

use crate::disable;
use crate::evaluation_result::EvaluationResult;
use crate::evaluation_result::StatisticsCounter;
use crate::integrands::HasIntegrand;
use crate::observables::Event;
use crate::utils;
use crate::utils::format_sample;
use crate::utils::F;
use crate::DiscreteGraphSamplingSettings;
use crate::Integrand;
use crate::IntegratorSettings;
use crate::SamplingSettings;
use crate::Settings;
use crate::INTERRUPTED;
use crate::{is_interrupted, set_interrupted};
use crate::{IntegratedPhase, IntegrationResult};
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use rayon::prelude::*;
use spenso::complex::Complex;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter;
use std::ops::MulAssign;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Duration;
use std::time::Instant;
use tabled::{Style, Table, Tabled};

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
    #[tabled(err = "err[%]")]
    err_perc: String,
    #[tabled(err = "PDF")]
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
    pub fn new() -> Self {
        Self {
            re: StatisticsAccumulator::new(),
            im: StatisticsAccumulator::new(),
        }
    }

    /// used to escalate the stability test on high evaluation values
    pub fn get_worst_case(&self) -> Complex<F<f64>> {
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
    pub fn add_sample(
        &mut self,
        result: Complex<F<f64>>,
        sample_weight: F<f64>,
        sample: Option<&Sample<F<f64>>>,
    ) {
        self.re.add_sample(result.re * sample_weight, sample);
        self.im.add_sample(result.im * sample_weight, sample);
    }

    pub fn merge(&mut self, other: &Self) {
        self.re.merge_samples_no_reset(&other.re);
        self.im.merge_samples_no_reset(&other.im);
    }

    pub fn update_iter(&mut self, use_weighted_average: bool) {
        self.re.update_iter(use_weighted_average);
        self.im.update_iter(use_weighted_average);
    }

    pub fn format_max_weights(&self, i_itg: usize) -> Vec<String> {
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
                #[allow(clippy::format_in_format_args)]
                Some(format!(
                    "|  {:<20} | {:<23} | {}",
                    format!(
                        "itg #{:-3} {} [{}] ",
                        format!("{:<3}", i_itg + 1),
                        format!("{:<2}", phase_str).blue(),
                        format!("{:<1}", sign_str).blue()
                    ),
                    format!("{:+.16e}", max_eval),
                    format!(
                        "( {} )",
                        if let Some(sample) = max_eval_sample {
                            format_sample(sample)
                        } else {
                            "N/A".to_string()
                        }
                    )
                ))
            } else {
                None
            }
        })
        .collect()
    }
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
    fn display_orientation_results(&self, settings: &Settings) {
        if settings.sampling.sample_orientations() {
            if let Grid::Discrete(graph_grids) = &self.grid {
                for (i_graph, graph_grid) in graph_grids.bins.iter().enumerate() {
                    info!("results for graph #{}", i_graph);
                    info!("-------------------------");
                    if let Grid::Discrete(orientation_grids) = graph_grid.sub_grid.as_ref().unwrap()
                    {
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
}

pub struct CoreResult {
    pub stats: StatisticsCounter,
    pub integral: ComplexAccumulator,
    pub grid: Grid<F<f64>>,
}

/// Integrate function used for local runs
pub fn havana_integrate<T>(
    settings: &Settings,
    user_data_generator: T,
    target: Option<Complex<F<f64>>>,
    state: Option<IntegrationState>,
    workspace: Option<PathBuf>,
) -> crate::IntegrationResult
where
    T: Fn(&Settings) -> UserData,
{
    let mut user_data = user_data_generator(settings);

    let mut integration_state = if let Some(integration_state) = state {
        integration_state
    } else {
        IntegrationState::new_from_settings(|| user_data.integrand[0].create_grid())
    };

    disable! {
        let grid_str = match &settings.sampling {
            SamplingSettings::MultiChanneling(_multi_channeling_settings) => {
                let cont_dimension = match &integration_state.grid {
                    Grid::Continuous(g) => g.continuous_dimensions.len(),
                    _ => unreachable!(),
                };

                // I don't specify the number of channels, because they are different for each graph
                format!(
                    "a continuous {}-dimensional grid with multi-channeling over lmbs",
                    cont_dimension
                )
            }
            SamplingSettings::Default => {
                let cont_dimension = match &integration_state.grid {
                    Grid::Continuous(g) => g.continuous_dimensions.len(),
                    _ => unreachable!(),
                };

                format!("a continuous {}-dimensional grid", cont_dimension)
            }
            SamplingSettings::DiscreteGraphs(discrete_graph_sampling_settings) => {
                let num_graphs = match &integration_state.grid {
                    Grid::Discrete(g) => g.bins.len(),
                    _ => unreachable!(),
                };


                let inner_settings_string = match momentum_sampling_settings {
                    DiscreteGraphSamplingSettings::Default => String::from(""),
                    DiscreteGraphSamplingSettings::DiscreteMultiChanneling(_) => {
                        String::from(" and a nested discrete grid over lmb-channels")
                    }
                    DiscreteGraphSamplingSettings::TropicalSampling(_) => {
                        format!(" and 🌴🥥 {} 🥥🌴", "tropical sampling".green().bold())
                    }
                    DiscreteGraphSamplingSettings::MultiChanneling(_) => {
                        String::from(" and multi-channeling over lmb-channels")
                    }
                };

                format!(
                    "a discrete grid with {} {}{}",
                    num_graphs,
                    if num_graphs > 1 { "graphs" } else { "graph" },
                    inner_settings_string
                )
            }
        };
    }

    let grid_str = "grid info currently unavailable";

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

    let mut n_samples_evaluated = 0;
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
        let n_points_per_core = repeatn(target_points_per_core, cores - 1).chain(
            rayon::iter::once(cur_points - target_points_per_core * (cores - 1)),
        );

        let current_max_evals = integration_state.integral.get_worst_case();

        let grids = repeatn(integration_state.grid.clone(), cores);

        let core_results: Vec<CoreResult> = user_data
            .integrand
            .par_iter_mut()
            .enumerate()
            .zip(grids)
            .zip(n_points_per_core)
            .map(|(((core_id, integrand), mut grid), n_points)| {
                // set the rng for the current core
                let mut rng = MonteCarloRng::new(
                    settings.integrator.seed,
                    cores * integration_state.iter + core_id,
                );

                let samples = (0..n_points)
                    .map(|_| {
                        let mut sample = Sample::new();
                        grid.sample(&mut rng, &mut sample);
                        sample
                    })
                    .collect_vec();

                let mut core_accumulator = ComplexAccumulator::new();

                let results = samples
                    .iter()
                    .map(|s| {
                        let result = integrand.evaluate_sample(
                            s,
                            s.get_weight(),
                            integration_state.iter,
                            false,
                            current_max_evals,
                        );

                        core_accumulator.add_sample(
                            result.integrand_result,
                            s.get_weight(),
                            Some(s),
                        );

                        let training_eval = match settings.integrator.integrated_phase {
                            IntegratedPhase::Real => result.integrand_result.re,
                            IntegratedPhase::Imag => result.integrand_result.im,
                            IntegratedPhase::Both => unimplemented!(),
                        };

                        if let Err(err) = grid.add_training_sample(s, training_eval) {
                            println!("Error adding training sample to grid: {}", err);
                        };

                        result
                    })
                    .take_while(|_| !INTERRUPTED.load(std::sync::atomic::Ordering::Relaxed)) // make sure ctrl+c does it's job
                    .collect_vec();

                let evaluation_statistics = StatisticsCounter::from_evaluation_results(&results);

                CoreResult {
                    stats: evaluation_statistics,
                    integral: core_accumulator,
                    grid,
                }
            })
            .collect();

        if is_interrupted() {
            warn!("{}", "Integration iterrupted by user".yellow());
            break 'integrateLoop;
        }

        for core_result in core_results.iter() {
            integration_state.stats = integration_state.stats.merged(&core_result.stats);
            integration_state.integral.merge(&core_result.integral);
            integration_state.all_integrals[0].merge(&core_result.integral);
        }

        // The grids need to merged by collecting them into the one of the first core, and then overwriting the grid in integration_state
        // This is because we clone the integration state for each core.

        let mut first_grid = core_results[0].grid.clone();
        for core_result in core_results[1..].iter() {
            first_grid
                .merge(&core_result.grid)
                .expect("could not merge grids");
        }

        integration_state.grid = first_grid;

        integration_state.grid.update(
            settings.integrator.discrete_dim_learning_rate,
            settings.integrator.continuous_dim_learning_rate,
        );

        integration_state.update_iter(false);

        integration_state.integral.update_iter(false);

        disable! {
            if settings.general.debug > 1 {
                // Currently not available when using havana from symbolica
                // if let Grid::Discrete(g) = &grid {
                //     g.bins[0]
                //         .plot(&format!("grid_disc_it_{}.svg", iter))
                //         .unwrap();
                // }
                let mut tabled_data = vec![];

                tabled_data.push(IntegralResult {
                    id: format!("Sum@it#{}", integration_state.integral.cur_iter),
                    n_samples: format!("{}", integration_state.integral.processed_samples),
                    n_samples_perc: format!("{:.3e}%", 100.),
                    integral: format!("{:.8e}", integration_state.integral.avg),
                    variance: format!(
                        "{:.8e}",
                        integration_state.integral.err
                            * F::<f64>::new_from_usize(
                                (integration_state.integral.processed_samples - 1).max(0)
                            )
                            .sqrt()
                    ),
                    err: format!("{:.8e}", integration_state.integral.err),
                    err_perc: format!(
                        "{:.3e}%",
                        (integration_state.integral.err
                            / (integration_state.integral.avg.abs()).max(F(1.0e-99)))
                        .abs()
                            * F(100.)
                    ),
                    pdf: String::from_str("N/A").unwrap(),
                });
                if let Grid::Discrete(g) = &integration_state.grid {
                    for (i, b) in g.bins.iter().enumerate() {
                        tabled_data.push(IntegralResult {
                            id: format!("chann#{}", i),
                            n_samples: format!("{}", b.accumulator.processed_samples),
                            n_samples_perc: format!(
                                "{:.3e}%",
                                ((b.accumulator.processed_samples as f64)
                                    / (integration_state.integral.processed_samples.max(1) as f64))
                                    * 100.
                            ),
                            integral: format!("{:.8e}", b.accumulator.avg),
                            variance: format!(
                                "{:.8e}",
                                b.accumulator.err
                                    * F::<f64>::new_from_usize(
                                        (b.accumulator.processed_samples - 1).max(0)
                                    )
                                    .sqrt()
                            ),
                            err: format!("{:.8e}", b.accumulator.err),
                            err_perc: format!(
                                "{:.3e}%",
                                (b.accumulator.err / (b.accumulator.avg.abs()).max(F(1.0e-99))).abs()
                                    * F(100.)
                            ),
                            pdf: format!("{:.8e}", b.pdf),
                        });
                    }
                }
                let mut f = BufWriter::new(
                    File::create(format!("results_it_{}.txt", integration_state.iter))
                        .expect("Could not create results file"),
                );
                writeln!(f, "{}", Table::new(tabled_data).with(Style::psql())).unwrap();
            }
        }

        integration_state.iter += 1;

        integration_state.num_points += cur_points;
        n_samples_evaluated += cur_points;

        // now merge all statistics and observables into the first
        let (first, others) = user_data.integrand.split_at_mut(1);
        for other_itg in others {
            first[0].merge_results(other_itg, integration_state.iter);
        }

        // now write the observables to disk
        if let Some(itg) = user_data.integrand.first_mut() {
            itg.update_results(integration_state.iter);
        }

        show_integration_status(
            &integration_state,
            cores,
            t_start.elapsed(),
            cur_points,
            n_samples_evaluated,
            &target,
            settings.integrator.show_max_wgt_info,
        );
        info!("");

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

    if integration_state.num_points > 0 {
        info!("");
        info!("{}", "Final integration results:".bold().green());
        info!("");
        show_integration_status(
            &integration_state,
            cores,
            t_start.elapsed(),
            0,
            n_samples_evaluated,
            &target,
            true,
        );
        info!("");
    } else {
        info!("");
        warn!(
            "{}",
            "No final integration results to display since no iteration completed.".yellow()
        );
        info!("");
    }

    integration_state.display_orientation_results(settings);

    disable! {
        IntegrationResult {
            neval: integration_state.integral.processed_samples as i64,
            fail: integration_state.integral.num_zero_evaluations as i32,
            result: integration_state
                .all_integrals
                .iter()
                .map(|res| res.avg)
                .collect::<Vec<_>>(),
            error: integration_state
                .all_integrals
                .iter()
                .map(|res| res.err)
                .collect::<Vec<_>>(),
            prob: integration_state
                .all_integrals
                .iter()
                .map(|res| res.chi_sq)
                .collect::<Vec<_>>(),
        }
    }
    todo!()
}

/// Batch integrate function used for distributed runs, used by the worker nodes.
/// Evaluates a batch of points and returns the results in a manner specified by the user.
pub fn batch_integrate(integrand: &mut Integrand, input: BatchIntegrateInput) -> BatchResult {
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
        input.num_cores,
        input.iter,
        input.max_eval,
    );

    let integrand_output = generate_integrand_output(
        integrand,
        &evaluation_results,
        &samples,
        input.integrand_output_settings,
        input.settings.integrator.integrated_phase,
    );

    let event_output = generate_event_output(
        evaluation_results,
        input.event_output_settings,
        input.settings,
    );

    BatchResult {
        statistics: metadata_statistics,
        integrand_data: integrand_output,
        event_data: event_output,
    }
}

/// Map the evaluation result on to the right output specified by the user.
fn generate_integrand_output(
    integrand: &Integrand,
    evaluation_results: &[EvaluationResult],
    samples: &[Sample<F<f64>>],
    integrand_output_settings: IntegralOutputSettings,
    integrated_phase: IntegratedPhase,
) -> BatchIntegrateOutput {
    match integrand_output_settings {
        IntegralOutputSettings::Default => {
            let integrand_values = evaluation_results
                .iter()
                .map(|result| result.integrand_result)
                .collect();

            BatchIntegrateOutput::Default(integrand_values, samples.to_vec())
        }
        IntegralOutputSettings::Accumulator => {
            let mut real_accumulator = StatisticsAccumulator::new();
            let mut imag_accumulator = StatisticsAccumulator::new();
            let mut grid = integrand.create_grid();

            for (result, sample) in evaluation_results.iter().zip(samples.iter()) {
                real_accumulator.add_sample(
                    result.integrand_result.re * sample.get_weight(),
                    Some(sample),
                );
                imag_accumulator.add_sample(
                    result.integrand_result.im * sample.get_weight(),
                    Some(sample),
                );

                match integrated_phase {
                    IntegratedPhase::Real => {
                        grid.add_training_sample(sample, result.integrand_result.re)
                            .unwrap();
                    }
                    IntegratedPhase::Imag => {
                        grid.add_training_sample(sample, result.integrand_result.im)
                            .unwrap();
                    }
                    IntegratedPhase::Both => {
                        unimplemented!()
                    }
                }
            }

            BatchIntegrateOutput::Accumulator((real_accumulator, imag_accumulator), grid)
        }
    }
}

/// Process events into histograms or event lists, as specified by the user.
fn generate_event_output(
    evaluation_results: Vec<EvaluationResult>,
    event_output_settings: EventOutputSettings,
    _settings: &Settings,
) -> EventOutput {
    match event_output_settings {
        EventOutputSettings::None => EventOutput::None,

        EventOutputSettings::EventList => {
            let event_list = evaluation_results
                .into_iter()
                .flat_map(|result| result.event_buffer) // the clone is necessary for the case where we return an accumulator +
                .collect();
            EventOutput::EventList { events: event_list }
        }
        EventOutputSettings::Histogram => todo!(), // process events using the observables
    }
}

/// This function actually evaluates the list of samples in parallel.
fn evaluate_sample_list(
    integrand: &mut Integrand,
    samples: &[Sample<F<f64>>],
    num_cores: usize,
    iter: usize,
    max_eval: Complex<F<f64>>,
) -> (Vec<EvaluationResult>, StatisticsCounter) {
    // todo!()
    let list_size = samples.len();
    let nvec_per_core = (list_size - 1) / num_cores + 1;

    let sample_chunks = samples.par_chunks(nvec_per_core);
    let integrands = (0..nvec_per_core)
        .map(|_| integrand.clone())
        .collect_vec()
        .into_par_iter();

    let mut evaluation_results_per_core = Vec::with_capacity(num_cores);

    sample_chunks
        .zip(integrands)
        .map(|(chunk, mut integrand)| {
            let cor_evals = chunk
                .iter()
                .map(|sample| {
                    integrand.evaluate_sample(sample, sample.get_weight(), iter, false, max_eval)
                })
                .collect_vec();

            cor_evals
        })
        .collect_into_vec(&mut evaluation_results_per_core);

    let evaluation_results = evaluation_results_per_core
        .into_iter()
        .flatten()
        .collect_vec();

    let meta_data_statistics = StatisticsCounter::from_evaluation_results(&evaluation_results);

    (evaluation_results, meta_data_statistics)
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
        grid: Grid<F<f64>>,
        num_points: usize,
        seed: u64,
        thread_id: usize,
    },
}

#[derive(Serialize, Deserialize)]
pub enum BatchIntegrateOutput {
    Default(Vec<Complex<F<f64>>>, Vec<Sample<F<f64>>>),
    Accumulator(
        (StatisticsAccumulator<F<f64>>, StatisticsAccumulator<F<f64>>),
        Grid<F<f64>>,
    ),
}

/// Different ways of processing events, EventList is a list of events, Histogram does accumulation of events on the worker nodes, so the
/// master node only has to merge the histograms.
#[derive(Serialize, Deserialize)]
pub enum EventOutput {
    None,
    EventList { events: Vec<Event> },
    Histogram { histograms: Vec<()> }, // placeholder for the actual histograms
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
    pub settings: &'a Settings,
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
    pub fn into_batch_integrate_input(self, settings: &Settings) -> BatchIntegrateInput<'_> {
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
    current_iter: usize,
}

impl MasterNode {
    pub fn new(grid: Grid<F<f64>>, integrator_settings: IntegratorSettings) -> Self {
        MasterNode {
            grid,
            integrator_settings,
            master_accumulator_im: StatisticsAccumulator::new(),
            master_accumulator_re: StatisticsAccumulator::new(),
            statistics: StatisticsCounter::new_empty(),
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
    pub fn update_accumulators_with_accumulators(
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
    pub fn update_iter(&mut self) {
        self.grid.update(
            self.integrator_settings.discrete_dim_learning_rate,
            self.integrator_settings.continuous_dim_learning_rate,
        );
        self.master_accumulator_re.update_iter(false);
        self.master_accumulator_im.update_iter(false);

        self.current_iter += 1;
    }

    /// Write the input for a batch job to a file.
    pub fn write_batch_input(
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
                grid: self.grid.clone(),
                num_points: num_samples,
                seed: self.integrator_settings.seed,
                thread_id: job_id,
            }
        } else {
            let mut rng = rand::thread_rng();
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
    pub fn process_batch_output(&mut self, output: BatchResult) -> Result<(), String> {
        self.update_metadata_statistics(output.statistics);

        match output.integrand_data {
            BatchIntegrateOutput::Default(results, samples) => {
                self.update_accumuators_with_samples(&samples, &results);
                self.update_grid_with_samples(&samples, &results)?;
            }
            BatchIntegrateOutput::Accumulator((real_accumulator, imag_accumulator), grid) => {
                self.update_accumulators_with_accumulators(real_accumulator, imag_accumulator);
                self.update_grid_with_grid(&grid)?;
            }
        }

        match output.event_data {
            EventOutput::None => {}
            _ => {
                todo!("event processing not yet implemented");
            }
        }

        Ok(())
    }

    /// Display the current status of the integration. Usually called after each iteration.
    pub fn display_status(&self) {
        print_integral_result(
            &self.master_accumulator_re,
            1,
            self.current_iter,
            "re",
            None,
        );

        print_integral_result(
            &self.master_accumulator_im,
            1,
            self.current_iter,
            "im",
            None,
        );

        self.statistics.display_status();
    }
}

#[allow(clippy::format_in_format_args)]
pub fn show_integration_status(
    integration_state: &IntegrationState,
    cores: usize,
    elapsed_time: Duration,
    cur_points: usize,
    n_samples_evaluated: usize,
    target: &Option<Complex<F<f64>>>,
    show_max_wgt_info: bool,
) {
    info!(
        "/  [ {} ] {}: n_pts={:-6.0}K {} {}",
        format!(
            "{:^7}",
            utils::format_wdhms(elapsed_time.as_secs() as usize)
        )
        .bold(),
        format!("Iteration #{:-4}", integration_state.iter)
            .bold()
            .green(),
        cur_points as f64 / 1000.,
        if integration_state.num_points >= 10_000_000 {
            format!(
                "n_tot={:-7.0}M",
                integration_state.num_points as f64 / 1_000_000.
            )
            .bold()
            .green()
        } else {
            format!(
                "n_tot={:-7.0}K",
                integration_state.num_points as f64 / 1000.
            )
            .bold()
            .green()
        },
        format!(
            "{:-22}",
            format!(
                "{:-9} /sample/core",
                if n_samples_evaluated == 0 {
                    "N/A".red()
                } else {
                    utils::format_evaluation_time_from_f64(
                        elapsed_time.as_secs_f64() / (n_samples_evaluated as f64) * (cores as f64),
                    )
                    .bold()
                    .blue()
                }
            )
        )
    );

    for (i_integrand, integrand) in integration_state.all_integrals.iter().enumerate() {
        print_integral_result(
            &integrand.re,
            i_integrand + 1,
            integration_state.iter,
            "re",
            if i_integrand == 0 {
                target.map(|o| o.re).or(None)
            } else {
                None
            },
        );
        print_integral_result(
            &integrand.im,
            i_integrand + 1,
            integration_state.iter,
            "im",
            if i_integrand == 0 {
                target.map(|o| o.im).or(None)
            } else {
                None
            },
        );
    }
    if show_max_wgt_info {
        info!("|  -------------------------------------------------------------------------------------------");
        info!(
            "|  {:<16} | {:<23} | {}",
            "Integrand", "Max Eval", "Max Eval xs",
        );

        for (i_itg, integral) in integration_state.all_integrals.iter().enumerate() {
            for max_weight_str in integral.format_max_weights(i_itg) {
                info!("{}", max_weight_str);
            }
        }

        info!("|  -------------------------------------------------------------------------------------------");
    }

    integration_state.stats.display_status();
}

pub fn print_integral_result(
    itg: &StatisticsAccumulator<F<f64>>,
    i_itg: usize,
    i_iter: usize,
    tag: &str,
    trgt: Option<F<f64>>,
) {
    info!(
        "|  itg #{:-3} {}: {} {} {} {} {}",
        format!("{:<3}", i_itg),
        format!("{:-2}", tag).blue().bold(),
        format!("{:-19}", utils::format_uncertainty(itg.avg, itg.err))
            .blue()
            .bold(),
        if itg.avg.is_non_zero() {
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
        },
        if itg.chi_sq / (F::<f64>::new_from_usize(i_iter)) > F(5.) {
            format!("{:-6.3} χ²/dof", itg.chi_sq.0 / (i_iter as f64)).red()
        } else {
            format!("{:-6.3} χ²/dof", itg.chi_sq.0 / (i_iter as f64)).normal()
        },
        if i_itg == 1 {
            if let Some(t) = trgt {
                if (t - itg.avg).abs().0 / itg.err.0 > 5.
                    || (t.abs().is_non_zero() && ((t - itg.avg).abs() / t.abs()).0 > 0.01)
                {
                    format!(
                        "Δ={:-7.3}σ, Δ={:-7.3}%",
                        (t - itg.avg).abs().0 / itg.err.0,
                        if t.abs() > F(0.) {
                            (t - itg.avg).abs().0 / t.abs().0 * 100.
                        } else {
                            0.
                        }
                    )
                    .red()
                } else {
                    format!(
                        "Δ={:-7.3}σ, Δ={:-7.3}%",
                        if t.is_non_zero() && itg.avg.is_non_zero() {
                            (t - itg.avg).abs().0 / itg.err.0
                        } else {
                            0.
                        },
                        if t.abs() > F(0.) {
                            (t - itg.avg).abs().0 / t.abs().0 * 100.
                        } else {
                            0.
                        }
                    )
                    .green()
                }
            } else {
                "".to_string().normal()
            }
        } else {
            "".to_string().normal()
        },
        if itg.avg.abs().0 != 0. {
            let mwi = itg.max_eval_negative.abs().max(itg.max_eval_positive.abs())
                / (itg.avg.abs() * (F::<f64>::new_from_usize(itg.processed_samples)));
            if mwi > F(1.) {
                format!("  mwi: {:<10.4e}", mwi.0).red()
            } else {
                format!("  mwi: {:<10.4e}", mwi.0).normal()
            }
        } else {
            format!("  mwi: {:<10.4e}", 0.).normal()
        }
    );
}
