//! This module steers the integration process.
//! It contains two main ways of integrating the integrand.
//! The havana_integrate function is mostly used for local runs.
//! The master node in combination with batch_integrate is for distributed runs.  

use color_eyre::Report;
use colored::Colorize;
use itertools::Itertools;
use serde::Deserialize;
use serde::Serialize;
use symbolica::numerical_integration::{Grid, MonteCarloRng, Sample, StatisticsAccumulator};

use crate::evaluation_result::EvaluationResult;
use crate::evaluation_result::StatisticsCounter;
use crate::integrands::HasIntegrand;
use crate::observables::Event;
use crate::observables::SerializableEvent;
use crate::utils;
use crate::DiscreteGraphSamplingSettings;
use crate::Integrand;
use crate::IntegratorSettings;
use crate::SamplingSettings;
use crate::Settings;
use crate::{IntegratedPhase, IntegrationResult};
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use num::Complex;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::str::FromStr;
use std::time::Instant;
use tabled::{Style, Table, Tabled};

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

/// Integrate function used for local runs
pub fn havana_integrate<F>(
    settings: &Settings,
    user_data_generator: F,
    target: Option<Complex<f64>>,
) -> crate::IntegrationResult
where
    F: Fn(&Settings) -> UserData,
{
    let mut num_points = 0;
    const N_INTEGRAND_ACCUMULATORS: usize = 2;

    let mut samples = vec![Sample::new(); settings.integrator.n_start];
    let mut f_evals = vec![vec![0.; N_INTEGRAND_ACCUMULATORS]; settings.integrator.n_start];
    let mut integral: StatisticsAccumulator<f64> = StatisticsAccumulator::new();
    let mut all_integrals = vec![StatisticsAccumulator::new(); N_INTEGRAND_ACCUMULATORS];
    let mut evaluation_results = vec![EvaluationResult::zero(); settings.integrator.n_start];
    let mut stats = StatisticsCounter::new_empty();

    let seed = settings.integrator.seed;
    let thread_id = 0; // Samples are generated on a single core in this function.
    let mut rng = MonteCarloRng::new(seed, thread_id);

    let mut user_data = user_data_generator(settings);

    let mut grid = user_data.integrand[0].create_grid();

    let grid_str = match &settings.sampling {
        SamplingSettings::MultiChanneling(_multi_channeling_settings) => {
            let cont_dimension = match &grid {
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
            let cont_dimension = match &grid {
                Grid::Continuous(g) => g.continuous_dimensions.len(),
                _ => unreachable!(),
            };

            format!("a continuous {}-dimensional grid", cont_dimension)
        }
        SamplingSettings::DiscreteGraphs(discrete_graph_sampling_settings) => {
            let num_graphs = match &grid {
                Grid::Discrete(g) => g.bins.len(),
                _ => unreachable!(),
            };

            let inner_settings_string = match discrete_graph_sampling_settings {
                DiscreteGraphSamplingSettings::Default => String::from(""),
                DiscreteGraphSamplingSettings::DiscreteMultiChanneling(_) => {
                    String::from(" and a nested discrete grid over lmb-channels")
                }
                DiscreteGraphSamplingSettings::TropicalSampling => {
                    String::from(" and tropical sampling")
                }
                DiscreteGraphSamplingSettings::MultiChanneling(_) => {
                    String::from(" and multi-channeling over lmb-channels")
                }
            };

            format!(
                "a discrete grid with {} graphs{}",
                num_graphs, inner_settings_string
            )
        }
    };

    let mut iter = 0;

    let cores = user_data.integrand.len();

    let t_start = Instant::now();

    info!("Integrating over a {} ...\n", grid_str);

    while num_points < settings.integrator.n_max {
        let cur_points = settings.integrator.n_start + settings.integrator.n_increase * iter;
        samples.resize(cur_points, Sample::new());
        f_evals.resize(cur_points, vec![0.; N_INTEGRAND_ACCUMULATORS]);
        evaluation_results.resize(cur_points, EvaluationResult::zero());

        for sample in &mut samples[..cur_points] {
            grid.sample(&mut rng, sample);
        }

        // the number of points per core for all cores but the last, which may have fewer
        let nvec_per_core = (cur_points - 1) / cores + 1;

        let current_max_eval = integral
            .max_eval_positive
            .abs()
            .max(integral.max_eval_negative.abs());

        user_data.integrand[..cores]
            .par_iter_mut()
            .zip(f_evals.par_chunks_mut(nvec_per_core))
            .zip(samples.par_chunks(nvec_per_core))
            .zip(evaluation_results.par_chunks_mut(nvec_per_core))
            .for_each(|(((integrand_f, ff), xi), result_list)| {
                for ((f_evals_i, s), result) in
                    ff.iter_mut().zip(xi.iter()).zip(result_list.iter_mut())
                {
                    let fc = integrand_f.evaluate_sample(
                        s,
                        s.get_weight(),
                        iter,
                        false,
                        current_max_eval,
                    );
                    f_evals_i[0] = fc.integrand_result.re;
                    f_evals_i[1] = fc.integrand_result.im;
                    *result = fc;
                }
            });

        for (s, f) in samples[..cur_points].iter().zip(&f_evals[..cur_points]) {
            let sel_f = match settings.integrator.integrated_phase {
                IntegratedPhase::Real => &f[0],
                IntegratedPhase::Imag => &f[1],
                IntegratedPhase::Both => unimplemented!(),
            };
            if let Err(err) = grid.add_training_sample(s, *sel_f) {
                warn!("WARNING: {}", err)
            } else {
                integral.add_sample(*sel_f * s.get_weight(), Some(s));
            }
        }

        let new_meta_data = StatisticsCounter::from_evaluation_results(&evaluation_results);
        stats = stats.merged(&new_meta_data);

        grid.update(settings.integrator.learning_rate);
        integral.update_iter();

        for i_integrand in 0..N_INTEGRAND_ACCUMULATORS {
            for (s, f) in samples[..cur_points].iter().zip(&f_evals[..cur_points]) {
                all_integrals[i_integrand].add_sample(f[i_integrand] * s.get_weight(), Some(s));
            }
            if !all_integrals[i_integrand].update_iter() {
                info!("WARNING: grid update failed, likely due to insufficient number of sample points to be considered.");
            }
        }

        if settings.general.debug > 1 {
            // Currently not available when using havana from symbolica
            // if let Grid::Discrete(g) = &grid {
            //     g.bins[0]
            //         .plot(&format!("grid_disc_it_{}.svg", iter))
            //         .unwrap();
            // }
            let mut tabled_data = vec![];

            tabled_data.push(IntegralResult {
                id: format!("Sum@it#{}", integral.cur_iter),
                n_samples: format!("{}", integral.processed_samples),
                n_samples_perc: format!("{:.3e}%", 100.),
                integral: format!("{:.8e}", integral.avg),
                variance: format!(
                    "{:.8e}",
                    integral.err * ((integral.processed_samples - 1).max(0) as f64).sqrt()
                ),
                err: format!("{:.8e}", integral.err),
                err_perc: format!(
                    "{:.3e}%",
                    (integral.err / (integral.avg.abs()).max(1.0e-99)).abs() * 100.
                ),
                pdf: String::from_str("N/A").unwrap(),
            });
            if let Grid::Discrete(g) = &grid {
                for (i, b) in g.bins.iter().enumerate() {
                    tabled_data.push(IntegralResult {
                        id: format!("chann#{}", i),
                        n_samples: format!("{}", b.accumulator.processed_samples),
                        n_samples_perc: format!(
                            "{:.3e}%",
                            ((b.accumulator.processed_samples as f64)
                                / (integral.processed_samples.max(1) as f64))
                                * 100.
                        ),
                        integral: format!("{:.8e}", b.accumulator.avg),
                        variance: format!(
                            "{:.8e}",
                            b.accumulator.err
                                * ((b.accumulator.processed_samples - 1).max(0) as f64).sqrt()
                        ),
                        err: format!("{:.8e}", b.accumulator.err),
                        err_perc: format!(
                            "{:.3e}%",
                            (b.accumulator.err / (b.accumulator.avg.abs()).max(1.0e-99)).abs()
                                * 100.
                        ),
                        pdf: format!("{:.8e}", b.pdf),
                    });
                }
            }
            let mut f = BufWriter::new(
                File::create(&format!("results_it_{}.txt", iter))
                    .expect("Could not create results file"),
            );
            writeln!(f, "{}", Table::new(tabled_data).with(Style::psql())).unwrap();
        }

        iter += 1;
        num_points += cur_points;

        info!(
            "/  [ {} ] {}: n_pts={:-6.0}K {} {} /sample/core ",
            format!(
                "{:^7}",
                utils::format_wdhms(t_start.elapsed().as_secs() as usize)
            )
            .bold(),
            format!("Iteration #{:-4}", iter).bold().green(),
            cur_points as f64 / 1000.,
            if num_points >= 10_000_000 {
                format!("n_tot={:-7.0}M", num_points as f64 / 1_000_000.)
                    .bold()
                    .green()
            } else {
                format!("n_tot={:-7.0}K", num_points as f64 / 1000.)
                    .bold()
                    .green()
            },
            format!(
                "{:-17.3} ms",
                (((t_start.elapsed().as_secs() as f64) * 1000.) / (num_points as f64))
                    * (cores as f64)
            )
            .bold()
            .blue()
        );

        for i_integrand in 0..(N_INTEGRAND_ACCUMULATORS / 2) {
            print_integral_result(
                &all_integrals[2 * i_integrand],
                i_integrand + 1,
                iter,
                "re",
                if i_integrand == 0 {
                    target.map(|o| o.re).or(None)
                } else {
                    None
                },
            );
            print_integral_result(
                &all_integrals[2 * i_integrand + 1],
                i_integrand + 1,
                iter,
                "im",
                if i_integrand == 0 {
                    target.map(|o| o.im).or(None)
                } else {
                    None
                },
            );
        }
        if settings.integrator.show_max_wgt_info {
            info!("|  -------------------------------------------------------------------------------------------");
            info!(
                "|  {:<16} | {:<23} | {}",
                "Integrand", "Max Eval", "Max Eval xs",
            );
            for i_integrand in 0..(N_INTEGRAND_ACCUMULATORS / 2) {
                for part in 0..=1 {
                    for sgn in 0..=1 {
                        if (if sgn == 0 {
                            all_integrals[2 * i_integrand + part].max_eval_positive
                        } else {
                            all_integrals[2 * i_integrand + part].max_eval_negative
                        }) == 0.
                        {
                            continue;
                        }

                        info!(
                            "|  {:<20} | {:<23} | {}",
                            format!(
                                "itg #{:-3} {} [{}] ",
                                format!("{:<3}", i_integrand + 1),
                                format!("{:<2}", if part == 0 { "re" } else { "im" }).blue(),
                                format!("{:<1}", if sgn == 0 { "+" } else { "-" }).blue()
                            ),
                            format!(
                                "{:+.16e}",
                                if sgn == 0 {
                                    all_integrals[2 * i_integrand + part].max_eval_positive
                                } else {
                                    all_integrals[2 * i_integrand + part].max_eval_negative
                                }
                            ),
                            format!(
                                "( {} )",
                                if let Some(sample) = if sgn == 0 {
                                    &all_integrals[2 * i_integrand + part].max_eval_positive_xs
                                } else {
                                    &all_integrals[2 * i_integrand + part].max_eval_negative_xs
                                } {
                                    match sample {
                                        Sample::Continuous(_w, v) => v
                                            .iter()
                                            .map(|&x| format!("{:.16}", x))
                                            .collect::<Vec<_>>()
                                            .join(", "),
                                        _ => "N/A".to_string(),
                                    }
                                } else {
                                    "N/A".to_string()
                                }
                            )
                        );
                    }
                }
            }
        }

        stats.display_status();

        // now merge all statistics and observables into the first
        let (first, others) = user_data.integrand[..cores].split_at_mut(1);
        for other_itg in others {
            first[0].merge_results(other_itg, iter);
        }

        // now write the observables to disk
        if let Some(itg) = user_data.integrand[..cores].first_mut() {
            itg.update_results(iter);
        }
        info!("");
    }

    IntegrationResult {
        neval: integral.processed_samples as i64,
        fail: integral.num_zero_evaluations as i32,
        result: all_integrals.iter().map(|res| res.avg).collect::<Vec<_>>(),
        error: all_integrals.iter().map(|res| res.err).collect::<Vec<_>>(),
        prob: all_integrals
            .iter()
            .map(|res| res.chi_sq)
            .collect::<Vec<_>>(),
    }
}

/// Batch integrate function used for distributed runs, used by the worker nodes.
/// Evaluates a batch of points and returns the results in a manner specified by the user.
pub fn batch_integrate(integrand: &Integrand, input: BatchIntegrateInput) -> BatchResult {
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
    samples: &[Sample<f64>],
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
    integrand: &Integrand,
    samples: &[Sample<f64>],
    num_cores: usize,
    iter: usize,
    max_eval: f64,
) -> (Vec<EvaluationResult>, StatisticsCounter) {
    let list_size = samples.len();
    let nvec_per_core = (list_size - 1) / num_cores + 1;

    let sample_chunks = samples.par_chunks(nvec_per_core);
    let mut evaluation_results_per_core = Vec::with_capacity(num_cores);

    sample_chunks
        .map(|chunk| {
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
        samples: Vec<Sample<f64>>,
    },
    Grid {
        grid: Grid<f64>,
        num_points: usize,
        seed: u64,
        thread_id: usize,
    },
}

pub enum BatchIntegrateOutput {
    Default(Vec<Complex<f64>>, Vec<Sample<f64>>),
    Accumulator(
        (StatisticsAccumulator<f64>, StatisticsAccumulator<f64>),
        Grid<f64>,
    ),
}

#[derive(Serialize, Deserialize)]
pub enum SerializableBatchIntegrateOutput {
    Default(Vec<(f64, f64)>, Vec<Sample<f64>>),
    Accumulator(
        (StatisticsAccumulator<f64>, StatisticsAccumulator<f64>),
        Grid<f64>,
    ),
}

impl SerializableBatchIntegrateOutput {
    pub fn from_batch_integrate_output(batch_integrate_output: &BatchIntegrateOutput) -> Self {
        match batch_integrate_output {
            BatchIntegrateOutput::Default(integrand_values, samples) => Self::Default(
                integrand_values.iter().map(|c| (c.re, c.im)).collect(),
                samples.clone(),
            ),
            BatchIntegrateOutput::Accumulator((real_accumulator, imag_accumulator), grid) => {
                Self::Accumulator(
                    (real_accumulator.clone(), imag_accumulator.clone()),
                    grid.clone(),
                )
            }
        }
    }

    pub fn into_batch_integrate_output(self) -> BatchIntegrateOutput {
        match self {
            Self::Default(integrand_values, samples) => BatchIntegrateOutput::Default(
                integrand_values
                    .into_iter()
                    .map(|(re, im)| Complex::new(re, im))
                    .collect(),
                samples,
            ),
            Self::Accumulator((real_accumulator, imag_accumulator), grid) => {
                BatchIntegrateOutput::Accumulator((real_accumulator, imag_accumulator), grid)
            }
        }
    }
}

/// Different ways of processing events, EventList is a list of events, Histogram does accumulation of events on the worker nodes, so the
/// master node only has to merge the histograms.
pub enum EventOutput {
    None,
    EventList { events: Vec<Event> },
    Histogram { histograms: Vec<()> }, // placeholder for the actual histograms
}

#[derive(Serialize, Deserialize)]
pub enum SerializableEventOutput {
    None,
    EventList { events: Vec<SerializableEvent> },
    Histogram { histograms: Vec<()> },
}

impl SerializableEventOutput {
    pub fn from_event_output(event_output: &EventOutput) -> Self {
        match event_output {
            EventOutput::None => Self::None,
            EventOutput::EventList { events } => Self::EventList {
                events: events.iter().map(SerializableEvent::from_event).collect(),
            },
            EventOutput::Histogram { histograms } => Self::Histogram {
                histograms: histograms.iter().map(|_| ()).collect(),
            },
        }
    }

    pub fn into_event_output(self) -> EventOutput {
        match self {
            Self::None => EventOutput::None,
            Self::EventList { events } => EventOutput::EventList {
                events: events
                    .into_iter()
                    .map(SerializableEvent::into_event)
                    .collect(),
            },
            Self::Histogram { histograms } => EventOutput::Histogram {
                histograms: histograms.into_iter().map(|_| ()).collect(),
            },
        }
    }
}

/// The result of evaluating a batch of points
pub struct BatchResult {
    pub statistics: StatisticsCounter,
    pub integrand_data: BatchIntegrateOutput,
    pub event_data: EventOutput,
}

#[derive(Serialize, Deserialize)]
pub struct SerializableBatchResult {
    pub statistics: StatisticsCounter,
    pub integrand_data: SerializableBatchIntegrateOutput,
    pub event_data: SerializableEventOutput,
}

impl SerializableBatchResult {
    pub fn from_batch_result(result: BatchResult) -> Self {
        Self {
            statistics: result.statistics,
            integrand_data: SerializableBatchIntegrateOutput::from_batch_integrate_output(
                &result.integrand_data,
            ),
            event_data: SerializableEventOutput::from_event_output(&result.event_data),
        }
    }

    pub fn into_batch_result(self) -> BatchResult {
        BatchResult {
            statistics: self.statistics,
            integrand_data: self.integrand_data.into_batch_integrate_output(),
            event_data: self.event_data.into_event_output(),
        }
    }
}

/// Input for the batch_integrate function, created by the master node
pub struct BatchIntegrateInput<'a> {
    // global run info:
    pub max_eval: f64,
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

#[derive(Serialize, Deserialize)]
pub struct SerializableBatchIntegrateInput {
    pub max_eval: f64,
    pub iter: usize,
    pub samples: SampleInput,
    pub integrand_output_settings: IntegralOutputSettings,
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
    grid: Grid<f64>,
    integrator_settings: IntegratorSettings,
    master_accumulator_re: StatisticsAccumulator<f64>,
    master_accumulator_im: StatisticsAccumulator<f64>,
    statistics: StatisticsCounter,
    current_iter: usize,
}

impl MasterNode {
    pub fn new(grid: Grid<f64>, integrator_settings: IntegratorSettings) -> Self {
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
    fn update_grid_with_grid(&mut self, other_grid: &Grid<f64>) -> Result<(), String> {
        self.grid.merge(other_grid)
    }

    /// Update the grid with the data from a set of samples.
    fn update_grid_with_samples(
        &mut self,
        samples_points: &[Sample<f64>],
        results: &[Complex<f64>],
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
        mut real_accumulator: StatisticsAccumulator<f64>,
        mut imaginary_accumulator: StatisticsAccumulator<f64>,
    ) {
        self.master_accumulator_re
            .merge_samples(&mut real_accumulator);

        self.master_accumulator_im
            .merge_samples(&mut imaginary_accumulator);
    }

    /// Update the accumulators with the data from a set of samples.
    fn update_accumuators_with_samples(
        &mut self,
        sample_points: &[Sample<f64>],
        results: &[Complex<f64>],
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
        self.grid.update(self.integrator_settings.learning_rate);
        self.master_accumulator_re.update_iter();
        self.master_accumulator_im.update_iter();

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
        let integrated_phase = self.integrator_settings.integrated_phase;

        let max_eval = match integrated_phase {
            IntegratedPhase::Real => self.master_accumulator_re.max_eval_positive,
            IntegratedPhase::Imag => self.master_accumulator_im.max_eval_positive,
            IntegratedPhase::Both => {
                unimplemented!("integrated phase both not yet implemented")
            }
        };

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

        let input_bytes = bincode::serialize(&input)?;
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

fn print_integral_result(
    itg: &StatisticsAccumulator<f64>,
    i_itg: usize,
    i_iter: usize,
    tag: &str,
    trgt: Option<f64>,
) {
    info!(
        "|  itg #{:-3} {}: {} {} {} {} {}",
        format!("{:<3}", i_itg),
        format!("{:-2}", tag).blue().bold(),
        format!("{:-19}", utils::format_uncertainty(itg.avg, itg.err))
            .blue()
            .bold(),
        if itg.avg != 0. {
            if (itg.err / itg.avg).abs() > 0.01 {
                format!(
                    "{:-8}",
                    format!("{:.3}%", (itg.err / itg.avg).abs() * 100.).red()
                )
            } else {
                format!(
                    "{:-8}",
                    format!("{:.3}%", (itg.err / itg.avg).abs() * 100.).green()
                )
            }
        } else {
            format!("{:-8}", "")
        },
        if itg.chi_sq / (i_iter as f64) > 5. {
            format!("{:-6.3} χ²/dof", itg.chi_sq / (i_iter as f64)).red()
        } else {
            format!("{:-6.3} χ²/dof", itg.chi_sq / (i_iter as f64)).normal()
        },
        if i_itg == 1 {
            if let Some(t) = trgt {
                if (t - itg.avg).abs() / itg.err > 5.
                    || (t.abs() != 0. && (t - itg.avg).abs() / t.abs() > 0.01)
                {
                    format!(
                        "Δ={:-7.3}σ, Δ={:-7.3}%",
                        (t - itg.avg).abs() / itg.err,
                        if t.abs() > 0. {
                            (t - itg.avg).abs() / t.abs() * 100.
                        } else {
                            0.
                        }
                    )
                    .red()
                } else {
                    format!(
                        "Δ={:-7.3}σ, Δ={:-7.3}%",
                        (t - itg.avg).abs() / itg.err,
                        if t.abs() > 0. {
                            (t - itg.avg).abs() / t.abs() * 100.
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
        if itg.avg.abs() != 0. {
            let mwi = itg.max_eval_negative.abs().max(itg.max_eval_positive.abs())
                / (itg.avg.abs() * (itg.processed_samples as f64));
            if mwi > 1. {
                format!("  mwi: {:-5.3}", mwi).red()
            } else {
                format!("  mwi: {:-5.3}", mwi).normal()
            }
        } else {
            format!("  mwi: {:-5.3}", 0.).normal()
        }
    );
}
