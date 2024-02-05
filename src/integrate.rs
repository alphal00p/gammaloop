use colored::Colorize;
use itertools::Itertools;
use symbolica::numerical_integration::{Grid, Sample, StatisticsAccumulator};

use crate::evaluation_result::EvaluationResult;
use crate::evaluation_result::MetaDataStatistics;
use crate::integrands::HasIntegrand;
use crate::observables::Event;
use crate::utils;
use crate::Integrand;
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

    let mut rng = rand::thread_rng();

    let mut user_data = user_data_generator(settings);

    let mut grid = user_data.integrand[0].create_grid();

    let grid_str = match &grid {
        Grid::Discrete(g) => format!(
            "top-level discrete {}-dimensional grid",
            format!("{}", g.bins.len()).bold().blue()
        ),
        Grid::Continuous(g) => {
            format!(
                "top-level continuous {}-dimensional grid",
                format!("{}", g.continuous_dimensions.len()).bold().blue()
            )
        }
    };

    let mut iter = 0;

    let cores = user_data.integrand.len();

    let t_start = Instant::now();
    info!(
        "gammaloop now integrates '{}' over a {} ...\n",
        format!("{}", settings.hard_coded_integrand).green(),
        grid_str
    );
    while num_points < settings.integrator.n_max {
        let cur_points = settings.integrator.n_start + settings.integrator.n_increase * iter;
        samples.resize(cur_points, Sample::new());
        f_evals.resize(cur_points, vec![0.; N_INTEGRAND_ACCUMULATORS]);

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
            .for_each(|((integrand_f, ff), xi)| {
                for (f_evals_i, s) in ff.iter_mut().zip(xi.iter()) {
                    let fc = integrand_f.evaluate_sample(
                        s,
                        s.get_weight(),
                        iter,
                        false,
                        current_max_eval,
                    );
                    f_evals_i[0] = fc.integrand_result.re;
                    f_evals_i[1] = fc.integrand_result.im;
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

pub fn batch_integrate(integrand: &Integrand, input: BatchIntegrateInput) -> BatchResult {
    let samples_generated_locally = matches!(
        &input.samples,
        BatchInput::Grid {
            grid: _,
            num_points: _
        }
    );

    match input.samples {
        BatchInput::SampleList { samples } => {
            let (evaluation_results, metadata_statistics) = evaluate_sample_list(
                integrand,
                samples,
                input.num_cores,
                input.iter,
                input.max_eval,
            );

            let integrand_output = generate_integrand_output(
                integrand,
                &evaluation_results,
                samples,
                input.integrand_output_settings,
                input.settings.integrator.integrated_phase,
                samples_generated_locally,
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
        BatchInput::Grid {
            mut grid,
            num_points,
        } => {
            // generate samples, then evalute them
            let mut rng = rand::thread_rng();
            let samples = (0..num_points)
                .map(|_| {
                    let mut sample = Sample::new();
                    grid.sample(&mut rng, &mut sample);
                    sample
                })
                .collect_vec();

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
                samples_generated_locally,
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
    }
}

fn generate_integrand_output(
    integrand: &Integrand,
    evaluation_results: &[EvaluationResult],
    samples: &[Sample<f64>],
    integrand_output_settings: IntegralOutputSettings,
    integrated_phase: IntegratedPhase,
    samples_generated_locally: bool,
) -> BatchIntegrateOutput {
    match integrand_output_settings {
        IntegralOutputSettings::Default => {
            let integrand_values = evaluation_results
                .iter()
                .map(|result| result.integrand_result)
                .collect();

            // if the points are generated by the grid, we need to pass them back to the master node if they are not being learned locally
            // In general, this combination of settings will probably not be used.
            let samples = if samples_generated_locally {
                Some(samples.to_vec())
            } else {
                None
            };

            BatchIntegrateOutput::Default(integrand_values, samples)
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

// core function of the batch_integrate function, used by all versions of IO
fn evaluate_sample_list(
    integrand: &Integrand,
    samples: &[Sample<f64>],
    num_cores: usize,
    iter: usize,
    max_eval: f64,
) -> (Vec<EvaluationResult>, MetaDataStatistics) {
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

    let meta_data_statistics = MetaDataStatistics::from_evaluation_results(&evaluation_results);

    (evaluation_results, meta_data_statistics)
}

pub enum BatchInput<'a> {
    SampleList { samples: &'a [Sample<f64>] },
    Grid { grid: Grid<f64>, num_points: usize },
}

pub enum BatchIntegrateOutput {
    Default(Vec<Complex<f64>>, Option<Vec<Sample<f64>>>),
    Accumulator(
        (StatisticsAccumulator<f64>, StatisticsAccumulator<f64>),
        Grid<f64>,
    ),
}

pub enum EventOutput {
    None,
    EventList { events: Vec<Event> },
    Histogram { histograms: Vec<()> }, // placeholder for the actual histograms
}

pub struct BatchResult {
    pub statistics: MetaDataStatistics,
    pub integrand_data: BatchIntegrateOutput,
    pub event_data: EventOutput,
}

pub struct BatchIntegrateInput<'a> {
    // global run info:
    pub max_eval: f64,
    pub iter: usize,
    pub settings: &'a Settings,
    // input data:
    pub samples: BatchInput<'a>,
    pub integrand_output_settings: IntegralOutputSettings,
    pub event_output_settings: EventOutputSettings,
    pub num_cores: usize,
}

pub enum IntegralOutputSettings {
    Default,
    Accumulator,
}

pub enum EventOutputSettings {
    None,
    EventList,
    Histogram,
}
