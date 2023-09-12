use colored::Colorize;
use symbolica::numerical_integration::{Grid, Sample, StatisticsAccumulator};

use crate::integrands::HasIntegrand;
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
    let mut integral = StatisticsAccumulator::new();
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

        user_data.integrand[..cores]
            .par_iter_mut()
            .zip(f_evals.par_chunks_mut(nvec_per_core))
            .zip(samples.par_chunks(nvec_per_core))
            .for_each(|((integrand_f, ff), xi)| {
                for (f_evals_i, s) in ff.iter_mut().zip(xi.iter()) {
                    let fc = integrand_f.evaluate_sample(s, s.get_weight(), iter, false);
                    f_evals_i[0] = fc.re;
                    f_evals_i[1] = fc.im;
                }
            });

        for (s, f) in samples[..cur_points].iter().zip(&f_evals[..cur_points]) {
            let sel_f = match settings.integrator.integrated_phase {
                IntegratedPhase::Real => &f[0],
                IntegratedPhase::Imag => &f[1],
                IntegratedPhase::Both => unimplemented!(),
            };
            grid.add_training_sample(s, *sel_f).unwrap();
            integral.add_sample(*sel_f * s.get_weight(), Some(s));
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
            writeln!(
                f,
                "{}",
                Table::new(tabled_data).with(Style::psql()).to_string()
            )
            .unwrap();
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
                        format!("{}", "").normal()
                    }
                } else {
                    format!("{}", "").normal()
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
                                    format!("{}", "N/A")
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
