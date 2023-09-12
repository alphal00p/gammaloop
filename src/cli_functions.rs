use clap::{App, Arg, SubCommand};
use color_eyre::Report;
use colored::Colorize;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use symbolica::numerical_integration::Sample;

use crate::{
    inspect::inspect,
    integrands::integrand_factory,
    integrands::HasIntegrand,
    integrate::{havana_integrate, UserData},
    utils::{print_banner, VERSION},
    Settings,
};
use num::Complex;
use std::str::FromStr;
use std::time::Instant;

pub fn cli(args: &Vec<String>) -> Result<(), Report> {
    let matches = App::new("gammaLoop")
        .version(VERSION)
        .about("New breed of Local Unitarity implementation")
        .arg(
            Arg::with_name("cores")
                .short("c")
                .long("cores")
                .value_name("NUMCORES")
                .help("Set the number of cores"),
        )
        .arg(
            Arg::with_name("config")
                .short("f")
                .long("config")
                .value_name("CONFIG_FILE")
                .default_value("./python/gammaloop/data/run_cards/rust_run_config.yaml")
                .help("Set the configuration file"),
        )
        .arg(
            Arg::with_name("target")
                .short("t")
                .long("target")
                .multiple(true)
                .allow_hyphen_values(true)
                .value_name("TARGET")
                .help("Specify the integration target a <real> <imag>"),
        )
        .arg(
            Arg::with_name("debug")
                .short("d")
                .long("debug")
                .value_name("LEVEL")
                .help("Set the debug level. Higher means more verbose."),
        )
        .arg(
            Arg::with_name("n_start")
                .long("n_start")
                .value_name("N_START")
                .help("Number of starting samples for the integrator"),
        )
        .arg(
            Arg::with_name("n_max")
                .long("n_max")
                .value_name("N_MAX")
                .help("Max number of starting samples to consider for integration"),
        )
        .arg(
            Arg::with_name("n_increase")
                .long("n_increase")
                .value_name("N_INCREASE")
                .help("Increase of number of sample points for each successive iteration"),
        )
        .subcommand(
            SubCommand::with_name("inspect")
                .about("Inspect a single input point")
                .arg(
                    Arg::with_name("point")
                        .short("p")
                        .required(true)
                        .min_values(2)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("use_f128")
                        .short("f128")
                        .long("use_f128")
                        .help("Use f128 evaluation"),
                )
                .arg(
                    Arg::with_name("force_radius")
                        .long("force_radius")
                        .help("force radius in parameterisation"),
                )
                .arg(
                    Arg::with_name("momentum_space")
                        .short("m")
                        .long("momentum_space")
                        .help("Set if the point is specified in momentum space"),
                )
                .arg(
                    Arg::with_name("debug")
                        .short("d")
                        .long("debug")
                        .value_name("LEVEL")
                        .help("Set the debug level. Higher means more verbose."),
                ),
        )
        .subcommand(
            SubCommand::with_name("bench")
                .about("Benchmark timing for individual evaluations of the integrand")
                .arg(
                    Arg::with_name("samples")
                        .required(true)
                        .long("samples")
                        .short("s")
                        .value_name("SAMPLES")
                        .help("Number of samples for benchmark"),
                ),
        )
        .get_matches_from(args);

    let mut settings: Settings = Settings::from_file(matches.value_of("config").unwrap())?;

    print_banner();
    if settings.general.debug > 0 {
        info!(
            "{}",
            format!("Debug mode enabled at level {}", settings.general.debug).red()
        );
        info!("");
    }

    if let Some(x) = matches.value_of("debug") {
        settings.general.debug = usize::from_str(x).unwrap();
    }
    if let Some(x) = matches.value_of("n_start") {
        settings.integrator.n_start = usize::from_str(x).unwrap();
    }
    if let Some(x) = matches.value_of("n_max") {
        settings.integrator.n_max = usize::from_str(x).unwrap();
    }
    if let Some(x) = matches.value_of("n_increase") {
        settings.integrator.n_increase = usize::from_str(x).unwrap();
    }

    let mut cores = 1;
    if let Some(x) = matches.value_of("cores") {
        cores = usize::from_str(x).unwrap();
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(cores)
        .build_global()
        .unwrap();

    let num_integrands = cores;

    let mut target = None;
    if let Some(t) = matches.values_of("target") {
        let tt: Vec<_> = t
            .map(|x| f64::from_str(x.trim_end_matches(',')).unwrap())
            .collect();
        if tt.len() != 2 {
            panic!("Expected two numbers for target");
        }
        target = Some(Complex::new(tt[0], tt[1]));
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        if let Some(x) = matches.value_of("debug") {
            settings.general.debug = usize::from_str(x).unwrap();
        }

        let mut integrand = integrand_factory(&settings);

        let pt = matches
            .values_of("point")
            .unwrap()
            .map(|x| f64::from_str(x.trim_end_matches(',')).unwrap())
            .collect::<Vec<_>>();
        let force_radius = matches.is_present("force_radius");

        let _result = inspect(
            &settings,
            &mut integrand,
            pt.clone(),
            force_radius,
            matches.is_present("momentum_space"),
            matches.is_present("use_f128"),
        );
    } else if let Some(matches) = matches.subcommand_matches("bench") {
        let n_samples: usize = matches.value_of("samples").unwrap().parse().unwrap();
        info!(
            "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
            format!("{}", settings.hard_coded_integrand).green(),
            format!("{}", n_samples).blue()
        );
        let mut integrand = integrand_factory(&settings);
        let now = Instant::now();
        for _i in 1..n_samples {
            integrand.evaluate_sample(
                &Sample::Continuous(
                    1.,
                    (0..integrand.get_n_dim())
                        .map(|_i| rand::random::<f64>())
                        .collect(),
                ),
                1.,
                1,
                false,
            );
        }
        let total_time = now.elapsed().as_secs_f64();
        info!(
            "\n> Total time: {} s for {} samples, {} ms per sample\n",
            format!("{:.1}", total_time).blue(),
            format!("{}", n_samples).blue(),
            format!("{:.5}", total_time * 1000. / (n_samples as f64)).green(),
        );
    } else {
        let user_data_generator = |settings: &Settings| UserData {
            integrand: (0..num_integrands)
                .map(|_i| integrand_factory(settings))
                .collect(),
        };
        let result = havana_integrate(&settings, user_data_generator, target);

        info!("");
        info!(
            "{}",
            format!(
                "Havana integration completed after {} sample evaluations.",
                format!("{:.2}M", (result.neval as f64) / (1000000. as f64))
                    .bold()
                    .blue()
            )
            .bold()
            .green()
        );
        info!("");
    }
    Ok(())
}
