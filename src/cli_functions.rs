use crate::{
    cross_section::Amplitude,
    inspect::inspect,
    integrands::{integrand_factory, HasIntegrand},
    integrate::{self, havana_integrate, SerializableBatchIntegrateInput, UserData},
    model::Model,
    utils::{print_banner, F, VERSION},
    Integrand, Settings,
};
use clap::{App, Arg, SubCommand};
use color_eyre::Report;
use colored::Colorize;
use eyre::eyre;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use spenso::complex::Complex;
use std::env;
use std::{fs, time::Instant};
use std::{path::PathBuf, str::FromStr};
use symbolica::numerical_integration::Sample;

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
                )
                .arg(
                    Arg::with_name("term")
                        .short("t")
                        .long("term")
                        .required(true)
                        .multiple(true)
                        .allow_hyphen_values(true)
                        .value_name("TERM")
                        .help("Specify the term to inspect"),
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
        .subcommand(
            SubCommand::with_name("batch")
                .about("Evaluate a batch of points, mostly for HPC use")
                .arg(
                    Arg::with_name("process_file")
                        .short("pf")
                        .long("process_file")
                        .required(true)
                        .value_name("PROCESS_FILE")
                        .help("Path to process output file"),
                )
                .arg(
                    Arg::with_name("batch_input_file")
                        .short("if")
                        .long("batch_input_file")
                        .required(true)
                        .value_name("BATCH_INPUT_FILE")
                        .help("Path to the batch input file"),
                )
                .arg(
                    Arg::with_name("name")
                        .short("n")
                        .long("name")
                        .required(true)
                        .value_name("NAME")
                        .help("Name of the amplitude to use"),
                )
                .arg(
                    Arg::with_name("output_name")
                        .short("on")
                        .long("output_name")
                        .required(true)
                        .value_name("NAME")
                        .help("Name of the output file"),
                ),
        )
        .get_matches_from(args);

    crate::set_interrupt_handler();

    if let Some(matches) = matches.subcommand_matches("batch") {
        let path_to_process_output = PathBuf::from_str(matches.value_of("process_file").unwrap())?;
        let path_to_batch_input = PathBuf::from_str(matches.value_of("batch_input_file").unwrap())?;
        let name = matches.value_of("name").unwrap();
        let output_name = matches.value_of("output_name").unwrap();

        return batch_branch(
            path_to_process_output,
            path_to_batch_input,
            name,
            output_name,
        );
    }

    let mut settings: Settings = Settings::from_file(matches.value_of("config").unwrap())?;

    if env::var("SYMBOLICA_LICENSE").is_err() {
        env::set_var("SYMBOLICA_LICENSE", "GAMMALOOP_USER");
    }
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
        target = Some(Complex::new(F(tt[0]), F(tt[1])));
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        if let Some(x) = matches.value_of("debug") {
            settings.general.debug = usize::from_str(x).unwrap();
        }

        let mut integrand = integrand_factory(&settings);

        let pt = matches
            .values_of("point")
            .unwrap()
            .map(|x| F(f64::from_str(x.trim_end_matches(',')).unwrap()))
            .collect::<Vec<_>>();
        let force_radius = matches.is_present("force_radius");
        let term = match matches.values_of("term") {
            Some(t) => t
                .map(|x| usize::from_str(x.trim_end_matches(',')).unwrap())
                .collect::<Vec<_>>(),
            None => vec![],
        };

        let _result = inspect(
            &settings,
            &mut integrand,
            pt.clone(),
            &term,
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
                    F(1.),
                    (0..integrand.get_n_dim())
                        .map(|_i| F(rand::random::<f64>()))
                        .collect(),
                ),
                F(1.),
                1,
                false,
                F(0.0),
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
        let result = havana_integrate(&settings, user_data_generator, target, None, None);

        info!("");
        info!(
            "{}",
            format!(
                "Havana integration completed after {} sample evaluations.",
                format!("{:.2}M", (result.neval as f64) / 1000000.)
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

fn batch_branch(
    process_output_file: PathBuf,
    batch_input_file: PathBuf,
    amplitude_name: &str,
    output_name: &str,
) -> Result<(), Report> {
    // much of this should be moved to the main cli function

    println!("settings passed by command line will be overwritten by configurations in the process output and batch input");

    // load the settings
    let path_to_settings = process_output_file.join("cards").join("run_card.yaml");
    let settings_string = std::fs::read_to_string(path_to_settings.clone())?;
    let settings: Settings = serde_yaml::from_str(&settings_string)?;

    // load the model, hardcoded to scalars.yaml for now
    let path_to_model = process_output_file
        .join("sources")
        .join("model")
        .join("scalars.yaml");

    let path_to_model_string = path_to_model
        .to_str()
        .ok_or_else(|| eyre!("could not convert path to string"))?
        .to_string();

    let model = Model::from_file(path_to_model_string)?;

    // load the amplitude
    let path_to_amplitude_yaml = process_output_file
        .join("sources")
        .join("amplitudes")
        .join(amplitude_name)
        .join("amplitude.yaml");

    // we should change all the file_path arguments to PathBuf or &Path
    let path_to_amplitude_yaml_as_string = path_to_amplitude_yaml.to_str().unwrap().to_string();

    // this is all very amplitude focused, will be generalized later when the structure is clearer
    let amplitude: Amplitude<_> = {
        let amp = Amplitude::from_file(&model, path_to_amplitude_yaml_as_string)?;

        let derived_data_path = process_output_file
            .join("sources")
            .join("amplitudes")
            .join(amplitude_name);

        amp.load_derived_data(&model, &derived_data_path, &settings)?
    };

    // load input data

    let batch_input_bytes = std::fs::read(batch_input_file)?;
    let serializable_batch_input =
        bincode::decode_from_slice::<SerializableBatchIntegrateInput, _>(
            &batch_input_bytes,
            bincode::config::standard(),
        )?
        .0;
    let batch_integrate_input = serializable_batch_input.into_batch_integrate_input(&settings);

    // construct integrand
    let mut integrand =
        Integrand::GammaLoopIntegrand(amplitude.generate_integrand(&path_to_settings)?);

    // integrate
    let batch_result = integrate::batch_integrate(&mut integrand, batch_integrate_input);

    // save result

    let batch_result_bytes = bincode::encode_to_vec(&batch_result, bincode::config::standard())?;
    fs::write(output_name, batch_result_bytes)?;

    Ok(())
}
