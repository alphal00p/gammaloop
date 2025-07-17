use crate::{
    cross_section::Amplitude,
    inspect::inspect,
    integrands::{integrand_factory, HasIntegrand},
    integrate::{self, havana_integrate, SerializableBatchIntegrateInput, UserData},
    model::Model,
    utils::{print_banner, F, VERSION},
    Integrand, Settings,
};
use clap::{Parser, Subcommand};
use color_eyre::Report;
use colored::Colorize;
use eyre::eyre;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use spenso::algebra::complex::Complex;
use std::env;
use std::path::PathBuf;
use std::{fs, time::Instant};
use symbolica::numerical_integration::Sample;
// pub mod repl;
/// Top‑level CLI definition using **clap**'s `derive` API (requires `features = ["derive"]`).
#[derive(Parser, Debug)]
#[command(
    name = "gammaLoop",
    version = VERSION,
    about = "New breed of Local Unitarity implementation",
)]
pub struct Cli {
    /// Number of Rayon worker threads (cores)
    #[arg(short = 'c', long, value_name = "NUMCORES", default_value_t = 1)]
    cores: usize,

    /// Path to the configuration file
    #[arg(
        short = 'f',
        long,
        value_name = "CONFIG_FILE",
        default_value = "./python/gammaloop/data/run_cards/rust_run_config.yaml"
    )]
    config: PathBuf,

    /// Target result to match (real, imag)
    #[arg(short = 't', long, value_name = "TARGET", num_args = 2)]
    target: Option<Vec<f64>>, // length must be 2

    /// Global debug verbosity level
    #[arg(short = 'd', long, value_name = "LEVEL")]
    debug: Option<usize>,

    /// Integrator starting samples
    #[arg(long, value_name = "N_START")]
    n_start: Option<usize>,

    /// Integrator maximum samples
    #[arg(long, value_name = "N_MAX")]
    n_max: Option<usize>,

    /// Integrator samples increase per iteration
    #[arg(long, value_name = "N_INCREASE")]
    n_increase: Option<usize>,

    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

impl Cli {
    pub fn run(self) -> Result<(), Report> {
        crate::set_interrupt_handler();

        // Load settings from YAML first so CLI flags can override.
        let mut settings: Settings = Settings::from_file(&self.config)?;

        // Override settings from CLI top‑level flags --------------------------
        if let Some(level) = self.debug {
            settings.general.debug = level;
        }
        if let Some(n) = self.n_start {
            settings.integrator.n_start = n;
        }
        if let Some(n) = self.n_max {
            settings.integrator.n_max = n;
        }
        if let Some(n) = self.n_increase {
            settings.integrator.n_increase = n;
        }

        // ---------------------------------------------------------------------
        // RAYON THREAD‑POOL SET‑UP -------------------------------------------
        // ---------------------------------------------------------------------
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.cores)
            .build_global()
            .unwrap();

        let num_integrands = self.cores; // keep old logic

        // Parse target if supplied -------------------------------------------
        let target = self.target.as_ref().map(|v| {
            assert_eq!(v.len(), 2, "--target expects exactly two numbers");
            Complex::new(F(v[0]), F(v[1]))
        });

        // Ensure SYMBOLICA licence variable is set before we do anything heavy.
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

        // =====================================================================
        // DISPATCH TO SUB‑COMMANDS OR DEFAULT INTEGRATION PATH  ===============
        // =====================================================================

        match self.command {
            Some(Commands::Batch {
                process_file,
                batch_input_file,
                name,
                output_name,
            }) => {
                return batch_branch(process_file, batch_input_file, &name, &output_name);
            }

            Some(Commands::Inspect {
                point,
                use_f128,
                force_radius,
                momentum_space,
                debug,
                term,
            }) => {
                if let Some(level) = debug {
                    settings.general.debug = level;
                }

                let mut integrand = integrand_factory(&settings);

                let pt = point.into_iter().map(F).collect::<Vec<_>>();

                let _ = inspect(
                    &settings,
                    &mut integrand,
                    pt,
                    &term,
                    force_radius,
                    momentum_space,
                    use_f128,
                );
                Ok(())
            }

            Some(Commands::Bench { samples }) => {
                info!(
                    "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
                    format!("{}", settings.hard_coded_integrand).green(),
                    format!("{}", samples).blue()
                );
                let mut integrand = integrand_factory(&settings);
                let now = Instant::now();
                for _ in 0..samples {
                    integrand.evaluate_sample(
                        &Sample::Continuous(
                            F(1.),
                            (0..integrand.get_n_dim())
                                .map(|_| F(rand::random::<f64>()))
                                .collect(),
                        ),
                        F(1.),
                        1,
                        false,
                        Complex::new_zero(),
                    );
                }
                let total_time = now.elapsed().as_secs_f64();
                info!(
                    "\n> Total time: {} s for {} samples, {} ms per sample\n",
                    format!("{:.1}", total_time).blue(),
                    format!("{}", samples).blue(),
                    format!("{:.5}", total_time * 1000. / (samples as f64)).green(),
                );
                Ok(())
            }

            None | _ => {
                let user_data_generator = |settings: &Settings| UserData {
                    integrand: (0..num_integrands)
                        .map(|_| integrand_factory(settings))
                        .collect(),
                };
                let result = havana_integrate(&settings, user_data_generator, target, None, None);

                info!("");
                info!(
                    "{}",
                    format!(
                        "Havana integration completed after {} sample evaluations.",
                        format!("{:.2}M", (result.neval as f64) / 1_000_000.)
                            .bold()
                            .blue()
                    )
                    .bold()
                    .green()
                );
                info!("");
                Ok(())
            }
        }
    }
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    Repl,

    /// Inspect a single phase‑space point / momentum configuration
    Inspect {
        /// The point to inspect (x y) or (p0 px ...)
        #[arg(short = 'p', num_args = 2.., value_name = "POINT")]
        // allow >2 for momentum‑space
        point: Vec<f64>,

        /// Evaluate in f128 precision
        #[arg(short = 'f', long = "use_f128")]
        use_f128: bool,

        /// Force the radius in the parameterisation
        #[arg(long)]
        force_radius: bool,

        /// Interpret point as momentum‑space coordinates
        #[arg(short = 'm', long)]
        momentum_space: bool,

        /// Override debug level for this inspection
        #[arg(short = 'd', long, value_name = "LEVEL")]
        debug: Option<usize>,

        /// Term(s) to inspect in the integrand
        #[arg(short = 't', long, value_name = "TERM", num_args = 1..)]
        term: Vec<usize>,
    },

    /// Benchmark raw integrand evaluation speed
    Bench {
        /// Number of random samples to evaluate
        #[arg(short = 's', long, value_name = "SAMPLES")]
        samples: usize,
    },

    /// HPC batch evaluation branch
    Batch {
        #[arg(value_name = "PROCESS_FILE")]
        process_file: PathBuf,
        #[arg(value_name = "BATCH_INPUT_FILE")]
        batch_input_file: PathBuf,
        #[arg(short = 'n', long, value_name = "NAME")]
        name: String,
        #[arg(value_name = "NAME")]
        output_name: String,
    },
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
