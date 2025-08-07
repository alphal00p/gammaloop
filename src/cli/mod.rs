use crate::{
    cross_section::Amplitude,
    feyngen::GenerationType,
    integrands::{integrand_factory, HasIntegrand},
    integrate::{self, SerializableBatchIntegrateInput},
    model::Model,
    new_cs::{ExportSettings, Process, ProcessDefinition},
    new_graph::Graph,
    utils::{print_banner, F, VERSION},
    Integrand, Settings,
};
use clap::{Parser, Subcommand, ValueEnum};
use clap_repl::{
    reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
    ClapEditor, ReadCommandOutput,
};
use color_eyre::Report;
use color_eyre::Result;
use colored::Colorize;
use console::style;
use dirs::home_dir;
use eyre::eyre;
use log::LevelFilter;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use spenso::algebra::complex::Complex;
use state::{State, LOG_FORMAT, LOG_LEVEL};
use std::{env, ops::ControlFlow, path::Path};
use std::{fs, time::Instant};
use std::{fs::File, path::PathBuf};
use symbolica::numerical_integration::Sample;

#[derive(Parser, Debug)]
#[command(
    name = "gammaLoop",
    version = VERSION,
    about = "New breed of Local Unitarity implementation",
)]
pub struct Cli {
    /// Path to the configuration file
    #[arg(short = 'f', long, default_value = "./gammaloop_state/settings.yaml")]
    config: PathBuf,

    /// Path to the state file
    #[arg(short = 's', long, default_value = "./gammaloop_state")]
    pub state_folder: PathBuf,

    /// Path to the model file
    #[arg(short = 'm', long)]
    pub model_file: Option<PathBuf>,

    /// Save state to file after each call
    #[arg(short = 'S', long, default_value_t = true)]
    save_state: bool,

    // /// Target result to match (real, imag)
    // #[arg(short = 't', long, value_name = "TARGET", num_args = 2)]
    // target: Option<Vec<f64>>, // length must be 2

    // /// Global debug verbosity level
    // #[arg(short = 'd', long, value_name = "LEVEL")]
    // debug: Option<usize>,

    // /// Integrator starting samples
    // #[arg(long, value_name = "N_START")]
    // n_start: Option<usize>,

    // /// Integrator maximum samples
    // #[arg(long, value_name = "N_MAX")]
    // n_max: Option<usize>,

    // /// Integrator samples increase per iteration
    // #[arg(long, value_name = "N_INCREASE")]
    // n_increase: Option<usize>,
    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

pub struct StateSaveSettings {
    pub override_state: bool,
}

impl From<()> for StateSaveSettings {
    fn from(_: ()) -> Self {
        StateSaveSettings {
            override_state: false,
        }
    }
}

impl From<bool> for StateSaveSettings {
    fn from(override_state: bool) -> Self {
        StateSaveSettings { override_state }
    }
}

impl Cli {
    fn override_settings(&mut self, other: Cli) {
        self.config = other.config;
        self.state_folder = other.state_folder;
        self.model_file = other.model_file;
        self.save_state = other.save_state;
    }

    fn run_command(
        &mut self,
        settings: &mut Settings,
        state: &mut State,
    ) -> Result<ControlFlow<StateSaveSettings>, Report> {
        match self.command.as_ref().ok_or(eyre!("missing command"))? {
            Commands::Quit { override_state } => {
                return Ok(ControlFlow::Break((*override_state).into()));
            }
            Commands::Batch {
                process_file,
                batch_input_file,
                name,
                output_name,
            } => {
                return Ok(ControlFlow::Continue(batch_branch(
                    process_file,
                    batch_input_file,
                    &name,
                    &output_name,
                )?));
            }
            Commands::Inspect(inspect) => inspect.run(state)?,
            Commands::Bench { samples } => {
                info!(
                    "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
                    format!("{}", settings.hard_coded_integrand).green(),
                    format!("{}", samples).blue()
                );
                let mut integrand = integrand_factory(&settings);
                let now = Instant::now();
                for _ in 0..*samples {
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
                    format!("{:.5}", total_time * 1000. / (*samples as f64)).green(),
                );
            }
            Commands::Import(s) => match s {
                Import::Amplitude { path } => {
                    let graphs = Graph::from_file(&path, &state.model)?;
                    let name = path.file_stem().unwrap().to_string_lossy().into_owned();
                    let process = Process::from_graph_list(
                        name,
                        graphs,
                        GenerationType::Amplitude,
                        ProcessDefinition::new_empty(),
                        None,
                    )?;

                    state.process_list.add_process(process);
                }
                Import::Model { path } => {
                    state.model = Model::from_file(path)?;
                }
            },
            Commands::Export(s) => match s {
                Export::Dot { path } => {
                    let exp_set = ExportSettings {
                        root_folder: path.clone(),
                    };
                    state.process_list.export_dot(&exp_set, &state.model)?
                }
                Export::Bin => {}
            },
            Commands::Set(s) => match s {
                Set::Log { level, format } => {
                    if let Some(level) = level {
                        *LOG_LEVEL.lock().unwrap() = *level;
                    }

                    *LOG_FORMAT.lock().unwrap() = *format;
                }
                Set::BaseDir { path } => {
                    self.state_folder = path.clone();
                }
            },

            Commands::Generate(g) => g.run(state, settings)?,

            Commands::List(l) => match l {
                List::Integrands { process_id } => {
                    let process = &state.process_list.processes[*process_id];
                    println!("Integrands for process {}:", process_id);
                    for integrand in process.get_integrand_names() {
                        println!("  {}", integrand);
                    }
                }
                List::Processes => {
                    println!("Processes:");
                    for (index, process) in state.process_list.processes.iter().enumerate() {
                        println!("  {}", process.definition.folder_name(&state.model, index));
                    }
                }
                List::Model => {
                    println!("{}", state.model.name)
                }
            },
            _ => {}
        }
        Ok(ControlFlow::Continue(()))
    }

    pub fn run(mut self) {
        let mut state = match State::load(&self.state_folder, self.model_file.clone()) {
            Ok(state) => state,
            Err(e) => {
                warn!("{e} loading default state");
                State::default()
            }
        };

        let mut settings = self.get_settings().unwrap();

        // let mut state = State::load(&cli.state_file,);

        let mut override_state_file = false;

        if let Some(_) = &self.command {
            let _ = self.run_command(&mut settings, &mut state).unwrap();
        } else {
            let prompt = DefaultPrompt {
                left_prompt: DefaultPromptSegment::Basic("γloop".to_owned()),
                ..DefaultPrompt::default()
            };

            // 2. Build the REPL – clap‑repl takes ownership and configures rustyline.
            let mut repl = ClapEditor::<Cli>::builder().with_prompt(Box::new(prompt));

            if let Some(home) = home_dir() {
                repl = repl.with_editor_hook(move |reed| {
                    // Do custom things with `Reedline` instance here
                    reed.with_history(Box::new(
                        FileBackedHistory::with_file(10000, home.join(".gammaLoop_history"))
                            .unwrap(),
                    ))
                })
            }
            let mut r = repl.build();

            loop {
                match r.read_command() {
                    ReadCommandOutput::Command(mut c) => {
                        match c.run_command(&mut settings, &mut state) {
                            Err(e) => eprintln!("{e}"),
                            Ok(ControlFlow::Break(StateSaveSettings { override_state })) => {
                                override_state_file = override_state;
                                break;
                            }
                            _ => {
                                self.override_settings(c);
                            }
                        }
                    }
                    ReadCommandOutput::EmptyLine => (),
                    ReadCommandOutput::ClapError(e) => {
                        e.print().unwrap();
                    }
                    ReadCommandOutput::ShlexError => {
                        println!(
                            "{} input was not valid and could not be processed",
                            style("Error:").red().bold()
                        );
                    }
                    ReadCommandOutput::ReedlineError(e) => {
                        panic!("{e}");
                    }
                    ReadCommandOutput::CtrlC => continue,
                    ReadCommandOutput::CtrlD => break,
                }
            }
        }

        if self.save_state {
            debug!("Saving State");
            match state.save(&self.state_folder, override_state_file, false) {
                Ok(()) => {}
                Err(e) => {
                    eprintln!("Failed to save state: {}", e);
                }
            };

            serde_yaml::to_writer(
                File::create(self.state_folder.join("settings.yaml")).unwrap(),
                &settings,
            )
            .unwrap();
        }
    }

    pub fn get_settings(&self) -> Result<Settings, Report> {
        crate::set_interrupt_handler();

        // Load settings from YAML first so CLI flags can override.
        let settings: Settings = Settings::from_file(&self.config).unwrap_or_default();

        // Override settings from CLI top‑level flags --------------------------
        // if let Some(level) = self.debug {
        //     settings.general.debug = level;
        // }
        // if let Some(n) = self.n_start {
        //     settings.integrator.n_start = n;
        // }
        // if let Some(n) = self.n_max {
        //     settings.integrator.n_max = n;
        // }
        // if let Some(n) = self.n_increase {
        //     settings.integrator.n_increase = n;
        // }

        *LOG_LEVEL.lock().unwrap() = match settings.general.debug {
            0 => LevelFilter::Off,
            1 => LevelFilter::Error,
            2 => LevelFilter::Warn,
            3 => LevelFilter::Info,
            4 => LevelFilter::Debug,
            _ => LevelFilter::Trace,
        };

        println!("LOG_LEVEL: {:?}", LOG_LEVEL.lock().unwrap());

        // Ensure SYMBOLICA licence variable is set before we do anything heavy.
        if env::var("SYMBOLICA_LICENSE").is_err() {
            env::set_var("SYMBOLICA_LICENSE", "GAMMALOOP_USER");
        }

        print_banner();
        Ok(settings)
    }
}

#[derive(Subcommand, Debug)]
pub enum Import {
    Model {
        #[arg(short = 'p')]
        path: PathBuf,
        // format: String,
    },
    Amplitude {
        #[arg(short = 'p')]
        path: PathBuf,
        // format: String,
    },
}

#[derive(Subcommand, Debug)]
pub enum Export {
    Dot {
        #[arg(short = 'p')]
        path: PathBuf,
    },
    Bin,
}

#[derive(Subcommand, Debug)]
pub enum Set {
    BaseDir {
        path: PathBuf,
    },
    Log {
        #[arg(short = 'l')]
        level: Option<LevelFilter>,
        // #[clap(subcommand)]
        #[arg(short, long, value_enum, default_value_t = LogFormat::Long)]
        format: LogFormat,
    },
}

#[derive(Subcommand, Debug)]
pub enum List {
    Model,
    Processes,
    Integrands {
        #[arg(short = 'p')]
        process_id: usize,
    },
}

pub enum Log {
    Level(LevelFilter),
    Format(LogFormat),
}

pub mod generate;
pub mod inspect;

#[derive(Subcommand, Debug)]
pub enum Commands {
    #[clap(subcommand)]
    List(List),
    #[clap(subcommand)]
    Set(Set),
    #[clap(subcommand)]
    Import(Import),
    #[clap(subcommand)]
    Export(Export),

    Generate(generate::Generate),

    /// Quit gammaloop
    Quit {
        #[arg(short = 'o', long, default_value_t = false)]
        override_state: bool,
    },
    /// Inspect a single phase‑space point / momentum configuration
    // #[clap(subcommand)]
    Inspect(inspect::Inspect),

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

pub mod settings;
pub mod state;

#[repr(usize)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, ValueEnum)]
pub enum LogFormat {
    Long,
    Short,
    Min,
    None,
}

// impl FromStr for LogFormat {
//     type Err = Report;
//     fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
//         Ok(match s {
//             "long" => LogFormat::Long,
//             "short" => LogFormat::Short,
//             "min" => LogFormat::Min,
//             "none" => LogFormat::None,
//             _ => Err(eyre!("Invalid log format"))?,
//         })
//     }
// }

impl Commands {}

fn batch_branch(
    process_output_file: impl AsRef<Path>,
    batch_input_file: impl AsRef<Path>,
    amplitude_name: &str,
    output_name: &str,
) -> Result<(), Report> {
    // much of this should be moved to the main cli function

    println!("settings passed by command line will be overwritten by configurations in the process output and batch input");

    // load the settings
    let path_to_settings = process_output_file
        .as_ref()
        .join("cards")
        .join("run_card.yaml");
    let settings_string = std::fs::read_to_string(path_to_settings.clone())?;
    let settings: Settings = serde_yaml::from_str(&settings_string)?;

    // load the model, hardcoded to scalars.yaml for now
    let path_to_model = process_output_file
        .as_ref()
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
        .as_ref()
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
            .as_ref()
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
