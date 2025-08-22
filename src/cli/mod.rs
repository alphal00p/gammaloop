use crate::{
    feyngen::GenerationType,
    graph::Graph,
    integrands::{integrand_factory, HasIntegrand},
    model::Model,
    processes::{ExportSettings, Process, ProcessDefinition},
    utils::{F, GIT_VERSION, GS},
};
use chrono::{Datelike, Local, Timelike};
use clap::{Args, Parser, Subcommand, ValueEnum};
use clap_repl::{
    reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
    ClapEditor, ReadCommandOutput,
};
use color_eyre::Result;
use color_eyre::{config::HookBuilder, Report};
use colored::Colorize;
use console::style;
use dirs::home_dir;
use integrate::Integrate;
use log::LevelFilter;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use state::{format_level, format_target, RunHistory, State, LOG_FORMAT, LOG_LEVEL};
use std::ops::ControlFlow;
use std::path::PathBuf;
use std::time::Instant;
use symbolica::numerical_integration::Sample;

pub mod generate;
pub mod inspect;
pub mod integrate;
pub mod settings;
pub mod state;

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
pub struct Cli {
    /// Path to the run-time settings file, by default is: `./gammaloop_state/runtime_settings.yaml`
    #[arg(short = 'r', long)]
    run_settings_path: Option<PathBuf>,

    /// Path to the gammaloop settings file, by default is: `./gammaloop_state/global_settings.yaml`
    #[arg(short = 'g', long)]
    global_settings_path: Option<PathBuf>,

    /// Path to the a run-history file
    #[arg(short = 'i', long)]
    run_history: Option<PathBuf>,

    /// Path to the state file
    #[arg(short = 's', long, default_value = "./gammaloop_state")]
    pub state_folder: PathBuf,

    /// Path to the model file
    #[arg(short = 'm', long)]
    pub model_file: Option<PathBuf>,

    /// Skip the
    #[arg(short = 'n', long, default_value_t = false, group = "saving")]
    no_save_state: bool,

    /// Save state to file after each call
    #[arg(short = 'o', long, default_value_t = false, group = "saving")]
    override_state: bool,

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

impl Cli {
    fn override_settings(&mut self, other: Cli) {
        self.run_settings_path = other.run_settings_path;
        self.global_settings_path = other.global_settings_path;
        self.state_folder = other.state_folder;
        self.model_file = other.model_file;
        self.override_state = other.override_state;
        self.no_save_state = other.no_save_state;
    }

    fn run_command(
        &mut self,
        command: Commands,
        run_history: &mut RunHistory,
        state: &mut State,
    ) -> Result<ControlFlow<()>, Report> {
        match command {
            Commands::Quit {} => {
                return Ok(ControlFlow::Break(()));
            }
            Commands::Inspect(inspect) => inspect.run(state)?,
            Commands::Bench { samples } => {
                info!(
                    "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
                    format!(
                        "{}",
                        run_history.default_runtime_settings.hard_coded_integrand
                    )
                    .green(),
                    format!("{}", samples).blue()
                );
                let mut integrand = integrand_factory(&run_history.default_runtime_settings);
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
            Commands::Save(s) => match s {
                Save::Dot { path } => {
                    let exp_set = ExportSettings {
                        root_folder: path.clone().unwrap_or(self.state_folder.clone()),
                    };
                    state.process_list.export_dot(&exp_set, &state.model)?
                }
                Save::State {} => {
                    debug!("Saving State, overriding: {}", self.override_state);
                    state.save(&self.state_folder, self.override_state, false)?;
                    run_history.save_yaml(&self.state_folder, self.override_state, false)?;
                }
            },
            Commands::Set(s) => match s {
                Set::Log { level, format } => {
                    if let Some(level) = level {
                        *LOG_LEVEL.lock().unwrap() = (level).into();
                    }

                    *LOG_FORMAT.lock().unwrap() = format;
                }
                Set::BaseDir { path } => {
                    self.state_folder = path.clone();
                }
            },

            Commands::Generate(g) => g.run(
                state,
                &self.state_folder,
                self.override_state,
                &run_history.global_settings,
                &run_history.default_runtime_settings,
            )?,

            Commands::Integrate(g) => {
                g.run(state)?;
            }

            Commands::Display(l) => match l {
                Display::Integrands { process_id } => {
                    let process = &state.process_list.processes[process_id];
                    println!("Integrands for process {}:", process_id);
                    for integrand in process.get_integrand_names() {
                        println!("  {}", integrand);
                    }
                }
                Display::Processes => {
                    println!("Processes:");
                    for (index, process) in state.process_list.processes.iter().enumerate() {
                        println!("  {}", process.definition.folder_name(&state.model, index));
                    }
                }
                Display::Model => {
                    println!("{}", state.model.name)
                }
            },
            Commands::Run(r) => {
                if let Some(r) = &self.run_history {
                    let mut new_run_history = RunHistory::from_file_yaml(r)?;
                    let res = new_run_history.run(self, state);

                    run_history.merge(new_run_history);
                    return res;
                }
            }
            _ => {}
        }
        Ok(ControlFlow::Continue(()))
    }

    pub fn run(mut self) -> Result<()> {
        setup_log().unwrap();
        self.initialize();
        let mut run_history = self.get_run_history()?;

        let mut state = match State::load(&self.state_folder, self.model_file.clone()) {
            Ok(state) => state,
            Err(e) => {
                warn!(
                    "{e} state folder:{} model_file:{} loading default state",
                    self.state_folder.display(),
                    self.model_file
                        .as_ref()
                        .map(|p| p.display().to_string())
                        .unwrap_or("None".to_string())
                );
                State::default()
            }
        };

        // let mut state = State::load(&cli.state_file,);

        if let Some(a) = self.command.take() {
            let _ = self.run_command(a.clone(), &mut run_history, &mut state)?;

            run_history.push(a);
        } else {
            print_banner();
            let prompt = DefaultPrompt {
                left_prompt: DefaultPromptSegment::Basic(format!(
                    "{}:γloop",
                    self.state_folder.display()
                )),
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
                    ReadCommandOutput::Command(mut cli) => {
                        if let Some(c) = cli.command.take() {
                            match cli.run_command(c.clone(), &mut run_history, &mut state) {
                                Err(e) => {
                                    eprintln!("{e:?}");
                                    self.override_settings(cli);
                                }
                                Ok(ControlFlow::Break(())) => {
                                    run_history.push(c);
                                    self.override_settings(cli);
                                    break;
                                }
                                _ => {
                                    run_history.push(c);
                                    self.override_settings(cli);
                                }
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

        if !self.no_save_state {
            debug!("Saving State, override: {}", self.override_state);
            state.save(&self.state_folder, self.override_state, false)?;
            run_history.save_yaml(&self.state_folder, self.override_state, false)?;
        }
        Ok(())
    }

    pub fn initialize(&self) {
        let (panic, eyre) = HookBuilder::default()
            .capture_span_trace_by_default(cfg!(debug_assertions))
            .into_hooks();
        panic.install();
        eyre.install().unwrap();

        let _ = GS.delta_vec;
        crate::set_interrupt_handler();
        // activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");

        crate::initialize_reps();
    }

    fn get_run_history(&self) -> Result<RunHistory> {
        let default_path = self.state_folder.join("run_history.yaml");
        let path = self.run_history.as_ref().unwrap_or(&default_path);
        RunHistory::new(
            path,
            self.run_settings_path.as_deref(),
            self.global_settings_path.as_deref(),
        )
    }
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone)]
pub enum Import {
    Model {
        // #[arg(short = 'p')]
        path: PathBuf,
        // format: String,
    },
    Amplitude {
        // #[arg(short = 'p')]
        path: PathBuf,
        // format: String,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone)]
pub enum Save {
    Dot {
        path: Option<PathBuf>,
    },
    State {
        // #[arg(short = 'o', long, default_value_t = false)]
        // override_state: bool,
    },
}

#[repr(usize)]
#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, ValueEnum, Serialize, Deserialize,
)]
pub enum LogLevel {
    /// A level lower than all log levels.
    Off,
    /// Corresponds to the `Error` log level.
    Error,
    /// Corresponds to the `Warn` log level.
    Warn,
    /// Corresponds to the `Info` log level.
    Info,
    /// Corresponds to the `Debug` log level.
    Debug,
    /// Corresponds to the `Trace` log level.
    Trace,
}

impl From<LogLevel> for LevelFilter {
    fn from(value: LogLevel) -> Self {
        match value {
            LogLevel::Debug => LevelFilter::Debug,
            LogLevel::Trace => LevelFilter::Trace,
            LogLevel::Warn => LevelFilter::Warn,
            LogLevel::Error => LevelFilter::Error,
            LogLevel::Off => LevelFilter::Off,
            LogLevel::Info => LevelFilter::Info,
        }
    }
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone)]
pub enum Set {
    BaseDir {
        path: PathBuf,
    },
    Log {
        #[arg(short = 'l')]
        level: Option<LogLevel>,
        // #[clap(subcommand)]
        #[arg(short, long, value_enum, default_value_t = LogFormat::Long)]
        format: LogFormat,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone)]
pub enum Display {
    Model,
    Processes,
    Integrands { process_id: usize },
}

#[derive(Debug, Args, Serialize, Deserialize, Clone)]
pub struct Run {}

pub enum Log {
    Level(LevelFilter),
    Format(LogFormat),
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone)]
pub enum Commands {
    #[clap(subcommand)]
    Display(Display),
    #[clap(subcommand)]
    Set(Set),
    #[clap(subcommand)]
    Import(Import),
    #[clap(subcommand)]
    Save(Save),

    Run(Run),

    Integrate(Integrate),

    Generate(generate::Generate),

    /// Quit gammaloop
    Quit {
        // #[arg(short = 'o', long, default_value_t = false)]
        // override_state: bool,
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

#[repr(usize)]
#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, ValueEnum, Serialize, Deserialize,
)]
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

pub(crate) fn print_banner() {
    println!(
        "\n{}{}\n",
        r"                                        _
                                       | |
   __ _  __ _ _ __ ___  _ __ ___   __ _| |     ___   ___  _ __
  / _` |/ _` | '_ ` _ \| '_ ` _ \ / _` | |    / _ \ / _ \| '_ \
 | (_| | (_| | | | | | | | | | | | (_| | |___| (_) | (_) | |_) |
  \__, |\__,_|_| |_| |_|_| |_| |_|\__,_|______\___/ \___/| .__/
   __/ |                                                 | |
"
        .to_string()
        .bold()
        .blue(),
        format!(
            r#"  |___/    version:{} {}          |_|    "#,
            format!("{:<15}", GIT_VERSION).green(),
            format!("{:>10}", LOG_LEVEL.lock().unwrap().as_str()).green(),
        )
        .bold()
        .blue(),
    );
}

pub(crate) fn setup_log() -> Result<()> {
    fern::Dispatch::new()
        .filter(|metadata| metadata.level() <= (*LOG_LEVEL.lock().unwrap()))
        // Perform allocation-free log formatting
        .format(|out, message, record| {
            let now = Local::now();
            match *LOG_FORMAT.lock().unwrap() {
                LogFormat::Long => out.finish(format_args!(
                    "[{}] @{} {}: {}",
                    format!(
                        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:03}",
                        now.year(),
                        now.month(),
                        now.day(),
                        now.hour(),
                        now.minute(),
                        now.second(),
                        now.timestamp_subsec_millis()
                    )
                    .bright_green(),
                    format_target(record.target().into(), record.level()),
                    format_level(record.level()),
                    message
                )),
                LogFormat::Short => out.finish(format_args!(
                    "[{}] {}: {}",
                    format!("{:02}:{:02}:{:02}", now.hour(), now.minute(), now.second())
                        .bright_green(),
                    format_level(record.level()),
                    message
                )),
                LogFormat::Min => out.finish(format_args!(
                    "{}: {}",
                    format_level(record.level()),
                    message
                )),
                LogFormat::None => out.finish(format_args!("{}", message)),
            }
        })
        // Output to stdout, files, and other Dispatch configurations
        .chain(std::io::stdout())
        .chain(fern::log_file("gammaloop_rust_output.log")?)
        // Apply globally
        .apply()?;

    Ok(())
}
