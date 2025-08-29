use crate::{
    feyngen::GenerationType,
    graph::Graph,
    initialisation::initialise,
    integrands::{integrand_factory, HasIntegrand},
    model::Model,
    processes::{ExportSettings, Process, ProcessDefinition},
    utils::{F, GIT_VERSION, GS},
};
use bincode_trait_derive::{Decode, Encode};
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
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use state::{current_log_spec, RunHistory, State};
use std::ops::ControlFlow;
use std::path::PathBuf;
use std::time::Instant;
use symbolica::numerical_integration::Sample;
use tracing::{init_tracing, LogLevel};

pub mod generate;
pub mod inspect;
pub mod integrate;
pub mod settings;
pub mod state;
pub mod tracing;

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
pub struct Cli {
    /// Path to the a run file as history, and as settings, by default is: `./gammaloop_state/run.yaml`
    #[arg(short = 'r', long)]
    run_history: Option<PathBuf>,

    /// Path to the state folder
    #[arg(short = 's', long, default_value = "./gammaloop_state")]
    pub state_folder: PathBuf,

    /// Path to the model file
    #[arg(short = 'm', long)]
    pub model_file: Option<PathBuf>,

    /// Skip saving state on exit
    #[arg(short = 'n', long, default_value_t = false, group = "saving")]
    no_save_state: bool,

    /// Save state to file after each call
    #[arg(short = 'o', long, default_value_t = false, group = "saving")]
    override_state: bool,

    // /// Debug level
    // #[arg(short = 'd', long, value_enum, default_value_t = LogLevel::Info)]
    // debug_level: LogLevel,
    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

impl Cli {
    pub fn new_test(state_folder: PathBuf) -> Self {
        Cli {
            run_history: None,
            state_folder,
            model_file: None,
            no_save_state: true,
            override_state: false,
            command: None,
        }
    }

    fn override_settings(&mut self, other: Cli) {
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
            Commands::Inspect(inspect) => {
                let _ = inspect.run(state)?;
            }
            Commands::Bench {
                samples,
                process_id,
                process_name,
                n_cores,
            } => {
                state.bench(samples, process_id, process_name, n_cores)?;
            }
            Commands::Import(s) => match s {
                Import::Amplitude { path } => {
                    state.import_amplitude(path, None)?;
                }
                Import::Model { path } => {
                    state.model = Model::from_file(path)?;
                }
            },
            Commands::Save(s) => match s {
                Save::Dot { path } => {
                    state.export_dots(path.unwrap_or(self.state_folder.clone()))?;
                }
                Save::State {} => {
                    state.save(&self.state_folder, self.override_state, false)?;
                    run_history.save_yaml(&self.state_folder, self.override_state, false)?;
                }
            },
            Commands::Set(s) => match s {
                Set::Log { level, .. } => {
                    if let Some(level) = level {
                        let spec = level.to_env_spec();
                        state.set_log_spec(spec)?;
                    }
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
            Commands::Run(Run { path }) => {
                let mut new_run_history = RunHistory::new(path)?;
                let res = new_run_history.run(self, state);

                run_history.merge(new_run_history);
                return res;
            }
            _ => {}
        }
        Ok(ControlFlow::Continue(()))
    }

    pub fn run(mut self) -> Result<()> {
        initialise();

        let mut state = match State::load(self.state_folder.clone(), self.model_file.clone()) {
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
                State::new(self.state_folder.clone())
            }
        };
        let mut run_history = self.get_run_history()?;

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

    // pub fn initialize(&self) {}

    fn get_run_history(&self) -> Result<RunHistory> {
        let default_path = self.state_folder.join("run.yaml");
        let path = self.run_history.as_ref().unwrap_or(&default_path);
        RunHistory::new(path)
    }
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema)]
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

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema)]
pub enum Save {
    Dot {
        path: Option<PathBuf>,
    },
    State {
        // #[arg(short = 'o', long, default_value_t = false)]
        // override_state: bool,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema)]
pub enum Set {
    BaseDir {
        path: PathBuf,
    },
    Log {
        #[arg(short = 'l')]
        level: Option<LogLevel>,
        // // #[clap(subcommand)]
        // #[arg(short, long, value_enum, default_value_t = LogFormat::Long)]
        // format: LogFormat,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema)]
pub enum Display {
    Model,
    Processes,
    Integrands { process_id: usize },
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema)]
pub struct Run {
    /// Path to a run file to execute
    path: PathBuf,
}

pub enum Log {
    Level(LevelFilter),
    Format(LogFormat),
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema)]
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
        /// The process id to inspect
        #[arg(short = 'i', long = "process-id", value_name = "ID")]
        process_id: usize,

        /// The name of the process to inspect
        #[arg(short = 'n', long = "name", value_name = "NAME")]
        process_name: String,
        /// Number of cores to parallelize over
        #[arg(short = 'c', long)]
        n_cores: usize,
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
    let spec = current_log_spec();
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
            r#"  |___/    version:{:<15} {:>10}          |_|    "#,
            GIT_VERSION.green(),
            spec.green(),
        )
        .bold()
        .blue(),
    );
}
