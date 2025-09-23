use std::{ops::ControlFlow, path::PathBuf, str::FromStr};

use clap::Subcommand;
use color_eyre::Report;
use gammalooprs::settings::RuntimeSettings;
use save::SaveState;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    state::{CommandHistory, RunHistory, State},
    GlobalCliSettings,
};
pub mod display;
pub use display::Display;
pub mod generate;
pub use generate::Generate;
pub mod import;
pub use import::Import;
pub mod inspect;
pub use inspect::Inspect;
pub mod integrate;
pub use integrate::Integrate;
pub mod reset;
pub use reset::Reset;
pub mod save;
pub use save::Save;
pub mod set;
pub mod shell;
pub use set::Set;
pub use shell::Shell;
pub mod run;
pub use run::Run;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
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

    #[clap(subcommand)]
    Reset(Reset),

    Integrate(Integrate),

    Generate(Generate),

    /// Quit gammaloop
    Quit(SaveState),
    /// Inspect a single phase‑space point / momentum configuration
    // #[clap(subcommand)]
    Inspect(Inspect),

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
        #[arg(value_name = "PROCESS_FILE", value_hint = clap::ValueHint::FilePath)]
        process_file: PathBuf,
        #[arg(value_name = "BATCH_INPUT_FILE", value_hint = clap::ValueHint::FilePath)]
        batch_input_file: PathBuf,
        #[arg(short = 'n', long, value_name = "NAME")]
        name: String,
        #[arg(value_name = "NAME")]
        output_name: String,
    },
    #[command(name = "!")]
    Shell(Shell),
}

impl FromStr for Commands {
    type Err = Report;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        Ok(CommandHistory::from_raw_string(s)?.command)
    }
}

impl Commands {
    pub fn run(
        self,
        state: &mut State,
        run_history: &mut RunHistory,
        global_settings: &mut GlobalCliSettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<SaveState>, Report> {
        if global_settings.dummy {
            if let Commands::Quit(s) = self {
                return Ok(ControlFlow::Break(s));
            }
            println!("Dummy mode: Command '{:?}' not executed.", self);
        } else {
            match self {
                Commands::Quit(s) => {
                    return Ok(ControlFlow::Break(s));
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
                Commands::Import(s) => s.run(state)?,
                Commands::Save(s) => s.run(
                    state,
                    run_history,
                    default_runtime_settings,
                    global_settings,
                )?,
                Commands::Set(s) => s.run(state, global_settings, default_runtime_settings)?,

                Commands::Generate(g) => g.run(
                    state,
                    &global_settings.state_folder,
                    global_settings.override_state,
                    &global_settings.global_settings,
                    &default_runtime_settings,
                )?,

                Commands::Integrate(g) => {
                    g.run(state)?;
                }

                Commands::Display(l) => l.run(state)?,
                Commands::Run(r) => {
                    return r.run(
                        state,
                        global_settings,
                        default_runtime_settings,
                        run_history,
                    )
                }
                Commands::Reset(r) => r.run(state)?,
                Commands::Batch {
                    process_file: _process_file,
                    batch_input_file: _batch_input_file,
                    name: _name,
                    output_name: _output_name,
                } => {
                    todo!("Batch command not implemented yet");
                }
                Commands::Shell(s) => {
                    s.run()?;
                }
            }
        }
        Ok(ControlFlow::Continue(()))
    }
}
