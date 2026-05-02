use std::{ops::ControlFlow, path::PathBuf, str::FromStr};

use clap::Subcommand;
use color_eyre::Report;
use gammalooprs::settings::RuntimeSettings;
use save::SaveState;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    completion::CompletionArgExt,
    state::{CommandHistory, ProcessRef, RunHistory, State},
    CLISettings,
};
use symbolica::atom::Atom;
pub mod commands_block;
pub use commands_block::StartCommandsBlock;
pub mod display;
pub use display::Display;
pub mod duplicate;
pub use duplicate::Duplicate;
pub mod generate;
pub use generate::Generate;
pub mod import;
pub use import::Import;
pub mod inspect;
pub use inspect::Inspect;
pub mod integrate;
pub use integrate::{Integrate, IntegrationOutput};
pub mod remove;
pub use remove::Remove;
pub mod save;
pub use save::Save;
pub mod set;
pub mod shell;
pub use set::Set;
pub use shell::Shell;
pub mod threedreps;
pub use threedreps::ThreeDRep;
pub mod run;
pub use run::Run;
pub mod evaluate;
pub mod evaluate_samples;
pub use evaluate::Evaluate;
pub mod renormalize;
pub use renormalize::Renormalize;
pub mod profile;
pub use profile::Profile;
pub(crate) mod process_settings;

#[derive(Debug, Clone, Default)]
pub enum CommandOutput {
    #[default]
    None,
    Evaluate(Atom),
    Integrate(IntegrationOutput),
}

#[derive(Debug, Clone)]
pub struct CommandExecution {
    pub flow: ControlFlow<SaveState>,
    pub output: CommandOutput,
}

impl CommandExecution {
    pub fn continue_with(output: CommandOutput) -> Self {
        Self {
            flow: ControlFlow::Continue(()),
            output,
        }
    }

    pub fn continue_without_output() -> Self {
        Self::continue_with(CommandOutput::None)
    }

    pub fn break_with(save_state: SaveState) -> Self {
        Self {
            flow: ControlFlow::Break(save_state),
            output: CommandOutput::None,
        }
    }
}

#[allow(clippy::large_enum_variant)]
#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Commands {
    #[clap(subcommand)]
    Display(Display),
    #[clap(subcommand)]
    Duplicate(Duplicate),
    #[clap(subcommand)]
    Set(Set),
    #[clap(subcommand)]
    Import(Import),
    #[clap(subcommand)]
    Save(Save),

    Run(Run),
    #[command(name = "start_commands_block")]
    StartCommandsBlock(StartCommandsBlock),
    #[command(name = "finish_commands_block")]
    FinishCommandsBlock,

    #[command(name = "remove")]
    #[clap(subcommand)]
    Remove(Remove),

    Integrate(Integrate),

    Generate(Generate),

    /// Quit gammaloop
    Quit(SaveState),
    /// Inspect a single phase‑space point / momentum configuration
    // #[clap(subcommand)]
    Inspect(Inspect),

    Evaluate(Evaluate),

    Renormalize(Renormalize),

    /// Benchmark raw integrand evaluation speed
    Bench {
        /// Number of random samples to evaluate
        #[arg(short = 's', long, value_name = "SAMPLES")]
        samples: usize,
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: ProcessRef,

        /// The integrand name to benchmark
        #[arg(
            short = 'i',
            long = "integrand-name",
            value_name = "NAME",
            completion_integrand_selector(crate::completion::SelectorKind::Any)
        )]
        integrand_name: String,
        /// Number of cores to parallelize over
        #[arg(short = 'c', long)]
        n_cores: usize,
    },
    #[clap(subcommand)]
    Profile(Profile),

    #[command(name = "3Drep")]
    #[clap(subcommand)]
    ThreeDRep(ThreeDRep),

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

    #[doc(hidden)]
    #[command(name = "__command_template", hide = true)]
    CommandTemplate,
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
        global_cli_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<CommandExecution, Report> {
        match self {
            Commands::Profile(p) => {
                p.run(state, global_cli_settings)?;
            }
            Commands::ThreeDRep(command) => {
                command.run(state, global_cli_settings, default_runtime_settings)?;
            }
            Commands::Quit(s) => {
                return Ok(CommandExecution::break_with(s));
            }
            Commands::Inspect(inspect) => {
                let _ = inspect.run(state)?;
            }
            Commands::Bench {
                samples,
                process,
                integrand_name,
                n_cores,
            } => {
                let process_id = process.resolve(&state.process_list)?;
                state.bench(samples, process_id, integrand_name, n_cores)?;
            }
            Commands::Import(s) => s.run(state, global_cli_settings)?,
            Commands::Save(s) => s.run(
                state,
                run_history,
                default_runtime_settings,
                global_cli_settings,
            )?,
            Commands::Set(s) => s.run(state, global_cli_settings, default_runtime_settings)?,
            Commands::Generate(g) => {
                let would_compile_into_active_state = global_cli_settings.session.read_only_state
                    && global_cli_settings.global.generation.evaluator.compile
                    && matches!(
                        g.mode.as_ref(),
                        None | Some(generate::GenerateCmd::Existing(_))
                    );
                if would_compile_into_active_state {
                    return Err(Report::msg(format!(
                        "Cannot compile generated integrands into '{}' because this session was started with --read-only-state. Disable `global.generation.evaluator.compile`, restart without --read-only-state, or save the state elsewhere first.",
                        global_cli_settings.state.folder.display()
                    )));
                }

                g.run(
                    state,
                    &global_cli_settings.state.folder,
                    global_cli_settings.override_state,
                    &global_cli_settings.global,
                    default_runtime_settings,
                )?
            }
            Commands::Integrate(g) => {
                return Ok(CommandExecution::continue_with(CommandOutput::Integrate(
                    g.run(state, global_cli_settings)?,
                )));
            }
            Commands::Evaluate(g) => {
                return Ok(CommandExecution::continue_with(CommandOutput::Evaluate(
                    g.run(state, global_cli_settings, default_runtime_settings)?,
                )));
            }
            Commands::Renormalize(r) => {
                _ = r.run(state, global_cli_settings)?;
            }
            Commands::Display(l) => {
                l.run(
                    state,
                    global_cli_settings,
                    default_runtime_settings,
                    run_history,
                )?;
            }
            Commands::Duplicate(command) => {
                command.run(state)?;
            }
            Commands::Run(r) => {
                return r.run(
                    state,
                    global_cli_settings,
                    default_runtime_settings,
                    run_history,
                );
            }
            Commands::StartCommandsBlock(_) | Commands::FinishCommandsBlock => {
                return Err(Report::msg(
                    "Command block definition commands must be handled by the CLI session",
                ));
            }
            Commands::Remove(r) => r.run(state)?,
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
            Commands::CommandTemplate => {
                return Err(Report::msg(
                    "Command templates can only be executed through `run` with `-D/--define` variables",
                ));
            }
        }
        Ok(CommandExecution::continue_without_output())
    }
}

#[cfg(test)]
mod tests {
    use super::Commands;
    use crate::{
        commands::generate::{Generate, GenerateCmd, ProcessArgs},
        state::{RunHistory, State},
        CLISettings,
    };
    use gammalooprs::settings::RuntimeSettings;

    #[test]
    fn generate_rejects_compilation_into_active_state_in_read_only_mode() {
        let mut state = State::new_test();
        let mut run_history = RunHistory::default();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();
        cli_settings.session.read_only_state = true;
        cli_settings.global.generation.evaluator.compile = true;

        let err = Commands::Generate(Generate {
            keep_sources: false,
            mode: Some(GenerateCmd::Existing(ProcessArgs {
                process: None,
                integrand_name: None,
            })),
        })
        .run(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut runtime_settings,
        )
        .unwrap_err();

        assert!(format!("{err:?}").contains("--read-only-state"));
    }
}
