use std::ops::ControlFlow;

use clap::Args;
use gammalooprs::settings::RuntimeSettings;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use color_eyre::{Result, eyre::eyre};

use crate::{
    CLISettings,
    command_parser::{split_command_line, split_command_list},
    state::{CommandHistory, RunHistory, State},
};

use super::{Commands, save::SaveState};

pub const MAX_RUN_DEPTH: usize = 100;

#[derive(Debug, Clone)]
pub enum PreparedCommand {
    Plain(CommandHistory),
    Run {
        command: CommandHistory,
        plan: PreparedRun,
    },
}

#[derive(Debug, Clone)]
pub struct PreparedCommandsBlock {
    pub name: String,
    pub commands: Vec<PreparedCommand>,
}

#[derive(Debug, Clone, Default)]
pub struct PreparedRun {
    pub blocks: Vec<PreparedCommandsBlock>,
    pub commands: Vec<PreparedCommand>,
}

impl PreparedRun {
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty() && self.commands.is_empty()
    }
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Run {
    /// Command block names to execute (in order)
    #[arg(value_name = "BLOCK_NAME")]
    pub(crate) block_names: Vec<String>,

    #[arg(short = 'c', long)]
    pub(crate) commands: Option<String>,
}

impl Run {
    pub fn selected_block_names(&self) -> &[String] {
        self.block_names.as_slice()
    }

    pub fn is_noop(&self) -> bool {
        self.block_names.is_empty()
            && self
                .commands
                .as_ref()
                .map(|commands| commands.trim().is_empty())
                .unwrap_or(true)
    }

    pub fn canonical_raw_string(&self) -> String {
        let mut parts = vec!["run".to_string()];
        parts.extend(self.block_names.iter().map(|name| shell_quote(name)));
        if let Some(commands) = self.commands.as_ref() {
            if !commands.trim().is_empty() {
                parts.push("-c".to_string());
                parts.push(shell_quote(commands));
            }
        }
        parts.join(" ")
    }

    pub fn parse_inline_commands(&self) -> Result<Vec<CommandHistory>> {
        let Some(commands) = self.commands.as_deref() else {
            return Ok(Vec::new());
        };

        split_command_list(commands)
            .map_err(|_| {
                eyre!("Could not parse run -c command list: unmatched quotes or trailing escape")
            })?
            .into_iter()
            .enumerate()
            .map(|(index, command)| {
                CommandHistory::from_raw_string(&command).map_err(|err| {
                    eyre!(
                        "Failed to parse run -c command #{} '{}': {}",
                        index + 1,
                        command,
                        err
                    )
                })
            })
            .collect()
    }

    pub fn prepare(&self, run_history: &RunHistory, depth: usize) -> Result<PreparedRun> {
        if depth > MAX_RUN_DEPTH {
            return Err(eyre!(
                "Maximum nested run depth of {} reached while preparing {}",
                MAX_RUN_DEPTH,
                self.canonical_raw_string()
            ));
        }

        let selected_blocks = run_history.select_command_blocks(self.selected_block_names())?;
        let commands = prepare_command_histories_with_context(
            &self.parse_inline_commands()?,
            run_history,
            depth + 1,
            "run -c",
        )?;
        let blocks = selected_blocks
            .into_iter()
            .map(|block| {
                let block_context = format!("command block '{}'", block.name);
                Ok(PreparedCommandsBlock {
                    name: block.name,
                    commands: prepare_command_histories_with_context(
                        &block.commands,
                        run_history,
                        depth + 1,
                        &block_context,
                    )?,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(PreparedRun { blocks, commands })
    }

    pub fn run(
        &self,
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
        run_history: &mut RunHistory,
    ) -> Result<ControlFlow<SaveState>> {
        let mut session_state = crate::session::CliSessionState::default();
        let mut session = crate::session::CliSession::new(
            state,
            run_history,
            global_settings,
            default_runtime_settings,
            &mut session_state,
        );
        session.execute_top_level(CommandHistory::new_with_raw(
            Commands::Run(self.clone()),
            self.canonical_raw_string(),
        ))
    }
}

pub fn prepare_command_histories(
    commands: &[CommandHistory],
    run_history: &RunHistory,
    depth: usize,
) -> Result<Vec<PreparedCommand>> {
    prepare_command_histories_with_context(commands, run_history, depth, "commands")
}

pub fn prepare_command_histories_with_context(
    commands: &[CommandHistory],
    run_history: &RunHistory,
    depth: usize,
    context: &str,
) -> Result<Vec<PreparedCommand>> {
    commands
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, command)| {
            PreparedCommand::prepare(command, run_history, depth).map_err(|err| {
                eyre!(
                    "Failed to validate {} command #{}: {}",
                    context,
                    index + 1,
                    err
                )
            })
        })
        .collect()
}

impl PreparedCommand {
    pub fn prepare(
        command: CommandHistory,
        run_history: &RunHistory,
        depth: usize,
    ) -> Result<Self> {
        match command.command.clone() {
            Commands::Run(run) => Ok(Self::Run {
                command,
                plan: run.prepare(run_history, depth)?,
            }),
            _ => Ok(Self::Plain(command)),
        }
    }
}

fn shell_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }

    if split_command_line(value)
        .map(|split| split.len() == 1 && split[0] == value)
        .unwrap_or(false)
        && !value.contains(';')
    {
        return value.to_string();
    }

    format!("'{}'", value.replace('\'', "'\"'\"'"))
}
