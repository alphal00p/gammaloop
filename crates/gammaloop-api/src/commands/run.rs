use std::str::FromStr;

use clap::Args;
use gammalooprs::settings::RuntimeSettings;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use color_eyre::{eyre::eyre, Result};

use crate::{
    command_parser::{normalize_clap_args, split_command_line, split_command_list},
    command_template::{
        contains_placeholder, expand, is_valid_placeholder_name, CommandEnvironment,
    },
    state::{CommandHistory, RunHistory, State},
    CLISettings,
};

use super::{CommandExecution, Commands};

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

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct RunVariable {
    pub key: String,
    pub value: String,
}

impl RunVariable {
    fn assignment(&self) -> String {
        format!("{}={}", self.key, self.value)
    }
}

impl FromStr for RunVariable {
    type Err = String;

    fn from_str(raw: &str) -> std::result::Result<Self, Self::Err> {
        let (key, value) = raw
            .split_once('=')
            .ok_or_else(|| "expected KEY=VALUE".to_string())?;
        if !is_valid_placeholder_name(key) {
            return Err(format!(
                "invalid variable name '{key}'; expected [A-Za-z_][A-Za-z0-9_]*"
            ));
        }
        Ok(Self {
            key: key.to_string(),
            value: value.to_string(),
        })
    }
}

fn parse_run_variable(raw: &str) -> std::result::Result<RunVariable, String> {
    raw.parse()
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Run {
    /// Command block names to execute (in order)
    #[arg(value_name = "BLOCK_NAME")]
    pub(crate) block_names: Vec<String>,

    /// Define a string variable available to selected command blocks; repeat as needed
    #[arg(
        short = 'D',
        long = "define",
        value_name = "KEY=VALUE",
        action = clap::ArgAction::Append,
        value_parser = parse_run_variable
    )]
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub(crate) defines: Vec<RunVariable>,

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
        for variable in &self.defines {
            parts.push("-D".to_string());
            parts.push(shell_quote(&variable.assignment()));
        }
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
        let mut active_blocks = Vec::new();
        self.prepare_with_context(
            run_history,
            depth,
            &CommandEnvironment::new(),
            &mut active_blocks,
        )
    }

    fn prepare_with_context(
        &self,
        run_history: &RunHistory,
        depth: usize,
        inherited_environment: &CommandEnvironment,
        active_blocks: &mut Vec<String>,
    ) -> Result<PreparedRun> {
        if depth > MAX_RUN_DEPTH {
            return Err(eyre!(
                "Maximum nested run depth of {} reached while preparing {}",
                MAX_RUN_DEPTH,
                self.canonical_raw_string()
            ));
        }

        let environment = self.environment(inherited_environment)?;
        let selected_block_names = self
            .selected_block_names()
            .iter()
            .map(|name| expand(name, &environment))
            .collect::<Result<Vec<_>>>()?;
        let selected_blocks = run_history.select_command_blocks(&selected_block_names)?;
        let commands =
            self.prepare_inline_commands(run_history, depth + 1, &environment, active_blocks)?;
        let mut blocks = Vec::with_capacity(selected_blocks.len());
        for block in selected_blocks {
            if let Some(first_recursive_index) =
                active_blocks.iter().position(|name| name == &block.name)
            {
                let mut cycle = active_blocks[first_recursive_index..].to_vec();
                cycle.push(block.name.clone());
                return Err(eyre!(
                    "Maximum nested run depth of {} reached while preparing {}. Command block recursion detected: {}",
                    MAX_RUN_DEPTH,
                    self.canonical_raw_string(),
                    cycle.join(" -> ")
                ));
            }

            active_blocks.push(block.name.clone());
            let block_context = format!("command block '{}'", block.name);
            let block_commands = prepare_command_histories_with_environment(
                &block.commands,
                run_history,
                depth + 1,
                &block_context,
                &environment,
                active_blocks,
            );
            active_blocks.pop();

            blocks.push(PreparedCommandsBlock {
                name: block.name,
                commands: block_commands?,
            });
        }

        Ok(PreparedRun { blocks, commands })
    }

    fn environment(&self, inherited: &CommandEnvironment) -> Result<CommandEnvironment> {
        let mut environment = inherited.clone();
        let mut defined_here = std::collections::BTreeSet::new();

        for variable in &self.defines {
            if !defined_here.insert(variable.key.clone()) {
                return Err(eyre!(
                    "Variable '{}' is defined more than once in {}",
                    variable.key,
                    self.canonical_raw_string()
                ));
            }
            let value = expand(&variable.value, &environment)?;
            environment.insert(variable.key.clone(), value);
        }

        Ok(environment)
    }

    fn prepare_inline_commands(
        &self,
        run_history: &RunHistory,
        depth: usize,
        environment: &CommandEnvironment,
        active_blocks: &mut Vec<String>,
    ) -> Result<Vec<PreparedCommand>> {
        let Some(commands) = self.commands.as_deref() else {
            return Ok(Vec::new());
        };

        let mut prepared_commands = Vec::new();
        for (index, command) in split_command_list(commands)
            .map_err(|_| {
                eyre!("Could not parse run -c command list: unmatched quotes or trailing escape")
            })?
            .into_iter()
            .enumerate()
        {
            let prepared = if contains_placeholder(&command) {
                prepare_raw_command_with_environment(
                    &command,
                    run_history,
                    depth,
                    &format!("run -c command #{}", index + 1),
                    environment,
                    active_blocks,
                )
                .map_err(|err| {
                    eyre!(
                        "Failed to validate run -c command #{} '{}': {}",
                        index + 1,
                        command,
                        err
                    )
                })?
            } else {
                let command_history = CommandHistory::from_raw_string(&command).map_err(|err| {
                    eyre!(
                        "Failed to parse run -c command #{} '{}': {}",
                        index + 1,
                        command,
                        err
                    )
                })?;
                PreparedCommand::prepare_with_environment(
                    command_history,
                    run_history,
                    depth,
                    environment,
                    active_blocks,
                )?
            };
            prepared_commands.push(prepared);
        }

        Ok(prepared_commands)
    }

    pub fn run(
        &self,
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
        run_history: &mut RunHistory,
    ) -> Result<CommandExecution> {
        let mut session_state = crate::session::CliSessionState::default();
        let mut session = crate::session::CliSession::new(
            state,
            run_history,
            global_settings,
            default_runtime_settings,
            &mut session_state,
        );
        session.execute_command(CommandHistory::new_with_raw(
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

fn prepare_command_histories_with_environment(
    commands: &[CommandHistory],
    run_history: &RunHistory,
    depth: usize,
    context: &str,
    environment: &CommandEnvironment,
    active_blocks: &mut Vec<String>,
) -> Result<Vec<PreparedCommand>> {
    let mut prepared_commands = Vec::with_capacity(commands.len());
    for (index, command) in commands.iter().cloned().enumerate() {
        let prepared = PreparedCommand::prepare_with_environment(
            command,
            run_history,
            depth,
            environment,
            active_blocks,
        )
        .map_err(|err| {
            eyre!(
                "Failed to validate {} command #{}: {}",
                context,
                index + 1,
                err
            )
        })?;
        prepared_commands.push(prepared);
    }
    Ok(prepared_commands)
}

fn prepare_raw_command_with_environment(
    raw: &str,
    run_history: &RunHistory,
    depth: usize,
    context: &str,
    environment: &CommandEnvironment,
    active_blocks: &mut Vec<String>,
) -> Result<PreparedCommand> {
    let command = command_history_from_template(raw, environment)
        .map_err(|err| eyre!("Failed to expand {context}: {err}"))?;
    PreparedCommand::prepare_with_environment(
        command,
        run_history,
        depth,
        environment,
        active_blocks,
    )
}

fn command_history_from_template(
    raw: &str,
    environment: &CommandEnvironment,
) -> Result<CommandHistory> {
    let args = split_command_line(raw)
        .map_err(|_| {
            eyre!("Could not parse command template: unmatched quotes or trailing escape")
        })?
        .into_iter()
        .map(|arg| expand(&arg, environment))
        .collect::<Result<Vec<_>>>()?;
    let raw_string = args
        .iter()
        .map(|arg| shell_quote(arg))
        .collect::<Vec<_>>()
        .join(" ");

    CommandHistory::from_args_and_raw(normalize_clap_args(args), raw_string)
        .map_err(|err| eyre!("Expanded command did not parse: {err}"))
}

fn expand_command_history(
    command: CommandHistory,
    environment: &CommandEnvironment,
) -> Result<CommandHistory> {
    let Some(raw) = command.raw_string() else {
        return Ok(command);
    };

    if command.is_template() || contains_placeholder(raw) {
        command_history_from_template(raw, environment)
    } else {
        Ok(command)
    }
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
            Commands::CommandTemplate => Ok(Self::Plain(command)),
            _ => Ok(Self::Plain(command)),
        }
    }

    fn prepare_with_environment(
        command: CommandHistory,
        run_history: &RunHistory,
        depth: usize,
        environment: &CommandEnvironment,
        active_blocks: &mut Vec<String>,
    ) -> Result<Self> {
        if let Commands::Run(run) = command.command.clone() {
            return Ok(Self::Run {
                command,
                plan: run.prepare_with_context(run_history, depth, environment, active_blocks)?,
            });
        }

        let command = expand_command_history(command, environment)?;
        match command.command.clone() {
            Commands::Run(run) => Ok(Self::Run {
                command,
                plan: run.prepare_with_context(run_history, depth, environment, active_blocks)?,
            }),
            Commands::CommandTemplate => Err(eyre!(
                "Command template '{}' could not be expanded",
                command.raw_string().unwrap_or("<unknown>")
            )),
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
