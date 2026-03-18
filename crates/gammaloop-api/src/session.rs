use std::{
    io::{self, IsTerminal, Write},
    ops::ControlFlow,
};

use color_eyre::{eyre::eyre, Result};
use colored::Colorize;
use gammalooprs::integrands::process::ProcessIntegrand;
use gammalooprs::processes::ProcessCollection;
use gammalooprs::settings::RuntimeSettings;
use tracing::{info, warn};

use crate::{
    commands::{
        process_settings::{serialize_runtime_named_settings, ProcessSettingsCompletionEntry},
        run::{prepare_command_histories_with_context, PreparedCommand, PreparedRun},
        save::SaveState,
        Commands, StartCommandsBlock,
    },
    repl::{IrProfileCompletionEntry, ProcessCompletionEntry, ProcessKind},
    state::{CommandHistory, CommandsBlock, RunHistory, State},
    CLISettings,
};

fn format_command_block_conflict_message(conflicting_blocks: &[String]) -> String {
    match conflicting_blocks {
        [] => String::new(),
        [single] => format!(
            "Run card command block '{}' redefines an existing block with different commands",
            single
        ),
        many => format!(
            "Run card command blocks {} redefine existing blocks with different commands",
            many.iter()
                .map(|name| format!("'{}'", name))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

fn prompt_command_block_conflict_override(conflicting_blocks: &[String]) -> Result<bool> {
    let prompt = format!(
        "{}. Proceed anyway? (it may break reproducibility of the state from run.toml) [y/n] > ",
        format_command_block_conflict_message(conflicting_blocks)
    );
    let mut stderr = io::stderr();
    let mut input = String::new();
    loop {
        eprint!("{prompt}");
        stderr.flush()?;
        input.clear();
        io::stdin().read_line(&mut input)?;
        match input.trim().to_ascii_lowercase().as_str() {
            "y" | "yes" => return Ok(true),
            "n" | "no" => return Ok(false),
            _ => {
                eprintln!("Please answer 'y' or 'n'.");
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct CliSessionState {
    pending_commands_block: Option<PendingCommandsBlock>,
}

#[derive(Debug, Clone)]
struct PendingCommandsBlock {
    name: String,
    commands: Vec<CommandHistory>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HistoryMode {
    Record,
    Suppress,
}

pub struct CliSession<'a> {
    state: &'a mut State,
    run_history: &'a mut RunHistory,
    cli_settings: &'a mut CLISettings,
    default_runtime_settings: &'a mut RuntimeSettings,
    session_state: &'a mut CliSessionState,
}

impl<'a> CliSession<'a> {
    pub fn new(
        state: &'a mut State,
        run_history: &'a mut RunHistory,
        cli_settings: &'a mut CLISettings,
        default_runtime_settings: &'a mut RuntimeSettings,
        session_state: &'a mut CliSessionState,
    ) -> Self {
        Self {
            state,
            run_history,
            cli_settings,
            default_runtime_settings,
            session_state,
        }
    }

    pub fn execute_top_level(&mut self, command: CommandHistory) -> Result<ControlFlow<SaveState>> {
        let prepared = PreparedCommand::prepare(command, self.run_history, 1)?;
        self.execute_prepared(prepared, HistoryMode::Record)
    }

    pub fn replay_run_history(&mut self) -> Result<ControlFlow<SaveState>> {
        self.run_history
            .apply_session_settings(self.cli_settings, self.default_runtime_settings)?;
        let prepared = prepare_command_histories_with_context(
            &self.run_history.commands.clone(),
            self.run_history,
            1,
            "run history",
        )?;
        self.execute_prepared_commands(prepared, HistoryMode::Suppress)
    }

    pub fn apply_boot_run_history(
        &mut self,
        boot_run_history: &RunHistory,
        effective_boot_run_history: &RunHistory,
        booted_existing_state: bool,
    ) -> Result<ControlFlow<SaveState>> {
        let interactive_prompt_available = io::stdin().is_terminal() && io::stdout().is_terminal();
        self.apply_boot_run_history_with_conflict_resolver(
            boot_run_history,
            effective_boot_run_history,
            booted_existing_state,
            |conflicting_blocks| {
                if interactive_prompt_available {
                    return prompt_command_block_conflict_override(conflicting_blocks);
                }

                Err(eyre!(
                    "{}. Cannot continue non-interactively; rerun in an interactive terminal or align the command block definitions.",
                    format_command_block_conflict_message(conflicting_blocks)
                ))
            },
        )
    }

    pub fn apply_boot_run_history_with_conflict_resolver<F>(
        &mut self,
        boot_run_history: &RunHistory,
        effective_boot_run_history: &RunHistory,
        booted_existing_state: bool,
        mut resolve_conflicts: F,
    ) -> Result<ControlFlow<SaveState>>
    where
        F: FnMut(&[String]) -> Result<bool>,
    {
        boot_run_history.validate()?;

        if booted_existing_state
            && !self
                .run_history
                .frozen_boot_settings_match(boot_run_history)
        {
            let warning = format!(
                "Boot card settings differ from the frozen settings stored in {}. This session is forced into --read-only-state. Line up cli_settings.global and default_runtime_settings in the boot card with the frozen settings in run.toml if you want to boot this state without read-only mode.",
                self.cli_settings.state.folder.join("run.toml").display()
            );
            self.cli_settings.session.read_only_state = true;
            self.cli_settings
                .session
                .startup_warnings
                .push(warning.clone());
            warn!("{warning}");
        }
        if !booted_existing_state {
            self.run_history.freeze_boot_settings_from(boot_run_history);
        }

        let mut merged_history = self.run_history.clone();
        let conflicting_blocks = merged_history
            .conflicting_command_block_names(&effective_boot_run_history.command_blocks);
        let overwrite_conflicting_blocks = if conflicting_blocks.is_empty() {
            false
        } else if self.cli_settings.session.read_only_state {
            return Err(eyre!(
                "{}. Cannot proceed while --read-only-state is enabled.",
                format_command_block_conflict_message(&conflicting_blocks)
            ));
        } else if resolve_conflicts(&conflicting_blocks)? {
            true
        } else {
            info!(
                "{}",
                "Boot run card application cancelled by user.".yellow()
            );
            return Ok(ControlFlow::Break(SaveState::default()));
        };
        merged_history.merge_command_blocks_with_overwrite(
            &effective_boot_run_history.command_blocks,
            overwrite_conflicting_blocks,
        )?;

        for block in &effective_boot_run_history.command_blocks {
            let block_context = format!("command block '{}'", block.name);
            let _ = prepare_command_histories_with_context(
                &block.commands,
                &merged_history,
                2,
                &block_context,
            )?;
        }

        let prepared = prepare_command_histories_with_context(
            &effective_boot_run_history.commands,
            &merged_history,
            1,
            "boot",
        )?;

        self.run_history.merge_command_blocks_with_overwrite(
            &effective_boot_run_history.command_blocks,
            overwrite_conflicting_blocks,
        )?;
        effective_boot_run_history
            .apply_session_settings(self.cli_settings, self.default_runtime_settings)?;
        self.execute_prepared_commands(prepared, HistoryMode::Record)
    }

    pub fn has_pending_commands_block(&self) -> bool {
        self.session_state.pending_commands_block.is_some()
    }

    pub fn pending_commands_block_name(&self) -> Option<String> {
        self.session_state
            .pending_commands_block
            .as_ref()
            .map(|pending| pending.name.clone())
    }

    pub fn prompt_state_label(&self) -> String {
        self.cli_settings.state.prompt_label()
    }

    pub fn current_commands_block_names(&self) -> Vec<String> {
        self.run_history
            .command_blocks
            .iter()
            .map(|block| block.name.clone())
            .collect()
    }

    pub fn current_process_entries(&self) -> Vec<ProcessCompletionEntry> {
        self.state
            .process_list
            .processes
            .iter()
            .enumerate()
            .map(|(id, process)| ProcessCompletionEntry {
                id,
                name: process.definition.folder_name.clone(),
                kind: match &process.collection {
                    ProcessCollection::Amplitudes(_) => ProcessKind::Amplitude,
                    ProcessCollection::CrossSections(_) => ProcessKind::CrossSection,
                },
                integrand_names: process
                    .collection
                    .get_integrand_names()
                    .into_iter()
                    .map(str::to_string)
                    .collect(),
            })
            .collect()
    }

    pub fn current_ir_profile_entries(&self) -> Vec<IrProfileCompletionEntry> {
        self.state
            .process_list
            .processes
            .iter()
            .filter_map(|process| {
                let ProcessCollection::CrossSections(cross_sections) = &process.collection else {
                    return None;
                };

                Some(
                    cross_sections
                        .iter()
                        .filter_map(|(integrand_name, cross_section)| {
                            let ProcessIntegrand::CrossSection(integrand) =
                                cross_section.integrand.as_ref()?
                            else {
                                return None;
                            };
                            let completion_entries = integrand.ir_profile_completion_entries();
                            let mut graph_names = completion_entries
                                .iter()
                                .map(|(graph_name, _)| graph_name.clone())
                                .collect::<Vec<_>>();
                            graph_names.sort();
                            graph_names.dedup();

                            let mut graph_limit_entries = completion_entries
                                .into_iter()
                                .flat_map(|(graph_name, limits)| {
                                    limits
                                        .into_iter()
                                        .map(move |limit| format!("{graph_name} {limit}"))
                                })
                                .collect::<Vec<_>>();
                            graph_limit_entries.sort();
                            graph_limit_entries.dedup();

                            Some(IrProfileCompletionEntry {
                                process_name: process.definition.folder_name.clone(),
                                integrand_name: integrand_name.clone(),
                                graph_names,
                                graph_limit_entries,
                            })
                        })
                        .collect::<Vec<_>>(),
                )
            })
            .flatten()
            .collect()
    }

    pub(crate) fn current_process_settings_entries(&self) -> Vec<ProcessSettingsCompletionEntry> {
        let mut entries = Vec::new();

        for (process_id, process) in self.state.process_list.processes.iter().enumerate() {
            for integrand_name in process.collection.get_integrand_names() {
                let Ok(integrand) = process.get_integrand(integrand_name) else {
                    continue;
                };
                let Ok(serialized) = serialize_runtime_named_settings(integrand.get_settings())
                else {
                    continue;
                };

                entries.push(ProcessSettingsCompletionEntry {
                    process_id,
                    process_name: process.definition.folder_name.clone(),
                    integrand_name: integrand_name.to_string(),
                    quantities: serialized.quantities,
                    observables: serialized.observables,
                    selectors: serialized.selectors,
                });
            }
        }

        entries.sort_by(|left, right| {
            left.process_id
                .cmp(&right.process_id)
                .then_with(|| left.integrand_name.cmp(&right.integrand_name))
        });
        entries
    }

    pub fn current_model_parameter_entries(
        &self,
    ) -> Vec<crate::repl::ModelParameterCompletionEntry> {
        let mut entries = self
            .state
            .model_parameters
            .keys()
            .filter_map(|name| {
                self.state
                    .model
                    .get_parameter_opt(name.to_string())
                    .map(|parameter| crate::repl::ModelParameterCompletionEntry {
                        name: name.to_string(),
                        parameter_type: parameter.parameter_type.clone(),
                    })
            })
            .collect::<Vec<_>>();
        entries.sort_by(|left, right| left.name.cmp(&right.name));
        entries
    }

    pub fn current_model_particle_names(&self) -> Vec<String> {
        let mut names = self
            .state
            .model
            .particles
            .iter()
            .flat_map(|particle| [particle.name.to_string(), particle.antiname.to_string()])
            .collect::<Vec<_>>();
        names.sort();
        names.dedup();
        names
    }

    pub fn current_model_coupling_names(&self) -> Vec<String> {
        let mut names = self
            .state
            .model
            .orders
            .iter()
            .map(|order| order.name.to_string())
            .collect::<Vec<_>>();
        names.sort();
        names.dedup();
        names
    }

    pub fn current_model_vertex_names(&self) -> Vec<String> {
        let mut names = self
            .state
            .model
            .vertex_rules
            .iter()
            .map(|vertex_rule| vertex_rule.0.name.to_string())
            .collect::<Vec<_>>();
        names.sort();
        names.dedup();
        names
    }

    pub fn dismiss_pending_commands_block(&mut self, trigger: &str) -> bool {
        let Some(pending) = self.session_state.pending_commands_block.take() else {
            return false;
        };

        info!(
            "{} {} {} {}",
            "Dismissing command block definition".blue(),
            pending.name.green(),
            "after".blue(),
            trigger.green()
        );
        true
    }

    fn execute_prepared_commands(
        &mut self,
        commands: Vec<PreparedCommand>,
        history_mode: HistoryMode,
    ) -> Result<ControlFlow<SaveState>> {
        for command in commands {
            if let ControlFlow::Break(save_state) = self.execute_prepared(command, history_mode)? {
                return Ok(ControlFlow::Break(save_state));
            }
        }
        Ok(ControlFlow::Continue(()))
    }

    fn execute_prepared(
        &mut self,
        command: PreparedCommand,
        history_mode: HistoryMode,
    ) -> Result<ControlFlow<SaveState>> {
        match command {
            PreparedCommand::Plain(command) => self.execute_plain(command, history_mode),
            PreparedCommand::Run { command, plan } => self.execute_run(command, plan, history_mode),
        }
    }

    fn execute_plain(
        &mut self,
        command: CommandHistory,
        history_mode: HistoryMode,
    ) -> Result<ControlFlow<SaveState>> {
        if self.session_state.pending_commands_block.is_some() {
            return self.handle_recording_mode(command);
        }

        let display_text = display_command(&command);
        info!("{} {}", "Running command".blue(), display_text.green());

        let record = match &command.command {
            Commands::StartCommandsBlock(start) => {
                self.start_commands_block(start.clone())?;
                false
            }
            Commands::FinishCommandsBlock => {
                self.finish_commands_block()?;
                false
            }
            _ => {
                if let ControlFlow::Break(save_state) = command.command.clone().run(
                    self.state,
                    self.run_history,
                    self.cli_settings,
                    self.default_runtime_settings,
                )? {
                    if history_mode == HistoryMode::Record {
                        self.record_command(&command);
                    }
                    return Ok(ControlFlow::Break(save_state));
                }
                history_mode == HistoryMode::Record
            }
        };

        if record {
            self.record_command(&command);
        }

        Ok(ControlFlow::Continue(()))
    }

    fn execute_run(
        &mut self,
        command: CommandHistory,
        plan: PreparedRun,
        history_mode: HistoryMode,
    ) -> Result<ControlFlow<SaveState>> {
        if self.session_state.pending_commands_block.is_some() {
            return self.handle_recording_mode(command);
        }

        if plan.is_empty() {
            return Ok(ControlFlow::Continue(()));
        }

        let display_text = display_command(&command);
        info!("{} {}", "Running command".blue(), display_text.green());

        for block in plan.blocks {
            info!("{} {}", "Starting command block".blue(), block.name.green());
            if let ControlFlow::Break(save_state) =
                self.execute_prepared_commands(block.commands, HistoryMode::Suppress)?
            {
                return Ok(ControlFlow::Break(save_state));
            }
        }

        if !plan.commands.is_empty() {
            info!(
                "{} {}",
                "Starting commands supplied with".blue(),
                "--commands/-c".green()
            );
        }

        if let ControlFlow::Break(save_state) =
            self.execute_prepared_commands(plan.commands, HistoryMode::Suppress)?
        {
            return Ok(ControlFlow::Break(save_state));
        }

        if history_mode == HistoryMode::Record {
            self.record_command(&command);
        }

        Ok(ControlFlow::Continue(()))
    }

    fn handle_recording_mode(&mut self, command: CommandHistory) -> Result<ControlFlow<SaveState>> {
        match command.command.clone() {
            Commands::StartCommandsBlock(_) => Err(eyre!(
                "Cannot start a new command block definition before finishing the current one"
            )),
            Commands::FinishCommandsBlock => {
                self.finish_commands_block()?;
                Ok(ControlFlow::Continue(()))
            }
            Commands::Quit(_) => {
                self.dismiss_pending_commands_block("quit");
                Ok(ControlFlow::Continue(()))
            }
            _ => {
                let stored = normalize_command_history(&command);
                let pending = self
                    .session_state
                    .pending_commands_block
                    .as_mut()
                    .expect("recording mode requires a pending block");
                pending.commands.push(stored);
                Ok(ControlFlow::Continue(()))
            }
        }
    }

    fn start_commands_block(&mut self, start: StartCommandsBlock) -> Result<()> {
        if self.session_state.pending_commands_block.is_some() {
            return Err(eyre!(
                "Cannot start a new command block definition before finishing the current one"
            ));
        }

        if self.run_history.command_block(&start.name).is_some() {
            info!(
                "{} {}",
                "Overwriting existing command block definition".blue(),
                start.name.green()
            );
        } else {
            info!(
                "{} {}",
                "Defining new command block definition".blue(),
                start.name.green()
            );
        }

        self.session_state.pending_commands_block = Some(PendingCommandsBlock {
            name: start.name,
            commands: Vec::new(),
        });
        Ok(())
    }

    fn finish_commands_block(&mut self) -> Result<()> {
        let Some(pending) = self.session_state.pending_commands_block.take() else {
            return Err(eyre!(
                "No command block definition is currently being recorded"
            ));
        };

        let block = CommandsBlock {
            name: pending.name.clone(),
            commands: pending.commands,
        };

        if let Some(existing_index) = self
            .run_history
            .command_blocks
            .iter()
            .position(|existing| existing.name == pending.name)
        {
            self.run_history.command_blocks[existing_index] = block;
        } else {
            self.run_history.command_blocks.push(block);
        }

        self.run_history.validate()?;
        Ok(())
    }

    fn record_command(&mut self, command: &CommandHistory) {
        let normalized = normalize_command_history(command);
        self.run_history
            .push_with_raw(normalized.command, normalized.raw_string);
    }
}

fn normalize_command_history(command: &CommandHistory) -> CommandHistory {
    let raw_string = command
        .raw_string
        .as_deref()
        .filter(|raw| raw_round_trips(raw, &command.command))
        .map(str::to_string)
        .or_else(|| match &command.command {
            Commands::Run(run) => Some(run.canonical_raw_string()),
            _ => None,
        });

    CommandHistory {
        command: command.command.clone(),
        raw_string,
    }
}

fn raw_round_trips(raw: &str, command: &Commands) -> bool {
    CommandHistory::from_raw_string(raw)
        .map(|parsed| parsed.command == *command)
        .unwrap_or(false)
}

pub(crate) fn display_command(command: &CommandHistory) -> String {
    if let Some(raw_string) = normalize_command_history(command).raw_string {
        raw_string
    } else {
        format!("{:?}", command.command)
    }
}
