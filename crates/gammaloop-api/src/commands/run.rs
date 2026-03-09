use std::{ops::ControlFlow, path::PathBuf};

use clap::Args;
use colored::Colorize;
use gammalooprs::settings::RuntimeSettings;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use color_eyre::Result;
use tracing::info;

use crate::{
    state::{CommandHistory, RunHistory, State},
    CLISettings,
};

use super::save::SaveState;

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Run {
    /// Path to a run file to execute
    #[arg(value_hint = clap::ValueHint::FilePath)]
    path: Option<PathBuf>,

    /// Command block names to execute (in order)
    #[arg(value_name = "BLOCK_NAME")]
    block_names: Vec<String>,

    #[arg(short = 'c', long)]
    commands: Option<String>,
}

impl Run {
    pub fn selected_block_names(&self) -> Option<&[String]> {
        (!self.block_names.is_empty()).then_some(self.block_names.as_slice())
    }

    pub fn load_run_history(&self) -> Result<Option<RunHistory>> {
        self.path
            .as_ref()
            .map(|path| RunHistory::load_with_block_selection(path, self.selected_block_names()))
            .transpose()
    }

    pub fn run(
        &self,
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
        run_history: &mut RunHistory,
    ) -> Result<ControlFlow<SaveState>> {
        let mut a: ControlFlow<SaveState> = ControlFlow::Continue(());
        if let Some(mut new_run_history) = self.load_run_history()? {
            let res = new_run_history.run(state, global_settings, default_runtime_settings);

            run_history.merge(new_run_history);
            a = res?;
        };

        if let Some(c) = self.commands.as_ref() {
            let commands: Vec<_> = c
                .split(';')
                .map(str::trim)
                .filter(|c| !c.is_empty())
                .collect();

            if !commands.is_empty() {
                info!(
                    "> {}{}{}",
                    "Starting commands supplied with the '".blue(),
                    "--comands/-c".green(),
                    "' CLI option.".blue()
                );
            }

            for c in commands {
                let cmd = CommandHistory::from_raw_string(c)?;
                info!(
                    "Running command: {}",
                    if let Some(rs) = cmd.raw_string.as_ref() {
                        rs.green()
                    } else {
                        format!("{:?}", cmd.command).normal()
                    }
                );
                let res = cmd.command.clone().run(
                    state,
                    run_history,
                    global_settings,
                    default_runtime_settings,
                )?;
                run_history.push_with_raw(cmd.command, cmd.raw_string);
                a = res;

                if let ControlFlow::Break(_) = &a {
                    break;
                }
            }
        }
        Ok(a)
    }
}
