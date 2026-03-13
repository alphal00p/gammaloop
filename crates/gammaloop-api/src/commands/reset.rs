use clap::Subcommand;
use gammalooprs::processes::ProcessList;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
};
use color_eyre::Result;
use colored::Colorize;
use eyre::eyre;
#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Reset {
    Processes {
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: Option<ProcessRef>,

        /// Restrict the reset to a single integrand within the selected process
        #[arg(
            short = 'i',
            long = "integrand-name",
            value_name = "NAME",
            completion_integrand_selector(crate::completion::SelectorKind::Any)
        )]
        integrand_name: Option<String>,
    },
}

impl Reset {
    pub fn run(&self, state: &mut State) -> Result<()> {
        match self {
            Self::Processes {
                process,
                integrand_name,
            } => {
                if process.is_none() && integrand_name.is_some() {
                    return Err(eyre!(
                        "{}",
                        "--integrand-name requires --process for `reset processes`"
                    )
                    .into());
                }

                if let Some(process) = process {
                    if let Some(integrand_name) = integrand_name {
                        let removed = state.remove_integrand(process, integrand_name)?;
                        if removed.removed_empty_process {
                            info!(
                                "{} {} {} {} {}",
                                "Removed integrand".blue(),
                                removed.integrand_name.green(),
                                "from process".blue(),
                                removed.process_name.green(),
                                "(process became empty and was removed)".blue()
                            );
                        } else {
                            info!(
                                "{} {} {} {}",
                                "Removed integrand".blue(),
                                removed.integrand_name.green(),
                                "from process".blue(),
                                removed.process_name.green()
                            );
                        }
                        return Ok(());
                    }

                    let removed = state.remove_process(Some(process))?;
                    info!(
                        "{} {}",
                        "Removed process".blue(),
                        removed.process_name.green()
                    );
                    return Ok(());
                }

                let n_processes = state.process_list.processes.len();
                state.process_list = ProcessList::default();
                info!(
                    "{} {} {}",
                    "Cleared all".blue(),
                    n_processes.to_string().green(),
                    "processes from the current state.".blue()
                );
            }
        }
        Ok(())
    }
}
