use clap::Subcommand;
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
pub enum Remove {
    Processes {
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: Option<ProcessRef>,

        /// Restrict the removal to a single integrand within the selected process
        #[arg(
            short = 'i',
            long = "integrand-name",
            value_name = "NAME",
            completion_integrand_selector(crate::completion::SelectorKind::Any)
        )]
        integrand_name: Option<String>,
    },
}

impl Remove {
    pub fn run(&self, state: &mut State) -> Result<()> {
        match self {
            Self::Processes {
                process,
                integrand_name,
            } => {
                if process.is_none() && integrand_name.is_some() {
                    return Err(eyre!(
                        "{}",
                        "--integrand-name requires --process for `remove processes`"
                    ));
                }

                let removed = state
                    .remove_selected_integrands(process.as_ref(), integrand_name.as_deref())?;
                for removed_integrand in removed {
                    info!(
                        "{} {} {} {}",
                        "Removed integrand".blue(),
                        removed_integrand.integrand_name.green(),
                        "from process".blue(),
                        removed_integrand.process_name.green()
                    );
                }
            }
        }
        Ok(())
    }
}
