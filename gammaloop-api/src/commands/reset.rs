use clap::Subcommand;
use gammalooprs::processes::ProcessList;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::state::State;
use color_eyre::Result;
#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Reset {
    Processes,
}

impl Reset {
    pub fn run(&self, state: &mut State) -> Result<()> {
        match self {
            Self::Processes => {
                let n_processes = state.process_list.processes.len();
                state.process_list = ProcessList::default();
                info!(
                    "All {} processes have been cleared from the current state.",
                    n_processes
                );
            }
        }
        Ok(())
    }
}
