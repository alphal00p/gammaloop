use clap::Subcommand;
use gammalooprs::{processes::ProcessList, status_info};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

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
                status_info!(
                    "All {} processes have been cleared from the current state.",
                    n_processes
                );
            }
        }
        Ok(())
    }
}
