use clap::Subcommand;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::state::State;
use color_eyre::Result;
use eyre::eyre;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Display {
    Model,
    Processes,
    Integrands { process_id: Option<usize> },
}

impl Display {
    pub fn run(&self, state: &State) -> Result<()> {
        match self {
            Display::Integrands { process_id } => {
                let process = if let Some(proc_id) = process_id {
                    if *proc_id >= state.process_list.processes.len() {
                        return Err(eyre!(
                            "Process ID {} invalid, only {} processes available",
                            proc_id,
                            state.process_list.processes.len()
                        ));
                    }
                    &state.process_list.processes[*proc_id]
                } else {
                    if state.process_list.processes.len() == 1 {
                        &state.process_list.processes[0]
                    } else {
                        return Err(eyre!(
                            "Multiple processes available, please specify the process."
                        ));
                    }
                };

                println!("Integrands for process {}:", process.definition.process_id);
                for integrand in process.get_integrand_names() {
                    println!("  {}", integrand);
                }
            }
            Display::Processes => {
                println!("Processes:");
                for process in state.process_list.processes.iter() {
                    println!(
                        "#{:-10}  {}",
                        process.definition.process_id, process.definition.folder_name
                    );
                }
            }
            Display::Model => {
                println!("{}", state.model.name)
            }
        }
        Ok(())
    }
}
