use clap::Subcommand;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::state::State;
use color_eyre::Result;
use eyre::eyre;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Display {
    Model {
        #[arg(short = 'c', long = "show-couplings", default_value_t = false)]
        show_couplings: bool,
        #[arg(short = 'v', long = "show-vertices", default_value_t = false)]
        show_vertices: bool,
        #[arg(short = 'r', long = "show-parameters", default_value_t = false)]
        show_parameters: bool,
        #[arg(short = 'p', long = "show-particles", default_value_t = false)]
        show_particles: bool,
        #[arg(short = 'a', long = "show-all", default_value_t = false)]
        show_all: bool,
    },
    Processes,
    Integrands {
        process_id: Option<usize>,
    },
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
                } else if state.process_list.processes.len() == 1 {
                    &state.process_list.processes[0]
                } else {
                    return Err(eyre!(
                        "Multiple processes available, please specify the process."
                    ));
                };

                info!("Integrands for process {}:", process.definition.process_id);
                for integrand in process.get_integrand_names() {
                    info!("  {}", integrand);
                }
            }
            Display::Processes => {
                info!("Processes:");
                for process in state.process_list.processes.iter() {
                    info!(
                        "#{:-10}  {}",
                        process.definition.process_id, process.definition.folder_name
                    );
                }
            }
            Display::Model {
                show_couplings,
                show_vertices,
                show_parameters,
                show_particles,
                show_all,
            } => {
                info!(
                    "\n{}",
                    state.model.get_description(
                        *show_particles || *show_all,
                        *show_parameters || *show_all,
                        *show_vertices || *show_all,
                        *show_couplings || *show_all,
                    )
                )
            }
        }
        Ok(())
    }
}
