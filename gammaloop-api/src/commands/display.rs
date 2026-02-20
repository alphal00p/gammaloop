use clap::Subcommand;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::state::{ProcessRef, State};
use color_eyre::Result;

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
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(short = 'p', long = "process", value_name = "PROCESS")]
        process: Option<ProcessRef>,
    },
}

impl Display {
    pub fn run(&self, state: &State) -> Result<()> {
        match self {
            Display::Integrands { process } => {
                let process_id = state.resolve_process_ref(process.as_ref())?;
                let process = &state.process_list.processes[process_id];

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
