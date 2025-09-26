use std::{fs, path::PathBuf};

use clap::Args;
use gammalooprs::{status_warn, utils::serde_utils::SmartSerde};

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;

use color_eyre::Result;
use colored::Colorize;
use gammalooprs::{
    integrate::{havana_integrate, print_integral_result, IntegrationState},
    settings::{runtime::IntegrationResult, RuntimeSettings},
    status_info,
    utils::F,
};
use tracing::info;

use crate::{state::State, CLISettings};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "IntegrationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Integrate {
    /// The process id to inspect
    #[arg(short = 'i', long = "process-id", value_name = "ID")]
    pub process_id: Option<usize>,

    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    pub integrand_name: Option<String>,

    /// The path to store results in
    #[arg(short = 'p', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,

    /// The path to run the integrationg within
    #[arg(short = 'w', long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    /// Specify the target integration result to compare against
    #[arg(short = 't', num_args = 2, long,
          value_delimiter = ',',          // allow --vals 1,-2,3
          allow_negative_numbers = true,  // treat -2 as a value, not a flag
          allow_hyphen_values=true
    )]
    pub target: Option<Vec<f64>>,

    /// Whether to restart the integration from scratch, or continue from a previous run if possible
    #[arg(short = 'r', long)]
    pub restart: bool,
}

impl Integrate {
    pub fn run(
        &self,
        state: &mut State,
        global_cli_settings: &CLISettings,
    ) -> Result<IntegrationResult> {
        let state_name: String = global_cli_settings
            .state_folder
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let workspace_path = if let Some(p) = self.workspace_path.clone() {
            p
        } else {
            PathBuf::from(format!("./{}_integration_workspace", state_name))
        };

        let result_path = if let Some(p) = self.result_path.clone() {
            p
        } else {
            PathBuf::from(format!(
                "./{}_integration_workspace/integration_result.json",
                state_name
            ))
        };

        let target = if let Some(t) = self.target.clone() {
            Some(Complex::new(F(t[0]), F(t[1])))
        } else {
            None
        };

        let process_id = if let Some(id) = self.process_id {
            if id >= state.process_list.processes.len() {
                return Err(color_eyre::eyre::eyre!(
                    "Invalid process id {}. Number of processes: {}",
                    id,
                    state.process_list.processes.len()
                ));
            }
            id
        } else {
            if state.process_list.processes.is_empty() {
                return Err(color_eyre::eyre::eyre!("No processes generated yet."));
            }
            if state.process_list.processes.is_empty() {
                return Err(color_eyre::eyre::eyre!("No processes generated yet."));
            } else if state.process_list.processes.len() > 1 {
                return Err(color_eyre::eyre::eyre!(
                    "There are {} processes available. Please specify a process id.",
                    state.process_list.processes.len()
                ));
            } else {
                0
            }
        };
        let all_integrand_names = state.process_list.processes[process_id]
            .collection
            .get_integrand_names();
        let integrand_name = if let Some(name) = self.integrand_name.clone() {
            if !all_integrand_names.contains(&name.as_str()) {
                return Err(color_eyre::eyre::eyre!(
                    "No integrand named '{}' in process id {}. Available integrands: {:?}",
                    name,
                    process_id,
                    all_integrand_names
                ));
            }
            name
        } else {
            if all_integrand_names.len() != 1 {
                return Err(color_eyre::eyre::eyre!(
                    "Multiple integrands in process id {}. Please specify one of: {:?}",
                    process_id,
                    all_integrand_names
                ));
            }
            all_integrand_names[0].to_string()
        };

        state.process_list.warm_up(&state.model)?;

        if self.restart && workspace_path.exists() {
            fs::remove_dir_all(&workspace_path)?;
        }

        let gloop_integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;

        status_info!("Gammaloop now integrates {}", integrand_name.green().bold());

        let path_to_state = workspace_path.join("integration_state");

        let integration_state = match fs::read(path_to_state) {
            Ok(state_bytes) => {
                status_info!(
                    "{}",
                    "Found integration state, result of previous integration:".yellow()
                );
                info!("");

                let state: IntegrationState = bincode::decode_from_slice::<IntegrationState, _>(
                    &state_bytes,
                    bincode::config::standard(),
                )
                .expect("Could not deserialize state")
                .0;

                let path_to_workspace_settings = workspace_path.join("settings.toml");

                let workspace_settings: RuntimeSettings =
                    RuntimeSettings::from_file(path_to_workspace_settings, "workspace settings")?;
                // force the settings to be the same as the ones used in the previous integration
                *gloop_integrand.get_mut_settings() = workspace_settings.clone();

                print_integral_result(
                    &state.integral.re,
                    1,
                    state.iter,
                    "re",
                    target.map(|c| c.re),
                );

                print_integral_result(
                    &state.integral.im,
                    2,
                    state.iter,
                    "im",
                    target.map(|c| c.im),
                );
                status_info!("");
                status_warn!("Any changes to the settings will be ignored, integrate with the {} option for changes to take effect","--restart".blue());
                status_info!("{}", "Resuming integration".yellow());

                Some(state)
            }

            Err(_) => {
                status_info!("No integration state found, starting new integration");
                None
            }
        };

        if !workspace_path.exists() {
            fs::create_dir_all(&workspace_path)?;
            info!(
                "Created workspace directory at {}",
                workspace_path.display()
            );
        }
        let settings = gloop_integrand.get_settings().clone();

        let result = havana_integrate(
            &settings,
            &state.model,
            |set| {
                gloop_integrand
                    .user_data_generator(global_cli_settings.global.n_cores.integrate, set)
            },
            target,
            integration_state,
            Some(workspace_path.clone()),
        );

        fs::write(&result_path, serde_json::to_string(&result)?)?;

        Ok(result)
    }
}
