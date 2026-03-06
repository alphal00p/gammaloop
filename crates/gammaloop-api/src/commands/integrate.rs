use std::{fs, path::PathBuf};

use clap::Args;
use gammalooprs::utils::serde_utils::SmartSerde;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;

use color_eyre::Result;
use colored::Colorize;
use gammalooprs::{
    integrate::{havana_integrate, print_integral_result, IntegrationState},
    settings::{runtime::IntegrationResult, RuntimeSettings},
    utils::F,
};
use tracing::{info, warn};

use crate::{
    state::{ProcessRef, State},
    CLISettings,
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "IntegrationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Integrate {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(short = 'p', long = "process", value_name = "PROCESS")]
    pub process: Option<ProcessRef>,

    /// The name of the process to inspect
    #[arg(short = 'i', long = "name", value_name = "NAME")]
    pub integrand_name: Option<String>,

    /// The path to store results in
    #[arg(short = 's', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,

    /// Number of cores to parallelize over
    #[arg(short = 'c', long)]
    pub n_cores: Option<usize>,

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
        let default_workspace_path = global_cli_settings
            .state_folder
            .join("integration_workspace");
        let workspace_path = if let Some(p) = self.workspace_path.clone() {
            p
        } else {
            default_workspace_path.clone()
        };

        let result_path = if let Some(p) = self.result_path.clone() {
            p
        } else {
            default_workspace_path.join("integration_result.json")
        };

        let target = self.target.clone().map(|t| Complex::new(F(t[0]), F(t[1])));

        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;

        if self.restart && workspace_path.exists() {
            fs::remove_dir_all(&workspace_path)?;
        }

        let gloop_integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;

        gloop_integrand.warm_up(&state.model)?;

        info!("Gammaloop now integrates {}", integrand_name.green().bold());

        let path_to_state = workspace_path.join("integration_state");

        let integration_state = match fs::read(path_to_state) {
            Ok(state_bytes) => {
                info!(
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

                let path_to_workspace_settings = workspace_path.join("settings.yaml");

                let workspace_settings: RuntimeSettings =
                    RuntimeSettings::from_file(path_to_workspace_settings, "workspace settings")?;
                // force the settings to be the same as the ones used in the previous integration
                if *gloop_integrand.get_mut_settings() != workspace_settings.clone() {
                    warn!("settings have changed with respect to workspace, reverting changes");
                    *gloop_integrand.get_mut_settings() = workspace_settings.clone();
                }

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
                info!("");
                warn!(
                    "Any changes to the settings will be ignored, integrate with the {} option for changes to take effect",
                    "--restart".blue()
                );
                info!("{}", "Resuming integration".yellow());

                Some(state)
            }

            Err(_) => {
                info!("No integration state found, starting new integration");
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

        let n_cores = self
            .n_cores
            .unwrap_or(global_cli_settings.global.n_cores.integrate);

        let result = havana_integrate(
            &settings,
            &state.model,
            |set| gloop_integrand.user_data_generator(n_cores, set),
            target,
            integration_state,
            Some(workspace_path.clone()),
        )?;

        fs::write(&result_path, serde_json::to_string(&result)?)?;

        Ok(result)
    }
}
