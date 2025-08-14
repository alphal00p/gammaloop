use std::{fs, path::PathBuf};

use clap::Args;
use clarabel::solver::traits::Info;
use itertools::Itertools;
use log::{info, warn};
use spenso::algebra::complex::Complex;

use crate::{
    integrate::{havana_integrate, print_integral_result, IntegrationState},
    utils::F,
    RuntimeSettings,
};
use color_eyre::Result;
use colored::{ColoredString, Colorize};

use super::state::State;

#[derive(Debug, Args)]
pub struct Integrate {
    /// The process id to inspect
    #[arg(short = 'i', long = "process-id", value_name = "ID")]
    process_id: usize,

    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    process_name: String,

    /// The path to store results in
    #[arg(short = 'p', long)]
    result_path: PathBuf,

    /// Number of cores to parallelize over
    #[arg(short = 'c', long)]
    n_cores: usize,

    /// The path to run the integrationg within
    #[arg(short = 'w', long)]
    workspace_path: PathBuf,

    /// Specify the target integration result to compare against
    #[arg(short = 't', num_args = 2, long)]
    target: Option<Vec<f64>>,
}

impl Integrate {
    pub fn run(&self, state: &mut State) -> Result<Vec<(f64, f64)>> {
        let target = if let Some(t) = self.target.clone() {
            Some(Complex::new(F(t[0]), F(t[1])))
        } else {
            None
        };

        let gloop_integrand = state
            .process_list
            .get_integrand_mut(self.process_id, &self.process_name)?;

        info!(
            "Gammaloop now integrates {}",
            self.process_name.green().bold()
        );

        let path_to_state = self.workspace_path.join("integration_state");

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

                let path_to_workspace_settings = self.workspace_path.join("settings.yaml");
                let workspace_settings_string = fs::read_to_string(path_to_workspace_settings)?;

                let workspace_settings: RuntimeSettings =
                    serde_yaml::from_str(&workspace_settings_string)?;

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
                info!("");
                warn!("Any changes to the settings will be ignored, integrate with the {} option for changes to take effect","--restart".blue());
                info!("{}", "Resuming integration".yellow());

                Some(state)
            }

            Err(_) => {
                info!("No integration state found, starting new integration");
                None
            }
        };

        if !self.workspace_path.exists() {
            fs::create_dir_all(&self.workspace_path)?;
            info!(
                "Created workspace directory at {}",
                self.workspace_path.display()
            );
        }
        let settings = gloop_integrand.get_settings().clone();

        let result = havana_integrate(
            &settings,
            |set| gloop_integrand.user_data_generator(self.n_cores, set),
            target,
            integration_state,
            Some(self.workspace_path.clone()),
        );

        fs::write(&self.result_path, serde_yaml::to_string(&result)?)?;

        Ok(result
            .result
            .iter()
            .tuple_windows()
            .map(|(re, im)| (re.0, im.0))
            .collect())
    }
}
