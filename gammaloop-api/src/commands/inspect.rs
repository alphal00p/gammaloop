use clap::Args;
use gammalooprs::utils::F;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;

use crate::{state::State, status_info};
use color_eyre::Result;

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Inspect {
    /// The process id to inspect
    #[arg(short = 'i', long = "process-id", value_name = "ID")]
    pub process_id: Option<usize>,
    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    pub integrand_name: Option<String>,
    /// The point to inspect (x y) or (p0 px ...)
    #[arg(short = 'p', num_args = 2.., value_name = "POINT",
          value_delimiter = ',',          // allow --vals 1,-2,3
          allow_negative_numbers = true,  // treat -2 as a value, not a flag
    )]
    // allow >2 for momentum‑space
    pub point: Vec<f64>,

    /// Evaluate in f128 precision
    #[arg(short = 'f', long = "use_f128")]
    pub use_f128: bool,

    /// Force the radius in the parameterisation
    #[arg(long)]
    pub force_radius: bool,

    /// Interpret point as momentum‑space coordinates
    #[arg(short = 'm', long)]
    pub momentum_space: bool,

    /// The discrete dimensions of the sample
    #[arg(short = 'd', long = "discrete-dim", value_name = "DIMS")]
    pub discrete_dim: Vec<usize>,
}

impl Default for Inspect {
    fn default() -> Self {
        Self {
            process_id: None,
            integrand_name: None,
            point: vec![],
            use_f128: false,
            force_radius: false,
            momentum_space: false,
            discrete_dim: vec![],
        }
    }
}

impl Inspect {
    pub fn run(&self, state: &mut State) -> Result<Complex<f64>> {
        state.process_list.warm_up(&state.model)?;

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

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;

        let pt = self.point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();

        let settings = integrand.get_settings().clone();

        let res = gammalooprs::inspect::inspect(
            &settings,
            integrand,
            &state.model,
            pt,
            &self.discrete_dim,
            self.force_radius,
            self.momentum_space,
            self.use_f128,
        );

        status_info!("Result: {}", res);

        Ok(res.map(|a| a.0))
    }
}
