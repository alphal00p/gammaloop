use clap::Args;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;

use crate::{status_info, utils::F};
use color_eyre::Result;

use super::state::State;

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Inspect {
    /// The process id to inspect
    #[arg(short = 'i', long = "process-id", value_name = "ID")]
    pub process_id: usize,
    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    pub process_name: String,
    /// The point to inspect (x y) or (p0 px ...)
    #[arg(short = 'p', num_args = 2.., value_name = "POINT")]
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

impl Inspect {
    pub fn run(&self, state: &mut State) -> Result<Complex<f64>> {
        state.process_list.warm_up()?;

        let integrand = state
            .process_list
            .get_integrand_mut(self.process_id, &self.process_name)?;

        let pt = self.point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();

        let settings = integrand.get_settings().clone();

        let res = crate::inspect::inspect(
            &settings,
            integrand,
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
