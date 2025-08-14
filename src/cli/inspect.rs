use clap::Args;
use clarabel::solver::traits::Info;
use log::info;

use crate::{utils::F, RuntimeSettings};
use color_eyre::Result;

use super::state::State;

#[derive(Debug, Args)]
pub struct Inspect {
    /// The process id to inspect
    #[arg(short = 'i', long = "process-id", value_name = "ID")]
    process_id: usize,
    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    process_name: String,
    /// The point to inspect (x y) or (p0 px ...)
    #[arg(short = 'p', num_args = 2.., value_name = "POINT")]
    // allow >2 for momentum‑space
    point: Vec<f64>,

    /// Evaluate in f128 precision
    #[arg(short = 'f', long = "use_f128")]
    use_f128: bool,

    /// Force the radius in the parameterisation
    #[arg(long)]
    force_radius: bool,

    /// Interpret point as momentum‑space coordinates
    #[arg(short = 'm', long)]
    momentum_space: bool,
}

impl Inspect {
    pub fn run(&self, state: &mut State) -> Result<()> {
        let integrand = state
            .process_list
            .get_integrand_mut(self.process_id, &self.process_name)?;

        let pt = self.point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();

        let settings = integrand.get_settings().clone();

        let res = crate::inspect::inspect(
            &settings,
            integrand,
            pt,
            &[],
            self.force_radius,
            self.momentum_space,
            self.use_f128,
        );

        info!("Result: {}", res);

        Ok(())
    }
}
