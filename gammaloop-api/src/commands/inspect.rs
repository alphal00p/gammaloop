use clap::Args;
use eyre::Ok;
use gammalooprs::utils::F;
use ndarray::{Array, ArrayBase, OwnedRepr, ViewRepr};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;

use crate::{state::State, status_info};
use color_eyre::{Result, Section};

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
    pub fn run(&self, state: &mut State) -> Result<(Option<f64>, Complex<f64>)> {
        state.process_list.warm_up(&state.model)?;

        let process_id = state.process_list.find_process(self.process_id)?;

        let integrand_name = state.process_list.processes[process_id]
            .collection
            .find_integrand(self.integrand_name.clone())
            .with_note(|| format!("in process id {process_id}"))?;

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;

        let pt = self.point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();

        let settings = integrand.get_settings().clone();

        let (inspect_res_jac, inspect_res_eval) = gammalooprs::inspect::inspect(
            &settings,
            integrand,
            &state.model,
            pt,
            &self.discrete_dim,
            self.force_radius,
            self.momentum_space,
            self.use_f128,
        );
        let res_to_return: Complex<f64> = if let Some(jac) = inspect_res_jac {
            status_info!("Jacobian for this point: {:+.6e}", jac);
            if jac == 0. {
                return Err(color_eyre::eyre::eyre!(
                    "Jacobian is zero at this point, cannot divide by zero."
                ));
            }
            let r = inspect_res_eval.map(|a| a.into());
            r / jac
        } else {
            inspect_res_eval.map(|a| a.into())
        };
        status_info!("Result: {}", inspect_res_eval);

        Ok((inspect_res_jac, res_to_return))

        //Ok(res.map(|a| a.0))
    }
}

pub struct BatchedInspect<'a> {
    pub process_id: Option<usize>,
    pub integrand_name: Option<String>,
    pub use_f128: bool,
    pub momentum_space: bool,
    pub points: ArrayBase<ViewRepr<&'a f64>, ndarray::Dim<[usize; 2]>>,
    pub discrete_dims: ArrayBase<ViewRepr<&'a usize>, ndarray::Dim<[usize; 2]>>,
}

impl<'a> BatchedInspect<'a> {
    // todo add jacobians to output
    pub fn run(
        &self,
        state: &mut State,
    ) -> Result<(
        ArrayBase<OwnedRepr<Complex<f64>>, ndarray::Dim<[usize; 1]>>,
        Option<ArrayBase<OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>>,
    )> {
        state.process_list.warm_up(&state.model)?;

        let process_id = state.process_list.find_process(self.process_id)?;

        let integrand_name = state.process_list.processes[process_id]
            .collection
            .find_integrand(self.integrand_name.clone())
            .with_note(|| format!("in process id {process_id}"))?;

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;

        let settings = integrand.get_settings().clone();

        let mut vec_res = vec![];
        let mut jac_res = vec![];

        for (point, discrete_dim) in self
            .points
            .outer_iter()
            .zip(self.discrete_dims.outer_iter())
        {
            let pt = point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();
            let discrete_dim = discrete_dim.iter().copied().collect::<Vec<usize>>();

            let (inspect_res_jac, inspect_res_eval) = gammalooprs::inspect::inspect(
                &settings,
                integrand,
                &state.model,
                pt,
                &discrete_dim,
                false,
                self.momentum_space,
                self.use_f128,
            );
            let res_to_return: Complex<f64> = if let Some(jac) = inspect_res_jac {
                jac_res.push(jac);
                if jac == 0. {
                    return Err(color_eyre::eyre::eyre!(
                        "Jacobian is zero at this point, cannot divide by zero."
                    ));
                }
                let r = inspect_res_eval.map(|a| a.into());
                r / jac
            } else {
                inspect_res_eval.map(|a| a.into())
            };

            vec_res.push(res_to_return);
        }

        let res_array = Array::from_vec(vec_res);
        let res_jac = if jac_res.is_empty() {
            None
        } else {
            Some(Array::from_vec(jac_res))
        };

        Ok((res_array, res_jac))
    }
}
