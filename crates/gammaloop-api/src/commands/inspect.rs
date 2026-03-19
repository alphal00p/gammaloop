use clap::Args;
use colored::Colorize;
use ndarray::Array2;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use tracing::info;

use crate::{
    commands::evaluate_samples::{evaluate_sample, EvaluateSamples},
    completion::CompletionArgExt,
    state::{ProcessRef, State},
};
use color_eyre::Result;
use eyre::eyre;
use gammalooprs::utils::F;

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Inspect {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,
    /// The integrand name to inspect
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,
    /// The point to inspect (x y) or (p0 px ...)
    #[arg(
        short = 'x',
        long = "point",
        num_args = 2..,
        value_name = "POINT",
        value_delimiter = ',',
        allow_negative_numbers = true,
    )]
    pub point: Vec<f64>,

    /// Evaluate in f128 precision
    #[arg(short = 'f', long = "use_arb_prec")]
    pub use_arb_prec: bool,

    /// Interpret point as momentum-space coordinates
    #[arg(short = 'm', long)]
    pub momentum_space: bool,

    /// The discrete dimensions of the sample
    #[arg(short = 'd', long = "discrete-dim", value_name = "DIMS")]
    pub discrete_dim: Vec<usize>,
}

impl Inspect {
    pub fn run(&self, state: &mut State) -> Result<(Option<f64>, Complex<f64>)> {
        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;
        let points = Array2::from_shape_vec((1, self.point.len()), self.point.clone())?;
        let discrete_dims =
            Array2::from_shape_vec((1, self.discrete_dim.len()), self.discrete_dim.clone())?;
        let evaluation = evaluate_sample(
            state,
            &EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(integrand_name.clone()),
                use_arb_prec: self.use_arb_prec,
                minimal_output: false,
                momentum_space: self.momentum_space,
                points: points.view(),
                discrete_dims: Some(discrete_dims.view()),
                graph_names: None,
                orientations: None,
            },
        )?
        .sample
        .evaluation;

        let raw_result = evaluation.integrand_result;
        let jacobian = evaluation.parameterization_jacobian.map(|jac| jac.0);
        let displayed_result = if self.momentum_space {
            raw_result
        } else {
            let jacobian = jacobian.ok_or_else(|| {
                eyre!("x-space inspect requires a parameterization jacobian, but none was returned")
            })?;
            raw_result.map(|entry| entry * F(jacobian))
        };

        let point_label = if self.momentum_space {
            "Input point in momentum space"
        } else {
            "Input point in unit hypercube xs"
        };
        info!(
            "\n{}:\n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
            point_label,
            format!(
                "( {} )",
                self.point
                    .iter()
                    .map(|x| format!("{x:.16}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
            .blue(),
            integrand_name.green(),
            format!(
                "( {:+.16e}, {:+.16e} i)",
                displayed_result.re, displayed_result.im
            )
            .blue(),
        );

        if let Some(jacobian) = jacobian {
            info!(
                "Parameterization jacobian for this point: {:+.16e}",
                jacobian
            );
        }

        Ok((jacobian, displayed_result.map(|entry| entry.0)))
    }
}
