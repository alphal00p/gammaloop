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
    #[arg(
        short = 'd',
        long = "discrete-dim",
        value_name = "DIMS",
        num_args = 1..,
        value_delimiter = ',',
        conflicts_with_all = ["graph_id", "orientation_id"],
    )]
    pub discrete_dim: Vec<usize>,

    /// Select a specific graph in momentum-space inspect
    #[arg(long = "graph-id", value_name = "GRAPH_ID")]
    pub graph_id: Option<usize>,

    /// Select a specific orientation of the selected graph in momentum-space inspect
    #[arg(
        long = "orientation-id",
        value_name = "ORIENTATION_ID",
        requires = "graph_id",
        conflicts_with = "discrete_dim"
    )]
    pub orientation_id: Option<usize>,
}

impl Inspect {
    fn validate_selector_mode(&self) -> Result<()> {
        if !self.momentum_space && (self.graph_id.is_some() || self.orientation_id.is_some()) {
            return Err(eyre!(
                "Graph and orientation selectors are only supported in momentum-space inspect."
            ));
        }
        Ok(())
    }

    fn resolve_graph_name(
        &self,
        state: &mut State,
        process_id: usize,
        integrand_name: &str,
    ) -> Result<Option<String>> {
        let Some(graph_id) = self.graph_id else {
            return Ok(None);
        };
        let integrand = state
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        let graph_name = integrand.graph_name_by_id(graph_id).ok_or_else(|| {
            eyre!(
                "Graph id {} is out of range for integrand '{}'; it has {} graphs.",
                graph_id,
                integrand_name,
                integrand.graph_count()
            )
        })?;
        Ok(Some(graph_name.to_string()))
    }

    fn validate_x_space_point(
        &self,
        state: &mut State,
        process_id: usize,
        integrand_name: &str,
    ) -> Result<()> {
        let integrand = state
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        let expected_dimension = integrand.expected_x_space_dimension(&self.discrete_dim)?;
        if self.point.len() != expected_dimension {
            return Err(eyre!(
                "Expected {} x-space coordinates for this integrand selection, got {}.",
                expected_dimension,
                self.point.len()
            ));
        }
        Ok(())
    }

    pub fn run(&self, state: &mut State) -> Result<(Option<f64>, Complex<f64>)> {
        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;
        self.validate_selector_mode()?;
        if !self.momentum_space {
            self.validate_x_space_point(state, process_id, &integrand_name)?;
        }
        let graph_name = self.resolve_graph_name(state, process_id, &integrand_name)?;
        let points = Array2::from_shape_vec((1, self.point.len()), self.point.clone())?;
        let discrete_dims =
            Array2::from_shape_vec((1, self.discrete_dim.len()), self.discrete_dim.clone())?;
        let result = evaluate_sample(
            state,
            &EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(integrand_name.clone()),
                use_arb_prec: self.use_arb_prec,
                minimal_output: false,
                return_generated_events: Some(true),
                momentum_space: self.momentum_space,
                points: points.view(),
                integrator_weights: None,
                discrete_dims: Some(discrete_dims.view()),
                graph_names: graph_name.map(|name| vec![Some(name)]),
                orientations: self
                    .orientation_id
                    .map(|orientation| vec![Some(orientation)]),
            },
        )?;
        let evaluation = &result.sample.evaluation;

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
        info!("\n{}", result);

        Ok((jacobian, displayed_result.map(|entry| entry.0)))
    }
}
