use std::{fs, path::PathBuf};

use clap::Args;
use gammalooprs::utils::GS;
use gammalooprs::{
    gammaloop_integrand::GLIntegrand,
    processes::{Amplitude, CrossSection, ProcessCollection},
    status_warn,
    utils::serde_utils::SmartSerde,
    uv::UltravioletGraph,
};

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
use symbolica::atom::Atom;
use tracing::info;

use crate::{state::State, CLISettings};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "EvaluationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Evaluate {
    /// The process id to inspect
    #[arg(short = 'i', long = "process-id", value_name = "ID")]
    pub process_id: Option<usize>,

    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    pub graphs_group_name: Option<String>,

    /// The path to store results in
    #[arg(short = 'p', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,
}

impl Evaluate {
    pub fn run(&self, state: &mut State, global_cli_settings: &CLISettings) -> Result<Atom> {
        let state_name: String = global_cli_settings
            .state_folder
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();

        let result_path = if let Some(p) = self.result_path.clone() {
            p
        } else {
            PathBuf::from(format!(
                "./{}_integration_workspace/integration_result.json",
                state_name
            ))
        };

        let (process_id, integrand_name) = state
            .process_list
            .find_integrand(self.process_id, self.graphs_group_name.as_ref())?;

        let amplitude: &Amplitude = match &state.process_list.processes[process_id].collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes.get(&integrand_name).unwrap(),
            ProcessCollection::CrossSections(_) => {
                return Err(color_eyre::eyre::eyre!(
                    "Evaluate command does not support cross-section graphs"
                ));
            }
        };

        for graph_term in amplitude.graphs.iter() {
            let g = &graph_term.graph;
            if g.n_externals() != 0 {
                return Err(color_eyre::eyre::eyre!(
                    "Graph named: {} has external legs. Analytical evaluation in gammaloop only supports vacuum graphs.",
                    graph_term.graph.name
                ));
            }
            for (i_gc, gc) in g.connected_components(&g.full_filter()).iter().enumerate() {
                info!(
                    "Evaluating connected component #{} of graph named: {}",
                    i_gc, graph_term.graph.name
                );
                let evaluated_expression_for_this_component =
                    graph_term.analytical_evaluation(gc)?;
            }
        }

        Ok(Atom::Zero)
    }
}
