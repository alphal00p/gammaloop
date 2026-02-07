use std::{fs, path::PathBuf};

use clap::Args;
use gammalooprs::processes::AmplitudeGraph;

use gammalooprs::processes::{Amplitude, ProcessCollection};

use gammalooprs::utils::vakint;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use color_eyre::Result;
use colored::Colorize;
use gammalooprs::settings::RuntimeSettings;
use symbolica::atom::{Atom, AtomCore};
use tracing::{info, warn};

use crate::{
    state::{ProcessRef, State},
    CLISettings,
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "EvaluationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Evaluate {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(short = 'p', long = "process", value_name = "PROCESS")]
    pub process: Option<ProcessRef>,

    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    pub graphs_group_name: Option<String>,

    /// The path to store results in
    #[arg(short = 'o', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,

    /// Whether to evaluate numerically or not the resulting analytical expression
    #[arg(short = 'm', long = "numerical")]
    pub numerical: bool,

    /// The number of terms in the epsilon expansion to compute
    /// Defaults is automatically inferred from input graph to get to the order O(epsilon^0)
    #[arg(short = 'e', long = "n-epsilon-terms")]
    pub number_of_terms_in_epsilon_expansion: Option<usize>,
}

impl Evaluate {
    pub fn run(
        &self,
        state: &mut State,
        global_cli_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<Atom> {
        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.graphs_group_name.as_ref())?;

        let amplitude: &Amplitude = match &state.process_list.processes[process_id].collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes.get(&integrand_name).unwrap(),
            ProcessCollection::CrossSections(_) => {
                return Err(color_eyre::eyre::eyre!(
                    "Evaluate command does not support cross-section graphs"
                ));
            }
        };

        let mut true_settings = global_cli_settings
            .global
            .generation
            .uv
            .vakint
            .true_settings();

        let vakint = vakint()?;

        if let Some(n_terms) = self.number_of_terms_in_epsilon_expansion {
            true_settings.number_of_terms_in_epsilon_expansion = n_terms as i64;
        }

        let mut full_evaluation = Atom::Zero;

        for graph_term in amplitude.graphs.iter() {
            let g = &graph_term.graph;
            let mut complete_evaluation_for_this_graph = Atom::num(1);
            if g.n_externals() != 0 {
                return Err(color_eyre::eyre::eyre!(
                    "Graph named: {} has external legs. Analytical evaluation in gammaloop only supports vacuum graphs.",
                    graph_term.graph.name
                ));
            }
            let connected_components = g.connected_components(&g.full_filter());
            if connected_components.len() > 1 && g.global_prefactor.num != Atom::num(1) {
                warn!(
                    "Graph named: {} has more than one ({}) connected components and its expression contains a global numerator ({}), which will only be applied to the first connected component. Make sure this is intended.",
                    graph_term.graph.name,
                    connected_components.len(),
                    g.global_prefactor.num
                );
            }
            for (i_gc, gc) in connected_components.iter().enumerate() {
                info!(
                    "Evaluating connected component #{} of graph named: {}",
                    i_gc + 1,
                    graph_term.graph.name.blue()
                );
                complete_evaluation_for_this_graph *= graph_term.analytical_evaluation(
                    &state.model,
                    gc,
                    self.numerical,
                    vakint,
                    &true_settings,
                    &global_cli_settings.global.generation.uv.vakint,
                    default_runtime_settings,
                    // Only include the overall global numerator on the first of the connected components
                    i_gc == 0,
                )?;
            }
            complete_evaluation_for_this_graph *= &g.global_prefactor.projector * &g.overall_factor;

            full_evaluation += complete_evaluation_for_this_graph;
        }

        if let Some(p) = self.result_path.clone() {
            info!(
                "Saving evaluation result for process {} to path: {}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                p.display()
            );
            fs::write(
                p,
                toml::to_string_pretty(&full_evaluation.to_canonical_string())?,
            )?;
        } else if self.numerical {
            let numerical_evaluation_result =
                AmplitudeGraph::to_numerical(full_evaluation.as_view(), &true_settings)?;
            info!(
                "Numerical evaluation of the analytical result for process {}:\n{}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                numerical_evaluation_result
            );
        } else {
            info!(
                "Analytical result for process {}:\n{}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                full_evaluation
            );
        }

        Ok(full_evaluation)
    }
}
