use std::fs::File;
use std::{fs, path::PathBuf};

use crate::state::{ProcessListExt, ProcessRef, State};
use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use gammalooprs::utils::symbolica_ext::TypstFormat;
use gammalooprs::uv::UVgenerationSettings;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::io::Write as _;
use symbolica::atom::{Atom, AtomCore};
use symbolica::printer::PrintOptions;
use tracing::info;

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "RenormalizeSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Renormalize {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(short = 'p', long = "process", value_name = "PROCESS")]
    pub process: Option<ProcessRef>,

    /// The name of the process to inspect
    #[arg(short = 'n', long = "name", value_name = "NAME")]
    pub integrand_name: Option<String>,

    /// The directory to store per-graph results in
    #[arg(short = 'o', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,
}

impl Renormalize {
    pub fn run(&self, state: &mut State) -> Result<Vec<Atom>> {
        let (process_id, _) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;
        let amplitude = state
            .process_list
            .get_amplitude_mut_ref(self.process.as_ref(), self.integrand_name.as_ref())?;

        let settings = UVgenerationSettings::default();

        let mut renormalization_part: Vec<Atom> = Vec::new();
        let output_dir = if let Some(path) = self.result_path.clone() {
            fs::create_dir_all(&path)?;
            Some(path)
        } else {
            None
        };

        const MAX_PRINT_CHARS: usize = 4000;
        let clip_expr = |expr: &Atom| -> (String, bool) {
            let expr_str = expr.to_string();
            if expr_str.len() > MAX_PRINT_CHARS {
                (format!("{}...", &expr_str[..MAX_PRINT_CHARS]), true)
            } else {
                (expr_str, false)
            }
        };

        for (index, graph_term) in amplitude.graphs.iter_mut().enumerate() {
            if graph_term.derived_data.cff_expression.is_none() {
                return Err(color_eyre::eyre::eyre!(
                    "Graph '{}' has no CFF expression. Generate integrands before renormalizing.",
                    graph_term.graph.name
                ));
            }

            let part = graph_term.renormalization_part(&settings)?;
            let (expr_to_print, clipped) = clip_expr(&part);
            if clipped {
                info!(
                    "Renormalization part for graph '{}' (clipped):\n{}",
                    graph_term.graph.name, expr_to_print
                );
            } else {
                info!(
                    "Renormalization part for graph '{}':\n{}",
                    graph_term.graph.name, expr_to_print
                );
            }

            if let Some(dir) = &output_dir {
                let sanitized_name = graph_term
                    .graph
                    .name
                    .chars()
                    .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
                    .collect::<String>();
                let file_name = format!("{:03}_{}.txt", index + 1, sanitized_name);
                let typst_name = format!("{:03}_{}.typ", index + 1, sanitized_name);
                let mut path = File::create(dir.join(file_name))?;
                write!(path, "{}", part.printer(PrintOptions::file()))?;
                let mut path = File::create(dir.join(typst_name))?;
                write!(path, "{}", part.typst_string())?;
            }
            renormalization_part.push(part);
        }

        if let Some(dir) = &output_dir {
            info!(
                "Saved renormalization part per graph for process {} in directory: {}",
                state.process_list.processes[process_id]
                    .definition
                    .folder_name
                    .green(),
                dir.display()
            );
        }

        Ok(renormalization_part)
    }
}
