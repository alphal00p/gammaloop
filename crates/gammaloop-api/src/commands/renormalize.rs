use std::fs::File;
use std::{fs, path::PathBuf};

use crate::completion::CompletionArgExt;
use crate::state::{ProcessListExt, ProcessRef, State};
use crate::CLISettings;
use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use gammalooprs::utils::symbolica_ext::TypstFormat;
use idenso::color::{ColorSimplifier, CS};
use idenso::metric::MetricSimplifier;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::io::Write as _;
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::parse;
use symbolica::printer::PrintOptions;
use tracing::info;

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, unsendable, name = "RenormalizeSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Renormalize {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Amplitude)
    )]
    pub process: Option<ProcessRef>,

    #[arg(short = 'a', long = "align-to-rqft", default_value_t = false)]
    pub align_to_rqft: bool,

    #[arg(short = 'f', long = "print-namespaces", default_value_t = false)]
    pub print_namespaces: bool,

    /// The name of the process to inspect
    /// The amplitude name to renormalize
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Amplitude)
    )]
    pub integrand_name: Option<String>,

    /// The directory to store per-graph results in
    #[arg(short = 'o', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,
}

impl Renormalize {
    pub fn run(&self, state: &mut State, global_cli_settings: &CLISettings) -> Result<Vec<Atom>> {
        let (process_id, _) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;
        let amplitude = state
            .process_list
            .get_amplitude_mut_ref(self.process.as_ref(), self.integrand_name.as_ref())?;

        let mut settings = global_cli_settings.global.generation.uv.clone();

        settings.only_integrated = true;
        settings.generate_integrated = true;

        let mut renormalization_part: Vec<Atom> = Vec::new();
        let output_dir = if let Some(path) = self.result_path.clone() {
            fs::create_dir_all(&path)?;
            Some(path)
        } else {
            None
        };

        const MAX_PRINT_CHARS: usize = 4000;
        let clip_expr = |expr: &Atom, print_namespaces: bool| -> (String, bool) {
            let expr_str = if print_namespaces {
                format!("{:#}", expr)
            } else {
                format!("{}", expr)
            };
            if expr_str.len() > MAX_PRINT_CHARS {
                (format!("{}...", &expr_str[..MAX_PRINT_CHARS]), true)
            } else {
                (expr_str, false)
            }
        };

        for (index, graph_term) in amplitude.graphs.iter_mut().enumerate() {
            let mut part = graph_term.renormalization_part(&settings)?.expression;

            part = state
                .model
                .apply_parameter_replacement_rules(&state.model.apply_coupling_replacement_rules(
                    &part.simplify_color().expand().simplify_metrics().to_dots(),
                ))
                .collect_factors();

            if self.align_to_rqft {
                part = (part
                    .replace(CS.tr)
                    .with(Atom::num((1, 2)))
                    .replace(CS.nc)
                    .with(CS.ca)
                    .replace(parse!("UFO::aS"))
                    .with(parse!("gs").pow(2) / (Atom::var(Symbol::PI) * 4))
                    / 8)
                .expand_num()
            }

            let (expr_to_print, clipped) = clip_expr(&part, self.print_namespaces);
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
