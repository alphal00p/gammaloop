use std::{
    env,
    path::{Path, PathBuf},
};

use clap::Subcommand;
use gammalooprs::graph::Graph;
use model::ImportModel;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
    CLISettings,
};
use color_eyre::Result;
use eyre::{eyre, Context};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Import {
    Model(ImportModel),
    Graphs {
        // #[arg(short = 'p')]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: PathBuf,

        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: Option<ProcessRef>,

        #[arg(short = 'i', completion_disable_special_value())]
        integrand_name: Option<String>,

        #[arg(short = 'o', default_value_t = false, conflicts_with = "append")]
        overwrite: bool,

        #[arg(short = 'a', default_value_t = false, conflicts_with = "overwrite")]
        append: bool,
    },
}

pub mod model;

impl Import {
    pub fn run(self, state: &mut State, cli_settings: &CLISettings) -> Result<()> {
        match self {
            Import::Graphs {
                path,
                process,
                integrand_name,
                overwrite,
                append,
            } => {
                let resolved_path =
                    Self::resolve_graph_import_path(&path, &cli_settings.state.folder)?;
                let default_process_name = resolved_path
                    .file_stem()
                    .ok_or_else(|| {
                        eyre!(
                            "Could not derive a process name from graph path '{}'",
                            resolved_path.display()
                        )
                    })?
                    .to_string_lossy()
                    .into_owned();
                let (process_name, process_id) = match process {
                    Some(ProcessRef::Id(id)) => (None, Some(id)),
                    Some(ProcessRef::Name(name)) => (Some(name), None),
                    Some(ProcessRef::Unqualified(value)) => {
                        let name_match = state
                            .process_list
                            .processes
                            .iter()
                            .position(|p| p.definition.folder_name == value);
                        if let Ok(id) = value.parse::<usize>() {
                            if name_match.is_some() {
                                return Err(eyre!(
                                    "Ambiguous process reference '{}'. Use '#{}' or 'name:{}' to disambiguate.",
                                    value,
                                    id,
                                    value
                                ));
                            }
                            (None, Some(id))
                        } else {
                            (Some(value), None)
                        }
                    }
                    None => (Some(default_process_name), None),
                };

                info!(
                    "Loading graphs from '{}'",
                    Self::display_graph_import_path(&resolved_path).display()
                );
                let graphs = Graph::from_path(&resolved_path, &state.model)?;
                state.import_graphs(
                    graphs,
                    process_name,
                    process_id,
                    integrand_name,
                    overwrite,
                    append,
                )
            }
            Import::Model(im) => im.run(state),
        }
    }

    fn resolve_graph_import_path(path: &Path, state_folder: &Path) -> Result<PathBuf> {
        let cwd = env::current_dir().wrap_err(
            "Failed to query the current working directory while resolving graph import path",
        )?;

        if path.is_absolute() {
            let normalized = Self::normalize_path_lexically(path);
            return normalized
                .exists()
                .then_some(normalized)
                .ok_or_else(|| eyre!("Graph file '{}' does not exist.", path.display()));
        }

        let state_root = state_folder.parent().unwrap_or(state_folder);
        let state_root_candidate;
        let absolute_state_root = Self::normalize_path_lexically(if state_root.is_absolute() {
            state_root
        } else {
            state_root_candidate = cwd.join(state_root);
            &state_root_candidate
        });
        let state_candidate = Self::normalize_path_lexically(&absolute_state_root.join(path));
        if state_candidate.exists() {
            return Ok(state_candidate);
        }

        let cwd_candidate = Self::normalize_path_lexically(&cwd.join(path));
        if cwd_candidate.exists() {
            return Ok(cwd_candidate);
        }

        Err(eyre!(
            "Could not find graph file '{}'. Tried '{}' (active state root) and '{}' (current working directory).",
            path.display(),
            state_candidate.display(),
            cwd_candidate.display()
        ))
    }

    fn display_graph_import_path(path: &Path) -> PathBuf {
        path.canonicalize()
            .unwrap_or_else(|_| Self::normalize_path_lexically(path))
    }

    fn normalize_path_lexically(path: &Path) -> PathBuf {
        use std::path::Component;

        let mut normalized = PathBuf::new();
        for component in path.components() {
            match component {
                Component::CurDir => {}
                Component::ParentDir => {
                    if normalized.components().next_back().is_some_and(|last| {
                        !matches!(last, Component::RootDir | Component::Prefix(_))
                    }) {
                        normalized.pop();
                    } else if !path.is_absolute() {
                        normalized.push(component.as_os_str());
                    }
                }
                Component::Normal(part) => normalized.push(part),
                Component::RootDir | Component::Prefix(_) => {
                    normalized.push(component.as_os_str());
                }
            }
        }
        normalized
    }
}
