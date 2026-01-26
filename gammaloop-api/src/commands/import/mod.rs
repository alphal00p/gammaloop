use std::path::PathBuf;

use clap::Subcommand;
use gammalooprs::graph::Graph;
use model::ImportModel;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::state::{ProcessRef, State};
use color_eyre::Result;
use eyre::eyre;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Import {
    Model(ImportModel),
    Graphs {
        // #[arg(short = 'p')]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: PathBuf,

        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(short = 'p', long = "process", value_name = "PROCESS")]
        process: Option<ProcessRef>,

        #[arg(short = 'i')]
        integrand_name: Option<String>,

        #[arg(short = 'o', default_value_t = false, conflicts_with = "append")]
        overwrite: bool,

        #[arg(short = 'a', default_value_t = false, conflicts_with = "overwrite")]
        append: bool,
    },
}

pub mod model;

impl Import {
    pub fn run(self, state: &mut State) -> Result<()> {
        match self {
            Import::Graphs {
                path,
                process,
                integrand_name,
                overwrite,
                append,
            } => {
                let default_process_name = path.file_stem().unwrap().to_string_lossy().into_owned();
                let process_count = state.process_list.processes.len();
                let (process_name, process_id) = match process {
                    Some(process_ref) => match process_ref {
                        ProcessRef::Id(id) => (None, Some(id)),
                        ProcessRef::Name(name) => {
                            let existing = state
                                .process_list
                                .processes
                                .iter()
                                .position(|p| p.definition.folder_name == name);
                            match existing {
                                Some(index) => (None, Some(index)),
                                None => (Some(name), None),
                            }
                        }
                        ProcessRef::Unqualified(value) => {
                            if let Ok(id) = value.parse::<usize>() {
                                if id >= process_count {
                                    return Err(eyre!(
                                        "Process ID {} invalid, only {} processes available",
                                        id,
                                        process_count
                                    ));
                                }
                                (None, Some(id))
                            } else {
                                let existing = state
                                    .process_list
                                    .processes
                                    .iter()
                                    .position(|p| p.definition.folder_name == value);
                                match existing {
                                    Some(index) => (None, Some(index)),
                                    None => (Some(value), None),
                                }
                            }
                        }
                    },
                    None => (Some(default_process_name), None),
                };

                let graphs = Graph::from_path(&path, &state.model)?;
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
}
