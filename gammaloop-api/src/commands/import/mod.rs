use std::path::PathBuf;

use clap::{arg, Subcommand};
use gammalooprs::graph::Graph;
use model::ImportModel;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::state::State;
use color_eyre::Result;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Import {
    Model(ImportModel),
    Graphs {
        // #[arg(short = 'p')]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: PathBuf,

        // TODO Add the ability to specify process name, and even perhaps remove that of specifying process id
        #[arg(short = 'i')]
        process_id: Option<usize>,

        #[arg(short = 'n')]
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
                process_id,
                integrand_name,
                overwrite,
                append,
            } => {
                let process_name = path.file_stem().unwrap().to_string_lossy().into_owned();

                let graphs = Graph::from_path(&path, &state.model)?;
                state.import_graphs(
                    graphs,
                    Some(process_name),
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
