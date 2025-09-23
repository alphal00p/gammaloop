use std::path::PathBuf;

use clap::{arg, Subcommand};
use model::ImportModel;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::state::State;
use color_eyre::Result;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Import {
    Model(ImportModel),
    Amplitude {
        // #[arg(short = 'p')]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: PathBuf,

        #[arg(short = 'i')]
        process_id: Option<usize>,

        #[arg(short = 'n')]
        integrand_name: Option<String>,
    },
}

pub mod model;

impl Import {
    pub fn run(self, state: &mut State) -> Result<()> {
        match self {
            Import::Amplitude {
                path,
                process_id,
                integrand_name,
            } => state.import_amplitude(path, None, process_id, integrand_name),
            Import::Model(im) => im.run(state),
        }
    }
}
