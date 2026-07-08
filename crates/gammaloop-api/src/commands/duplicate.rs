use clap::{Args, Subcommand};
use color_eyre::Result;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Duplicate {
    Integrand(DuplicateIntegrand),
}

impl Duplicate {
    pub fn run(self, state: &mut State) -> Result<()> {
        match self {
            Duplicate::Integrand(command) => command.run(state),
        }
    }
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct DuplicateIntegrand {
    /// Source process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// Source integrand name
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// Destination process name
    #[arg(long = "output_process_name", value_name = "NAME")]
    pub output_process_name: String,

    /// Destination integrand name
    #[arg(long = "output_integrand_name", value_name = "NAME")]
    pub output_integrand_name: String,
}

impl DuplicateIntegrand {
    pub fn run(self, state: &mut State) -> Result<()> {
        state.duplicate_integrand(
            self.process.as_ref(),
            self.integrand_name.as_ref(),
            &self.output_process_name,
            &self.output_integrand_name,
        )
    }
}
