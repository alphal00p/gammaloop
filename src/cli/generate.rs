use clap::Args;
use color_eyre::Result;
use log::debug;

use crate::{ProcessSettings, Settings};

use super::state::State;

#[derive(Debug, Args)]
/// Generate integrands
pub struct Generate {}

impl Generate {
    pub fn run(&self, state: &mut State, settings: &mut Settings) -> Result<()> {
        debug!("Preprocessing");
        state
            .process_list
            .preprocess(&state.model, ProcessSettings::default())?;
        debug!("Generating integrands");
        state
            .process_list
            .generate_integrands(&settings, &state.model)
    }
}
