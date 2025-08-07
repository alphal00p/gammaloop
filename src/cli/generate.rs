use clap::Args;
use color_eyre::Result;

use crate::{ProcessSettings, Settings};

use super::state::State;

#[derive(Debug, Args)]
/// Generate integrands
pub struct Generate {}

impl Generate {
    pub fn run(state: &mut State, settings: &mut Settings) -> Result<()> {
        state
            .process_list
            .preprocess(&state.model, ProcessSettings::default())?;

        state
            .process_list
            .generate_integrands(&settings, &state.model)
    }
}
