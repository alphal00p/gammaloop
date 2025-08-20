use clap::Args;
use color_eyre::Result;
use log::debug;

use crate::settings::{runtime, GlobalSettings, RuntimeSettings};

use super::state::State;

#[derive(Debug, Args)]
/// Generate integrands
pub struct Generate {}

impl Generate {
    pub fn run(
        &self,
        state: &mut State,
        generation_settings: &GlobalSettings,
        runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        debug!("Preprocessing");

        state
            .process_list
            .preprocess(&state.model, generation_settings)?;
        debug!("Generating integrands");
        state
            .process_list
            .generate_integrands(&state.model, runtime_settings.into())
    }
}
