use std::path::Path;

use clap::Args;
use color_eyre::Result;
use log::{debug, info};
use serde::{Deserialize, Serialize};

use crate::{
    settings::{runtime, GlobalSettings, RuntimeSettings},
    status_info,
};

use super::state::State;

#[derive(Debug, Args, Serialize, Deserialize, Clone)]
/// Generate integrands
pub struct Generate {}

impl Generate {
    pub fn run(
        &self,
        state: &mut State,
        compile_folder: impl AsRef<Path>,
        override_existing_compiled: bool,
        generation_settings: &GlobalSettings,
        runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        status_info!("Preprocessing");

        state
            .process_list
            .preprocess(&state.model, generation_settings)?;
        status_info!("Generating integrands");
        state.process_list.generate_integrands(
            &state.model,
            generation_settings,
            runtime_settings.into(),
        )?;

        if generation_settings.generation.evaluator_settings.compile {
            status_info!("Compiling integrands");
            state.process_list.compile(
                compile_folder,
                override_existing_compiled,
                generation_settings,
                &state.model,
            )?
        }
        Ok(())
    }
}
