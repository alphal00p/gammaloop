use std::path::Path;

use clap::Args;
use color_eyre::Result;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::settings::{GlobalSettings, RuntimeSettings};

use super::state::State;

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema)]
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
        state.generate_integrands(generation_settings, runtime_settings.into())?;

        if generation_settings.generation.evaluator_settings.compile {
            state.compile_integrands(
                compile_folder,
                override_existing_compiled,
                generation_settings,
            )?
        }
        Ok(())
    }
}
