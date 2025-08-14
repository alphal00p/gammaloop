use clap::Args;
use color_eyre::Result;
use log::debug;

use crate::{GenerationSettings, RuntimeSettings};

use super::state::State;

#[derive(Debug, Args)]
/// Generate integrands
pub struct Generate {
    /// Enable thresholds
    #[arg(long, default_value_t = false)]
    pub disable_thresholds: bool,
}

impl Generate {
    pub fn run(
        &self,
        state: &mut State,
        generation_settings: &mut GenerationSettings,
        runtime_settings: &mut RuntimeSettings,
    ) -> Result<()> {
        debug!("Preprocessing");

        let mut generation_settings_for_this_generation = generation_settings.clone();
        if self.disable_thresholds {
            debug!("Thresholds are disabled");
            generation_settings_for_this_generation.enable_thresholds = false;
        }
        state
            .process_list
            .preprocess(&state.model, generation_settings_for_this_generation)?;
        debug!("Generating integrands");
        state
            .process_list
            .generate_integrands(runtime_settings, &state.model)
    }
}
