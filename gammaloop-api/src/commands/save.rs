use std::{fs, path::PathBuf};

use clap::{arg, Args, Subcommand};
use color_eyre::{eyre::eyre, owo_colors::OwoColorize, Result};
use gammalooprs::{
    settings::RuntimeSettings,
    status_info,
    utils::serde_utils::{SmartSerde, SHOWDEFAULTS},
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::warn;

use crate::{
    state::{set_serialize_commands_as_strings, RunHistory, State},
    templates::Templates,
    write_schemas, CLISettings,
};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Save {
    Dot {
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: Option<PathBuf>,
    },
    State(SaveState),
    /// regenerate the schema files
    Schema {},
}

impl Save {
    pub fn run(
        self,
        state: &mut State,
        run_history: &RunHistory,
        default_runtime_settings: &RuntimeSettings,
        global_settings: &CLISettings,
    ) -> Result<()> {
        match self {
            Save::Dot { path } => {
                // Use original default location (state folder) or custom path if provided
                let target_dir = path.unwrap_or(global_settings.state_folder.clone());

                // Extract embedded templates to build/templates relative to target directory
                if let Err(e) = Templates::extract_to_build_dir(&target_dir) {
                    println!(
                        "Warning: Could not extract templates to drawings/templates: {}",
                        e
                    );
                }

                // Generate dynamic edge styles based on the model
                let template_path = target_dir.join("drawings/templates/edge-style.typ");
                if let Err(e) = state.model.generate_edge_style_template(&template_path) {
                    warn!("Warning: Could not generate dynamic edge styles: {}", e);
                }

                // Export dot files to original location
                state.export_dots(&target_dir)?;

                // Create Justfile with draw recipe
                let justfile_path = target_dir.join("justfile");
                let justfile_content =
                    "# Generate drawings from dot files\ndraw *INPUTS:\n    linnet --build-dir drawings --input {{INPUTS}} . -o drawings.pdf\n";
                if let Err(e) = fs::write(&justfile_path, justfile_content) {
                    warn!(
                        "Warning: Could not create justfile at {}: {}",
                        justfile_path.display(),
                        e
                    );
                }

                Ok(())
            }
            Save::State(s) => s.save(
                state,
                run_history,
                default_runtime_settings,
                global_settings,
            ),

            Save::Schema {} => write_schemas(),
        }
    }
}

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct SaveState {
    /// Path to save the state to, by default is the current state folder
    #[arg(short = 'p', long, value_hint = clap::ValueHint::FilePath)]
    pub path: Option<PathBuf>,

    /// Save state to file after each call
    #[arg(short = 'o', long,num_args(0..=1), default_missing_value = "true",
           value_parser = clap::builder::BoolishValueParser::new(),)]
    pub override_state: Option<bool>,

    #[arg(short = 'n', long, default_value_t = false)]
    pub no_save_state: bool,

    /// Try to serialize using strings when saving run history
    #[arg(long, num_args(0..=1),default_missing_value = "true",
           value_parser = clap::builder::BoolishValueParser::new(),)]
    pub try_strings: Option<bool>,

    #[arg(long, num_args(0..=1),default_missing_value = "true",
           value_parser = clap::builder::BoolishValueParser::new())]
    pub strict: Option<bool>,
}

impl SaveState {
    pub fn save(
        &self,
        state: &mut State,
        run_history: &RunHistory,
        default_runtime_settings: &RuntimeSettings,
        global_settings: &CLISettings,
    ) -> Result<()> {
        if self.no_save_state {
            // status_info!("Skipping saving state as per user request");
            return Ok(());
        }
        println!(
            "Saving state to {}..",
            global_settings.state_folder.display()
        );
        // let root_folder = root_folder.join("gammaloop_state");

        // check if the export root exists, if not create it, if it does return error
        let mut selected_root_folder = self
            .path
            .clone()
            .unwrap_or(global_settings.state_folder.clone());
        if !selected_root_folder.exists() {
            fs::create_dir_all(&selected_root_folder)?;
        } else {
            if self.strict.unwrap_or(false) {
                return Err(eyre!(
                    "Export root already exists, please choose a different path or remove the existing directory",
                ));
            }

            if !self
                .override_state
                .unwrap_or(global_settings.override_state)
            {
                while selected_root_folder.clone().exists() {
                    eprint!(
                        "Gammaloop export root {} already exists. Specify '{}' for overwriting, '{}' for not saving, or '{}' to specify where to save current state to:\n > ",
                        selected_root_folder.display().to_string().green(),
                        "o".red().bold(),
                        "n".blue().bold(),
                        "<NEW_PATH>".green().bold()
                    );
                    let mut user_input = String::new();
                    std::io::stdin()
                        .read_line(&mut user_input)
                        .expect("Could not read user-specified gammaloop state export destination");
                    //user_input = user_input.trim().into();
                    match user_input.trim() {
                        "o" => {
                            status_info!(
                                "Overwriting existing gammaloop state at {}",
                                selected_root_folder.display().to_string().green()
                            );
                            break;
                        }
                        "n" => {
                            return Ok(());
                        }
                        new_path => {
                            selected_root_folder = new_path.into();
                            continue;
                        }
                    }
                }
            }
        }

        state.save(&selected_root_folder, true, false)?;

        set_serialize_commands_as_strings(self.try_strings.unwrap_or(global_settings.try_strings));
        run_history.save_toml(&selected_root_folder, true, false)?;
        set_serialize_commands_as_strings(false);

        SHOWDEFAULTS.store(true, std::sync::atomic::Ordering::Relaxed);
        default_runtime_settings.to_file(
            &selected_root_folder.join("default_runtime_settings.toml"),
            true,
        )?;
        global_settings.to_file(&selected_root_folder.join("cli_settings.toml"), true)?;

        SHOWDEFAULTS.store(false, std::sync::atomic::Ordering::Relaxed);

        Ok(())
    }
}
