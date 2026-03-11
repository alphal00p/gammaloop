use std::{fs, path::PathBuf};

use clap::{Args, Subcommand};
use color_eyre::{eyre::eyre, owo_colors::OwoColorize, Result};
use gammalooprs::{
    processes::{
        DotExportSettings, StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings,
    },
    settings::RuntimeSettings,
    utils::serde_utils::{SmartSerde, SHOWDEFAULTS},
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tracing::{info, warn};

use crate::{
    state::{set_serialize_commands_as_strings, RunHistory, State},
    templates::Assets,
    write_schemas, CLISettings, DEFAULT_RUNTIME_SETTINGS_FILENAME, GLOBAL_SETTINGS_FILENAME,
};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Save {
    Dot {
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: Option<PathBuf>,
        #[arg(short = 'c', default_value_t = false)]
        combine_diagrams: bool,
        #[arg(short = 'n', long,num_args(0..=1), default_missing_value = "true",
               value_parser = clap::builder::BoolishValueParser::new(),)]
        with_uv: Option<bool>,
        #[arg(short = 'u', long,num_args(0..=1), default_missing_value = "true",
               value_parser = clap::builder::BoolishValueParser::new(),)]
        output_full_numerator: Option<bool>,
        #[arg(short = 'g', long,num_args(0..=1), default_missing_value = "true",
               value_parser = clap::builder::BoolishValueParser::new(),)]
        do_gamma_algebra: Option<bool>,
        #[arg(long,num_args(0..=1), default_missing_value = "true",
               value_parser = clap::builder::BoolishValueParser::new(),)]
        do_color_algebra: Option<bool>,
    },
    Standalone {
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: Option<PathBuf>,
        #[arg(long, default_value_t = false, conflicts_with = "rust")]
        python: bool,
        #[arg(long, default_value_t = false, conflicts_with = "python")]
        rust: bool,
        #[arg(long, default_value_t = false, conflicts_with = "json")]
        binary: bool,
        #[arg(long, default_value_t = false, conflicts_with = "binary")]
        json: bool,
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
            Save::Dot {
                path,
                combine_diagrams,
                output_full_numerator,
                do_color_algebra,
                do_gamma_algebra,
                with_uv,
            } => {
                // Use original default location (state folder) or custom path if provided
                let target_dir = path.unwrap_or(global_settings.state.folder.clone());
                info!("Saving dot files to {}", target_dir.display());

                // Extract embedded templates to drawings/templates relative to target directory
                if let Err(e) = Assets::extract_templates(&target_dir) {
                    warn!(
                        "Warning: Could not extract templates to drawings/templates: {}",
                        e
                    );
                }

                // Generate dynamic edge styles based on the model
                let template_path = target_dir.join("drawings/templates/edge-style.typ");
                if let Err(e) = state.model.generate_edge_style_template(&template_path) {
                    warn!("Warning: Could not generate dynamic edge styles: {}", e);
                }

                let settings = DotExportSettings {
                    do_color_algebra: do_color_algebra.unwrap_or(false),
                    do_gamma_algebra: do_gamma_algebra.unwrap_or(false),
                    output_full_numerator: output_full_numerator.unwrap_or(false),
                    split_xs_by_initial_states: true,
                    with_uv: with_uv.unwrap_or(false),
                    combine_diagrams,
                };

                // Export dot files to original location
                state.export_dots(&target_dir, &settings)?;

                // Create Justfile with draw recipe from embedded template
                if let Err(e) = Assets::extract_justfile(&target_dir) {
                    warn!(
                        "Warning: Could not create justfile at {}: {}",
                        target_dir.join("justfile").display(),
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
            Save::Standalone {
                path,
                python,
                rust,
                json,
                binary,
            } => {
                let target_dir = path.unwrap_or(global_settings.state.folder.clone());
                let mode = match (python, rust) {
                    (true, false) => StandaloneExportMode::Python,
                    (false, true) | (false, false) => StandaloneExportMode::Rust,
                    (true, true) => unreachable!("clap enforces mutual exclusivity"),
                };

                let format = match (json, binary) {
                    (true, false) => StandaloneDataFormat::Json,
                    (false, true) | (false, false) => StandaloneDataFormat::Binary,
                    (true, true) => unreachable!("clap enforces mutual exclusivity"),
                };
                let settings = StandaloneExportSettings { mode, format };
                state.process_list.export_standalone(&target_dir, &settings)
            }

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
            // info!("Skipping saving state as per user request");
            return Ok(());
        }
        // println!(
        //     "Saving state to {}..",
        //     global_settings.state_folder.display()
        // );
        // let root_folder = root_folder.join("gammaloop_state");

        // check if the export root exists, if not create it, if it does return error
        let mut selected_root_folder = self
            .path
            .clone()
            .unwrap_or(global_settings.state.folder.clone());
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
                            info!(
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
            selected_root_folder.join(DEFAULT_RUNTIME_SETTINGS_FILENAME),
            true,
        )?;
        global_settings.to_file(selected_root_folder.join(GLOBAL_SETTINGS_FILENAME), true)?;

        SHOWDEFAULTS.store(false, std::sync::atomic::Ordering::Relaxed);

        Ok(())
    }
}
