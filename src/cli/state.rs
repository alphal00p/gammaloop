use std::{
    collections::HashMap,
    fs::{self, File},
    path::{Path, PathBuf},
    sync::{LazyLock, Mutex},
};

use chrono::{Datelike, Local, Timelike};
use color_eyre::Result;
use colored::{ColoredString, Colorize};
use eyre::eyre;
use log::{debug, LevelFilter};

use crate::{integrands::Integrand, model::Model, new_cs::ProcessList, GammaLoopContextContainer};

use super::LogFormat;

pub struct State {
    pub model: Model,
    pub process_list: ProcessList,
    pub model_path: Option<PathBuf>,
}

impl Default for State {
    fn default() -> Self {
        let a = Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
        };
        a.setup_log().unwrap();
        a
    }
}

pub static LOG_LEVEL: LazyLock<Mutex<LevelFilter>> =
    LazyLock::new(|| Mutex::new(LevelFilter::Debug));
pub static LOG_FORMAT: LazyLock<Mutex<LogFormat>> = LazyLock::new(|| Mutex::new(LogFormat::Long));

impl State {
    pub fn load(root_folder: &Path, model_path: Option<PathBuf>) -> Result<Self> {
        // let root_folder = root_folder.join("gammaloop_state");

        let model = if let Some(model_path) = &model_path {
            debug!("Loading model from {}", model_path.display());
            Model::from_file(model_path)?
        } else {
            let model_dir = root_folder.join("model.yaml");
            debug!(
                "Loading model from default location: {}",
                model_dir.display()
            );

            Model::from_serializable_model(serde_yaml::from_reader(File::open(model_dir)?)?)
        };

        let state = symbolica::state::State::import(
            &mut fs::File::open(root_folder.join("symbolica_state.bin"))?,
            None,
        )?;

        let context = GammaLoopContextContainer {
            state_map: &state,
            model: &model,
        };

        let process_list_data = fs::read(root_folder.join("process_list.bin"))?;

        let (process_list, _) = bincode::decode_from_slice_with_context(
            &process_list_data,
            bincode::config::standard(),
            context,
        )?;

        Ok(State {
            model,
            model_path,
            process_list,
        })
    }

    pub fn save(&self, root_folder: &Path, override_state_file: bool, strict: bool) -> Result<()> {
        // let root_folder = root_folder.join("gammaloop_state");

        // check if the export root exists, if not create it, if it does return error
        if !root_folder.exists() {
            fs::create_dir_all(root_folder)?;
        } else {
            if strict {
                return Err(eyre!(
                    "Export root already exists, please choose a different path or remove the existing directory",
                ));
            }

            if !override_state_file {
                return Err(eyre!(
                    "Export root already exists, please choose a different path or remove the existing directory",
                ));
            }
        }

        let mut state_file =
        // info!("Hi");
            fs::File::create(root_folder.join("symbolica_state.bin"))?;

        symbolica::state::State::export(&mut state_file)?;

        let binary = bincode::encode_to_vec(&self.process_list, bincode::config::standard())?;
        fs::write(root_folder.join("process_list.bin"), binary)?;

        // let binary = bincode::encode_to_vec(&self.integrands, bincode::config::standard())?;
        // fs::write(root_folder.join("process_list.bin"), binary)?;?
        let model_yaml = serde_yaml::to_string(&self.model.to_serializable())?;

        fs::write(root_folder.join("model.yaml"), model_yaml)?;
        Ok(())
    }

    pub(crate) fn setup_log(&self) -> Result<()> {
        fern::Dispatch::new()
            .filter(|metadata| metadata.level() <= (*LOG_LEVEL.lock().unwrap()))
            // Perform allocation-free log formatting
            .format(|out, message, record| {
                let now = Local::now();
                match *LOG_FORMAT.lock().unwrap() {
                    LogFormat::Long => out.finish(format_args!(
                        "[{}] @{} {}: {}",
                        format!(
                            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:03}",
                            now.year(),
                            now.month(),
                            now.day(),
                            now.hour(),
                            now.minute(),
                            now.second(),
                            now.timestamp_subsec_millis()
                        )
                        .bright_green(),
                        format_target(record.target().into(), record.level()),
                        format_level(record.level()),
                        message
                    )),
                    LogFormat::Short => out.finish(format_args!(
                        "[{}] {}: {}",
                        format!("{:02}:{:02}:{:02}", now.hour(), now.minute(), now.second())
                            .bright_green(),
                        format_level(record.level()),
                        message
                    )),
                    LogFormat::Min => out.finish(format_args!(
                        "{}: {}",
                        format_level(record.level()),
                        message
                    )),
                    LogFormat::None => out.finish(format_args!("{}", message)),
                }
            })
            // Output to stdout, files, and other Dispatch configurations
            .chain(std::io::stdout())
            .chain(fern::log_file("gammaloop_rust_output.log")?)
            // Apply globally
            .apply()?;

        Ok(())
    }
}

pub(crate) fn format_level(level: log::Level) -> ColoredString {
    match level {
        log::Level::Error => format!("{:<8}", "ERROR").red(),
        log::Level::Warn => format!("{:<8}", "WARNING").yellow(),
        log::Level::Info => format!("{:<8}", "INFO").into(),
        log::Level::Debug => format!("{:<8}", "DEBUG").bright_black(),
        log::Level::Trace => format!("{:<8}", "TRACE").into(),
    }
}

pub(crate) fn format_target(target: String, level: log::Level) -> ColoredString {
    let split_targets = target.split("::").collect::<Vec<_>>();
    //[-2..].iter().join("::");
    let start = split_targets.len().saturating_sub(2);
    let mut shortened_path = split_targets[start..].join("::");
    if level < log::Level::Debug && shortened_path.len() > 20 {
        shortened_path = format!("{}...", shortened_path.chars().take(17).collect::<String>());
    }
    format!("{:<20}", shortened_path).bright_blue()
}
