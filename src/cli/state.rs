use std::{
    fs::{self, File},
    io::{Read, Write},
    ops::ControlFlow,
    path::{Path, PathBuf},
    sync::{LazyLock, Mutex},
};

use color_eyre::{Result, Section};
use colored::{ColoredString, Colorize};
use eyre::{eyre, Context};
use log::{debug, LevelFilter};
use serde::{Deserialize, Serialize};

use crate::{
    model::{Model, SerializableModel},
    processes::ProcessList,
    settings::{GlobalSettings, RuntimeSettings},
    GammaLoopContextContainer,
};

use super::{Cli, Commands, LogFormat};

#[derive(Debug, Serialize, Deserialize, Default)]
pub struct RunHistory {
    pub default_runtime_settings: RuntimeSettings,
    pub global_settings: GlobalSettings,
    #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    pub commands: Vec<Commands>,
}

impl RunHistory {
    pub(crate) fn run(&mut self, cli: &mut Cli, state: &mut State) -> Result<ControlFlow<()>> {
        for command in self.commands.clone() {
            if let ControlFlow::Break(_) = cli.run_command(command, self, state)? {
                return Ok(ControlFlow::Break(()));
            }
        }
        Ok(ControlFlow::Continue(()))
    }
    pub fn merge(&mut self, other: Self) {
        self.commands.extend(other.commands);
        self.default_runtime_settings = other.default_runtime_settings;
        self.global_settings = other.global_settings;
    }

    pub(crate) fn from_file_toml(filename: impl AsRef<Path>) -> Result<Self> {
        let filename = filename.as_ref();
        let mut f = File::open(filename)
            .wrap_err_with(|| format!("Could not open run history file {}", filename.display()))
            .suggestion("Does the path exist?")?;
        let mut buf = String::new();
        f.read_to_string(&mut buf)
            .wrap_err("Could not load string history file")?;

        Ok(toml::from_str(&buf)?)
    }

    pub(crate) fn from_file_yaml(filename: impl AsRef<Path>) -> Result<Self> {
        let filename = filename.as_ref();
        let mut f = File::open(filename)
            .wrap_err_with(|| format!("Could not open run history file {}", filename.display()))
            .suggestion("Does the path exist?")?;
        let mut buf = String::new();
        f.read_to_string(&mut buf)
            .wrap_err("Could not load string history file")?;

        Ok(serde_yaml::from_str(&buf)?)
    }
    pub fn new(
        state_folder: impl AsRef<Path>,
        run_settings_path: Option<&Path>,
        global_settings_path: Option<&Path>,
    ) -> Result<Self> {
        // Load settings from YAML first so CLI flags can override.
        let runtime_settings: Option<RuntimeSettings> = run_settings_path
            .map(RuntimeSettings::from_file)
            .transpose()
            .with_context(|| "Error trying to read runtime settings from file:")?;

        let global_settings: Option<GlobalSettings> = global_settings_path
            .map(GlobalSettings::from_file)
            .transpose()
            .with_context(|| "Error trying to read global settings from file:")?;

        let mut runhistory: Self =
            Self::from_file_yaml(state_folder.as_ref().join("run_history.yaml"))
                .unwrap_or_default();

        if let Some(runtime_settings) = runtime_settings {
            debug!("Overriding runtime settings from file");
            runhistory.default_runtime_settings = runtime_settings;
        }

        if let Some(global_settings) = global_settings {
            debug!("Overriding global settings from file");
            runhistory.global_settings = global_settings;
        }

        *LOG_LEVEL.lock().unwrap() = match runhistory.default_runtime_settings.general.debug {
            0 => LevelFilter::Off,
            1 => LevelFilter::Error,
            2 => LevelFilter::Warn,
            3 => LevelFilter::Info,
            4 => LevelFilter::Debug,
            _ => LevelFilter::Trace,
        };

        println!("LOG_LEVEL: {:?}", LOG_LEVEL.lock().unwrap());

        Ok(runhistory)
    }

    pub fn save_toml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        strict: bool,
    ) -> Result<()> {
        File::create(root_folder.join("run_history.toml"))?
            .write(toml::to_string_pretty(self)?.as_bytes())?;
        Ok(())
    }

    pub fn save_yaml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        strict: bool,
    ) -> Result<()> {
        File::create(root_folder.join("run_history.yaml"))?
            .write(serde_yaml::to_string(self)?.as_bytes())?;
        Ok(())
    }
}

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

            Model::from_serializable_model(SerializableModel::from_file(model_dir)?)
        };

        debug!("Loaded model: {}", model.name);

        let state = symbolica::state::State::import(
            &mut fs::File::open(root_folder.join("symbolica_state.bin"))
                .context("Trying to open symbolica state binary")?,
            None,
        )?;

        let context = GammaLoopContextContainer {
            state_map: &state,
            model: &model,
        };

        let process_list = ProcessList::load(root_folder, context)
            .context("Trying to load processList")
            .unwrap();

        Ok(State {
            model,
            model_path,
            process_list,
        })
    }

    pub fn save(
        &mut self,
        root_folder: &Path,
        override_state_file: bool,
        strict: bool,
    ) -> Result<()> {
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
        self.process_list
            .save(root_folder, override_state_file, &self.model)?;

        // let binary = bincode::encode_to_vec(&self.integrands, bincode::config::standard())?;
        // fs::write(root_folder.join("process_list.bin"), binary)?;?
        let model_yaml = serde_yaml::to_string(&self.model.to_serializable())?;

        fs::write(root_folder.join("model.yaml"), model_yaml)?;
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
