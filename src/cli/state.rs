use std::{
    fmt,
    fs::{self, File},
    io::{Read, Write},
    ops::ControlFlow,
    path::{Path, PathBuf},
    sync::{LazyLock, Mutex, OnceLock},
};

use color_eyre::{Result, Section};
use colored::{ColoredString, Colorize};
use eyre::{eyre, Context};
use log::{debug, warn, LevelFilter};
use serde::{Deserialize, Serialize};
use tracing_appender::rolling;
use tracing_subscriber::{
    layer::SubscriberExt, reload, util::SubscriberInitExt, EnvFilter, Layer, Registry,
};

use crate::{
    model::{Model, SerializableModel},
    processes::ProcessList,
    settings::{GlobalSettings, RuntimeSettings},
    status_debug, status_warn, GammaLoopContextContainer,
};

use super::{tracing::FILTER_HANDLE, Cli, Commands, LogFormat};

#[derive(Debug, Serialize, Deserialize, Default)]
pub struct RunHistory {
    pub default_runtime_settings: RuntimeSettings,
    pub global_settings: GlobalSettings,
    #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    pub commands: Vec<Commands>,
}

impl RunHistory {
    pub fn push(&mut self, command: Commands) {
        if !matches!(&command, Commands::Quit { .. }) {
            self.commands.push(command);
        }
    }
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
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let runhistory: Self = match Self::from_file_yaml(path.as_ref()) {
            Ok(r) => {
                status_debug!(
                    "Loaded run history from YAML file {}",
                    path.as_ref().display()
                );
                r
            }
            Err(e) => {
                status_warn!(
                    "Could not load run history from YAML at {}: {}, loading default",
                    path.as_ref().display(),
                    e
                );

                Default::default()
            }
        };

        let spec = runhistory.global_settings.debug_level.to_env_spec();
        let _ = set_log_spec(spec)?;

        Ok(runhistory)
    }

    pub fn save_toml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        strict: bool,
    ) -> Result<()> {
        File::create(root_folder.join("run.toml"))?
            .write(toml::to_string_pretty(self)?.as_bytes())?;
        Ok(())
    }

    pub fn save_yaml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        strict: bool,
    ) -> Result<()> {
        File::create(root_folder.join("run.yaml"))?
            .write(serde_yaml::to_string(self)?.as_bytes())?;
        Ok(())
    }
}

pub struct State {
    pub model: Model,
    pub process_list: ProcessList,
    pub model_path: Option<PathBuf>,
    save_path: PathBuf,
    log_filter: reload::Handle<EnvFilter, Registry>,
}

impl State {
    pub fn new(save_path: PathBuf) -> Self {
        let handle =
            super::tracing::init_tracing("info,symbolica::poly::gcd=off", &save_path.join("logs"));

        let a = Self {
            save_path,
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
        };
        a
    }
}

// just so you can show it in the banner etc.
pub(super) static LOG_SPEC: OnceLock<Mutex<String>> = OnceLock::new();

pub fn set_log_spec(spec: &str) -> color_eyre::Result<()> {
    let handle = FILTER_HANDLE.get().expect("tracing not initialized");
    handle.modify(|f| *f = EnvFilter::new(spec))?;
    if let Some(s) = LOG_SPEC.get() {
        *s.lock().unwrap() = spec.to_string();
    }
    Ok(())
}

pub fn current_log_spec() -> String {
    LOG_SPEC
        .get()
        .map(|s| s.lock().unwrap().clone())
        .unwrap_or_else(|| "info".into())
}

impl State {
    pub fn set_log_spec(&self, spec: &str) -> Result<()> {
        set_log_spec(spec)
    }

    pub fn load(save_path: PathBuf, model_path: Option<PathBuf>) -> Result<Self> {
        // let root_folder = root_folder.join("gammaloop_state");

        let model = if let Some(model_path) = &model_path {
            debug!("Loading model from {}", model_path.display());
            Model::from_file(model_path)?
        } else {
            let model_dir = save_path.join("model.yaml");
            debug!(
                "Loading model from default location: {}",
                model_dir.display()
            );

            Model::from_serializable_model(SerializableModel::from_file(model_dir)?)
        };

        debug!("Loaded model: {}", model.name);

        let state = symbolica::state::State::import(
            &mut fs::File::open(save_path.join("symbolica_state.bin"))
                .context("Trying to open symbolica state binary")?,
            None,
        )?;

        let context = GammaLoopContextContainer {
            state_map: &state,
            model: &model,
        };

        let process_list = ProcessList::load(&save_path, context)
            .context("Trying to load processList")
            .unwrap();

        let mut state = State::new(save_path);

        state.process_list = process_list;
        state.model = model;
        state.model_path = model_path;
        Ok(state)
    }

    pub fn compile(
        &mut self,
        root_folder: &Path,
        override_compiled: bool,
        settings: &GlobalSettings,
    ) -> Result<()> {
        fs::create_dir_all(root_folder)?;

        self.process_list
            .compile(root_folder, override_compiled, settings, &self.model)?;
        Ok(())
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

pub(crate) fn format_level(level: tracing::Level) -> ColoredString {
    match level {
        tracing::Level::ERROR => format!("{:<8}", "ERROR").red(),
        tracing::Level::WARN => format!("{:<8}", "WARNING").yellow(),
        tracing::Level::INFO => format!("{:<8}", "INFO").into(),
        tracing::Level::DEBUG => format!("{:<8}", "DEBUG").bright_black(),
        tracing::Level::TRACE => format!("{:<8}", "TRACE").into(),
    }
}

pub(crate) fn format_target(target: String, level: tracing::Level) -> ColoredString {
    // println!("Target: {}", target);
    let split_targets = target.split("::").collect::<Vec<_>>();
    //[-2..].iter().join("::");
    let start = split_targets.len().saturating_sub(2);
    let mut shortened_path = split_targets[start..].join("::");
    if level < tracing::Level::DEBUG && shortened_path.len() > 20 {
        shortened_path = format!("{}...", shortened_path.chars().take(17).collect::<String>());
    }
    format!("{:<20}", shortened_path).bright_blue()
}
