use std::{
    fs::{self, File},
    io::{self, Read, Write},
    ops::ControlFlow,
    path::{Path, PathBuf},
    sync::{Mutex, OnceLock},
    time::Instant,
};

use color_eyre::{Result, Section};
use colored::{ColoredString, Colorize};
use eyre::{eyre, Context};
use log::debug;
use schemars::{generate::SchemaSettings, schema_for, JsonSchema, Schema, SchemaGenerator};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;
use tracing::info;
use tracing_subscriber::{reload, EnvFilter, Registry};

use crate::{
    feyngen::GenerationType,
    graph::Graph,
    integrands::HasIntegrand,
    model::{Model, SerializableModel},
    processes::{ExportSettings, Process, ProcessDefinition, ProcessList},
    settings::{runtime::LockedRuntimeSettings, GlobalSettings, RuntimeSettings},
    status_debug, status_warn,
    utils::F,
    GammaLoopContextContainer,
};

use super::{tracing::FILTER_HANDLE, Cli, Commands};

#[derive(Debug, Serialize, Deserialize, Default, JsonSchema)]
pub struct RunHistory {
    pub default_runtime_settings: RuntimeSettings,
    pub global_settings: GlobalSettings,
    #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    #[schemars(with = "Vec<Commands>")]
    pub commands: Vec<Commands>,
}

impl RunHistory {
    pub fn push(&mut self, command: Commands) {
        if !matches!(&command, Commands::Quit { .. }) {
            self.commands.push(command);
        }
    }

    pub fn schema() -> Schema {
        schema_for!(RunHistory)
    }
    pub fn run(&mut self, cli: &mut Cli, state: &mut State) -> Result<ControlFlow<()>> {
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

    pub fn from_file_yaml(filename: impl AsRef<Path>) -> Result<Self> {
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

        File::create(root_folder.join("run_schema.json"))?
            .write(serde_json::to_string_pretty(&Self::schema())?.as_bytes())?;
        Ok(())
    }
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "GammaLoopState")
)]
#[derive(Clone)]
pub struct State {
    pub model: Model,
    pub process_list: ProcessList,
    pub model_path: Option<PathBuf>,
    pub save_path: PathBuf,
    pub log_filter: reload::Handle<EnvFilter, Registry>,
}

impl State {
    pub fn import_model(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.model = Model::from_file(path)?;
        Ok(())
    }

    pub fn generate_integrands(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<()> {
        self.process_list.preprocess(&self.model, global_settings)?;
        self.process_list
            .generate_integrands(&self.model, global_settings, runtime_default)?;
        Ok(())
    }

    pub fn compile_integrands(
        &mut self,
        folder: impl AsRef<Path>,
        override_existing: bool,
        global_settings: &GlobalSettings,
    ) -> Result<()> {
        self.process_list
            .compile(folder, override_existing, global_settings, &self.model)?;
        Ok(())
    }

    pub fn export_dots(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let exp_set = ExportSettings {
            root_folder: path.as_ref().to_path_buf(),
        };
        self.process_list.export_dot(&exp_set, &self.model)?;
        Ok(())
    }

    pub fn import_amplitude(&mut self, path: impl AsRef<Path>, name: Option<String>) -> Result<()> {
        let graphs = Graph::from_file(&path, &self.model)?;
        let name = name.unwrap_or(
            path.as_ref()
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .into_owned(),
        );
        let process = Process::from_graph_list(
            name,
            graphs,
            GenerationType::Amplitude,
            ProcessDefinition::new_empty(),
            None,
        )?;

        self.process_list.add_process(process);
        Ok(())
    }

    pub fn bench(
        &mut self,
        samples: usize,
        process_id: usize,
        process_name: String,
        n_cores: usize,
    ) -> Result<()> {
        let integrand = self
            .process_list
            .get_integrand_mut(process_id, process_name)?;
        let name = integrand.name();

        info!(
            "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
            name.green(),
            samples.to_string().blue()
        );

        let now = Instant::now();
        for _ in 0..samples {
            integrand.evaluate_sample(
                &Sample::Continuous(
                    F(1.),
                    (0..integrand.get_n_dim())
                        .map(|_| F(rand::random::<f64>()))
                        .collect(),
                ),
                F(1.),
                1,
                false,
                Complex::new_zero(),
            );
        }
        let total_time = now.elapsed().as_secs_f64();
        info!(
            "\n> Total time: {} s for {} samples, {} ms per sample\n",
            format!("{:.1}", total_time).blue(),
            format!("{}", samples).blue(),
            format!("{:.5}", total_time * 1000. / (samples as f64)).green(),
        );

        Ok(())
    }

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

    pub fn new_test(state_folder: PathBuf) -> Self {
        let handle = super::tracing::init_test_tracing();

        let a = Self {
            save_path: state_folder,
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
        };
        a
    }

    pub fn new_bench(state_folder: PathBuf) -> Self {
        let handle = super::tracing::init_bench_tracing();

        let a = Self {
            save_path: state_folder,
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
        };
        a
    }

    pub fn new_test_cli(&self) -> Cli {
        Cli {
            run_history: None,
            state_folder: self.save_path.clone(),
            model_file: None,
            no_save_state: true,
            override_state: false,
            command: None,
        }
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
        let mut selected_root_folder = PathBuf::from(root_folder);
        let mut user_input = String::new();
        if !root_folder.exists() {
            fs::create_dir_all(root_folder)?;
        } else {
            if strict {
                return Err(eyre!(
                    "Export root already exists, please choose a different path or remove the existing directory",
                ));
            }

            if !override_state_file {
                while selected_root_folder.exists() {
                    println!(
                        "Gammaloop export root {} already exists. Specify 'o' for overwriting, 'n' for not saving, or '<NEW_PATH>' to specify where to save current state to:",
                        selected_root_folder.display()
                    );
                    user_input.clear();
                    io::stdin()
                        .read_line(&mut user_input)
                        .expect("Could not read user-specified gammaloop state export destination");
                    //user_input = user_input.trim().into();
                    match user_input.trim() {
                        "o" => break,
                        "n" => {
                            return Ok(());
                        }
                        new_path => {
                            selected_root_folder = PathBuf::from(new_path);
                            continue;
                        }
                    }
                }
            }
        }

        let mut state_file =
        // info!("Hi");
            fs::File::create(root_folder.join("symbolica_state.bin"))?;

        symbolica::state::State::export(&mut state_file)?;
        self.process_list
            .save(&selected_root_folder, override_state_file, &self.model)?;

        // let binary = bincode::encode_to_vec(&self.integrands, bincode::config::standard())?;
        // fs::write(root_folder.join("process_list.bin"), binary)?;?
        let model_yaml = serde_yaml::to_string(&self.model.to_serializable())?;

        fs::write(selected_root_folder.join("model.yaml"), model_yaml)?;
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
