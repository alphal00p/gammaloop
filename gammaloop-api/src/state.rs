use std::{
    fs::{self},
    io::{self},
    ops::ControlFlow,
    path::{Path, PathBuf},
    sync::{atomic::AtomicBool, Mutex, OnceLock},
    time::Instant,
};

use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use schemars::{schema_for, JsonSchema, Schema};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;
use tracing::debug;
use tracing::{info, warn};
use tracing_subscriber::{reload, EnvFilter, Registry};

use gammalooprs::{
    feyngen::GenerationType,
    graph::Graph,
    integrands::HasIntegrand,
    model::{InputParamCard, Model},
    processes::{ExportSettings, Process, ProcessCollection, ProcessDefinition, ProcessList},
    settings::{runtime::LockedRuntimeSettings, GlobalSettings, RuntimeSettings},
    status_debug, status_info, status_warn,
    utils::{
        serde_utils::{get_schema_folder, SmartSerde},
        tracing::{init_bench_tracing, init_test_tracing, FILTER_HANDLE},
        F,
    },
    GammaLoopContextContainer,
};

use crate::generate::ProcessArgs;

use super::{Cli, Commands};

// Static flag to control serialization behavior
static SERIALIZE_COMMANDS_AS_STRINGS: AtomicBool = AtomicBool::new(false);

/// Set whether CommandHistory should serialize as strings when the raw_string is available
pub fn set_serialize_commands_as_strings(value: bool) {
    SERIALIZE_COMMANDS_AS_STRINGS.store(value, std::sync::atomic::Ordering::Relaxed);
}

/// Get the current setting for CommandHistory serialization behavior
pub fn get_serialize_commands_as_strings() -> bool {
    SERIALIZE_COMMANDS_AS_STRINGS.load(std::sync::atomic::Ordering::Relaxed)
}

/// Represents a command with optional raw string representation
///
/// This struct stores both the parsed command and optionally the original
/// string that was used to create it. This allows for preserving the exact
/// user input while still having access to the structured command data.
#[derive(Debug, Clone, JsonSchema, PartialEq)]
pub struct CommandHistory {
    /// The parsed command
    pub command: Commands,
    /// The original string representation of the command, if available
    pub raw_string: Option<String>,
}

impl Serialize for CommandHistory {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        if get_serialize_commands_as_strings() {
            if let Some(ref raw_string) = self.raw_string {
                raw_string.serialize(serializer)
            } else {
                self.command.serialize(serializer)
            }
        } else {
            self.command.serialize(serializer)
        }
    }
}

impl<'de> Deserialize<'de> for CommandHistory {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::Error;

        // First try to deserialize as a string
        let value = serde_yaml::Value::deserialize(deserializer)?;

        if let serde_yaml::Value::String(s) = &value {
            if let Ok(a) = Self::from_raw_string(s) {
                return Ok(a);
            }
        }

        // Fall back to deserializing as Commands directly
        match Commands::deserialize(value) {
            Ok(command) => Ok(CommandHistory {
                command,
                raw_string: None,
            }),
            Err(e) => Err(D::Error::custom(format!(
                "Failed to deserialize as both string and Commands: {}",
                e
            ))),
        }
    }
}

impl CommandHistory {
    /// Create a new CommandHistory with just a command (no raw string)
    pub fn new(command: Commands) -> Self {
        Self {
            command,
            raw_string: None,
        }
    }

    /// Create a new CommandHistory with both command and raw string
    pub fn new_with_raw(command: Commands, raw_string: String) -> Self {
        Self {
            command,
            raw_string: Some(raw_string),
        }
    }

    /// Create a CommandHistory from a command (alias for new)
    pub fn from_command(command: Commands) -> Self {
        Self::new(command)
    }

    /// Parse a raw string into a CommandHistory
    ///
    /// This function attempts to parse the raw string using clap, and if successful,
    /// creates a CommandHistory with both the parsed command and the original string.
    pub fn from_raw_string(raw_string: &str) -> Result<Self, clap::Error> {
        use crate::Cli;
        use clap::Parser;

        let args: Vec<&str> = raw_string.split_whitespace().collect();
        let cli = Cli::try_parse_from(std::iter::once("gammaloop").chain(args.iter().copied()))?;

        if let Some(command) = cli.command {
            Ok(Self::new_with_raw(command, raw_string.into()))
        } else {
            Err(clap::Error::new(clap::error::ErrorKind::MissingSubcommand))
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Default, JsonSchema, PartialEq)]
pub struct RunHistory {
    pub default_runtime_settings: RuntimeSettings,
    pub global_settings: GlobalSettings,
    #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    #[schemars(with = "Vec<CommandHistory>")]
    pub commands: Vec<CommandHistory>,
}

impl SmartSerde for RunHistory {
    fn has_schema_path(&self) -> Option<Result<PathBuf>> {
        Some(get_schema_folder().map(|f| f.join("runhistory.json")))
    }
}

impl RunHistory {
    /// Add a command to the run history
    pub fn push(&mut self, command: Commands) {
        self.push_with_raw(command, None);
    }

    /// Add a command with optional raw string to the run history
    ///
    /// If raw_string is provided, it will be stored alongside the command
    /// for potential later serialization as a string.
    pub fn push_with_raw(&mut self, command: Commands, raw_string: Option<String>) {
        if !matches!(&command, Commands::Quit { .. }) {
            self.commands.push(CommandHistory {
                command,
                raw_string,
            });
        }
    }

    pub fn schema() -> Schema {
        schema_for!(RunHistory)
    }
    pub fn run(
        &mut self,
        cli: &mut Cli,
        state: &mut State,
        global_settings: &mut GlobalSettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<()>> {
        for command_history in self.commands.clone() {
            status_info!("Running command: {:?}", command_history.command);
            if let ControlFlow::Break(_) = cli.run_command(
                command_history.command,
                state,
                self,
                global_settings,
                default_runtime_settings,
            )? {
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

    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let runhistory: Self = match Self::from_file(path.as_ref(), "run history") {
            Ok(r) => {
                status_debug!(
                    "Loaded run history from YAML file {}",
                    path.as_ref().display()
                );
                r
            }
            Err(e) => {
                status_warn!(
                    "Could not load run history at {}: {}, loading default",
                    path.as_ref().display(),
                    e
                );

                Default::default()
            }
        };

        let spec = runhistory.global_settings.debug_level.as_str();
        let _ = set_log_spec(spec)?;

        Ok(runhistory)
    }

    pub fn save_toml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        _strict: bool,
    ) -> Result<()> {
        self.to_file(root_folder.join("run.toml"), override_state_file)?;

        //Self::schema().to_file(root_folder.join("run_schema.json"))?;
        Ok(())
    }

    pub fn save_yaml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        _strict: bool,
    ) -> Result<()> {
        self.to_file(root_folder.join("run.yaml"), override_state_file)?;
        //Self::schema().to_file(root_folder.join("run_schema.json"))?;
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
    pub model_parameters: InputParamCard<F<f64>>,
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

    pub fn generate_integrand(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        process: &ProcessArgs,
    ) -> Result<()> {
        let p = &mut self.process_list.processes[process.process_id];
        if let Some(name) = &process.name {
            match &mut p.collection {
                ProcessCollection::Amplitudes(a) => {
                    if let Some(a) = a.get_mut(name) {
                        a.preprocess(&self.model, &global_settings.generation)?;
                        a.build_integrand(&self.model, global_settings, runtime_default)?;
                    } else {
                        return Err(eyre!(
                            "No amplitude named '{}' in process id {}",
                            name,
                            process.process_id
                        ));
                    }
                }
                ProcessCollection::CrossSections(a) => {
                    // a[name].preprocess(&self.model, &global_settings.generation)?;
                }
            }
        } else {
            p.preprocess(&self.model, global_settings)?;
            p.generate_integrands(&self.model, global_settings, runtime_default)?;
        }

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
                &self.model,
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

    pub fn new(save_path: PathBuf, log_file_name: Option<String>) -> Self {
        let handle = super::tracing::init_tracing("info", &save_path.join("logs"), log_file_name);

        let a = Self {
            save_path,
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
            model_parameters: InputParamCard::default(),
        };
        a
    }

    pub fn new_test(state_folder: PathBuf) -> Self {
        let handle = init_test_tracing();

        let a = Self {
            save_path: state_folder,
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
            model_parameters: InputParamCard::default(),
        };
        a
    }

    pub fn new_bench(state_folder: PathBuf) -> Self {
        let handle = init_bench_tracing();

        let a = Self {
            save_path: state_folder,
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_path: None,
            model_parameters: InputParamCard::default(),
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
            level: None,
            debug: false,
            trace_logs_filename: None,
            no_skip_default: false,
            try_strings: false,
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

    pub fn load(
        save_path: PathBuf,
        model_path: Option<PathBuf>,
        trace_logs_filename: Option<String>,
    ) -> Result<Self> {
        // let root_folder = root_folder.join("gammaloop_state");

        let mut model = if let Some(model_path) = &model_path {
            warn!("Loading model from {}", model_path.display());
            Model::from_file(model_path)?
        } else {
            let model_dir = save_path.join("model.json");
            warn!(
                "Loading model from default location: {}",
                model_dir.display()
            );
            Model::from_file(model_dir)?
        };

        debug!("Loaded model: {}", model.name);

        let input_param_card = if save_path.join("model_parameters.json").exists() {
            let a = InputParamCard::from_file(save_path.join("model_parameters.json"))?;

            model.apply_param_card(&a);
            a
        } else {
            InputParamCard::default_from_model(&model)
        };

        let state = symbolica::state::State::import(
            &mut fs::File::open(save_path.join("symbolica_state.bin"))
                .context("Trying to open symbolica state binary")?,
            None,
        )?;

        let context: GammaLoopContextContainer<'_> = GammaLoopContextContainer {
            state_map: &state,
            model: &model,
        };

        let process_list = ProcessList::load(&save_path, context)
            .context("Trying to load processList")
            .unwrap();

        let mut state = State::new(save_path, trace_logs_filename);

        state.process_list = process_list;
        state.model = model;
        state.model_parameters = input_param_card;
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
        self.model
            .to_serializable()
            .to_file(selected_root_folder.join("model.json"), override_state_file)?;
        self.model_parameters.to_file(
            selected_root_folder.join("model_parameters.json"),
            override_state_file,
        )?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use gammalooprs::{
        momentum::{Dep, ExternalMomenta, SignOrZero},
        settings::{runtime::kinematic::Externals, KinematicsSettings},
        utils::serde_utils::SHOWDEFAULTS,
    };

    use super::*;

    #[test]
    fn test_run_history() {
        use crate::state::RunHistory;
        //SHOWDEFAULTS.store(true, std::sync::atomic::Ordering::Relaxed);
        let mut run_history: RunHistory = Default::default();
        let kinematics_settings = KinematicsSettings {
            e_cm: 100.0,
            externals: Externals::Constant {
                momenta: vec![
                    ExternalMomenta::Independent([F(1.), F(2.), F(3.), F(4.)]),
                    ExternalMomenta::Dependent(Dep::Dep),
                ],
                helicities: vec![SignOrZero::Plus, SignOrZero::Minus],
            },
        };

        run_history.default_runtime_settings.kinematics = kinematics_settings;
        let toml = toml::to_string_pretty(&run_history).unwrap();
        println!("{}", toml);
        let deserialized: RunHistory = toml::from_str(&toml).unwrap();
        assert_eq!(run_history, deserialized);
        SHOWDEFAULTS.store(false, std::sync::atomic::Ordering::Relaxed);

        run_history.to_file("test_path.toml", true).unwrap();
        let deserialized_from_file = RunHistory::from_file("test_path.toml", " ").unwrap();
        assert_eq!(run_history, deserialized_from_file);
    }

    #[test]
    fn test_command_history_serialization() {
        use super::{set_serialize_commands_as_strings, CommandHistory};
        use crate::Commands;

        // Test basic construction
        let cmd_history = CommandHistory::new(Commands::Quit {});
        assert_eq!(cmd_history.raw_string, None);

        // Test with raw string
        let cmd_history_with_raw =
            CommandHistory::new_with_raw(Commands::Quit {}, "quit".to_string());
        assert_eq!(cmd_history_with_raw.raw_string, Some("quit".to_string()));

        // Test serialization as Commands (default behavior)
        set_serialize_commands_as_strings(false);
        let yaml = serde_yaml::to_string(&cmd_history).unwrap();
        let deserialized_yaml: CommandHistory = serde_yaml::from_str(&yaml).unwrap();
        let json = serde_json::to_string(&cmd_history).unwrap();
        let deserialized_json: CommandHistory = serde_json::from_str(&json).unwrap();
        let toml = toml::to_string(&cmd_history).unwrap();
        let deserialized_toml: CommandHistory = toml::from_str(&toml).unwrap();
        assert_eq!(cmd_history, deserialized_yaml);
        assert_eq!(cmd_history, deserialized_json);
        assert_eq!(cmd_history, deserialized_toml);

        // Test serialization as string
        set_serialize_commands_as_strings(true);
        let yaml_string = serde_yaml::to_string(&cmd_history_with_raw).unwrap();
        let json_string = serde_json::to_string(&cmd_history_with_raw).unwrap();
        let toml_string = toml::to_string(&cmd_history_with_raw).unwrap();
        assert!(yaml_string.contains("quit"));
        assert!(json_string.contains("quit"));
        assert!(toml_string.contains("quit"));

        // Reset flag
        set_serialize_commands_as_strings(false);
    }
}

#[derive(Args, Debug, Clone)]
pub struct ExistingArgs {
    pub process_id: u32,
    pub name: Option<String>,
}
