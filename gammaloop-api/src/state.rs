use std::{
    fs::{self},
    io::{self},
    ops::ControlFlow,
    path::{Path, PathBuf},
    sync::atomic::AtomicBool,
    time::Instant,
};

use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::{
    processes::{Amplitude, CrossSection},
    utils::serde_utils::IsDefault,
};
use linnet::half_edge::subgraph::SubGraph;
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
    initialisation::{initialise, initialise_with_settings},
    integrands::HasIntegrand,
    model::{InputParamCard, Model},
    processes::{ExportSettings, Process, ProcessCollection, ProcessDefinition, ProcessList},
    settings::{runtime::LockedRuntimeSettings, GlobalSettings, RuntimeSettings},
    status_debug, status_info,
    utils::{
        serde_utils::{get_schema_folder, SmartSerde},
        tracing::{init_bench_tracing, init_test_tracing},
        F,
    },
    GammaLoopContextContainer,
};

use crate::{
    commands::{save::SaveState, Commands},
    tracing::{set_file_log_filter, set_stderr_log_filter},
    CLISettings,
};

pub trait SyncSettings {
    fn sync_settings(&self) -> Result<()>;
}

impl SyncSettings for CLISettings {
    fn sync_settings(&self) -> Result<()> {
        set_file_log_filter(&self.global.logfile_directive)?;
        set_stderr_log_filter(&self.global.display_directive)?;
        initialise_with_settings(Some(&self.global))?;
        Ok(())
    }
}

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
        use serde::de::{self, Visitor};
        use std::fmt;

        struct CommandHistoryVisitor;

        impl<'de> Visitor<'de> for CommandHistoryVisitor {
            type Value = CommandHistory;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a string or a Commands structure")
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                CommandHistory::from_raw_string(value)
                    .map_err(|_| E::custom(format!("Failed to parse command string: '{}'", value)))
            }

            fn visit_string<E>(self, value: String) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                self.visit_str(&value)
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: de::SeqAccess<'de>,
            {
                // Handle TOML array format for enums like [Quit]
                let command =
                    Commands::deserialize(de::value::SeqAccessDeserializer::new(&mut seq))?;
                Ok(CommandHistory {
                    command,
                    raw_string: None,
                })
            }

            fn visit_map<A>(self, map: A) -> Result<Self::Value, A::Error>
            where
                A: de::MapAccess<'de>,
            {
                // Handle map format for enums
                let command = Commands::deserialize(de::value::MapAccessDeserializer::new(map))?;
                Ok(CommandHistory {
                    command,
                    raw_string: None,
                })
            }
        }

        deserializer.deserialize_any(CommandHistoryVisitor)
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
        use crate::Repl;
        use clap::Parser;

        let args: Vec<&str> = raw_string.split_whitespace().collect();
        let cli = Repl::try_parse_from(std::iter::once("gammaloop").chain(args.iter().copied()))?;

        Ok(Self::new_with_raw(cli.command, raw_string.into()))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, JsonSchema, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct RunHistory {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub default_runtime_settings: RuntimeSettings,

    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub cli_settings: CLISettings,
    // #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    // #[schemars(with = "Vec<CommandHistory>")]
    pub commands: Vec<CommandHistory>,
}

impl SmartSerde for RunHistory {
    fn has_schema_path(&self, online: bool) -> Option<Result<PathBuf>> {
        Some(get_schema_folder(online).map(|f| f.join("runhistory.json")))
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
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<SaveState>> {
        for command_history in self.commands.clone() {
            status_info!("Running command: {:?}", command_history.command);
            if let ControlFlow::Break(a) = command_history.command.run(
                state,
                self,
                global_settings,
                default_runtime_settings,
            )? {
                return Ok(ControlFlow::Break(a));
            }
        }
        Ok(ControlFlow::Continue(()))
    }
    pub fn merge(&mut self, other: Self) {
        self.commands.extend(other.commands);
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        status_debug!("Loaded run history from file {}", path.as_ref().display());

        let runhistory = Self::from_file(path.as_ref(), "run history")?;
        runhistory.cli_settings.sync_settings()?;

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
        let generation_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(global_settings.n_cores.generate)
            .build()?;

        self.process_list
            .preprocess(&self.model, global_settings, &generation_pool)?;
        self.process_list.generate_integrands(
            &self.model,
            global_settings,
            runtime_default,
            &generation_pool,
        )?;
        Ok(())
    }

    pub fn generate_integrand(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        process_id: usize,
        integrand_name: Option<String>,
    ) -> Result<()> {
        let generation_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(global_settings.n_cores.generate)
            .build()?;

        let p = &mut self.process_list.processes[process_id];
        if let Some(name) = &integrand_name {
            match &mut p.collection {
                ProcessCollection::Amplitudes(a) => {
                    if let Some(a) = a.get_mut(name) {
                        a.preprocess(&self.model, &global_settings.generation, &generation_pool)?;
                        a.build_integrand(
                            &self.model,
                            global_settings,
                            runtime_default,
                            &generation_pool,
                        )?;
                    } else {
                        return Err(eyre!(
                            "No amplitude named '{}' in process id {}",
                            name,
                            process_id
                        ));
                    }
                }
                ProcessCollection::CrossSections(cs) => {
                    if let Some(cs) = cs.get_mut(name) {
                        cs.preprocess(
                            &self.model,
                            &p.definition,
                            &global_settings.generation,
                            &generation_pool,
                        )?;
                        cs.build_integrand(
                            &self.model,
                            global_settings,
                            runtime_default,
                            &generation_pool,
                        )?;
                    } else {
                        return Err(eyre!(
                            "No cross section named '{}' in process id {}",
                            name,
                            process_id
                        ));
                    }
                }
            }
        } else {
            p.preprocess(&self.model, global_settings, &generation_pool)?;
            p.generate_integrands(
                &self.model,
                global_settings,
                runtime_default,
                &generation_pool,
            )?;
        }

        Ok(())
    }

    pub fn compile_integrands(
        &mut self,
        folder: impl AsRef<Path>,
        override_existing: bool,
        global_settings: &GlobalSettings,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> Result<()> {
        let compile_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(global_settings.n_cores.compile)
            .build()?;
        self.process_list.compile(
            folder,
            override_existing,
            global_settings,
            process_id,
            integrand_name,
            &compile_pool,
        )?;
        Ok(())
    }

    pub fn export_dots(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let exp_set = ExportSettings {
            root_folder: path.as_ref().to_path_buf(),
        };
        self.process_list.export_dot(&exp_set)?;
        Ok(())
    }

    pub fn import_graphs(
        &mut self,
        graphs: Vec<Graph>,
        process_name: Option<String>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        overwrite: bool,
        append: bool,
    ) -> Result<()> {
        let generation_type = if graphs.iter().all(|g| g.initial_state_cut.nedges(g) == 0) {
            GenerationType::Amplitude
        } else if graphs.iter().all(|g| g.initial_state_cut.nedges(g) > 0) {
            GenerationType::CrossSection
        } else {
            return Err(eyre!(
                "Mix of amplitude and cross section graphs in the same file is not supported"
            ));
        };

        let integrand_base_name = integrand_name.clone().unwrap_or("default".to_string());
        let process = if let Some(proc_id) = process_id {
            if proc_id >= self.process_list.processes.len() {
                return Err(eyre!(
                    "Process ID {} invalid, only {} processes available",
                    proc_id,
                    self.process_list.processes.len()
                ));
            }
            Some(&mut self.process_list.processes[proc_id])
        } else {
            let p_name = match process_name {
                Some(n) => n,
                None => {
                    return Err(eyre!(
                        "Either process ID or process name must be provided when importing graphs"
                    ))
                }
            };
            if let Some(existing_proc) = self
                .process_list
                .processes
                .iter_mut()
                .find(|p| p.definition.folder_name == p_name)
            {
                Some(existing_proc)
            } else {
                let process_defintion =
                    ProcessDefinition::from_graph_list(&graphs, generation_type, &self.model)?;
                let process = Process::from_graph_list(
                    p_name,
                    integrand_base_name.clone(),
                    // TODO: avoid clone here
                    graphs.clone(),
                    generation_type,
                    Some(process_defintion),
                    None,
                    &self.model,
                )?;

                self.process_list.add_process(process);
                None
            }
        };
        if let Some(p) = process {
            let existing_names = p.get_integrand_names();
            let integrand_name = if existing_names.contains(&integrand_base_name.as_str()) {
                if append {
                    let mut integrand_i = 0;
                    while existing_names
                        .iter()
                        .any(|ce| *ce == format!("{}_{}", integrand_base_name, integrand_i))
                    {
                        integrand_i += 1;
                    }
                    format!("{}_{}", integrand_base_name, integrand_i)
                } else if overwrite {
                    p.collection.remove_integrand(&integrand_base_name)?;
                    integrand_base_name.clone()
                } else {
                    return Err(eyre!(
                        "Integrand name '{}' already exists in process '{}', use either --overwrite or --append flag when loading graphs",
                        integrand_base_name,
                        p.definition.folder_name
                    ));
                }
            } else {
                integrand_base_name.clone()
            };

            match generation_type {
                GenerationType::Amplitude => p
                    .collection
                    .add_amplitude(Amplitude::from_graph_list(integrand_name.clone(), graphs)?),
                GenerationType::CrossSection => {
                    p.collection
                        .add_cross_section(CrossSection::from_graph_list(
                            integrand_name.clone(),
                            graphs,
                            &self.model,
                        )?)
                }
            }
        }

        Ok(())
    }

    pub fn bench(
        &mut self,
        samples: usize,
        process_id: usize,
        process_name: String,
        _n_cores: usize,
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

    pub fn new(log_dir: impl AsRef<Path>, log_file_name: Option<String>) -> Self {
        let _ = initialise();
        let handle = super::tracing::init_tracing(&log_dir.as_ref().join("logs"), log_file_name);

        let a = Self {
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
        };
        a
    }

    pub fn new_test() -> Self {
        let handle = init_test_tracing();

        let a = Self {
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
        };
        a
    }

    pub fn new_bench() -> Self {
        let handle = init_bench_tracing();

        let a = Self {
            log_filter: handle,
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
        };
        a
    }
}

impl State {
    pub fn load(
        save_path: PathBuf,
        model_path: Option<PathBuf>,
        trace_logs_filename: Option<String>,
    ) -> Result<Self> {
        // let root_folder = root_folder.join("gammaloop_state");

        let mut model = if let Some(model_path) = &model_path {
            info!("Loading model from {}", model_path.display());
            Model::from_file(model_path)?
        } else {
            let model_dir = save_path.join("model.json");
            info!(
                "Loading model from default location: {}",
                model_dir.display()
            );
            Model::from_file(model_dir)?
        };

        debug!("Loaded model: {}", model.name);

        let input_param_card = if save_path.join("model_parameters.json").exists() {
            let a = InputParamCard::from_file(save_path.join("model_parameters.json"))?;

            let _ = model.apply_param_card(&a);
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
        Ok(state)
    }

    pub fn compile(
        &mut self,
        root_folder: &Path,
        override_compiled: bool,
        settings: &GlobalSettings,
    ) -> Result<()> {
        fs::create_dir_all(root_folder)?;

        let compile_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(settings.n_cores.compile)
            .build()?;
        self.process_list.compile(
            root_folder,
            override_compiled,
            settings,
            None,
            None,
            &compile_pool,
        )?;
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
            .save(&selected_root_folder, override_state_file)?;

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
        improve_ps::PhaseSpaceImprovementSettings,
        momentum::{Dep, ExternalMomenta, SignOrZero},
        settings::{runtime::kinematic::Externals, KinematicsSettings},
        utils::serde_utils::SHOWDEFAULTS,
    };

    use crate::commands::{
        save::SaveState,
        set::{Set, SetArgs},
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
                improvement_settings: PhaseSpaceImprovementSettings::default(),
                helicities: vec![SignOrZero::Plus, SignOrZero::Minus],
                f_64_cache: None,
                f_128_cache: None,
            },
        };

        run_history.push(Commands::Set(Set::Global {
            input: SetArgs::Stored,
        }));

        run_history.default_runtime_settings.kinematics = kinematics_settings;
        set_serialize_commands_as_strings(true);
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
        use crate::commands::Commands;

        // Test basic construction
        let cmd_history = CommandHistory::new(Commands::Quit(SaveState::default()));
        assert_eq!(cmd_history.raw_string, None);

        // Test with raw string
        let cmd_history_with_raw =
            CommandHistory::new_with_raw(Commands::Quit(SaveState::default()), "quit".to_string());
        assert_eq!(cmd_history_with_raw.raw_string, Some("quit".to_string()));

        // Test serialization as Commands (default behavior)
        set_serialize_commands_as_strings(false);

        let json = serde_json::to_string(&cmd_history).unwrap();
        let deserialized_json: CommandHistory = serde_json::from_str(&json).unwrap();
        let toml = toml::to_string(&cmd_history).unwrap();
        let deserialized_toml: CommandHistory = toml::from_str(&toml).unwrap();
        assert_eq!(cmd_history, deserialized_json);
        assert_eq!(cmd_history, deserialized_toml);

        // Test serialization as string
        set_serialize_commands_as_strings(true);

        let json_string = serde_json::to_string_pretty(&cmd_history_with_raw).unwrap();
        assert!(json_string.contains("quit"));

        // Reset flag
        set_serialize_commands_as_strings(false);
    }

    #[test]
    fn test_command_history_toml_and_json_formats() {
        use super::{set_serialize_commands_as_strings, CommandHistory};
        use crate::commands::Commands;

        // Test different command types
        let quit_cmd = CommandHistory::new(Commands::Quit(SaveState::default()));
        let quit_with_raw =
            CommandHistory::new_with_raw(Commands::Quit(SaveState::default()), "quit".to_string());

        // Test JSON serialization/deserialization
        {
            // Test Commands format in JSON
            set_serialize_commands_as_strings(false);
            let json = serde_json::to_string_pretty(&quit_cmd).unwrap();
            let deserialized: CommandHistory = serde_json::from_str(&json).unwrap();
            assert_eq!(quit_cmd, deserialized);

            // Test string format in JSON
            set_serialize_commands_as_strings(true);
            let json_string = serde_json::to_string_pretty(&quit_with_raw).unwrap();
            assert!(json_string.contains("quit"));
            let deserialized_string: CommandHistory = serde_json::from_str(&json_string).unwrap();
            assert_eq!(quit_with_raw, deserialized_string);
        }

        // Test TOML serialization/deserialization with wrapper struct
        {
            #[derive(serde::Serialize, serde::Deserialize)]
            struct CommandWrapper {
                command: CommandHistory,
            }

            // Test Commands format in TOML
            set_serialize_commands_as_strings(false);
            let wrapper = CommandWrapper {
                command: quit_cmd.clone(),
            };
            let toml = toml::to_string_pretty(&wrapper).unwrap();
            let deserialized_wrapper: CommandWrapper = toml::from_str(&toml).unwrap();
            assert_eq!(quit_cmd, deserialized_wrapper.command);

            // Test string format in TOML
            set_serialize_commands_as_strings(true);
            let wrapper_string = CommandWrapper {
                command: quit_with_raw.clone(),
            };
            let toml_string = toml::to_string_pretty(&wrapper_string).unwrap();
            assert!(toml_string.contains("quit"));
            let deserialized_string_wrapper: CommandWrapper = toml::from_str(&toml_string).unwrap();
            assert_eq!(quit_with_raw, deserialized_string_wrapper.command);
        }

        // Test cross-format compatibility: serialize in one format, deserialize in another
        {
            set_serialize_commands_as_strings(false);

            // Serialize as JSON, deserialize the Commands directly from JSON Value
            let json = serde_json::to_string(&quit_cmd).unwrap();
            let json_value: serde_json::Value = serde_json::from_str(&json).unwrap();
            let from_json: CommandHistory = serde_json::from_value(json_value).unwrap();
            assert_eq!(quit_cmd, from_json);
        }
        // Reset flag
        set_serialize_commands_as_strings(false);
    }
}

#[derive(Args, Debug, Clone)]
pub struct ExistingArgs {
    pub process_id: u32,
    pub name: Option<String>,
}
