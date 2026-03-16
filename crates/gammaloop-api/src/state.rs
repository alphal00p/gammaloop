use std::{
    collections::HashSet,
    fs::{self},
    io::{self},
    ops::ControlFlow,
    path::{Path, PathBuf},
    str::FromStr,
    sync::atomic::AtomicBool,
    time::Instant,
};

use clap::Args;
use color_eyre::{Result, Section};
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::{
    processes::{Amplitude, CrossSection},
    utils::serde_utils::IsDefault,
};
use linnet::half_edge::subgraph::SubGraphLike;
use schemars::{schema_for, JsonSchema, Schema};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;
use toml::Value as TomlValue;
use tracing::debug;
use tracing::info;

use gammalooprs::{
    feyngen::GenerationType,
    graph::Graph,
    initialisation::initialise,
    integrands::HasIntegrand,
    model::{InputParamCard, Model},
    processes::{DotExportSettings, Process, ProcessCollection, ProcessDefinition, ProcessList},
    settings::{runtime::LockedRuntimeSettings, GlobalSettings, RuntimeSettings},
    utils::{
        serde_utils::{get_schema_folder, SmartSerde},
        tracing::{init_bench_tracing, init_test_tracing},
        F,
    },
    GammaLoopContextContainer,
};

use crate::{
    command_parser::split_command_line,
    commands::{save::SaveState, Commands},
    tracing::{set_file_log_filter, set_log_style, set_stderr_log_filter},
    CLISettings,
};

#[derive(Debug, Clone, PartialEq, Eq, JsonSchema)]
pub enum ProcessRef {
    Id(usize),
    Name(String),
    Unqualified(String),
}

impl Serialize for ProcessRef {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match self {
            ProcessRef::Id(id) => serializer.serialize_u64(*id as u64),
            ProcessRef::Name(name) => serializer.serialize_str(&format!("name:{name}")),
            ProcessRef::Unqualified(value) => serializer.serialize_str(value),
        }
    }
}

impl<'de> Deserialize<'de> for ProcessRef {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, Visitor};
        use std::fmt;

        struct ProcessRefVisitor;

        impl<'de> Visitor<'de> for ProcessRefVisitor {
            type Value = ProcessRef;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a process reference string or numeric id")
            }

            fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Ok(ProcessRef::Id(value as usize))
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                ProcessRef::from_str(value).map_err(E::custom)
            }

            fn visit_string<E>(self, value: String) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                self.visit_str(&value)
            }
        }

        deserializer.deserialize_any(ProcessRefVisitor)
    }
}

impl FromStr for ProcessRef {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        if let Some(rest) = value.strip_prefix('#') {
            let id = rest
                .parse::<usize>()
                .map_err(|_| format!("Invalid process id in '{value}'"))?;
            return Ok(ProcessRef::Id(id));
        }
        if let Some(rest) = value.strip_prefix("name:") {
            if rest.is_empty() {
                return Err("Process name cannot be empty".to_string());
            }
            return Ok(ProcessRef::Name(rest.to_string()));
        }
        if value.is_empty() {
            return Err("Process reference cannot be empty".to_string());
        }
        Ok(ProcessRef::Unqualified(value.to_string()))
    }
}

#[test]
fn try_complicated() {
    "GGHHH3loop_no_iterative_optimization_3L"
        .parse::<ProcessRef>()
        .unwrap();
    // ProcessRef::
}

impl ProcessRef {
    pub fn resolve(&self, process_list: &ProcessList) -> Result<usize> {
        let process_count = process_list.processes.len();
        match self {
            ProcessRef::Id(id) => {
                if *id >= process_count {
                    return Err(eyre!(
                        "Process ID {} invalid, only {} processes available",
                        id,
                        process_count
                    ));
                }
                Ok(*id)
            }
            ProcessRef::Name(name) => process_list
                .processes
                .iter()
                .position(|p| p.definition.folder_name == *name)
                .ok_or_else(|| {
                    eyre!(
                        "No process named '{}'. Use 'display processes' to list available processes.",
                        name
                    )
                }),
            ProcessRef::Unqualified(value) => {
                let name_match = process_list
                    .processes
                    .iter()
                    .position(|p| p.definition.folder_name == *value);
                if let Ok(id) = value.parse::<usize>() {
                    let id_valid = id < process_count;
                    match (id_valid, name_match) {
                        (true, Some(_)) => Err(eyre!(
                            "Ambiguous process reference '{}'. Use '#{}' or 'name:{}' to disambiguate.",
                            value,
                            id,
                            value
                        )),
                        (true, None) => Ok(id),
                        (false, Some(index)) => Ok(index),
                        (false, None) => Err(eyre!(
                            "No process named '{}'. Use 'display processes' to list available processes.",
                            value
                        )),
                    }
                } else if let Some(index) = name_match {
                    Ok(index)
                } else {
                    Err(eyre!(
                        "No process named '{}'. Use 'display processes' to list available processes.",
                        value
                    ))
                }
            }
        }
    }
}

pub trait ProcessListExt {
    fn find_integrand_ref(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)>;
    fn get_amplitude_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut Amplitude>;
    fn get_cross_section_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut CrossSection>;
}

impl ProcessListExt for ProcessList {
    fn find_integrand_ref(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)> {
        let process_id = match process {
            Some(process_ref) => process_ref.resolve(self)?,
            None => self.find_process(None)?,
        };
        let integrand_name = self.processes[process_id]
            .collection
            .find_integrand(integrand_name.cloned())
            .with_note(|| format!("in process id {process_id}"))?;
        Ok((process_id, integrand_name))
    }

    fn get_amplitude_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut Amplitude> {
        let (process_id, integrand_name) = self.find_integrand_ref(process, integrand_name)?;
        let process = &mut self.processes[process_id];
        match &mut process.collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                amplitudes.get_mut(&integrand_name).ok_or_else(|| {
                    eyre!(
                        "No amplitude named '{}' in process '{}'",
                        integrand_name,
                        process.definition.folder_name
                    )
                })
            }
            ProcessCollection::CrossSections(_) => Err(eyre!(
                "Process '{}' does not contain amplitudes",
                process.definition.folder_name
            )),
        }
    }

    fn get_cross_section_mut_ref(
        &mut self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<&mut CrossSection> {
        let (process_id, integrand_name) = self.find_integrand_ref(process, integrand_name)?;
        let process = &mut self.processes[process_id];
        match &mut process.collection {
            ProcessCollection::CrossSections(crosssections) => {
                crosssections.get_mut(&integrand_name).ok_or_else(|| {
                    eyre!(
                        "No cross section named '{}' in process '{}'",
                        integrand_name,
                        process.definition.folder_name
                    )
                })
            }
            ProcessCollection::Amplitudes(_) => Err(eyre!(
                "Process '{}' does not contain crosssections",
                process.definition.folder_name
            )),
        }
    }
}

pub trait SyncSettings {
    fn sync_settings(&self) -> Result<()>;
}

impl SyncSettings for CLISettings {
    fn sync_settings(&self) -> Result<()> {
        // println!("Syncing settings {}", self.global.logfile_directive);
        set_file_log_filter(&self.global.logfile_directive)?;
        set_stderr_log_filter(&self.global.display_directive)?;
        set_log_style(self.global.log_style.clone());
        Ok(())
    }
}

// Static flag to control serialization behavior
static SERIALIZE_COMMANDS_AS_STRINGS: AtomicBool = AtomicBool::new(false);

fn is_commands_blocks_empty(commands_blocks: &Vec<CommandsBlock>) -> bool {
    commands_blocks.is_empty()
}

fn should_persist_command(command: &Commands) -> bool {
    !matches!(
        command,
        Commands::Quit(_) | Commands::StartCommandsBlock(_) | Commands::FinishCommandsBlock
    )
}

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
                CommandHistory::from_raw_string(value).map_err(|err| {
                    E::custom(format!(
                        "Failed to parse command string '{}': {}",
                        value, err
                    ))
                })
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
        use clap::error::ErrorKind;
        use clap::Parser;

        let args = split_command_line(raw_string).map_err(|_| {
            clap::Error::raw(
                ErrorKind::InvalidValue,
                "Could not parse command: unmatched quotes or trailing escape",
            )
        })?;
        let cli = Repl::try_parse_from(
            std::iter::once("gammaloop").chain(args.iter().map(String::as_str)),
        )?;

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

    #[serde(skip_serializing_if = "is_commands_blocks_empty")]
    pub commands_blocks: Vec<CommandsBlock>,
    // #[serde(with = "serde_yaml::with::singleton_map_recursive")]
    // #[schemars(with = "Vec<CommandHistory>")]
    pub commands: Vec<CommandHistory>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, JsonSchema, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct CommandsBlock {
    pub name: String,
    pub commands: Vec<CommandHistory>,
}

#[derive(Debug, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct RawRunHistoryToml {
    default_runtime_settings: RuntimeSettings,
    cli_settings: CLISettings,
    commands_blocks: Vec<RawCommandsBlockToml>,
    commands: Vec<TomlValue>,
}

#[derive(Debug, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct RawCommandsBlockToml {
    name: String,
    commands: Vec<TomlValue>,
}

impl CommandsBlock {
    pub fn semantically_eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.commands.len() == other.commands.len()
            && self
                .commands
                .iter()
                .zip(other.commands.iter())
                .all(|(left, right)| left.command == right.command)
    }
}

impl SmartSerde for RunHistory {
    fn has_schema_path(&self, online: bool) -> Option<Result<PathBuf>> {
        Some(get_schema_folder(online).map(|f| f.join("runhistory.json")))
    }
}

impl RunHistory {
    pub fn freeze_boot_settings_from(&mut self, boot_run_history: &RunHistory) {
        self.cli_settings.global = boot_run_history.cli_settings.global.clone();
        self.default_runtime_settings = boot_run_history.default_runtime_settings.clone();
    }

    pub fn frozen_boot_settings_match(&self, boot_run_history: &RunHistory) -> bool {
        self.cli_settings.global == boot_run_history.cli_settings.global
            && self.default_runtime_settings == boot_run_history.default_runtime_settings
    }

    /// Add a command to the run history
    pub fn push(&mut self, command: Commands) {
        self.push_with_raw(command, None);
    }

    /// Add a command with optional raw string to the run history
    ///
    /// If raw_string is provided, it will be stored alongside the command
    /// for potential later serialization as a string.
    pub fn push_with_raw(&mut self, command: Commands, raw_string: Option<String>) {
        if should_persist_command(&command) {
            self.commands.push(CommandHistory {
                command,
                raw_string,
            });
        }
    }

    pub fn schema() -> Schema {
        schema_for!(RunHistory)
    }

    pub fn validate(&self) -> Result<()> {
        let mut seen_names = HashSet::with_capacity(self.commands_blocks.len());
        for block in &self.commands_blocks {
            if block.name.trim().is_empty() {
                return Err(eyre!(
                    "Run card `commands_blocks` contains a block with an empty name"
                ));
            }
            if !seen_names.insert(block.name.clone()) {
                return Err(eyre!(
                    "Run card `commands_blocks` contains duplicate block name '{}'",
                    block.name
                ));
            }
        }
        Ok(())
    }

    pub fn command_block(&self, name: &str) -> Option<&CommandsBlock> {
        self.commands_blocks.iter().find(|block| block.name == name)
    }

    pub fn select_commands_blocks(
        &self,
        selected_block_names: &[String],
    ) -> Result<Vec<CommandsBlock>> {
        let mut selected = Vec::with_capacity(selected_block_names.len());
        for name in selected_block_names {
            let block = self.command_block(name).ok_or_else(|| {
                eyre!(
                    "Unknown command block '{}'. Available command blocks: {}",
                    name,
                    self.commands_blocks
                        .iter()
                        .map(|block| block.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            })?;
            selected.push(block.clone());
        }
        Ok(selected)
    }

    pub fn merge_commands_blocks(&mut self, commands_blocks: &[CommandsBlock]) -> Result<()> {
        for new_block in commands_blocks {
            match self.command_block(&new_block.name) {
                Some(existing_block) if existing_block.semantically_eq(new_block) => {}
                Some(_) => {
                    return Err(eyre!(
                        "Run card command block '{}' redefines an existing block with different commands",
                        new_block.name
                    ));
                }
                None => self.commands_blocks.push(new_block.clone()),
            }
        }
        self.validate()
    }

    pub fn apply_session_settings(
        &self,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<()> {
        if self.cli_settings.global != GlobalSettings::default() {
            global_settings.global = self.cli_settings.global.clone();
            global_settings.sync_settings()?;
        }
        if self.default_runtime_settings != RuntimeSettings::default() {
            *default_runtime_settings = self.default_runtime_settings.clone();
        }
        Ok(())
    }

    pub fn run(
        &mut self,
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<SaveState>> {
        let mut session_state = crate::session::CliSessionState::default();
        let mut session = crate::session::CliSession::new(
            state,
            self,
            global_settings,
            default_runtime_settings,
            &mut session_state,
        );
        session.replay_run_history()
    }

    fn filtered_for_save(&self) -> Self {
        let mut filtered = self.clone();
        filtered
            .commands
            .retain(|command_history| should_persist_command(&command_history.command));
        filtered
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        debug!("Loaded run history from file {}", path.display());

        let runhistory = match path.extension().and_then(|ext| ext.to_str()) {
            Some("toml") => Self::load_toml(path)?,
            _ => Self::from_file(path, "run history")?,
        };
        runhistory.validate()?;
        Ok(runhistory)
    }

    pub fn save_toml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        _strict: bool,
    ) -> Result<()> {
        self.filtered_for_save()
            .to_file(root_folder.join("run.toml"), override_state_file)?;

        //Self::schema().to_file(root_folder.join("run_schema.json"))?;
        Ok(())
    }

    pub fn save_yaml(
        &self,
        root_folder: &Path,
        override_state_file: bool,
        _strict: bool,
    ) -> Result<()> {
        self.filtered_for_save()
            .to_file(root_folder.join("run.yaml"), override_state_file)?;
        //Self::schema().to_file(root_folder.join("run_schema.json"))?;
        Ok(())
    }

    fn load_toml(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("Could not open run history file {}", path.display()))?;
        let raw_run_history: RawRunHistoryToml =
            toml::from_str(&raw).wrap_err("Error parsing run history toml")?;

        let commands = raw_run_history
            .commands
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                parse_toml_command_history(value, &format!("top-level command #{}", index + 1))
            })
            .collect::<Result<Vec<_>>>()?;

        let commands_blocks = raw_run_history
            .commands_blocks
            .into_iter()
            .map(|block| {
                let RawCommandsBlockToml { name, commands } = block;
                let commands = commands
                    .into_iter()
                    .enumerate()
                    .map(|(index, value)| {
                        parse_toml_command_history(
                            value,
                            &format!("command block '{}' command #{}", name, index + 1),
                        )
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(CommandsBlock { name, commands })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            default_runtime_settings: raw_run_history.default_runtime_settings,
            cli_settings: raw_run_history.cli_settings,
            commands_blocks,
            commands,
        })
    }
}

fn parse_toml_command_history(value: TomlValue, context: &str) -> Result<CommandHistory> {
    match value {
        TomlValue::String(raw) => CommandHistory::from_raw_string(&raw)
            .map_err(|err| eyre!("Failed to parse {} '{}': {}", context, raw, err)),
        other => other
            .try_into::<CommandHistory>()
            .map_err(|err| eyre!("Failed to parse {}: {}", context, err)),
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
}

const STATE_MANIFEST_FILE: &str = "state_manifest.toml";
const CURRENT_STATE_MANIFEST_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
struct StateManifest {
    version: u32,
}

impl Default for StateManifest {
    fn default() -> Self {
        Self {
            version: CURRENT_STATE_MANIFEST_VERSION,
        }
    }
}

fn run_state_migration_checks(manifest: &StateManifest, save_path: &Path) -> Result<()> {
    if manifest.version > CURRENT_STATE_MANIFEST_VERSION {
        return Err(eyre!(
            "State version {} is newer than this binary supports (max {}). Please upgrade gammaloop.",
            manifest.version,
            CURRENT_STATE_MANIFEST_VERSION
        ));
    }

    match manifest.version {
        1 => {
            if !save_path.join("symbolica_state.bin").exists() {
                return Err(eyre!(
                    "Saved state at '{}' is missing required file symbolica_state.bin",
                    save_path.display()
                ));
            }
            if !save_path.join("processes").exists() {
                return Err(eyre!(
                    "Saved state at '{}' is missing required folder processes/",
                    save_path.display()
                ));
            }
            if !save_path.join("model.json").exists() {
                return Err(eyre!(
                    "Saved state at '{}' is missing required file model.json",
                    save_path.display()
                ));
            }
            Ok(())
        }
        _ => Err(eyre!(
            "State version {} is not supported by this binary.",
            manifest.version
        )),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StateFolderKind {
    Missing,
    Scratch,
    Saved,
    Invalid(String),
}

fn is_scratch_state_entry(entry: &fs::DirEntry) -> bool {
    entry.file_type().map(|ft| ft.is_dir()).unwrap_or(false)
        && entry.file_name().to_string_lossy() == "logs"
}

pub fn classify_state_folder(save_path: &Path) -> Result<StateFolderKind> {
    if !save_path.exists() {
        return Ok(StateFolderKind::Missing);
    }
    if !save_path.is_dir() {
        return Ok(StateFolderKind::Invalid(format!(
            "'{}' exists but is not a directory",
            save_path.display()
        )));
    }

    let manifest_path = save_path.join(STATE_MANIFEST_FILE);
    if manifest_path.exists() {
        let manifest = load_state_manifest(save_path)?;
        return Ok(match run_state_migration_checks(&manifest, save_path) {
            Ok(()) => StateFolderKind::Saved,
            Err(err) => StateFolderKind::Invalid(err.to_string()),
        });
    }

    let mut entries = fs::read_dir(save_path)
        .with_context(|| format!("Trying to read state folder '{}'", save_path.display()))?;
    if entries.by_ref().all(|entry| {
        entry
            .map(|entry| is_scratch_state_entry(&entry))
            .unwrap_or(false)
    }) {
        return Ok(StateFolderKind::Scratch);
    }

    Ok(StateFolderKind::Invalid(format!(
        "State folder '{}' is missing required file {}",
        save_path.display(),
        STATE_MANIFEST_FILE
    )))
}

fn load_state_manifest(save_path: &Path) -> Result<StateManifest> {
    let manifest_path = save_path.join(STATE_MANIFEST_FILE);
    let raw_manifest = fs::read_to_string(&manifest_path).with_context(|| {
        format!(
            "Trying to read state manifest file {}",
            manifest_path.display()
        )
    })?;
    let manifest = toml::from_str::<StateManifest>(&raw_manifest).with_context(|| {
        format!(
            "Trying to parse state manifest file {}",
            manifest_path.display()
        )
    })?;

    run_state_migration_checks(&manifest, save_path)?;
    Ok(manifest)
}

fn save_state_manifest(save_path: &Path) -> Result<()> {
    let manifest = StateManifest::default();
    let raw_manifest =
        toml::to_string_pretty(&manifest).context("Trying to serialize state manifest to TOML")?;
    fs::write(save_path.join(STATE_MANIFEST_FILE), raw_manifest).with_context(|| {
        format!(
            "Trying to write state manifest file {}",
            save_path.join(STATE_MANIFEST_FILE).display()
        )
    })?;
    Ok(())
}

impl State {
    pub fn import_model(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.model = Model::from_file(path)?;
        Ok(())
    }

    pub fn remove_process(&mut self, process: Option<&ProcessRef>) -> Result<RemovedProcess> {
        let process_id = self.resolve_process_ref(process)?;
        let removed = self.process_list.processes.remove(process_id);
        Ok(RemovedProcess {
            process_id,
            process_name: removed.definition.folder_name,
        })
    }

    pub fn remove_integrand(
        &mut self,
        process: &ProcessRef,
        integrand_name: &str,
    ) -> Result<RemovedIntegrand> {
        let process_id = process.resolve(&self.process_list)?;
        let process_entry = &mut self.process_list.processes[process_id];
        let process_name = process_entry.definition.folder_name.clone();
        let canonical_integrand_name = process_entry
            .collection
            .find_integrand(Some(integrand_name.to_string()))?;
        process_entry
            .collection
            .remove_integrand(&canonical_integrand_name)?;

        let removed_empty_process = process_entry.collection.get_integrand_names().is_empty();
        if removed_empty_process {
            self.process_list.processes.remove(process_id);
        }

        Ok(RemovedIntegrand {
            process_id,
            process_name,
            integrand_name: canonical_integrand_name,
            removed_empty_process,
        })
    }

    pub fn resolve_process_ref(&self, process: Option<&ProcessRef>) -> Result<usize> {
        match process {
            Some(process_ref) => process_ref.resolve(&self.process_list),
            None => self.process_list.find_process(None),
        }
    }

    pub fn find_integrand_ref(
        &self,
        process: Option<&ProcessRef>,
        integrand_name: Option<&String>,
    ) -> Result<(usize, String)> {
        let process_id = self.resolve_process_ref(process)?;
        let integrand_name = self.process_list.processes[process_id]
            .collection
            .find_integrand(integrand_name.cloned())
            .with_note(|| format!("in process id {process_id}"))?;
        Ok((process_id, integrand_name))
    }

    pub fn generate_integrands(
        &mut self,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<()> {
        let generation_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(global_settings.n_cores.generate)
            .build()?;

        self.process_list.preprocess(
            &self.model,
            global_settings,
            &runtime_default,
            &generation_pool,
        )?;
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
                        a.preprocess(
                            &self.model,
                            &global_settings.generation,
                            &runtime_default,
                            &generation_pool,
                        )?;
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
            p.preprocess(
                &self.model,
                global_settings,
                &runtime_default,
                &generation_pool,
            )?;
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

    pub fn export_dots(
        &mut self,
        path: impl AsRef<Path>,
        settings: &DotExportSettings,
    ) -> Result<()> {
        self.process_list.export_dot(path, settings)?;
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
                    ));
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
        integrand_name: String,
        _n_cores: usize,
    ) -> Result<()> {
        let integrand = self
            .process_list
            .get_integrand_mut(process_id, integrand_name)?;
        let name = integrand.name();

        info!(
            "\nBenchmarking runtime of integrand '{}' over {} samples...\n",
            name.green(),
            samples.to_string().blue()
        );

        let now = Instant::now();
        for _ in 0..samples {
            let _ = integrand.evaluate_sample(
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
        super::tracing::init_tracing(log_dir.as_ref().join("logs"), log_file_name);

        Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
        }
    }

    pub fn new_test() -> Self {
        let _ = init_test_tracing();

        Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
        }
    }

    pub fn new_bench() -> Self {
        let _ = init_bench_tracing();

        Self {
            model: Model::default(),
            process_list: ProcessList::default(),
            model_parameters: InputParamCard::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedProcess {
    pub process_id: usize,
    pub process_name: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedIntegrand {
    pub process_id: usize,
    pub process_name: String,
    pub integrand_name: String,
    pub removed_empty_process: bool,
}

impl State {
    pub fn load(
        save_path: PathBuf,
        model_path: Option<PathBuf>,
        trace_logs_filename: Option<String>,
    ) -> Result<Self> {
        // let root_folder = root_folder.join("gammaloop_state");
        let manifest = load_state_manifest(&save_path)?;
        debug!("Loading state manifest version {}", manifest.version);

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

        fs::create_dir_all(&selected_root_folder)?;

        let mut state_file =
        // info!("Hi");
            fs::File::create(selected_root_folder.join("symbolica_state.bin"))?;

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
        save_state_manifest(&selected_root_folder)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use gammalooprs::{
        momentum::{Dep, ExternalMomenta, SignOrZero},
        settings::runtime::kinematic::improvement::PhaseSpaceImprovementSettings,
        settings::{runtime::kinematic::Externals, KinematicsSettings, RuntimeSettings},
        utils::serde_utils::SHOWDEFAULTS,
    };
    use tempfile::tempdir;

    use crate::commands::{
        display::Display,
        save::SaveState,
        set::{ProcessSetArgs, Set, SetArgs},
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
    fn run_history_applies_default_runtime_before_commands() {
        let mut run_history = RunHistory::default();
        run_history.default_runtime_settings = toml::from_str(
            r#"
[sampling]
type = "discrete_graph_sampling"
[sampling.sampling_type]
subtype = "tropical"
"#,
        )
        .unwrap();

        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        let mut default_runtime_settings = RuntimeSettings::default();

        let _ = run_history
            .run(&mut state, &mut cli_settings, &mut default_runtime_settings)
            .unwrap();

        assert_eq!(
            default_runtime_settings.sampling,
            run_history.default_runtime_settings.sampling
        );
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

    #[test]
    fn run_history_push_with_raw_skips_quit_and_definition_commands() {
        use super::RunHistory;
        use crate::commands::{run::Run, StartCommandsBlock};

        let mut run_history = RunHistory::default();
        run_history.push_with_raw(
            Commands::Run(Run {
                block_names: vec!["block_a".to_string()],
                commands: None,
            }),
            Some("run block_a".to_string()),
        );
        run_history.push_with_raw(
            Commands::Quit(SaveState::default()),
            Some("quit".to_string()),
        );
        run_history.push_with_raw(
            Commands::StartCommandsBlock(StartCommandsBlock {
                name: "block_a".to_string(),
            }),
            Some("start_commands_block block_a".to_string()),
        );
        run_history.push_with_raw(
            Commands::FinishCommandsBlock,
            Some("finish_commands_block".to_string()),
        );

        assert_eq!(run_history.commands.len(), 1);
        assert_eq!(
            run_history.commands[0].raw_string.as_deref(),
            Some("run block_a")
        );
    }

    #[test]
    fn run_history_filtered_for_save_preserves_commands_blocks() {
        let mut run_history = RunHistory::default();
        run_history.commands_blocks = vec![
            CommandsBlock {
                name: "generate".to_string(),
                commands: vec![CommandHistory::new_with_raw(
                    Commands::Display(Display::Processes),
                    "display processes".to_string(),
                )],
            },
            CommandsBlock {
                name: "integrate".to_string(),
                commands: vec![CommandHistory::new_with_raw(
                    CommandHistory::from_raw_string("quit -o").unwrap().command,
                    "quit -o".to_string(),
                )],
            },
        ];
        run_history.push_with_raw(
            CommandHistory::from_raw_string("display processes")
                .unwrap()
                .command,
            Some("display processes".to_string()),
        );

        let filtered = run_history.filtered_for_save();
        set_serialize_commands_as_strings(true);
        let toml = toml::to_string_pretty(&filtered).unwrap();
        set_serialize_commands_as_strings(false);

        assert_eq!(filtered.commands_blocks.len(), 2);
        assert!(toml.contains("commands = ["));
        assert!(toml.contains("[[commands_blocks]]"));
        assert!(toml.contains("quit -o"));
    }

    #[test]
    fn command_history_parses_multiline_set_string() {
        let raw = "set process -p epem_a_tth -i LO string '[integrator]\nn_start = 1000\n'";
        let cmd = CommandHistory::from_raw_string(raw).unwrap();
        assert_eq!(cmd.raw_string.as_deref(), Some(raw));

        match cmd.command {
            Commands::Set(Set::Process { input, .. }) => match input {
                ProcessSetArgs::String { string } => {
                    assert_eq!(string, "[integrator]\nn_start = 1000\n");
                }
                other => panic!("Expected string set input, got {other:?}"),
            },
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn command_history_parses_hash_process_refs() {
        let cmd = CommandHistory::from_raw_string("display integrands -p #12").unwrap();
        match cmd.command {
            Commands::Display(Display::Integrands { process }) => {
                assert_eq!(process, Some(ProcessRef::Id(12)));
            }
            other => panic!("Expected display integrands command, got {other:?}"),
        }
    }

    #[test]
    fn run_history_parses_triple_quoted_set_kv_command() {
        let toml = r#"
commands = [
    """set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[1.0,2.0,3.0,4.0],[5.0,6.0,-7.0,-8.0]],"helicities":[1,1]}}'""",
]

[default_runtime_settings.general]
use_picobarns = true
"#;

        let run_history: RunHistory = toml::from_str(toml).unwrap();
        assert_eq!(run_history.commands.len(), 1);

        let expected_cmd = r#"set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[1.0,2.0,3.0,4.0],[5.0,6.0,-7.0,-8.0]],"helicities":[1,1]}}'"#;
        let command_history = &run_history.commands[0];
        assert_eq!(command_history.raw_string.as_deref(), Some(expected_cmd));

        match &command_history.command {
            Commands::Set(Set::DefaultRuntime {
                input: SetArgs::Kv { pairs },
            }) => {
                assert_eq!(pairs.len(), 1);
                assert_eq!(pairs[0].key, "kinematics.externals");
                assert_eq!(
                    pairs[0].value,
                    r#"{"type":"constant","data":{"momenta":[[1.0,2.0,3.0,4.0],[5.0,6.0,-7.0,-8.0]],"helicities":[1,1]}}"#
                );
            }
            other => panic!("Expected set default-runtime kv command, got {other:?}"),
        }

        assert!(run_history.default_runtime_settings.general.use_picobarns);
    }

    #[test]
    fn run_history_load_preserves_commands_blocks() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[commands_blocks]]
name = "zeta"
commands = ["quit -n"]

[[commands_blocks]]
name = "alpha"
commands = ["quit -o"]
"#,
        )
        .unwrap();

        let run_history = RunHistory::load(&run_path).unwrap();
        assert!(run_history.commands.is_empty());
        assert_eq!(run_history.commands_blocks.len(), 2);
    }

    #[test]
    fn run_history_selects_named_commands_blocks_in_order() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[commands_blocks]]
name = "first"
commands = ["quit -n"]

[[commands_blocks]]
name = "second"
commands = ["quit -o"]
"#,
        )
        .unwrap();

        let requested = vec!["second".to_string(), "first".to_string()];
        let run_history = RunHistory::load(&run_path).unwrap();
        let selected = run_history
            .select_commands_blocks(requested.as_slice())
            .unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(
            selected[0].commands[0].raw_string.as_deref(),
            Some("quit -o")
        );
        assert_eq!(
            selected[1].commands[0].raw_string.as_deref(),
            Some("quit -n")
        );
    }

    #[test]
    fn run_history_selection_rejects_unknown_command_block() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[commands_blocks]]
name = "first"
commands = ["quit -n"]
"#,
        )
        .unwrap();

        let requested = vec!["missing".to_string()];
        let run_history = RunHistory::load(&run_path).unwrap();
        let err = run_history
            .select_commands_blocks(requested.as_slice())
            .unwrap_err();
        let message = format!("{err}");
        assert!(message.contains("Unknown command block 'missing'"));
        assert!(message.contains("first"));
    }

    #[test]
    fn run_history_selection_rejects_missing_command_block_when_only_legacy_commands_exist() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
commands = ["quit -o"]
"#,
        )
        .unwrap();

        let requested = vec!["first".to_string()];
        let run_history = RunHistory::load(&run_path).unwrap();
        let err = run_history
            .select_commands_blocks(requested.as_slice())
            .unwrap_err();
        assert!(format!("{err}").contains("Unknown command block"));
    }

    #[test]
    fn run_history_load_accepts_commands_and_commands_blocks_together() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
commands = ["quit -o"]

[[commands_blocks]]
name = "first"
commands = ["quit -n"]
"#,
        )
        .unwrap();

        let run_history = RunHistory::load(&run_path).unwrap();
        assert_eq!(run_history.commands.len(), 1);
        assert_eq!(run_history.commands_blocks.len(), 1);
        assert_eq!(
            run_history.commands[0].raw_string.as_deref(),
            Some("quit -o")
        );
    }

    #[test]
    fn run_history_load_rejects_duplicate_command_block_names() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        fs::write(
            &run_path,
            r#"
[[commands_blocks]]
name = "first"
commands = ["quit -o"]

[[commands_blocks]]
name = "first"
commands = ["quit -n"]
"#,
        )
        .unwrap();

        let err = RunHistory::load(&run_path).unwrap_err();
        assert!(format!("{err}").contains("duplicate block name"));
    }

    #[test]
    fn state_manifest_roundtrip_current_version() {
        let temp = tempdir().unwrap();
        save_state_manifest(temp.path()).unwrap();

        let manifest = load_state_manifest(temp.path()).unwrap();
        assert_eq!(manifest.version, CURRENT_STATE_MANIFEST_VERSION);
    }

    #[test]
    fn state_manifest_rejects_future_versions() {
        let temp = tempdir().unwrap();
        let future_manifest = StateManifest {
            version: CURRENT_STATE_MANIFEST_VERSION + 1,
        };
        fs::write(
            temp.path().join(STATE_MANIFEST_FILE),
            toml::to_string_pretty(&future_manifest).unwrap(),
        )
        .unwrap();

        let err = load_state_manifest(temp.path()).unwrap_err();
        assert!(format!("{err}").contains("newer than this binary supports"));
    }

    #[test]
    fn state_folder_classifies_saved_layout() {
        let temp = tempdir().unwrap();
        save_state_manifest(temp.path()).unwrap();
        fs::write(temp.path().join("model.json"), "{}").unwrap();
        fs::write(temp.path().join("symbolica_state.bin"), []).unwrap();
        fs::create_dir_all(temp.path().join("processes")).unwrap();

        assert_eq!(
            classify_state_folder(temp.path()).unwrap(),
            StateFolderKind::Saved
        );
    }

    #[test]
    fn state_folder_classifies_logs_only_folder_as_scratch() {
        let temp = tempdir().unwrap();
        fs::create_dir_all(temp.path().join("logs")).unwrap();
        fs::write(temp.path().join("logs").join("gammalog.jsonl"), "").unwrap();

        assert_eq!(
            classify_state_folder(temp.path()).unwrap(),
            StateFolderKind::Scratch
        );
    }
}

#[derive(Args, Debug, Clone)]
pub struct ExistingArgs {
    pub process_id: u32,
    pub name: Option<String>,
}
