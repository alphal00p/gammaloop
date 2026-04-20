#[cfg(all(
    feature = "no_pyo3",
    any(
        feature = "python_api",
        feature = "python_abi",
        feature = "python_stubgen",
        feature = "pyo3-extension-module",
        feature = "ufo_support",
    )
))]
compile_error!(
    "feature `no_pyo3` is incompatible with python/pyo3 features (python_api, python_abi, \
python_stubgen, pyo3-extension-module, ufo_support). Use --no-default-features --features \
cli,no_pyo3."
);

use ::tracing::debug;
use ::tracing::info;
use ::tracing::level_filters::LevelFilter;
use ::tracing::warn;
use clap::parser::ValueSource;
use clap::{CommandFactory, FromArgMatches, Parser, ValueEnum};
use clap_complete::shells::{Bash, Elvish, Fish, PowerShell, Zsh};
use clap_complete_nushell::Nushell;
use commands::save::SaveState;
use commands::Commands;

use gammalooprs::utils::serde_utils::{SerdeFileError, SHOWDEFAULTS};
use reedline::FileBackedHistory;
use repl::ClapEditor;
use repl::ReadCommandOutput;
use session::{CliSession, CliSessionState};
use tracing::{
    configure_file_log_boot_mode, get_stderr_log_filter_label, set_file_log_filter,
    set_file_log_filter_override, set_log_format_override, set_log_style, set_stderr_log_filter,
    set_stderr_log_filter_override,
};

// use clap_repl::{
//     reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
//     ClapEditor, ReadCommandOutput,
// };

use color_eyre::Result;
use colored::Colorize;
use console::{measure_text_width, style};
use dirs::home_dir;
use eyre::{eyre, Context};
use gammalooprs::{
    initialisation::initialise,
    processes::ProcessCollection,
    settings::{GlobalSettings, RuntimeSettings},
    utils::serde_utils::IsDefault,
    utils::{
        serde_utils::{get_schema_folder, is_false, is_true, SmartSerde},
        tracing::{LogFormat, LogLevel},
        GIT_VERSION,
    },
};
use reedline::{Prompt, PromptEditMode, PromptHistorySearch};
use schemars::{schema_for, JsonSchema};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use state::{
    classify_state_folder, CommandHistory, RunHistory, State, StateFolderKind, SyncSettings,
};
use std::{
    borrow::Cow, ffi::OsString, io::IsTerminal, path::Path, path::PathBuf, sync::atomic::Ordering,
};
use std::{fs, fs::File, ops::ControlFlow, time::Duration, time::Instant};
use walkdir::WalkDir;

// use tracing::LogLevel;
mod command_parser;
pub(crate) mod completion;
pub mod integrand_info;
pub(crate) mod model_parameters;
#[cfg(feature = "python_api")]
pub mod python;
pub mod repl;
pub mod session;
pub(crate) mod settings_tree;

pub mod state;
pub mod templates;
pub mod tracing;

#[cfg(test)]
pub(crate) static LOG_TEST_MUTEX: std::sync::Mutex<()> = std::sync::Mutex::new(());

pub(crate) const GLOBAL_SETTINGS_FILENAME: &str = "global_settings.toml";
pub(crate) const DEFAULT_RUNTIME_SETTINGS_FILENAME: &str = "default_runtime_settings.toml";
const BANNER_ART: &str = r"              ██         ▄████████▄  ▄████████▄  ██████████▄
  ██          ▀▀         ▀▀      ▀▀  ▀▀      ▀▀  ▀▀       ██
    ██  ▄██████████████████████████████████████████████████▀
     ▀██▀     ▄▄         ▄▄      ▄▄  ▄▄      ▄▄  ▄▄
    ██  ██    █████████  ▀████████▀  ▀████████▀  ██
   ██    ██
   ██    ██
";

pub(crate) fn render_smart_toml<T: SmartSerde>(value: &T) -> Result<String> {
    let mut toml_string = if let Some(schema_path) = value.has_schema_path(true) {
        let schema_path = schema_path?;
        format!("#:schema {}\n", schema_path.display())
    } else {
        String::new()
    };
    toml_string.push_str(&toml::to_string_pretty(value)?);
    Ok(toml_string)
}

fn lexically_normalized_absolute_path(path: &Path) -> Result<PathBuf> {
    let absolute_path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };

    let mut normalized = PathBuf::new();
    for component in absolute_path.components() {
        match component {
            std::path::Component::Prefix(prefix) => normalized.push(prefix.as_os_str()),
            std::path::Component::RootDir => normalized.push(component.as_os_str()),
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                normalized.pop();
            }
            std::path::Component::Normal(part) => normalized.push(part),
        }
    }

    Ok(normalized)
}

fn path_lies_within(base: &Path, candidate: &Path) -> Result<bool> {
    let normalized_base = lexically_normalized_absolute_path(base)?;
    let normalized_candidate = lexically_normalized_absolute_path(candidate)?;
    Ok(normalized_candidate == normalized_base
        || normalized_candidate.strip_prefix(&normalized_base).is_ok())
}

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
pub struct Repl {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Parser, Debug)]
#[command(
    name = "gammaLoop",
    version,
    about,
    subcommand_precedence_over_arg = true
)]
#[command(next_line_help = true)]
pub struct OneShot {
    /// Remove the resolved state folder before startup so the session starts from a blank state
    #[arg(long, default_value_t = false)]
    pub clean_state: bool,

    /// Optional TOML card to load at boot time
    #[arg(value_hint = clap::ValueHint::FilePath)]
    pub boot_commands_path: Option<PathBuf>,

    /// Path to the state folder
    #[arg(short = 's', long, default_value = "./gammaloop_state", value_hint = clap::ValueHint::DirPath)]
    pub state_folder: PathBuf,

    /// Internal flag indicating whether `state_folder` was explicitly set on CLI.
    #[arg(skip = false)]
    state_folder_explicitly_set: bool,

    /// Path to the model file
    #[arg(short = 'm', long, value_hint = clap::ValueHint::FilePath)]
    pub model_file: Option<PathBuf>,

    /// Skip saving state on exit
    #[arg(short = 'n', long, default_value_t = false, group = "saving")]
    no_save_state: bool,

    /// Save state to file after each call
    #[arg(short = 'o', long, default_value_t = false, group = "saving")]
    override_state: bool,

    /// Set the name of the file containing all traces from gammaloop (logs) for current session
    #[arg(short = 't', long = "trace-logs-filename")]
    trace_logs_filename: Option<String>,

    /// Set log level for current session
    #[arg(short = 'l')]
    level: Option<LogLevel>,

    /// Set logfile log level for current session
    #[arg(short = 'L', long = "logfile-level")]
    logfile_level: Option<LogLevel>,

    /// Type of prefix for the logging format
    #[arg(short = 'p', long = "logging-prefix")]
    logging_prefix: Option<LogFormat>,

    /// Prevent writes into the state folder for the lifetime of this session
    #[arg(long, default_value_t = false)]
    read_only_state: bool,

    /// Path to a global settings TOML file.
    #[arg(short = 'g', long = "settings-global", value_hint = clap::ValueHint::FilePath)]
    settings_global_path: Option<PathBuf>,

    /// Path to a default runtime settings TOML file.
    #[arg(
        short = 'r',
        long = "settings-runtime-defaults",
        value_hint = clap::ValueHint::FilePath
    )]
    settings_runtime_defaults_path: Option<PathBuf>,

    /// Try to serialize using strings when saving run history
    #[arg(long)]
    no_try_strings: bool,

    /// Generate a shell completion script for the gammaloop executable
    #[arg(long = "completions", value_enum)]
    completions: Option<CompletionShell>,

    // /// Debug level
    // #[arg(short = 'd', long, value_enum, default_value_t = LogLevel::Info)]
    // debug_level: LogLevel,
    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
enum CompletionShell {
    Bash,
    Elvish,
    Fish,
    PowerShell,
    Zsh,
    Nushell,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct CLISettings {
    #[serde(skip_serializing_if = "is_true")]
    pub try_strings: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub override_state: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub state: StateSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub global: GlobalSettings,
    #[serde(skip)]
    #[schemars(skip)]
    pub session: SessionSettings,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct StateSettings {
    pub folder: PathBuf,
    #[serde(
        default,
        skip_serializing_if = "skip_optional_nonempty_string",
        serialize_with = "serialize_optional_nonempty_string",
        deserialize_with = "deserialize_optional_nonempty_string"
    )]
    pub name: Option<String>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct SessionSettings {
    #[serde(skip)]
    #[schemars(skip)]
    pub read_only_state: bool,
    #[serde(skip)]
    #[schemars(skip)]
    pub(crate) read_only_state_origin: Option<ReadOnlyStateOrigin>,
    #[serde(skip)]
    #[schemars(skip)]
    pub startup_warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReadOnlyStateOrigin {
    UserRequested,
    BootSettingsMismatch,
}

impl SessionSettings {
    pub(crate) fn set_user_requested_read_only_state(&mut self, read_only_state: bool) {
        self.read_only_state = read_only_state;
        self.read_only_state_origin = read_only_state.then_some(ReadOnlyStateOrigin::UserRequested);
    }

    pub(crate) fn force_read_only_state(&mut self, origin: ReadOnlyStateOrigin) {
        if !self.read_only_state {
            self.read_only_state_origin = Some(origin);
        }
        self.read_only_state = true;
    }

    pub(crate) fn is_read_only_due_to_boot_settings_mismatch(&self) -> bool {
        self.read_only_state
            && self.read_only_state_origin == Some(ReadOnlyStateOrigin::BootSettingsMismatch)
    }
}

impl Default for StateSettings {
    fn default() -> Self {
        Self {
            folder: "./gammaloop_state".into(),
            name: None,
        }
    }
}

impl StateSettings {
    pub fn prompt_label(&self) -> String {
        self.name
            .as_deref()
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .map(str::to_string)
            .unwrap_or_else(|| self.folder.display().to_string())
    }
}

fn deserialize_optional_nonempty_string<'de, D>(
    deserializer: D,
) -> std::result::Result<Option<String>, D::Error>
where
    D: Deserializer<'de>,
{
    let value = Option::<String>::deserialize(deserializer)?;
    Ok(value.and_then(|name| {
        let trimmed = name.trim();
        (!trimmed.is_empty()).then(|| trimmed.to_string())
    }))
}

fn skip_optional_nonempty_string(value: &Option<String>) -> bool {
    value.is_none() && !SHOWDEFAULTS.load(Ordering::Relaxed)
}

fn serialize_optional_nonempty_string<S>(
    value: &Option<String>,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match value {
        Some(name) => serializer.serialize_str(name),
        None if SHOWDEFAULTS.load(Ordering::Relaxed) => serializer.serialize_str(""),
        None => serializer.serialize_none(),
    }
}

impl Default for CLISettings {
    fn default() -> Self {
        CLISettings {
            try_strings: true,
            override_state: false,
            state: StateSettings::default(),
            global: GlobalSettings::default(),
            session: SessionSettings::default(),
        }
    }
}

impl PartialEq for CLISettings {
    fn eq(&self, other: &Self) -> bool {
        self.try_strings == other.try_strings
            && self.override_state == other.override_state
            && self.state == other.state
            && self.global == other.global
    }
}

impl CLISettings {
    pub fn override_with(&mut self, cli: &OneShot) {
        self.try_strings = !cli.no_try_strings;

        self.override_state = cli.override_state;

        self.state.folder = cli.state_folder.clone();
        self.session
            .set_user_requested_read_only_state(cli.read_only_state);
    }

    pub(crate) fn ensure_write_target_outside_active_state(
        &self,
        target: &Path,
        operation: &str,
    ) -> Result<()> {
        if !self.session.read_only_state {
            return Ok(());
        }

        if !path_lies_within(&self.state.folder, target)? {
            return Ok(());
        }

        Err(eyre!(
            "Cannot {operation} at '{}' because this session was started with --read-only-state and the target lies inside the active state folder '{}'. Choose a path outside the active state folder or restart without --read-only-state.",
            target.display(),
            self.state.folder.display()
        ))
    }
}

impl SmartSerde for CLISettings {}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct StateLoadOption {
    pub clean_state: bool,
    pub boot_commands_path: Option<PathBuf>,
    pub state_folder: Option<PathBuf>,
    pub model_file: Option<PathBuf>,
    pub trace_logs_filename: Option<String>,
    pub level: Option<LogLevel>,
    pub logfile_level: Option<LogLevel>,
    pub logging_prefix: Option<LogFormat>,
    pub read_only_state: bool,
    pub settings_global_path: Option<PathBuf>,
    pub settings_runtime_defaults_path: Option<PathBuf>,
}

pub struct LoadedState {
    pub state: State,
    pub run_history: RunHistory,
    pub cli_settings: CLISettings,
    pub default_runtime_settings: RuntimeSettings,
    pub session_state: CliSessionState,
    pub state_load_summary: Option<StateLoadSummary>,
}

impl LoadedState {
    pub fn cli_session(&mut self) -> CliSession<'_> {
        CliSession::new(
            &mut self.state,
            &mut self.run_history,
            &mut self.cli_settings,
            &mut self.default_runtime_settings,
            &mut self.session_state,
        )
    }
}

pub struct Parsed {
    pub cli: OneShot,
    pub input_string: String,
    pub matches: clap::ArgMatches,
}

impl StateLoadOption {
    fn into_oneshot(self) -> OneShot {
        let state_folder_explicitly_set = self.state_folder.is_some();
        OneShot {
            clean_state: self.clean_state,
            boot_commands_path: self.boot_commands_path,
            state_folder: self
                .state_folder
                .unwrap_or_else(|| PathBuf::from("./gammaloop_state")),
            state_folder_explicitly_set,
            model_file: self.model_file,
            no_save_state: true,
            override_state: false,
            trace_logs_filename: self.trace_logs_filename,
            level: self.level,
            logfile_level: self.logfile_level,
            logging_prefix: self.logging_prefix,
            read_only_state: self.read_only_state,
            settings_global_path: self.settings_global_path,
            settings_runtime_defaults_path: self.settings_runtime_defaults_path,
            no_try_strings: false,
            completions: None,
            command: None,
        }
    }

    pub fn load(self) -> Result<LoadedState> {
        initialise()?;
        let mut one_shot = self.into_oneshot();
        let (loaded_state, boot_exit) = one_shot.bootstrap_session()?;
        if boot_exit.is_some() {
            return Err(eyre::eyre!(
                "Boot run history requested to exit, which is not supported by the API state-load entry point"
            ));
        }
        Ok(loaded_state)
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
struct SettingsFileOverrides {
    global: Option<GlobalSettings>,
    default_runtime_settings: Option<RuntimeSettings>,
}

impl SettingsFileOverrides {
    fn apply(
        &self,
        cli_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<()> {
        if let Some(global) = &self.global {
            cli_settings.global = global.clone();
            cli_settings.sync_settings()?;
        }
        if let Some(runtime_settings) = &self.default_runtime_settings {
            *default_runtime_settings = runtime_settings.clone();
        }
        Ok(())
    }

    fn apply_to_run_history(&self, run_history: &RunHistory) -> RunHistory {
        let mut overridden = run_history.clone();
        if let Some(global) = &self.global {
            overridden.cli_settings.global = global.clone();
        }
        if let Some(runtime_settings) = &self.default_runtime_settings {
            overridden.default_runtime_settings = runtime_settings.clone();
        }
        overridden
    }
}

#[derive(Debug, Clone)]
pub struct StateLoadSummary {
    pub elapsed: Duration,
    pub serialized_size_bytes: Option<u64>,
    pub total_graphs: usize,
}

#[derive(Clone)]
struct SessionPrompt {
    left_prompt: String,
}

impl Prompt for SessionPrompt {
    fn render_prompt_left(&self) -> Cow<'_, str> {
        Cow::Borrowed(&self.left_prompt)
    }

    fn render_prompt_right(&self) -> Cow<'_, str> {
        Cow::Borrowed("")
    }

    fn render_prompt_indicator(&self, _prompt_mode: PromptEditMode) -> Cow<'_, str> {
        Cow::Borrowed("> ")
    }

    fn render_prompt_multiline_indicator(&self) -> Cow<'_, str> {
        Cow::Borrowed("::: ")
    }

    fn render_prompt_history_search_indicator(
        &self,
        history_search: PromptHistorySearch,
    ) -> Cow<'_, str> {
        Cow::Owned(format!("(reverse-search: {}) ", history_search.term))
    }
}

fn make_repl_prompt(state_label: &str, pending_block_name: Option<&str>) -> Box<dyn Prompt> {
    let left_prompt = match pending_block_name {
        Some(name) => format!(
            "{} | γloop {} ",
            state_label,
            format!("[defining: {name}]").blue()
        ),
        None => format!("{} | γloop ", state_label),
    };
    Box::new(SessionPrompt { left_prompt })
}

fn format_duration_human(duration: Duration) -> String {
    if duration.as_secs() >= 60 {
        let minutes = duration.as_secs() / 60;
        let seconds = duration.as_secs_f64() - (minutes * 60) as f64;
        format!("{minutes}m {seconds:.1}s")
    } else if duration.as_secs_f64() >= 1.0 {
        format!("{:.2}s", duration.as_secs_f64())
    } else {
        format!("{}ms", duration.as_millis())
    }
}

fn format_size_human(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];

    let mut value = bytes as f64;
    let mut unit = 0usize;
    while value >= 1024.0 && unit < UNITS.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }

    if unit == 0 {
        format!("{bytes} {}", UNITS[unit])
    } else {
        format!("{value:.2} {}", UNITS[unit])
    }
}

fn serialized_size_of_path(path: &Path) -> Option<u64> {
    let mut total = 0u64;
    for entry in WalkDir::new(path) {
        let entry = entry.ok()?;
        if entry.file_type().is_file() {
            total = total.checked_add(entry.metadata().ok()?.len())?;
        }
    }
    Some(total)
}

fn total_graph_count(state: &State) -> usize {
    state
        .process_list
        .processes
        .iter()
        .map(|process| match &process.collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes
                .values()
                .map(|amplitude| amplitude.graphs.len())
                .sum::<usize>(),
            ProcessCollection::CrossSections(cross_sections) => cross_sections
                .values()
                .map(|cross_section| cross_section.supergraphs.len())
                .sum::<usize>(),
        })
        .sum::<usize>()
}

fn banner_footer_line(spec: &str) -> String {
    format!(
        r#"   ▀██████▀   version:{:<15} log level:{}         "#,
        GIT_VERSION, spec,
    )
}

fn banner_width(spec: &str) -> usize {
    BANNER_ART
        .lines()
        .map(measure_text_width)
        .chain(std::iter::once(measure_text_width(&banner_footer_line(
            spec,
        ))))
        .max()
        .unwrap_or_default()
}

fn print_state_load_summary(summary: &StateLoadSummary) {
    let spec_label = get_stderr_log_filter_label();
    let banner_width = banner_width(&spec_label);
    let plain_summary = format!(
        "State: load {} | disk {} | #graphs {}",
        format_duration_human(summary.elapsed),
        summary
            .serialized_size_bytes
            .map(format_size_human)
            .unwrap_or_else(|| "unknown".to_string()),
        summary.total_graphs
    );
    let left_padding =
        " ".repeat(banner_width.saturating_sub(measure_text_width(&plain_summary)) / 2);
    let state_label = "State:".blue();
    let load_label = "load".blue();
    let disk_label = "disk".blue();
    let graphs_label = "#graphs".blue();
    let load_value = format_duration_human(summary.elapsed).green();
    let disk_value = summary
        .serialized_size_bytes
        .map(format_size_human)
        .unwrap_or_else(|| "unknown".to_string())
        .green();
    let graph_value = summary.total_graphs.to_string().green();

    println!(
        "{}{} {} {} | {} {} | {} {}\n",
        left_padding,
        state_label,
        load_label,
        load_value,
        disk_label,
        disk_value,
        graphs_label,
        graph_value
    );
}

impl OneShot {
    pub fn new_cli_settings(&self, global: GlobalSettings) -> CLISettings {
        let mut session = SessionSettings::default();
        session.set_user_requested_read_only_state(self.read_only_state);
        CLISettings {
            try_strings: !self.no_try_strings,
            override_state: self.override_state,
            state: StateSettings {
                folder: self.state_folder.clone(),
                ..StateSettings::default()
            },
            global,
            session,
        }
    }

    pub fn new_test(state_folder: PathBuf) -> Self {
        OneShot {
            state_folder,
            state_folder_explicitly_set: false,
            boot_commands_path: None,
            model_file: None,
            no_save_state: true,
            override_state: false,
            command: None,
            level: None,
            logfile_level: None,
            logging_prefix: None,
            read_only_state: false,
            settings_global_path: None,
            settings_runtime_defaults_path: None,
            no_try_strings: false,
            completions: None,
            clean_state: false,
            trace_logs_filename: None,
        }
    }

    fn current_state_folder_kind(&self) -> Result<StateFolderKind> {
        classify_state_folder(&self.state_folder)
    }

    fn clean_resolved_state_folder(&self) -> Result<()> {
        if !self.clean_state {
            return Ok(());
        }

        if self.read_only_state {
            return Err(eyre!(
                "Cannot remove the active state folder '{}' because this session was started with --read-only-state. Restart without --read-only-state to use --clean-state.",
                self.state_folder.display()
            ));
        }

        let metadata = match fs::symlink_metadata(&self.state_folder) {
            Ok(metadata) => metadata,
            Err(err) if err.kind() == std::io::ErrorKind::NotFound => return Ok(()),
            Err(err) => {
                return Err(err).wrap_err_with(|| {
                    format!(
                        "Failed to inspect state path '{}' before cleaning",
                        self.state_folder.display()
                    )
                })
            }
        };

        if metadata.file_type().is_dir() && !metadata.file_type().is_symlink() {
            fs::remove_dir_all(&self.state_folder).wrap_err_with(|| {
                format!(
                    "Failed to remove state folder '{}' before startup",
                    self.state_folder.display()
                )
            })?;
        } else {
            fs::remove_file(&self.state_folder).wrap_err_with(|| {
                format!(
                    "Failed to remove state path '{}' before startup",
                    self.state_folder.display()
                )
            })?;
        }

        Ok(())
    }

    fn initial_cli_settings_for_startup(
        &self,
        state_folder_kind: &StateFolderKind,
        boot_run_history: Option<&RunHistory>,
        settings_file_overrides: &SettingsFileOverrides,
    ) -> Result<CLISettings> {
        let mut cli_settings = match state_folder_kind {
            StateFolderKind::Saved => Self::load_global_settings_file(&self.state_folder)?,
            StateFolderKind::Missing | StateFolderKind::Scratch | StateFolderKind::Unmanifested => {
                boot_run_history
                    .map(|run_history| run_history.cli_settings.clone())
                    .unwrap_or_else(|| self.new_cli_settings(GlobalSettings::default()))
            }
            StateFolderKind::Invalid(reason) => {
                return Err(eyre::eyre!(reason.clone()));
            }
        };
        cli_settings.override_with(self);
        if let Some(global) = &settings_file_overrides.global {
            cli_settings.global = global.clone();
        }
        Ok(cli_settings)
    }

    fn configure_startup_tracing(&self, cli_settings: &CLISettings) -> Result<()> {
        if self.read_only_state && !matches!(self.logfile_level, None | Some(LogLevel::Off)) {
            return Err(eyre::eyre!(
                "--read-only-state is incompatible with enabling logfile output"
            ));
        }

        set_file_log_filter(&cli_settings.global.logfile_directive)?;
        set_stderr_log_filter(&cli_settings.global.display_directive)?;
        set_log_style(cli_settings.global.log_style.clone());
        set_log_format_override(self.logging_prefix);
        set_stderr_log_filter_override(
            self.level
                .map(|level| level.to_cli_display_directive_spec().to_string()),
        )?;

        let (file_override, hard_disable_file_logs, hard_disable_reason) = if self.read_only_state {
            (
                Some(LogLevel::Off.to_cli_logfile_directive_spec().to_string()),
                true,
                Some("--read-only-state"),
            )
        } else if matches!(self.logfile_level, Some(LogLevel::Off)) {
            (
                Some(LogLevel::Off.to_cli_logfile_directive_spec().to_string()),
                true,
                Some("--logfile-level off"),
            )
        } else {
            (
                self.logfile_level
                    .map(|level| level.to_cli_logfile_directive_spec().to_string()),
                false,
                None,
            )
        };

        set_file_log_filter_override(file_override)?;
        configure_file_log_boot_mode(hard_disable_file_logs, hard_disable_reason)?;
        Ok(())
    }

    fn load_boot_run_history(&self) -> Result<Option<RunHistory>> {
        self.boot_commands_path
            .as_ref()
            .map(RunHistory::load)
            .transpose()
    }

    fn load_settings_file_overrides(&self) -> Result<SettingsFileOverrides> {
        let global = self
            .settings_global_path
            .as_ref()
            .map(|path| CLISettings::from_file_typed(path).map(|settings| settings.global))
            .transpose()?;
        let default_runtime_settings = self
            .settings_runtime_defaults_path
            .as_ref()
            .map(RuntimeSettings::from_file_typed)
            .transpose()?;

        Ok(SettingsFileOverrides {
            global,
            default_runtime_settings,
        })
    }

    fn load_global_settings_file(state_folder: &std::path::Path) -> Result<CLISettings> {
        let global_settings_path = state_folder.join(GLOBAL_SETTINGS_FILENAME);
        match CLISettings::from_file_typed(&global_settings_path) {
            Ok(settings) => Ok(settings),
            Err(SerdeFileError::FileError(_)) => Ok(CLISettings::default()),
            Err(err) => Err(err.into()),
        }
    }

    fn subcmd_input_string(
        argv: &[OsString],
        cmd: &clap::Command,
        matches: &clap::ArgMatches,
    ) -> Option<String> {
        let (sc_name, _) = matches.subcommand()?; // no subcommand -> None

        // Find the Command that corresponds to the canonical subcommand name
        let sc = cmd.get_subcommands().find(|c| c.get_name() == sc_name)?;

        // Tokens that could have been used to invoke it: canonical + all aliases
        let mut names: Vec<&str> = vec![sc_name];
        names.extend(sc.get_all_aliases());

        // Locate the index in argv where the subcommand token appears
        let idx = argv.iter().position(|t| {
            let s = t.to_string_lossy();
            names.iter().any(|&n| s == n)
        })?;

        // Join from the subcommand onward
        Some(
            argv[idx..]
                .iter()
                .map(|s| s.to_string_lossy().into_owned())
                .collect::<Vec<_>>()
                .join(" "),
        )
    }

    /// Parse from env args *and* capture ArgMatches (explicit vs defaults).
    pub fn parse_env_with_capture() -> Result<Parsed, clap::Error> {
        let argv: Vec<OsString> = command_parser::normalize_clap_args(
            std::env::args_os()
                .map(|arg| arg.to_string_lossy().into_owned())
                .collect(),
        )
        .into_iter()
        .map(OsString::from)
        .collect();

        // Build a Command (same as derive(Parser)) and get matches
        let mut cmd = <OneShot as CommandFactory>::command();
        let matches = cmd.clone().try_get_matches_from(&argv)?;

        let cli = <OneShot as FromArgMatches>::from_arg_matches(&matches)
            .map_err(|e| e.format(&mut cmd))?;
        let mut cli = cli;
        cli.state_folder_explicitly_set =
            matches.value_source("state_folder") == Some(ValueSource::CommandLine);

        let input_string = OneShot::subcmd_input_string(&argv, &cmd, &matches).unwrap_or_default();
        Ok(Parsed {
            input_string,
            matches,
            cli,
        })
    }

    fn resolve_initial_state_folder(
        &self,
        boot_run_history: Option<&RunHistory>,
    ) -> Result<PathBuf> {
        if self.state_folder_explicitly_set {
            return Ok(self.state_folder.clone());
        }

        if let Some(run_history) = boot_run_history {
            return Ok(run_history.cli_settings.state.folder.clone());
        }

        Ok(self.state_folder.clone())
    }

    fn load_with_boot_context(
        &mut self,
        state_folder_kind: StateFolderKind,
        boot_run_history: Option<&RunHistory>,
        settings_file_overrides: &SettingsFileOverrides,
    ) -> Result<(
        State,
        RunHistory,
        CLISettings,
        RuntimeSettings,
        Option<StateLoadSummary>,
    )> {
        let startup_cli_settings = self.initial_cli_settings_for_startup(
            &state_folder_kind,
            boot_run_history,
            settings_file_overrides,
        )?;
        self.configure_startup_tracing(&startup_cli_settings)?;

        let (state, run_history, cli_settings, default_runtime_settings, load_summary) =
            match state_folder_kind {
                StateFolderKind::Saved => {
                    let load_started = Instant::now();
                    let mut state = State::load(
                        self.state_folder.clone(),
                        self.model_file.clone(),
                        self.trace_logs_filename.clone(),
                    )
                    .wrap_err_with(|| {
                        format!(
                            "Failed to load existing state from {}",
                            self.state_folder.display()
                        )
                    })?;
                    let allow_symjit_fallback =
                        startup_cli_settings.global.generation.evaluator.compile
                            && matches!(
                                startup_cli_settings
                                    .global
                                    .generation
                                    .compile
                                    .compilation_mode,
                                gammalooprs::settings::global::CompilationMode::Symjit
                            );
                    state.activate_loaded_integrand_backends(allow_symjit_fallback)?;

                    let default_runtime = match RuntimeSettings::from_file_typed(
                        self.state_folder.join(DEFAULT_RUNTIME_SETTINGS_FILENAME),
                    ) {
                        Ok(a) => a,
                        Err(SerdeFileError::FileError(_)) => RuntimeSettings::default(),
                        Err(e) => return Err(e.into()),
                    };
                    let run_path = self.state_folder.join("run.toml");
                    let run_history = if run_path.exists() {
                        RunHistory::load(run_path)?
                    } else {
                        RunHistory::default()
                    };

                    let load_summary = StateLoadSummary {
                        elapsed: load_started.elapsed(),
                        serialized_size_bytes: serialized_size_of_path(&self.state_folder),
                        total_graphs: total_graph_count(&state),
                    };

                    (
                        state,
                        run_history,
                        startup_cli_settings,
                        default_runtime,
                        Some(load_summary),
                    )
                }
                StateFolderKind::Missing | StateFolderKind::Scratch => {
                    info!(
                        "{} {}",
                        "Initializing new state in".blue(),
                        self.state_folder.display().to_string().green()
                    );

                    let mut run_history = RunHistory::default();
                    if let Some(boot_run_history) = boot_run_history {
                        run_history.freeze_boot_settings_from(boot_run_history);
                    }

                    (
                        State::new(self.state_folder.clone(), self.trace_logs_filename.clone()),
                        run_history,
                        startup_cli_settings,
                        RuntimeSettings::default(),
                        None,
                    )
                }
                StateFolderKind::Unmanifested => {
                    warn!(
                        "State folder '{}' has no {}; treating leftover contents as blank scratch state.",
                        self.state_folder.display(),
                        "state_manifest.toml",
                    );
                    info!(
                        "{} {}",
                        "Initializing new state in".blue(),
                        self.state_folder.display().to_string().green()
                    );

                    let mut run_history = RunHistory::default();
                    if let Some(boot_run_history) = boot_run_history {
                        run_history.freeze_boot_settings_from(boot_run_history);
                    }

                    (
                        State::new(self.state_folder.clone(), self.trace_logs_filename.clone()),
                        run_history,
                        startup_cli_settings,
                        RuntimeSettings::default(),
                        None,
                    )
                }
                StateFolderKind::Invalid(_) => unreachable!(),
            };

        cli_settings.sync_settings()?;

        Ok((
            state,
            run_history,
            cli_settings,
            default_runtime_settings,
            load_summary,
        ))
    }

    pub fn load(
        &mut self,
    ) -> Result<(
        State,
        RunHistory,
        CLISettings,
        RuntimeSettings,
        Option<StateLoadSummary>,
    )> {
        let boot_run_history = self.load_boot_run_history()?;
        let settings_file_overrides = self.load_settings_file_overrides()?;
        self.state_folder = self.resolve_initial_state_folder(boot_run_history.as_ref())?;
        self.clean_resolved_state_folder()?;
        let state_folder_kind = self.current_state_folder_kind()?;
        if let StateFolderKind::Invalid(reason) = &state_folder_kind {
            return Err(eyre::eyre!(reason.clone()));
        }
        self.load_with_boot_context(
            state_folder_kind,
            boot_run_history.as_ref(),
            &settings_file_overrides,
        )
    }

    fn bootstrap_session(&mut self) -> Result<(LoadedState, Option<SaveState>)> {
        let boot_run_history = self.load_boot_run_history()?;
        let settings_file_overrides = self.load_settings_file_overrides()?;
        self.state_folder = self.resolve_initial_state_folder(boot_run_history.as_ref())?;
        self.clean_resolved_state_folder()?;
        let state_folder_kind = self.current_state_folder_kind()?;
        if let StateFolderKind::Invalid(reason) = &state_folder_kind {
            return Err(eyre::eyre!(reason.clone()));
        }
        let booted_existing_state = matches!(state_folder_kind, StateFolderKind::Saved);

        let (
            mut state,
            mut run_history,
            mut cli_settings,
            mut default_runtime_settings,
            state_load_summary,
        ) = self.load_with_boot_context(
            state_folder_kind,
            boot_run_history.as_ref(),
            &settings_file_overrides,
        )?;
        settings_file_overrides.apply(&mut cli_settings, &mut default_runtime_settings)?;

        let mut session_state = CliSessionState::default();
        let mut session = CliSession::new(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut default_runtime_settings,
            &mut session_state,
        );

        let mut boot_exit = None;
        if let Some(boot_run_history) = boot_run_history.as_ref() {
            let effective_boot_run_history =
                settings_file_overrides.apply_to_run_history(boot_run_history);
            if let ControlFlow::Break(save_state) = session.apply_boot_run_history(
                boot_run_history,
                &effective_boot_run_history,
                booted_existing_state,
            )? {
                boot_exit = Some(save_state);
            }
        }

        Ok((
            LoadedState {
                state,
                run_history,
                cli_settings,
                default_runtime_settings,
                session_state,
                state_load_summary,
            },
            boot_exit,
        ))
    }

    pub fn run(mut self, raw: String) -> Result<()> {
        if let Some(shell) = self.completions {
            print!("{}", generate_completion_script(shell));
            return Ok(());
        }

        initialise()?;

        let (
            LoadedState {
                mut state,
                mut run_history,
                mut cli_settings,
                mut default_runtime_settings,
                mut session_state,
                state_load_summary,
            },
            boot_exit,
        ) = self.bootstrap_session()?;

        let boot_requested_exit = boot_exit.is_some();
        let mut save_state = boot_exit.unwrap_or_default();
        let mut session = CliSession::new(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut default_runtime_settings,
            &mut session_state,
        );

        if !boot_requested_exit {
            let had_initial_command = self.command.is_some();
            let mut enter_repl = !had_initial_command;

            if let Some(a) = self.command.take() {
                let command = if raw.trim().is_empty() {
                    CommandHistory::new(a)
                } else {
                    CommandHistory::new_with_raw(a, raw)
                };
                match session.execute_command(command)?.flow {
                    ControlFlow::Break(a) => save_state = a,
                    ControlFlow::Continue(()) => {
                        enter_repl =
                            std::io::stdin().is_terminal() && std::io::stdout().is_terminal();
                    }
                }
            }

            if enter_repl {
                if !had_initial_command {
                    print_banner();
                    if let Some(summary) = state_load_summary.as_ref() {
                        print_state_load_summary(summary);
                    }
                }
                run_repl_session(&mut session, &mut save_state);
            }
        }

        let implicit_read_only_exit =
            cli_settings.session.read_only_state && save_state == SaveState::default();
        let requested_active_state_save_after_auto_read_only = !self.no_save_state
            && !implicit_read_only_exit
            && !save_state.no_save_state
            && cli_settings
                .session
                .is_read_only_due_to_boot_settings_mismatch()
            && path_lies_within(
                &cli_settings.state.folder,
                &save_state
                    .path
                    .clone()
                    .unwrap_or_else(|| cli_settings.state.folder.clone()),
            )?;
        if requested_active_state_save_after_auto_read_only {
            warn!(
                "Skipping save state to {} because boot card settings differ from the frozen settings stored in {} and this session was forced into --read-only-state. The active state was left unchanged; save to a path outside the active state folder if you want a separate snapshot.",
                save_state
                    .path
                    .as_deref()
                    .unwrap_or(&cli_settings.state.folder)
                    .display(),
                cli_settings.state.folder.join("run.toml").display()
            );
        }
        if !self.no_save_state
            && !implicit_read_only_exit
            && !requested_active_state_save_after_auto_read_only
        {
            debug!("Saving State, override: {}", self.override_state);
            save_state.save(
                &mut state,
                &run_history,
                &default_runtime_settings,
                &cli_settings,
            )?
        }
        Ok(())
    }

    // pub fn initialize(&self) {}
}

fn run_repl_session(session: &mut CliSession<'_>, save_state: &mut SaveState) {
    let completion_state = repl::SharedCompletionState::new();
    completion_state.update_from_session(session);
    let mut repl = ClapEditor::<Repl>::builder()
        .with_prompt(make_repl_prompt(&session.prompt_state_label(), None))
        .with_completion_state(completion_state.clone());

    if let Some(home) = home_dir() {
        repl = repl.with_editor_hook(move |reed| {
            reed.with_history(Box::new(
                FileBackedHistory::with_file(10000, home.join(".gammaLoop_history")).unwrap(),
            ))
        })
    }
    let mut editor = repl.build();
    let refresh_repl_state = |editor: &mut repl::ClapEditor<Repl>,
                              session: &session::CliSession<'_>| {
        completion_state.update_from_session(session);
        let prompt_state_label = session.prompt_state_label();
        let pending_block_name = session.pending_commands_block_name();
        editor.set_prompt(make_repl_prompt(
            &prompt_state_label,
            pending_block_name.as_deref(),
        ));
    };

    loop {
        match editor.read_command() {
            ReadCommandOutput::Command(command, raw_input) => {
                match session
                    .execute_command(CommandHistory::new_with_raw(command.command, raw_input))
                {
                    Err(e) => {
                        eprintln!("{e:?}");
                    }
                    Ok(execution) => match execution.flow {
                        ControlFlow::Break(a) => {
                            refresh_repl_state(&mut editor, session);
                            *save_state = a;
                            break;
                        }
                        ControlFlow::Continue(()) => {
                            refresh_repl_state(&mut editor, session);
                        }
                    },
                }
            }
            ReadCommandOutput::EmptyLine => (),
            ReadCommandOutput::ClapError(e) => {
                e.print().unwrap();
            }
            ReadCommandOutput::ShlexError => {
                println!(
                    "{} input was not valid and could not be processed",
                    style("Error:").red().bold()
                );
            }
            ReadCommandOutput::ReedlineError(e) => {
                panic!("{e}");
            }
            ReadCommandOutput::CtrlC => {
                if session.dismiss_pending_commands_block("Ctrl-C") {
                    editor.set_prompt(make_repl_prompt(&session.prompt_state_label(), None));
                    continue;
                }
                continue;
            }
            ReadCommandOutput::CtrlD => {
                if session.dismiss_pending_commands_block("Ctrl-D") {
                    editor.set_prompt(make_repl_prompt(&session.prompt_state_label(), None));
                    continue;
                }
                break;
            }
        }
    }
}

fn generate_completion_script(shell: CompletionShell) -> String {
    let mut command = OneShot::command();
    let mut output = Vec::new();
    match shell {
        CompletionShell::Bash => {
            clap_complete::generate(Bash, &mut command, "gammaloop", &mut output)
        }
        CompletionShell::Elvish => {
            clap_complete::generate(Elvish, &mut command, "gammaloop", &mut output)
        }
        CompletionShell::Fish => {
            clap_complete::generate(Fish, &mut command, "gammaloop", &mut output)
        }
        CompletionShell::PowerShell => {
            clap_complete::generate(PowerShell, &mut command, "gammaloop", &mut output)
        }
        CompletionShell::Zsh => {
            clap_complete::generate(Zsh, &mut command, "gammaloop", &mut output)
        }
        CompletionShell::Nushell => {
            clap_complete::generate(Nushell, &mut command, "gammaloop", &mut output)
        }
    }
    let script = String::from_utf8(output).expect("clap completion script must be valid UTF-8");
    match shell {
        CompletionShell::Bash => patch_bash_completion_for_repo_wrapper(script),
        CompletionShell::Fish => {
            patch_fish_completion_for_repo_wrapper(normalize_fish_completion(script))
        }
        _ => script,
    }
}

fn patch_bash_completion_for_repo_wrapper(mut script: String) -> String {
    script.push_str(
        "\n# Support the repository wrapper script path as well.\n\
if [[ $(type -t _gammaloop) == function ]]; then\n\
    complete -F _gammaloop -o bashdefault -o default ./gammaloop\n\
fi\n",
    );
    script
}

fn patch_fish_completion_for_repo_wrapper(script: String) -> String {
    let mut wrapper_lines = Vec::new();
    for line in script.lines() {
        if line.starts_with("complete -c gammaloop") {
            wrapper_lines.push(line.replacen(
                "complete -c gammaloop",
                "complete -c './gammaloop'",
                1,
            ));
        }
    }
    if wrapper_lines.is_empty() {
        script
    } else {
        format!(
            "{script}\n# Support the repository wrapper script path as well.\n{}\n",
            wrapper_lines.join("\n")
        )
    }
}

fn normalize_fish_completion(script: String) -> String {
    let mut normalized = String::with_capacity(script.len());
    let bytes = script.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i..].starts_with(b"-a \"") {
            normalized.push_str("-a \"");
            i += 4;
            while i < bytes.len() {
                let b = bytes[i];
                if b == b'"' {
                    normalized.push('"');
                    i += 1;
                    break;
                }
                if b == b'\n' {
                    normalized.push(' ');
                } else {
                    normalized.push(b as char);
                }
                i += 1;
            }
        } else {
            normalized.push(bytes[i] as char);
            i += 1;
        }
    }
    normalized
}

pub mod commands;

pub enum Log {
    Level(LevelFilter),
    Format(LogFormat),
}

// impl FromStr for LogFormat {
//     type Err = Report;
//     fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
//         Ok(match s {
//             "long" => LogFormat::Long,
//             "short" => LogFormat::Short,
//             "min" => LogFormat::Min,
//             "none" => LogFormat::None,
//             _ => Err(eyre!("Invalid log format"))?,
//         })
//     }
// }

pub(crate) fn print_banner() {
    let spec = get_stderr_log_filter_label();
    println!(
        "\n{}{}\n",
        BANNER_ART.to_string().bold().blue(),
        banner_footer_line(&spec)
            .replace(GIT_VERSION, &GIT_VERSION.green().to_string())
            .replace(&spec, &spec.green().to_string())
            .bold()
            .blue(),
    );
}

pub fn write_schemas() -> Result<()> {
    let global_schema = schema_for!(GlobalSettings);
    let runtime_schema = schema_for!(RuntimeSettings);
    let runhistory_schema = schema_for!(RunHistory);
    let folder = get_schema_folder(false)?;

    let mut global_file = File::create(folder.join("global.json"))?;
    let mut runtime_file = File::create(folder.join("runtime.json"))?;
    let mut runhistory_file = File::create(folder.join("runhistory.json"))?;

    serde_json::to_writer_pretty(&mut global_file, &global_schema)
        .wrap_err("Could not write global schema")?;
    serde_json::to_writer_pretty(&mut runtime_file, &runtime_schema)
        .wrap_err("Could not write runtime schema")?;
    serde_json::to_writer_pretty(&mut runhistory_file, &runhistory_schema)
        .wrap_err("Could not write runhistory schema")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::{ffi::OsString, fs, path::PathBuf};

    use clap::Parser;
    use serde_json::Value as JsonValue;
    use tempfile::tempdir;

    use gammalooprs::{
        settings::{GlobalSettings, RuntimeSettings},
        utils::serde_utils::{ShowDefaultsGuard, SmartSerde},
    };

    use crate::commands::{save::SaveState, Duplicate, Remove};
    use crate::settings_tree::serialize_settings_with_defaults;
    use crate::state::ProcessRef;
    use crate::tracing::{
        configure_file_log_boot_mode, get_file_log_filter, get_stderr_log_filter,
        set_file_log_filter, set_file_log_filter_override, set_log_format_override, set_log_style,
        set_stderr_log_filter, set_stderr_log_filter_override,
    };

    use super::{
        generate_completion_script, CLISettings, CommandHistory, Commands, CompletionShell,
        LogFormat, LogLevel, OneShot, Repl, RunHistory, State, StateSettings,
        GLOBAL_SETTINGS_FILENAME, LOG_TEST_MUTEX,
    };
    use crate::state::CommandsBlock;

    const ENV_FILE_LOG_FILTER: &str = "GL_LOGFILE_FILTER";
    const ENV_DISPLAY_LOG_FILTER: &str = "GL_DISPLAY_FILTER";
    const ENV_ALL_LOG_FILTER: &str = "GL_ALL_LOG_FILTER";

    struct LogEnvGuard {
        file: Option<OsString>,
        display: Option<OsString>,
        all: Option<OsString>,
    }

    impl LogEnvGuard {
        fn capture() -> Self {
            Self {
                file: std::env::var_os(ENV_FILE_LOG_FILTER),
                display: std::env::var_os(ENV_DISPLAY_LOG_FILTER),
                all: std::env::var_os(ENV_ALL_LOG_FILTER),
            }
        }

        fn clear() {
            unsafe {
                std::env::remove_var(ENV_FILE_LOG_FILTER);
                std::env::remove_var(ENV_DISPLAY_LOG_FILTER);
                std::env::remove_var(ENV_ALL_LOG_FILTER);
            }
        }

        fn set(var: &str, value: &str) {
            unsafe {
                std::env::set_var(var, value);
            }
        }
    }

    impl Drop for LogEnvGuard {
        fn drop(&mut self) {
            fn restore_var(name: &str, value: &Option<OsString>) {
                unsafe {
                    if let Some(value) = value {
                        std::env::set_var(name, value);
                    } else {
                        std::env::remove_var(name);
                    }
                }
            }

            restore_var(ENV_FILE_LOG_FILTER, &self.file);
            restore_var(ENV_DISPLAY_LOG_FILTER, &self.display);
            restore_var(ENV_ALL_LOG_FILTER, &self.all);
        }
    }

    fn reset_tracing_state() {
        LogEnvGuard::clear();
        configure_file_log_boot_mode(false, None).unwrap();
        set_stderr_log_filter_override(None).unwrap();
        set_file_log_filter_override(None).unwrap();
        set_stderr_log_filter("info").unwrap();
        set_file_log_filter("off").unwrap();
        set_log_style(Default::default());
        set_log_format_override(None);
    }

    fn write_run_card(path: &PathBuf, state_folder: &str) {
        let run_history = RunHistory {
            cli_settings: CLISettings {
                state: StateSettings {
                    folder: state_folder.into(),
                    ..StateSettings::default()
                },
                ..CLISettings::default()
            },
            ..RunHistory::default()
        };
        fs::write(path, toml::to_string_pretty(&run_history).unwrap()).unwrap();
    }

    #[test]
    fn resolve_initial_state_folder_prefers_explicit_cli_value() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        write_run_card(&run_path, "./from_run_card");

        let mut one_shot =
            OneShot::try_parse_from(["gammaloop", run_path.to_string_lossy().as_ref()]).unwrap();
        one_shot.state_folder = "./from_cli".into();
        one_shot.state_folder_explicitly_set = true;

        assert_eq!(
            one_shot.resolve_initial_state_folder(None).unwrap(),
            PathBuf::from("./from_cli")
        );
    }

    #[test]
    fn resolve_initial_state_folder_uses_run_card_when_cli_not_explicit() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        write_run_card(&run_path, "./from_run_card");

        let one_shot =
            OneShot::try_parse_from(["gammaloop", run_path.to_string_lossy().as_ref()]).unwrap();
        let boot_run_history = one_shot.load_boot_run_history().unwrap();

        assert_eq!(
            one_shot
                .resolve_initial_state_folder(boot_run_history.as_ref())
                .unwrap(),
            PathBuf::from("./from_run_card")
        );
    }

    #[test]
    fn oneshot_accepts_clean_state_flag() {
        let parsed = OneShot::try_parse_from(["gammaloop", "--clean-state"]).unwrap();
        assert!(parsed.clean_state);
    }

    #[test]
    fn clean_state_removes_resolved_run_card_state_before_validation() {
        let temp = tempdir().unwrap();
        let state_path = temp.path().join("from_run_card");
        let run_path = temp.path().join("run.toml");
        write_run_card(&run_path, state_path.to_string_lossy().as_ref());

        fs::create_dir_all(&state_path).unwrap();
        fs::write(state_path.join("stale.txt"), "stale").unwrap();

        let mut one_shot = OneShot::try_parse_from([
            "gammaloop",
            "--clean-state",
            run_path.to_string_lossy().as_ref(),
        ])
        .unwrap();

        let (_state, run_history, cli_settings, runtime_settings, summary) =
            one_shot.load().unwrap();

        assert!(summary.is_none());
        assert!(run_history.commands.is_empty());
        assert_eq!(cli_settings.state.folder, state_path);
        assert_eq!(runtime_settings, RuntimeSettings::default());
        assert!(!cli_settings.state.folder.exists());
    }

    #[test]
    fn clean_state_rejects_read_only_state_mode() {
        let temp = tempdir().unwrap();
        let state_path = temp.path().join("read_only_state");
        fs::create_dir_all(&state_path).unwrap();
        fs::write(state_path.join("stale.txt"), "stale").unwrap();

        let mut one_shot = OneShot::try_parse_from([
            "gammaloop",
            "--clean-state",
            "--read-only-state",
            "-s",
            state_path.to_string_lossy().as_ref(),
        ])
        .unwrap();

        let err = match one_shot.load() {
            Ok(_) => panic!("read-only clean-state startup should fail"),
            Err(err) => err,
        };
        assert!(format!("{err:?}").contains("--read-only-state"));
        assert!(state_path.join("stale.txt").exists());
    }

    #[test]
    fn one_shot_run_history_persists_all_executed_commands() {
        let temp = tempdir().unwrap();
        let state_path = temp.path().join("state");
        let run_path = temp.path().join("boot.toml");
        let run_history = RunHistory {
            cli_settings: CLISettings {
                state: StateSettings {
                    folder: state_path.clone(),
                    ..StateSettings::default()
                },
                ..CLISettings::default()
            },
            command_blocks: vec![
                CommandsBlock {
                    name: "cmdBlockA".to_string(),
                    commands: vec![CommandHistory::from_raw_string(
                        "set global kv global.logfile_directive=error",
                    )
                    .unwrap()],
                },
                CommandsBlock {
                    name: "cmdBlockB".to_string(),
                    commands: vec![CommandHistory::from_raw_string(
                        "set default-runtime kv general.mu_r=12.0",
                    )
                    .unwrap()],
                },
            ],
            commands: vec![CommandHistory::from_raw_string(
                "set global kv global.display_directive=warn",
            )
            .unwrap()],
            ..RunHistory::default()
        };
        fs::write(&run_path, run_history.to_toml_string(true).unwrap()).unwrap();

        let one_shot = OneShot::try_parse_from([
            "gammaloop",
            run_path.to_string_lossy().as_ref(),
            "run",
            "cmdBlockA",
            "cmdBlockB",
            "-c",
            "set default-runtime kv general.m_uv=7.0; quit -o",
        ])
        .unwrap();

        one_shot
            .run(
                "run cmdBlockA cmdBlockB -c \"set default-runtime kv general.m_uv=7.0; quit -o\""
                    .to_string(),
            )
            .unwrap();

        let persisted = RunHistory::load(state_path.join("run.toml")).unwrap();
        let commands = persisted
            .commands
            .iter()
            .map(crate::session::display_command)
            .collect::<Vec<_>>();
        assert_eq!(
            commands,
            vec![
                "set global kv global.display_directive=warn",
                "run cmdBlockA cmdBlockB -c 'set default-runtime kv general.m_uv=7.0'",
            ]
        );
    }

    #[test]
    fn auto_read_only_boot_settings_mismatch_ignores_quit_override_for_active_state() {
        let temp = tempdir().unwrap();
        let state_path = temp.path().join("state");
        let boot_path = temp.path().join("boot.toml");

        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        cli_settings.state.folder = state_path.clone();
        let mut saved_run_history = RunHistory::default();
        saved_run_history.default_runtime_settings.general.mu_r = 11.0;
        SaveState {
            override_state: Some(true),
            ..Default::default()
        }
        .save(
            &mut state,
            &saved_run_history,
            &RuntimeSettings::default(),
            &cli_settings,
        )
        .unwrap();
        let saved_run_toml = fs::read_to_string(state_path.join("run.toml")).unwrap();

        let mut boot_runtime = RuntimeSettings::default();
        boot_runtime.general.mu_r = 29.0;
        let boot_run_history = RunHistory {
            default_runtime_settings: boot_runtime,
            ..RunHistory::default()
        };
        fs::write(
            &boot_path,
            boot_run_history
                .to_toml_string(cli_settings.try_strings)
                .unwrap(),
        )
        .unwrap();

        let mut one_shot = OneShot::try_parse_from([
            "gammaloop",
            "-s",
            state_path.to_string_lossy().as_ref(),
            boot_path.to_string_lossy().as_ref(),
            "quit",
            "-o",
        ])
        .unwrap();
        one_shot.state_folder_explicitly_set = true;

        one_shot.run("quit -o".to_string()).unwrap();

        assert_eq!(
            fs::read_to_string(state_path.join("run.toml")).unwrap(),
            saved_run_toml
        );
    }

    #[test]
    fn oneshot_rejects_removed_fresh_state_flag() {
        assert!(OneShot::try_parse_from(["gammaloop", "--fresh-state"]).is_err());
        assert!(OneShot::try_parse_from(["gammaloop", "-f"]).is_err());
    }

    #[test]
    fn oneshot_parses_run_subcommand_block_names() {
        let parsed =
            OneShot::try_parse_from(["gammaloop", "run", "generation", "integration"]).unwrap();
        let Some(Commands::Run(run)) = parsed.command else {
            panic!("expected run command");
        };
        assert_eq!(
            run.selected_block_names(),
            &["generation".to_string(), "integration".to_string()][..]
        );
    }

    #[test]
    fn oneshot_treats_subcommand_names_as_run_block_names() {
        let parsed = OneShot::try_parse_from([
            "gammaloop",
            "run",
            "generate",
            "integrate_euclidean",
            "quit",
        ])
        .unwrap();
        let Some(Commands::Run(run)) = parsed.command else {
            panic!("expected run command");
        };
        assert_eq!(
            run.selected_block_names(),
            &[
                "generate".to_string(),
                "integrate_euclidean".to_string(),
                "quit".to_string()
            ][..]
        );
    }

    #[test]
    fn oneshot_requires_explicit_generate_mode() {
        assert!(OneShot::try_parse_from(["gammaloop", "generate", "g", "g", ">", "h"]).is_err());
    }

    #[test]
    fn oneshot_parses_subcommand_without_run_card() {
        let parsed = OneShot::try_parse_from(["gammaloop", "display", "integrand"]).unwrap();
        assert!(matches!(parsed.command, Some(Commands::Display(_))));
    }

    #[test]
    fn oneshot_accepts_boot_card_positional() {
        let parsed = OneShot::try_parse_from(["gammaloop", "card.toml"]).unwrap();
        assert_eq!(parsed.boot_commands_path, Some(PathBuf::from("card.toml")));
    }

    #[test]
    fn oneshot_allows_boot_card_positional_with_run_subcommand() {
        let parsed = OneShot::try_parse_from([
            "gammaloop",
            "card.toml",
            "-s",
            "./GL_OUTPUT/triangle",
            "run",
            "-c",
            "quit -o",
        ])
        .unwrap();

        assert_eq!(parsed.boot_commands_path, Some(PathBuf::from("card.toml")));
        assert_eq!(parsed.state_folder, PathBuf::from("./GL_OUTPUT/triangle"));
        assert!(matches!(parsed.command, Some(Commands::Run(_))));
    }

    #[test]
    fn oneshot_rejects_removed_boot_card_flag() {
        assert!(OneShot::try_parse_from(["gammaloop", "-c", "card.toml"]).is_err());
    }

    #[test]
    fn oneshot_accepts_trace_logs_short_flag() {
        let parsed =
            OneShot::try_parse_from(["gammaloop", "-t", "trace.log", "card.toml"]).unwrap();
        assert_eq!(parsed.trace_logs_filename.as_deref(), Some("trace.log"));
        assert_eq!(parsed.boot_commands_path, Some(PathBuf::from("card.toml")));
    }

    #[test]
    fn oneshot_leaves_logging_prefix_unspecified_by_default() {
        let parsed = OneShot::try_parse_from(["gammaloop"]).unwrap();
        assert_eq!(parsed.logging_prefix, None);
    }

    #[test]
    fn oneshot_accepts_logging_prefix_long_flag_override() {
        let parsed = OneShot::try_parse_from(["gammaloop", "--logging-prefix", "long"]).unwrap();
        assert_eq!(parsed.logging_prefix, Some(LogFormat::Long));
    }

    #[test]
    fn oneshot_accepts_logging_prefix_full_flag_override() {
        let parsed = OneShot::try_parse_from(["gammaloop", "--logging-prefix", "full"]).unwrap();
        assert_eq!(parsed.logging_prefix, Some(LogFormat::Full));
    }

    #[test]
    fn oneshot_rejects_removed_logging_prefix_underscore_flag() {
        assert!(OneShot::try_parse_from(["gammaloop", "--logging_prefix", "long"]).is_err());
    }

    #[test]
    fn oneshot_accepts_settings_global_short_flag() {
        let parsed = OneShot::try_parse_from(["gammaloop", "-g", "global.toml"]).unwrap();
        assert_eq!(
            parsed.settings_global_path,
            Some(PathBuf::from("global.toml"))
        );
    }

    #[test]
    fn oneshot_accepts_settings_runtime_defaults_short_flag() {
        let parsed = OneShot::try_parse_from(["gammaloop", "-r", "runtime.toml"]).unwrap();
        assert_eq!(
            parsed.settings_runtime_defaults_path,
            Some(PathBuf::from("runtime.toml"))
        );
    }

    #[test]
    fn oneshot_accepts_logfile_level_and_read_only_state_flags() {
        let parsed =
            OneShot::try_parse_from(["gammaloop", "--read-only-state", "--logfile-level", "off"])
                .unwrap();
        assert!(parsed.read_only_state);
        assert_eq!(parsed.logfile_level, Some(LogLevel::Off));
    }

    #[test]
    fn configure_startup_tracing_uses_global_defaults_without_overrides() {
        let _guard = LOG_TEST_MUTEX.lock().unwrap_or_else(|err| err.into_inner());
        let _env = LogEnvGuard::capture();
        reset_tracing_state();

        let cli = OneShot::try_parse_from(["gammaloop"]).unwrap();
        let cli_settings = CLISettings::default();

        cli.configure_startup_tracing(&cli_settings).unwrap();

        assert_eq!(get_stderr_log_filter(), "info");
        assert_eq!(get_file_log_filter(), "off");
    }

    #[test]
    fn configure_startup_tracing_prefers_all_env_override_over_specific_envs() {
        let _guard = LOG_TEST_MUTEX.lock().unwrap_or_else(|err| err.into_inner());
        let _env = LogEnvGuard::capture();
        reset_tracing_state();

        LogEnvGuard::set(ENV_ALL_LOG_FILTER, "trace");
        LogEnvGuard::set(ENV_DISPLAY_LOG_FILTER, "warn");
        LogEnvGuard::set(ENV_FILE_LOG_FILTER, "error");

        let cli = OneShot::try_parse_from(["gammaloop"]).unwrap();
        let mut cli_settings = CLISettings::default();
        cli_settings.global.display_directive = "info".into();
        cli_settings.global.logfile_directive = "off".into();

        cli.configure_startup_tracing(&cli_settings).unwrap();

        assert_eq!(get_stderr_log_filter(), "trace");
        assert_eq!(get_file_log_filter(), "trace");
    }

    #[test]
    fn configure_startup_tracing_cli_overrides_supersede_settings_and_env() {
        let _guard = LOG_TEST_MUTEX.lock().unwrap_or_else(|err| err.into_inner());
        let _env = LogEnvGuard::capture();
        reset_tracing_state();

        LogEnvGuard::set(ENV_DISPLAY_LOG_FILTER, "warn");
        LogEnvGuard::set(ENV_FILE_LOG_FILTER, "error");

        let cli = OneShot::try_parse_from(["gammaloop", "-l", "debug", "--logfile-level", "trace"])
            .unwrap();
        let mut cli_settings = CLISettings::default();
        cli_settings.global.display_directive = "info".into();
        cli_settings.global.logfile_directive = "off".into();

        cli.configure_startup_tracing(&cli_settings).unwrap();

        assert_eq!(
            get_stderr_log_filter(),
            "gammaloop_api=debug,gammalooprs=debug"
        );
        assert_eq!(
            get_file_log_filter(),
            "gammaloop_api=trace,gammalooprs=trace"
        );
    }

    #[test]
    fn configure_startup_tracing_read_only_state_forces_logfile_off() {
        let _guard = LOG_TEST_MUTEX.lock().unwrap_or_else(|err| err.into_inner());
        let _env = LogEnvGuard::capture();
        reset_tracing_state();

        let cli = OneShot::try_parse_from(["gammaloop", "--read-only-state"]).unwrap();
        let mut cli_settings = CLISettings::default();
        cli_settings.global.logfile_directive = "debug".into();

        cli.configure_startup_tracing(&cli_settings).unwrap();

        assert_eq!(get_stderr_log_filter(), "info");
        assert_eq!(get_file_log_filter(), "gammaloop_api=off,gammalooprs=off");
    }

    #[test]
    fn oneshot_accepts_completions_flag() {
        let parsed = OneShot::try_parse_from(["gammaloop", "--completions", "bash"]).unwrap();
        assert_eq!(parsed.completions, Some(CompletionShell::Bash));
    }

    #[test]
    fn bash_completion_script_binds_repo_wrapper_path() {
        let script = generate_completion_script(CompletionShell::Bash);
        assert!(script.contains("complete -F _gammaloop -o bashdefault -o default gammaloop"));
        assert!(script.contains("complete -F _gammaloop -o bashdefault -o default ./gammaloop"));
    }

    #[test]
    fn fish_completion_script_binds_repo_wrapper_path() {
        let script = generate_completion_script(CompletionShell::Fish);
        assert!(script.contains("complete -c gammaloop"));
        assert!(script.contains("complete -c './gammaloop'"));
        assert!(!script.contains("-a \"true\\t''\nfalse\\t''\""));
    }

    #[test]
    fn oneshot_accepts_nushell_completions_flag() {
        let parsed = OneShot::try_parse_from(["gammaloop", "--completions", "nushell"]).unwrap();
        assert_eq!(parsed.completions, Some(CompletionShell::Nushell));
    }

    #[test]
    fn nushell_completion_script_is_exportable_module() {
        let script = generate_completion_script(CompletionShell::Nushell);
        assert!(script.contains("module completions"));
        assert!(script.contains("export use completions *"));
    }

    #[test]
    fn repl_parses_normalized_integrand_selector_flags() {
        let parsed =
            Repl::try_parse_from(["gammaloop", "integrate", "-p", "triangle", "-i", "LO"]).unwrap();

        let Commands::Integrate(integrate) = parsed.command else {
            panic!("expected integrate command");
        };
        assert_eq!(
            integrate.process,
            vec![ProcessRef::Unqualified("triangle".to_string())]
        );
        assert_eq!(integrate.integrand_name, vec!["LO".to_string()]);
    }

    #[test]
    fn repl_parses_duplicate_integrand_output_selectors() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "duplicate",
            "integrand",
            "-p",
            "box",
            "-i",
            "scalar_box",
            "--output_process_name",
            "box_copy",
            "--output_integrand_name",
            "scalar_box_copy",
        ])
        .unwrap();

        let Commands::Duplicate(Duplicate::Integrand(command)) = parsed.command else {
            panic!("expected duplicate integrand command");
        };
        assert_eq!(
            command.process,
            Some(ProcessRef::Unqualified("box".to_string()))
        );
        assert_eq!(command.integrand_name, Some("scalar_box".to_string()));
        assert_eq!(command.output_process_name, "box_copy");
        assert_eq!(command.output_integrand_name, "scalar_box_copy");
    }

    #[test]
    fn repl_parses_inspect_process_and_point_with_normalized_flags() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "inspect",
            "-p",
            "triangle",
            "-i",
            "LO",
            "-x",
            "0.1",
            "0.2",
            "0.3",
        ])
        .unwrap();

        let Commands::Inspect(inspect) = parsed.command else {
            panic!("expected inspect command");
        };
        assert_eq!(
            inspect.process,
            Some(ProcessRef::Unqualified("triangle".to_string()))
        );
        assert_eq!(inspect.integrand_name, Some("LO".to_string()));
        assert_eq!(inspect.point, vec![0.1, 0.2, 0.3]);
    }

    #[test]
    fn repl_parses_inspect_discrete_dims_from_one_flag_occurrence() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "inspect",
            "-p",
            "triangle",
            "-i",
            "LO",
            "-m",
            "-x",
            "0.1",
            "0.2",
            "0.3",
            "-d",
            "0",
            "0",
        ])
        .unwrap();

        let Commands::Inspect(inspect) = parsed.command else {
            panic!("expected inspect command");
        };
        assert_eq!(inspect.discrete_dim, vec![0, 0]);
    }

    #[test]
    fn repl_parses_inspect_discrete_dims_from_repeated_flags() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "inspect",
            "-p",
            "triangle",
            "-i",
            "LO",
            "-m",
            "-x",
            "0.1",
            "0.2",
            "0.3",
            "-d",
            "0",
            "-d",
            "0",
        ])
        .unwrap();

        let Commands::Inspect(inspect) = parsed.command else {
            panic!("expected inspect command");
        };
        assert_eq!(inspect.discrete_dim, vec![0, 0]);
    }

    #[test]
    fn repl_parses_inspect_graph_and_orientation_ids() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "inspect",
            "-p",
            "triangle",
            "-i",
            "LO",
            "-m",
            "-x",
            "0.1",
            "0.2",
            "0.3",
            "--graph-id",
            "2",
            "--orientation-id",
            "1",
        ])
        .unwrap();

        let Commands::Inspect(inspect) = parsed.command else {
            panic!("expected inspect command");
        };
        assert_eq!(inspect.graph_id, Some(2));
        assert_eq!(inspect.orientation_id, Some(1));
        assert!(inspect.discrete_dim.is_empty());
    }

    #[test]
    fn repl_parses_remove_processes_with_normalized_selectors() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "remove",
            "processes",
            "-p",
            "triangle",
            "-i",
            "LO",
        ])
        .unwrap();

        let Commands::Remove(Remove::Processes {
            process,
            integrand_name,
        }) = parsed.command
        else {
            panic!("expected remove processes command");
        };
        assert_eq!(
            process,
            Some(ProcessRef::Unqualified("triangle".to_string()))
        );
        assert_eq!(integrand_name, Some("LO".to_string()));
    }

    #[test]
    fn settings_file_overrides_apply_without_touching_state_folder() {
        let temp = tempdir().unwrap();
        let global_path = temp.path().join("global.toml");
        let runtime_path = temp.path().join("runtime.toml");

        let mut global_file_settings = CLISettings::default();
        global_file_settings.state.folder = "./ignored".into();
        global_file_settings.global.display_directive = "warn".into();
        fs::write(
            &global_path,
            toml::to_string_pretty(&global_file_settings).unwrap(),
        )
        .unwrap();

        let mut runtime_file_settings = RuntimeSettings::default();
        runtime_file_settings.general.mu_r = 37.0;
        fs::write(
            &runtime_path,
            toml::to_string_pretty(&runtime_file_settings).unwrap(),
        )
        .unwrap();

        let cli = OneShot::try_parse_from([
            "gammaloop",
            "--settings-global",
            global_path.to_string_lossy().as_ref(),
            "--settings-runtime-defaults",
            runtime_path.to_string_lossy().as_ref(),
        ])
        .unwrap();

        let overrides = cli.load_settings_file_overrides().unwrap();
        let mut cli_settings = CLISettings::default();
        cli_settings.state.folder = "./keep".into();
        cli_settings.global = GlobalSettings::default();
        let mut runtime_settings = RuntimeSettings::default();

        overrides
            .apply(&mut cli_settings, &mut runtime_settings)
            .unwrap();

        assert_eq!(cli_settings.state.folder, PathBuf::from("./keep"));
        assert_eq!(cli_settings.global.display_directive, "warn");
        assert_eq!(runtime_settings.general.mu_r, 37.0);
    }

    #[test]
    fn state_name_empty_string_deserializes_to_none() {
        let settings: CLISettings = toml::from_str(
            r#"
                [state]
                folder = "./named_state"
                name = ""
            "#,
        )
        .unwrap();

        assert_eq!(settings.state.folder, PathBuf::from("./named_state"));
        assert_eq!(settings.state.name, None);
    }

    #[test]
    fn state_prompt_label_prefers_name_and_falls_back_to_folder() {
        let named = StateSettings {
            folder: "./named_state".into(),
            name: Some("demo".into()),
        };
        assert_eq!(named.prompt_label(), "demo");

        let unnamed = StateSettings {
            folder: "./named_state".into(),
            name: Some("   ".into()),
        };
        assert_eq!(unnamed.prompt_label(), "./named_state");
    }

    #[test]
    fn state_name_is_present_in_completion_serialization() {
        let serialized =
            serialize_settings_with_defaults(&CLISettings::default(), "CLI settings").unwrap();
        let state = serialized
            .get("state")
            .and_then(JsonValue::as_object)
            .expect("state object must be present");

        assert_eq!(
            state.get("name").and_then(JsonValue::as_str),
            Some(""),
            "state.name should stay visible to settings completion even when unset"
        );
    }

    #[test]
    fn global_settings_defaults_stay_visible_in_completion_serialization() {
        let serialized =
            serialize_settings_with_defaults(&CLISettings::default(), "CLI settings").unwrap();
        let global = serialized
            .get("global")
            .and_then(JsonValue::as_object)
            .expect("global object must be present");

        assert_eq!(
            global.get("display_directive").and_then(JsonValue::as_str),
            Some("info")
        );
        assert_eq!(
            global.get("logfile_directive").and_then(JsonValue::as_str),
            Some("off")
        );
        assert_eq!(
            global
                .get("generation")
                .and_then(JsonValue::as_object)
                .and_then(|generation| generation.get("compile"))
                .and_then(JsonValue::as_object)
                .and_then(|compile| compile.get("compiler"))
                .and_then(JsonValue::as_str),
            Some(gammalooprs::settings::global::default_external_compiler())
        );
    }

    #[test]
    fn runtime_settings_custom_defaults_stay_visible_in_completion_serialization() {
        let serialized = serialize_settings_with_defaults(
            &RuntimeSettings::default(),
            "default runtime settings",
        )
        .unwrap();
        let serialized_text = serde_json::to_string(&serialized).unwrap();
        assert!(serialized_text.contains("large_deformation_check"));
    }

    #[test]
    fn saved_global_settings_file_keeps_default_directives_visible() {
        let temp = tempdir().unwrap();
        let _show_defaults_guard = ShowDefaultsGuard::new(true);
        CLISettings::default()
            .to_file(temp.path().join(GLOBAL_SETTINGS_FILENAME), true)
            .unwrap();
        let saved = fs::read_to_string(temp.path().join(GLOBAL_SETTINGS_FILENAME)).unwrap();

        assert!(saved.contains("display_directive = \"info\""), "{saved}");
        assert!(saved.contains("logfile_directive = \"off\""), "{saved}");
        assert!(saved.contains("log_format = \"Long\""), "{saved}");
        assert!(
            saved.contains(&format!(
                "compiler = \"{}\"",
                gammalooprs::settings::global::default_external_compiler()
            )),
            "{saved}"
        );
    }

    #[test]
    fn settings_file_overrides_replace_boot_card_settings() {
        let mut run_history = RunHistory::default();
        run_history.cli_settings.global.display_directive = "error".into();
        run_history.default_runtime_settings.general.mu_r = 11.0;

        let overridden_global = GlobalSettings {
            display_directive: "warn".into(),
            ..Default::default()
        };
        let mut overridden_runtime = RuntimeSettings::default();
        overridden_runtime.general.mu_r = 29.0;

        let overrides = super::SettingsFileOverrides {
            global: Some(overridden_global),
            default_runtime_settings: Some(overridden_runtime),
        };

        let overridden = overrides.apply_to_run_history(&run_history);

        assert_eq!(overridden.cli_settings.global.display_directive, "warn");
        assert_eq!(overridden.default_runtime_settings.general.mu_r, 29.0);
        assert_eq!(run_history.cli_settings.global.display_directive, "error");
        assert_eq!(run_history.default_runtime_settings.general.mu_r, 11.0);
    }

    #[test]
    fn load_global_settings_file_uses_default_when_missing() {
        let temp = tempdir().unwrap();

        let mut new_settings = CLISettings::default();
        new_settings.global.display_directive = "warn".into();
        fs::write(
            temp.path().join(GLOBAL_SETTINGS_FILENAME),
            toml::to_string_pretty(&new_settings).unwrap(),
        )
        .unwrap();

        let loaded = OneShot::load_global_settings_file(temp.path()).unwrap();
        assert_eq!(loaded.global.display_directive, "warn");

        fs::remove_file(temp.path().join(GLOBAL_SETTINGS_FILENAME)).unwrap();
        let loaded = OneShot::load_global_settings_file(temp.path()).unwrap();
        assert_eq!(loaded, CLISettings::default());
        assert_eq!(loaded.global.log_style.log_format, LogFormat::Long);
    }

    #[test]
    fn load_treats_logs_only_folder_as_blank_state() {
        let temp = tempdir().unwrap();
        std::fs::create_dir_all(temp.path().join("logs")).unwrap();
        std::fs::write(temp.path().join("logs").join("gammalog.jsonl"), "").unwrap();

        let mut cli = OneShot::new_test(temp.path().to_path_buf());
        let (_state, run_history, cli_settings, runtime_settings, summary) = cli.load().unwrap();

        assert!(summary.is_none());
        assert!(run_history.commands.is_empty());
        assert_eq!(cli_settings.state.folder, temp.path().to_path_buf());
        assert_eq!(runtime_settings, RuntimeSettings::default());
    }

    #[test]
    fn load_treats_unmanifested_state_folder_as_blank_state() {
        let temp = tempdir().unwrap();
        std::fs::create_dir_all(temp.path().join("processes").join("amplitudes")).unwrap();

        let mut cli = OneShot::new_test(temp.path().to_path_buf());
        let (_state, run_history, cli_settings, runtime_settings, summary) = cli.load().unwrap();

        assert!(summary.is_none());
        assert!(run_history.commands.is_empty());
        assert_eq!(cli_settings.state.folder, temp.path().to_path_buf());
        assert_eq!(runtime_settings, RuntimeSettings::default());
    }
}
