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
use clap::parser::ValueSource;
use clap::{CommandFactory, FromArgMatches, Parser};
use commands::save::SaveState;
use commands::Commands;

use gammalooprs::utils::serde_utils::{SerdeFileError, SHOWDEFAULTS};
use reedline::FileBackedHistory;
use repl::ClapEditor;
use repl::ReadCommandOutput;
use session::{CliSession, CliSessionState};
use tracing::{get_stderr_log_filter, set_log_format_override};

// use clap_repl::{
//     reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
//     ClapEditor, ReadCommandOutput,
// };

use color_eyre::Result;
use colored::Colorize;
use console::{measure_text_width, style};
use dirs::home_dir;
use eyre::Context;
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

use state::{CommandHistory, RunHistory, State, SyncSettings};
use std::{borrow::Cow, ffi::OsString, path::Path, path::PathBuf, sync::atomic::Ordering};
use std::{fs::File, ops::ControlFlow, time::Duration, time::Instant};
use symbolica::activate_oem_license;
use walkdir::WalkDir;

// use tracing::LogLevel;
mod command_parser;
pub(crate) mod completion;
#[cfg(feature = "python_api")]
pub mod python;
pub mod repl;
pub mod session;
pub(crate) mod settings_tree;

pub mod state;
pub mod templates;
pub mod tracing;

pub(crate) const GLOBAL_SETTINGS_FILENAME: &str = "global_settings.toml";
pub(crate) const DEFAULT_RUNTIME_SETTINGS_FILENAME: &str = "default_runtime_settings.toml";
const SETTINGS_GLOBAL_SHORTCUT: &str = "-sg";
const SETTINGS_RUNTIME_DEFAULTS_SHORTCUT: &str = "-sr";
const BANNER_ART: &str = r"              ██         ▄████████▄  ▄████████▄  ██████████▄
  ██          ▀▀         ▀▀      ▀▀  ▀▀      ▀▀  ▀▀       ██
    ██  ▄██████████████████████████████████████████████████▀
     ▀██▀     ▄▄         ▄▄      ▄▄  ▄▄      ▄▄  ▄▄
    ██  ██    █████████  ▀████████▀  ▀████████▀  ██
   ██    ██
   ██    ██
";

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
    /// Don't try to load state, just start with a new one
    #[arg(short = 'f', long, default_value_t = false)]
    pub fresh_state: bool,

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

    /// Type of prefix for the logging format
    #[arg(short = 'p', long = "logging_prefix")]
    logging_prefix: Option<LogFormat>,

    /// Path to a global settings TOML file. Also accepts `-sg`.
    #[arg(long = "settings-global", value_hint = clap::ValueHint::FilePath)]
    settings_global_path: Option<PathBuf>,

    /// Path to a default runtime settings TOML file. Also accepts `-sr`.
    #[arg(
        long = "settings-runtime-defaults",
        value_hint = clap::ValueHint::FilePath
    )]
    settings_runtime_defaults_path: Option<PathBuf>,

    /// Try to serialize using strings when saving run history
    #[arg(long)]
    no_try_strings: bool,

    // /// Debug level
    // #[arg(short = 'd', long, value_enum, default_value_t = LogLevel::Info)]
    // debug_level: LogLevel,
    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, PartialEq)]
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
        }
    }
}

impl CLISettings {
    pub fn override_with(&mut self, cli: &OneShot) {
        self.try_strings = !cli.no_try_strings;

        self.override_state = cli.override_state;

        self.state.folder = cli.state_folder.clone();
    }
}

impl SmartSerde for CLISettings {}

pub struct Parsed {
    pub cli: OneShot,
    pub input_string: String,
    pub matches: clap::ArgMatches,
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
    elapsed: Duration,
    serialized_size_bytes: Option<u64>,
    total_graphs: usize,
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
    let spec = get_stderr_log_filter();
    let banner_width = banner_width(&spec);
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
        CLISettings {
            try_strings: !self.no_try_strings,
            override_state: self.override_state,
            state: StateSettings {
                folder: self.state_folder.clone(),
                ..StateSettings::default()
            },
            global,
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
            level: Some(LogLevel::Info),
            logging_prefix: None,
            settings_global_path: None,
            settings_runtime_defaults_path: None,
            no_try_strings: false,
            fresh_state: false,
            trace_logs_filename: None,
        }
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

    fn normalize_cli_shortcuts(argv: &[OsString]) -> Vec<OsString> {
        argv.iter()
            .map(|arg| {
                if arg == std::ffi::OsStr::new(SETTINGS_GLOBAL_SHORTCUT) {
                    OsString::from("--settings-global")
                } else if arg == std::ffi::OsStr::new(SETTINGS_RUNTIME_DEFAULTS_SHORTCUT) {
                    OsString::from("--settings-runtime-defaults")
                } else if let Some(arg) = arg.to_str() {
                    if let Some(value) = arg.strip_prefix("-sg=") {
                        OsString::from(format!("--settings-global={value}"))
                    } else if let Some(value) = arg.strip_prefix("-sr=") {
                        OsString::from(format!("--settings-runtime-defaults={value}"))
                    } else {
                        arg.into()
                    }
                } else {
                    arg.clone()
                }
            })
            .collect()
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
        let argv: Vec<OsString> = std::env::args_os().collect();
        let normalized_argv = Self::normalize_cli_shortcuts(&argv);

        // Build a Command (same as derive(Parser)) and get matches
        let mut cmd = <OneShot as CommandFactory>::command();
        let matches = cmd.clone().try_get_matches_from(&normalized_argv)?;

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

    pub fn load(
        &mut self,
    ) -> Result<(
        State,
        RunHistory,
        CLISettings,
        RuntimeSettings,
        Option<StateLoadSummary>,
    )> {
        let state_exists = self.state_folder.exists();

        let (state, run_history, cli_settings, default_runtime_settings, load_summary) =
            if !self.fresh_state && state_exists {
                let load_started = Instant::now();
                let state = State::load(
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

                let mut global = match Self::load_global_settings_file(&self.state_folder) {
                    Ok(a) => a,
                    Err(e) => return Err(e),
                };

                global.override_with(self);

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
                    global,
                    default_runtime,
                    Some(load_summary),
                )
            } else {
                if self.fresh_state && state_exists {
                    info!(
                        "{} {}",
                        "Starting fresh state in".blue(),
                        self.state_folder.display().to_string().green()
                    );
                } else {
                    info!(
                        "{} {}",
                        "Initializing new state in".blue(),
                        self.state_folder.display().to_string().green()
                    );
                }

                (
                    State::new(self.state_folder.clone(), self.trace_logs_filename.clone()),
                    RunHistory::default(),
                    self.new_cli_settings(GlobalSettings::default()),
                    RuntimeSettings::default(),
                    None,
                )
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

    pub fn run(mut self, raw: String) -> Result<()> {
        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        };
        initialise()?;
        set_log_format_override(self.logging_prefix);

        let boot_run_history = self.load_boot_run_history()?;
        let settings_file_overrides = self.load_settings_file_overrides()?;
        self.state_folder = self.resolve_initial_state_folder(boot_run_history.as_ref())?;

        let (
            mut state,
            mut run_history,
            mut cli_settings,
            mut default_runtime_settings,
            state_load_summary,
        ) = self.load()?;
        settings_file_overrides.apply(&mut cli_settings, &mut default_runtime_settings)?;

        let mut save_state = SaveState::default();
        let mut session_state = CliSessionState::default();
        let mut session = CliSession::new(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut default_runtime_settings,
            &mut session_state,
        );

        let mut boot_requested_exit = false;
        if let Some(boot_run_history) = boot_run_history.as_ref() {
            let effective_boot_run_history =
                settings_file_overrides.apply_to_run_history(boot_run_history);
            if let ControlFlow::Break(a) =
                session.apply_boot_run_history(&effective_boot_run_history)?
            {
                save_state = a;
                boot_requested_exit = true;
            }
        }

        if !boot_requested_exit {
            if let Some(a) = self.command.take() {
                let command = if raw.trim().is_empty() {
                    CommandHistory::new(a)
                } else {
                    CommandHistory::new_with_raw(a, raw)
                };
                let flow = session.execute_top_level(command)?;

                if let ControlFlow::Break(a) = flow {
                    save_state = a;
                }
            } else {
                print_banner();
                if let Some(summary) = state_load_summary.as_ref() {
                    print_state_load_summary(summary);
                }
                // 2. Build the REPL – clap‑repl takes ownership and configures rustyline.
                let completion_state = repl::new_completion_state();
                repl::set_commands_block_names(
                    &completion_state,
                    session.current_commands_block_names(),
                );
                repl::set_process_entries(&completion_state, session.current_process_entries());
                repl::set_model_parameter_names(
                    &completion_state,
                    session.current_model_parameter_names(),
                );
                repl::set_model_particle_names(
                    &completion_state,
                    session.current_model_particle_names(),
                );
                repl::set_model_coupling_names(
                    &completion_state,
                    session.current_model_coupling_names(),
                );
                repl::set_model_vertex_names(
                    &completion_state,
                    session.current_model_vertex_names(),
                );
                let mut repl = ClapEditor::<Repl>::builder()
                    .with_prompt(make_repl_prompt(&session.prompt_state_label(), None))
                    .with_completion_state(completion_state.clone());

                if let Some(home) = home_dir() {
                    repl = repl.with_editor_hook(move |reed| {
                        reed.with_history(Box::new(
                            FileBackedHistory::with_file(10000, home.join(".gammaLoop_history"))
                                .unwrap(),
                        ))
                    })
                }
                let mut r = repl.build();

                loop {
                    match r.read_command() {
                        ReadCommandOutput::Command(command, raw_input) => {
                            match session.execute_top_level(CommandHistory::new_with_raw(
                                command.command,
                                raw_input,
                            )) {
                                Err(e) => {
                                    eprintln!("{e:?}");
                                }
                                Ok(ControlFlow::Break(a)) => {
                                    repl::set_commands_block_names(
                                        &completion_state,
                                        session.current_commands_block_names(),
                                    );
                                    repl::set_process_entries(
                                        &completion_state,
                                        session.current_process_entries(),
                                    );
                                    repl::set_model_parameter_names(
                                        &completion_state,
                                        session.current_model_parameter_names(),
                                    );
                                    repl::set_model_particle_names(
                                        &completion_state,
                                        session.current_model_particle_names(),
                                    );
                                    repl::set_model_coupling_names(
                                        &completion_state,
                                        session.current_model_coupling_names(),
                                    );
                                    repl::set_model_vertex_names(
                                        &completion_state,
                                        session.current_model_vertex_names(),
                                    );
                                    let prompt_state_label = session.prompt_state_label();
                                    let pending_block_name = session.pending_commands_block_name();
                                    r.set_prompt(make_repl_prompt(
                                        &prompt_state_label,
                                        pending_block_name.as_deref(),
                                    ));
                                    save_state = a;
                                    break;
                                }
                                Ok(ControlFlow::Continue(())) => {
                                    repl::set_commands_block_names(
                                        &completion_state,
                                        session.current_commands_block_names(),
                                    );
                                    repl::set_process_entries(
                                        &completion_state,
                                        session.current_process_entries(),
                                    );
                                    repl::set_model_parameter_names(
                                        &completion_state,
                                        session.current_model_parameter_names(),
                                    );
                                    repl::set_model_particle_names(
                                        &completion_state,
                                        session.current_model_particle_names(),
                                    );
                                    repl::set_model_coupling_names(
                                        &completion_state,
                                        session.current_model_coupling_names(),
                                    );
                                    repl::set_model_vertex_names(
                                        &completion_state,
                                        session.current_model_vertex_names(),
                                    );
                                    let prompt_state_label = session.prompt_state_label();
                                    let pending_block_name = session.pending_commands_block_name();
                                    r.set_prompt(make_repl_prompt(
                                        &prompt_state_label,
                                        pending_block_name.as_deref(),
                                    ));
                                }
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
                                r.set_prompt(make_repl_prompt(&session.prompt_state_label(), None));
                                continue;
                            }
                            continue;
                        }
                        ReadCommandOutput::CtrlD => {
                            if session.dismiss_pending_commands_block("Ctrl-D") {
                                r.set_prompt(make_repl_prompt(&session.prompt_state_label(), None));
                                continue;
                            }
                            break;
                        }
                    }
                }
            }
        }

        drop(session);

        if !self.no_save_state {
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
    let spec = get_stderr_log_filter();
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

    use gammalooprs::settings::{GlobalSettings, RuntimeSettings};

    use crate::commands::Reset;
    use crate::settings_tree::serialize_settings_with_defaults;
    use crate::state::ProcessRef;

    use super::{
        CLISettings, Commands, LogFormat, OneShot, Repl, RunHistory, StateSettings,
        GLOBAL_SETTINGS_FILENAME,
    };

    fn parse_with_shortcuts(args: &[&str]) -> OneShot {
        let argv: Vec<OsString> = args.iter().map(OsString::from).collect();
        OneShot::try_parse_from(OneShot::normalize_cli_shortcuts(&argv)).unwrap()
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
    fn oneshot_parses_subcommand_without_run_card() {
        let parsed = OneShot::try_parse_from(["gammaloop", "display", "integrands"]).unwrap();
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
    fn oneshot_accepts_logging_prefix_override() {
        let parsed = OneShot::try_parse_from(["gammaloop", "-p", "long"]).unwrap();
        assert_eq!(parsed.logging_prefix, Some(LogFormat::Long));
    }

    #[test]
    fn oneshot_accepts_settings_global_shortcut() {
        let parsed = parse_with_shortcuts(&["gammaloop", "-sg", "global.toml"]);
        assert_eq!(
            parsed.settings_global_path,
            Some(PathBuf::from("global.toml"))
        );
    }

    #[test]
    fn oneshot_accepts_settings_runtime_defaults_shortcut() {
        let parsed = parse_with_shortcuts(&["gammaloop", "-sr", "runtime.toml"]);
        assert_eq!(
            parsed.settings_runtime_defaults_path,
            Some(PathBuf::from("runtime.toml"))
        );
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
            Some(ProcessRef::Unqualified("triangle".to_string()))
        );
        assert_eq!(integrate.integrand_name, Some("LO".to_string()));
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
    fn repl_parses_reset_processes_with_normalized_selectors() {
        let parsed = Repl::try_parse_from([
            "gammaloop",
            "reset",
            "processes",
            "-p",
            "triangle",
            "-i",
            "LO",
        ])
        .unwrap();

        let Commands::Reset(Reset::Processes {
            process,
            integrand_name,
        }) = parsed.command
        else {
            panic!("expected reset processes command");
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
        runtime_file_settings.general.mu_r_sq = 37.0;
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
        assert_eq!(runtime_settings.general.mu_r_sq, 37.0);
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
    fn settings_file_overrides_replace_boot_card_settings() {
        let mut run_history = RunHistory::default();
        run_history.cli_settings.global.display_directive = "error".into();
        run_history.default_runtime_settings.general.mu_r_sq = 11.0;

        let mut overridden_global = GlobalSettings::default();
        overridden_global.display_directive = "warn".into();
        let mut overridden_runtime = RuntimeSettings::default();
        overridden_runtime.general.mu_r_sq = 29.0;

        let overrides = super::SettingsFileOverrides {
            global: Some(overridden_global),
            default_runtime_settings: Some(overridden_runtime),
        };

        let overridden = overrides.apply_to_run_history(&run_history);

        assert_eq!(overridden.cli_settings.global.display_directive, "warn");
        assert_eq!(overridden.default_runtime_settings.general.mu_r_sq, 29.0);
        assert_eq!(run_history.cli_settings.global.display_directive, "error");
        assert_eq!(run_history.default_runtime_settings.general.mu_r_sq, 11.0);
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
    }
}
