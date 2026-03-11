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

use gammalooprs::utils::serde_utils::SerdeFileError;
use reedline::DefaultPrompt;
use reedline::DefaultPromptSegment;
use reedline::FileBackedHistory;
use repl::ClapEditor;
use repl::ReadCommandOutput;
use session::{CliSession, CliSessionState};
use tracing::get_stderr_log_filter;

// use clap_repl::{
//     reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
//     ClapEditor, ReadCommandOutput,
// };

use color_eyre::Result;
use colored::Colorize;
use console::style;
use dirs::home_dir;
use eyre::Context;
use gammalooprs::{
    initialisation::initialise,
    settings::{GlobalSettings, RuntimeSettings},
    utils::serde_utils::IsDefault,
    utils::{
        serde_utils::{get_schema_folder, is_false, is_true, SmartSerde},
        tracing::{LogFormat, LogLevel},
        GIT_VERSION,
    },
};
use schemars::{schema_for, JsonSchema};
use serde::{Deserialize, Serialize};

use state::{CommandHistory, RunHistory, State, SyncSettings};
use std::{ffi::OsString, path::PathBuf};
use std::{fs::File, ops::ControlFlow};
use symbolica::activate_oem_license;

// use tracing::LogLevel;
mod command_parser;
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
    pub state_folder: PathBuf,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub global: GlobalSettings,
}

impl Default for CLISettings {
    fn default() -> Self {
        CLISettings {
            try_strings: true,
            override_state: false,
            state_folder: "./gammaloop_state".into(),
            global: GlobalSettings::default(),
        }
    }
}

impl CLISettings {
    pub fn override_with(&mut self, cli: &OneShot) {
        self.try_strings = !cli.no_try_strings;

        self.override_state = cli.override_state;

        self.state_folder = cli.state_folder.clone();
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

impl OneShot {
    pub fn new_cli_settings(&self, global: GlobalSettings) -> CLISettings {
        CLISettings {
            try_strings: !self.no_try_strings,
            override_state: self.override_state,
            state_folder: self.state_folder.clone(),
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
            return Ok(run_history.cli_settings.state_folder.clone());
        }

        Ok(self.state_folder.clone())
    }

    pub fn load(&mut self) -> Result<(State, RunHistory, CLISettings, RuntimeSettings)> {
        let state_exists = self.state_folder.exists();

        let (state, run_history, cli_settings, default_runtime_settings) =
            if !self.fresh_state && state_exists {
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

                (state, run_history, global, default_runtime)
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
                )
            };

        cli_settings.sync_settings()?;

        Ok((state, run_history, cli_settings, default_runtime_settings))
    }

    pub fn run(mut self, raw: String) -> Result<()> {
        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        };
        initialise()?;

        let boot_run_history = self.load_boot_run_history()?;
        let settings_file_overrides = self.load_settings_file_overrides()?;
        self.state_folder = self.resolve_initial_state_folder(boot_run_history.as_ref())?;

        let (mut state, mut run_history, mut cli_settings, mut default_runtime_settings) =
            self.load()?;
        settings_file_overrides.apply(&mut cli_settings, &mut default_runtime_settings)?;

        println!(":{}", cli_settings.global.display_directive);

        let mut save_state = SaveState::default();
        let prompt_state_folder = cli_settings.state_folder.clone();
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
                let prompt = DefaultPrompt {
                    left_prompt: DefaultPromptSegment::Basic(format!(
                        "{} | γloop ",
                        prompt_state_folder.display()
                    )),
                    ..DefaultPrompt::default()
                };

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
                    .with_prompt(Box::new(prompt))
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
                                continue;
                            }
                            continue;
                        }
                        ReadCommandOutput::CtrlD => {
                            if session.dismiss_pending_commands_block("Ctrl-D") {
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
        r"              ██         ▄████████▄  ▄████████▄  ██████████▄
  ██          ▀▀         ▀▀      ▀▀  ▀▀      ▀▀  ▀▀       ██
    ██  ▄██████████████████████████████████████████████████▀
     ▀██▀     ▄▄         ▄▄      ▄▄  ▄▄      ▄▄  ▄▄
    ██  ██    █████████  ▀████████▀  ▀████████▀  ██
   ██    ██
   ██    ██
"
        .to_string()
        .bold()
        .blue(),
        format!(
            r#"   ▀██████▀   version:{:<15} log level:{}         "#,
            GIT_VERSION.green(),
            spec.green(),
        )
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
    use tempfile::tempdir;

    use gammalooprs::settings::{GlobalSettings, RuntimeSettings};

    use crate::state::ProcessRef;

    use super::{CLISettings, Commands, OneShot, Repl, RunHistory, GLOBAL_SETTINGS_FILENAME};

    fn parse_with_shortcuts(args: &[&str]) -> OneShot {
        let argv: Vec<OsString> = args.iter().map(OsString::from).collect();
        OneShot::try_parse_from(OneShot::normalize_cli_shortcuts(&argv)).unwrap()
    }

    fn write_run_card(path: &PathBuf, state_folder: &str) {
        let run_history = RunHistory {
            cli_settings: CLISettings {
                state_folder: state_folder.into(),
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
    fn settings_file_overrides_apply_without_touching_state_folder() {
        let temp = tempdir().unwrap();
        let global_path = temp.path().join("global.toml");
        let runtime_path = temp.path().join("runtime.toml");

        let mut global_file_settings = CLISettings::default();
        global_file_settings.state_folder = "./ignored".into();
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
        cli_settings.state_folder = "./keep".into();
        cli_settings.global = GlobalSettings::default();
        let mut runtime_settings = RuntimeSettings::default();

        overrides
            .apply(&mut cli_settings, &mut runtime_settings)
            .unwrap();

        assert_eq!(cli_settings.state_folder, PathBuf::from("./keep"));
        assert_eq!(cli_settings.global.display_directive, "warn");
        assert_eq!(runtime_settings.general.mu_r_sq, 37.0);
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
