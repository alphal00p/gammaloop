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

use eyre::eyre;
use gammalooprs::utils::serde_utils::SerdeFileError;
use reedline::DefaultPrompt;
use reedline::DefaultPromptSegment;
use reedline::FileBackedHistory;
use repl::ClapEditor;
use repl::ReadCommandOutput;
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

use state::{RunHistory, State, SyncSettings};
use std::{ffi::OsString, path::PathBuf};
use std::{fs::File, ops::ControlFlow};
use symbolica::activate_oem_license;

// use tracing::LogLevel;
mod command_parser;
#[cfg(feature = "python_api")]
pub mod python;
pub mod repl;

pub mod state;
pub mod templates;
pub mod tracing;

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
pub struct Repl {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
#[command(args_conflicts_with_subcommands = true)]
pub struct OneShot {
    /// Don't actually run anything, just build up run card
    #[arg(short = 'd', long, default_value_t = false)]
    pub dry_run: bool,

    /// Don't try to load state, just start with a new one
    #[arg(short = 'f', long, default_value_t = false)]
    pub fresh_state: bool,

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
    #[arg(short = 't')]
    trace_logs_filename: Option<String>,

    /// Set log level for current session
    #[arg(short = 'l')]
    level: Option<LogLevel>,

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
    #[serde(skip_serializing_if = "is_false")]
    pub dummy: bool,
    pub state_folder: PathBuf,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub global: GlobalSettings,
}

impl Default for CLISettings {
    fn default() -> Self {
        CLISettings {
            try_strings: true,
            override_state: false,
            dummy: false,
            state_folder: "./gammaloop_state".into(),
            global: GlobalSettings::default(),
        }
    }
}

impl CLISettings {
    pub fn override_with(&mut self, cli: &OneShot) {
        self.try_strings = !cli.no_try_strings;
        self.dummy = cli.dry_run;

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
impl OneShot {
    pub fn new_cli_settings(&self, global: GlobalSettings) -> CLISettings {
        CLISettings {
            dummy: self.dry_run,
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
            dry_run: false,
            model_file: None,
            no_save_state: true,
            override_state: false,
            command: None,
            level: Some(LogLevel::Info),
            no_try_strings: false,
            fresh_state: false,
            trace_logs_filename: None,
        }
    }

    fn run_command(&self) -> Option<&commands::run::Run> {
        match self.command.as_ref() {
            Some(Commands::Run(run)) => Some(run),
            _ => None,
        }
    }

    fn load_initial_run_history(&self) -> Result<Option<RunHistory>> {
        self.run_command()
            .map(commands::run::Run::load_run_history)
            .transpose()
            .map(Option::flatten)
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

    fn resolve_initial_state_folder(&self) -> Result<PathBuf> {
        if self.state_folder_explicitly_set {
            return Ok(self.state_folder.clone());
        }

        if let Some(run_history) = self.load_initial_run_history()? {
            return Ok(run_history.cli_settings.state_folder.clone());
        }

        Ok(self.state_folder.clone())
    }

    pub fn load(&mut self) -> Result<(State, RunHistory, CLISettings, RuntimeSettings)> {
        let (state, run_history, cli_settings, default_runtime_settings) = match if self.fresh_state
        {
            Err(eyre!("Fresh state requested, skipping load"))
        } else {
            State::load(
                self.state_folder.clone(),
                self.model_file.clone(),
                self.trace_logs_filename.clone(),
            )
        } {
            Ok(state) => {
                let mut global =
                    match CLISettings::from_file_typed(self.state_folder.join("cli_settings.toml"))
                    {
                        Ok(a) => a,
                        Err(SerdeFileError::FileError(_)) => {
                            println!(
                                "File error for {} loading default ",
                                &self.state_folder.join("cli_settings.toml").display()
                            );
                            CLISettings::default()
                        }

                        Err(e) => return Err(e.into()),
                    };

                global.override_with(self);

                let default_runtime = match RuntimeSettings::from_file_typed(
                    self.state_folder.join("default_runtime_settings.toml"),
                ) {
                    Ok(a) => a,
                    Err(SerdeFileError::FileError(_)) => RuntimeSettings::default(),
                    Err(e) => return Err(e.into()),
                };
                let run_path = self.state_folder.join("run.toml");
                let run_history = if run_path.exists() {
                    RunHistory::load(self.state_folder.join("run.toml"))?
                } else {
                    RunHistory::default()
                };

                (state, run_history, global, default_runtime)
            }
            Err(e) => {
                info!(
                    "No valid state folder found ({e}) :{} model_file:{}. Creating new default state.",
                    self.state_folder.display(),
                    self.model_file
                        .as_ref()
                        .map(|p| p.display().to_string())
                        .unwrap_or("None".to_string())
                );
                let state = State::new(self.state_folder.clone(), self.trace_logs_filename.clone());

                if let Some(mut run_history) = self.load_initial_run_history()? {
                    // Freeze state folder for the whole session.
                    run_history.cli_settings.state_folder = self.state_folder.clone();
                    let mut global = run_history.cli_settings.clone();
                    global.override_with(self);
                    let default_runtime = run_history.default_runtime_settings.clone();
                    (state, run_history, global, default_runtime)
                } else {
                    (
                        state,
                        RunHistory::default(),
                        self.new_cli_settings(GlobalSettings::default()),
                        RuntimeSettings::default(),
                    )
                }
            }
        };

        cli_settings.sync_settings()?;

        Ok((state, run_history, cli_settings, default_runtime_settings))
    }

    pub fn run(mut self, raw: String) -> Result<()> {
        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        };
        initialise()?;

        self.state_folder = self.resolve_initial_state_folder()?;

        let (mut state, mut run_history, mut cli_settings, mut default_runtime_settings) =
            self.load()?;

        println!(":{}", cli_settings.global.display_directive);

        let mut save_state = SaveState::default();

        if let Some(a) = self.command.take() {
            let flow = a.clone().run(
                &mut state,
                &mut run_history,
                &mut cli_settings,
                &mut default_runtime_settings,
            )?;

            if let ControlFlow::Break(a) = flow {
                save_state = a;
            }
            run_history.push_with_raw(a, Some(raw));
        } else {
            print_banner();
            let prompt = if cli_settings.dummy {
                DefaultPrompt {
                    left_prompt: DefaultPromptSegment::Basic(format!(
                        "{} | γloop DUMMY ",
                        cli_settings.state_folder.display()
                    )),
                    ..DefaultPrompt::default()
                }
            } else {
                DefaultPrompt {
                    left_prompt: DefaultPromptSegment::Basic(format!(
                        "{} | γloop ",
                        cli_settings.state_folder.display()
                    )),
                    ..DefaultPrompt::default()
                }
            };

            // 2. Build the REPL – clap‑repl takes ownership and configures rustyline.
            let mut repl = ClapEditor::<Repl>::builder().with_prompt(Box::new(prompt));

            if let Some(home) = home_dir() {
                repl = repl.with_editor_hook(move |reed| {
                    // Do custom things with `Reedline` instance here
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
                        match command.command.clone().run(
                            &mut state,
                            &mut run_history,
                            &mut cli_settings,
                            &mut default_runtime_settings,
                        ) {
                            Err(e) => {
                                eprintln!("{e:?}");
                            }
                            Ok(ControlFlow::Break(a)) => {
                                save_state = a;
                                run_history.push_with_raw(command.command, Some(raw_input.clone()));
                                break;
                            }
                            _ => {
                                run_history.push_with_raw(command.command, Some(raw_input));
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
                    ReadCommandOutput::CtrlC => continue,
                    ReadCommandOutput::CtrlD => break,
                }
            }
        }

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
    use std::{fs, path::PathBuf};

    use clap::Parser;
    use tempfile::tempdir;

    use super::{CLISettings, Commands, OneShot, RunHistory};

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
            OneShot::try_parse_from(["gammaloop", "run", run_path.to_string_lossy().as_ref()])
                .unwrap();
        one_shot.state_folder = "./from_cli".into();
        one_shot.state_folder_explicitly_set = true;

        assert_eq!(
            one_shot.resolve_initial_state_folder().unwrap(),
            PathBuf::from("./from_cli")
        );
    }

    #[test]
    fn resolve_initial_state_folder_uses_run_card_when_cli_not_explicit() {
        let temp = tempdir().unwrap();
        let run_path = temp.path().join("run.toml");
        write_run_card(&run_path, "./from_run_card");

        let one_shot =
            OneShot::try_parse_from(["gammaloop", "run", run_path.to_string_lossy().as_ref()])
                .unwrap();

        assert_eq!(
            one_shot.resolve_initial_state_folder().unwrap(),
            PathBuf::from("./from_run_card")
        );
    }

    #[test]
    fn oneshot_parses_run_subcommand_block_names() {
        let parsed =
            OneShot::try_parse_from(["gammaloop", "run", "card.toml", "generation", "integration"])
                .unwrap();
        let run = parsed.run_command().unwrap();
        assert_eq!(
            run.selected_block_names(),
            Some(&["generation".to_string(), "integration".to_string()][..])
        );
    }

    #[test]
    fn oneshot_treats_subcommand_names_as_run_block_names() {
        let parsed = OneShot::try_parse_from([
            "gammaloop",
            "run",
            "card.toml",
            "generate",
            "integrate_euclidean",
            "quit",
        ])
        .unwrap();
        let run = parsed.run_command().unwrap();
        assert_eq!(
            run.selected_block_names(),
            Some(
                &[
                    "generate".to_string(),
                    "integrate_euclidean".to_string(),
                    "quit".to_string()
                ][..]
            )
        );
    }

    #[test]
    fn oneshot_parses_subcommand_without_run_card() {
        let parsed = OneShot::try_parse_from(["gammaloop", "display", "integrands"]).unwrap();
        assert!(matches!(parsed.command, Some(Commands::Display(_))));
    }

    #[test]
    fn oneshot_rejects_run_card_without_run_subcommand() {
        assert!(OneShot::try_parse_from(["gammaloop", "card.toml"]).is_err());
    }
}
