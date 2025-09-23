use ::tracing::debug;
use ::tracing::info;
use ::tracing::level_filters::LevelFilter;
use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand};
use gammalooprs::processes::ProcessList;
use reedline::DefaultPrompt;
use reedline::DefaultPromptSegment;
use reedline::FileBackedHistory;
use repl::ClapEditor;
use repl::ReadCommandOutput;
use state::set_serialize_commands_as_strings;
use tracing::get_stderr_log_filter;
use tracing::set_stderr_log_filter;

// use clap_repl::{
//     reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
//     ClapEditor, ReadCommandOutput,
// };
use crate::state::CommandHistory;
use color_eyre::Report;
use color_eyre::Result;
use colored::Colorize;
use console::style;
use dirs::home_dir;
use eyre::{eyre, Context};
use gammalooprs::{
    initialisation::initialise,
    settings::{GlobalSettings, RuntimeSettings},
    status_info,
    utils::{
        serde_utils::{get_schema_folder, SmartSerde, SHOWDEFAULTS},
        tracing::{LogFormat, LogLevel},
        GIT_VERSION,
    },
};
use integrate::Integrate;
use schemars::{schema_for, JsonSchema};
use serde::{Deserialize, Serialize};
use set::Set;
use state::{RunHistory, State};
use std::str::FromStr;
use std::{ffi::OsString, fs, path::PathBuf};
use std::{fs::File, ops::ControlFlow};
use symbolica::activate_oem_license;
// use tracing::LogLevel;

pub mod generate;
pub mod import_model;
pub mod inspect;
pub mod integrate;
#[cfg(feature = "python_api")]
pub mod python;
pub mod repl;
pub mod set;
pub mod state;
pub mod tracing;

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
pub struct Cli {
    /// Path to the a run file as history, and as settings, by default is: `./gammaloop_state/run.yaml`
    #[arg(short = 'r', long, value_hint = clap::ValueHint::FilePath)]
    pub run_history: Option<PathBuf>,

    /// Path to the state folder
    #[arg(short = 's', long, default_value = "./gammaloop_state", value_hint = clap::ValueHint::DirPath)]
    pub state_folder: PathBuf,

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

    #[arg(short = 'd', default_value_t = false)]
    debug: bool,

    /// Do not skip fields that are default when saving state
    #[arg(long)]
    no_skip_default: bool,

    /// Try to serialize using strings when saving run history
    #[arg(long)]
    try_strings: bool,

    // /// Debug level
    // #[arg(short = 'd', long, value_enum, default_value_t = LogLevel::Info)]
    // debug_level: LogLevel,
    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

impl FromStr for Commands {
    type Err = Report;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        Ok(CommandHistory::from_raw_string(s)?.command)
    }
}

pub struct Parsed {
    pub cli: Cli,
    pub input_string: String,
    pub matches: clap::ArgMatches,
}
impl Cli {
    pub fn new_test(state_folder: PathBuf) -> Self {
        Cli {
            run_history: None,
            state_folder,
            model_file: None,
            no_save_state: true,
            override_state: false,
            command: None,
            level: Some(LogLevel::Info),
            debug: false,
            no_skip_default: false,
            try_strings: false,
            trace_logs_filename: None,
        }
    }

    /// Parse from env args *and* capture ArgMatches (explicit vs defaults).
    pub fn parse_env_with_capture() -> Result<Parsed, clap::Error> {
        let argv: Vec<OsString> = std::env::args_os().collect();

        // Build a Command (same as derive(Parser)) and get matches
        let mut cmd = <Cli as CommandFactory>::command();
        let matches = cmd.clone().try_get_matches_from(&argv)?;

        // Recreate the typed struct from matches (don’t lose info)
        // Recreate the typed struct from matches (don't lose info)
        let cli =
            <Cli as FromArgMatches>::from_arg_matches(&matches).map_err(|e| e.format(&mut cmd))?;

        Ok(Parsed {
            input_string: argv
                .into_iter()
                .map(|s| s.to_string_lossy().into_owned()) // handles non-UTF8 gracefully
                .collect::<Vec<_>>()
                .join(" "),
            matches,
            cli,
        })
    }

    fn override_settings(&mut self, other: Cli) {
        self.state_folder = other.state_folder;
        self.model_file = other.model_file;
        self.override_state = other.override_state;
        self.no_save_state = other.no_save_state;
        self.try_strings = other.try_strings;
        self.no_skip_default = other.no_skip_default;
    }

    pub fn run_command(
        &mut self,
        command: Commands,
        state: &mut State,
        run_history: &mut RunHistory,
        global_settings: &mut GlobalSettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<()>, Report> {
        if let Some(level) = self.level {
            let spec = level.to_env_spec();
            set_stderr_log_filter(spec)?;
        }
        if self.debug {
            set_stderr_log_filter(LogLevel::Debug.to_env_spec())?;
        }
        match command {
            Commands::Quit {} => {
                return Ok(ControlFlow::Break(()));
            }
            Commands::Inspect(inspect) => {
                let _ = inspect.run(state)?;
            }
            Commands::Bench {
                samples,
                process_id,
                process_name,
                n_cores,
            } => {
                state.bench(samples, process_id, process_name, n_cores)?;
            }
            Commands::Import(s) => match s {
                Import::Amplitude {
                    path,
                    process_id,
                    integrand_name,
                } => {
                    state.import_amplitude(path, None, process_id, integrand_name)?;
                }
                Import::Model(im) => {
                    im.run(state)?;
                }
            },
            Commands::Save(s) => match s {
                Save::Dot { path } => {
                    state.export_dots(path.unwrap_or(self.state_folder.clone()))?;
                }
                Save::State { path, .. } => {
                    self.save(
                        state,
                        run_history,
                        default_runtime_settings,
                        global_settings,
                        path,
                        false,
                    )?;
                }
                Save::Schema {} => {
                    write_schemas()?;
                }
            },
            Commands::Set(s) => s.run(
                state,
                global_settings,
                default_runtime_settings,
                &mut self.state_folder,
            )?,

            Commands::Generate(g) => g.run(
                state,
                &self.state_folder,
                self.override_state,
                &global_settings,
                &default_runtime_settings,
            )?,

            Commands::Integrate(g) => {
                g.run(state)?;
            }

            Commands::Display(l) => match l {
                Display::Integrands { process_id } => {
                    let process = if let Some(proc_id) = process_id {
                        if proc_id >= state.process_list.processes.len() {
                            return Err(eyre!(
                                "Process ID {} invalid, only {} processes available",
                                proc_id,
                                state.process_list.processes.len()
                            ));
                        }
                        &state.process_list.processes[proc_id]
                    } else {
                        if state.process_list.processes.len() == 1 {
                            &state.process_list.processes[0]
                        } else {
                            return Err(eyre!(
                                "Multiple processes available, please specify the process."
                            ));
                        }
                    };

                    println!("Integrands for process {}:", process.definition.process_id);
                    for integrand in process.get_integrand_names() {
                        println!("  {}", integrand);
                    }
                }
                Display::Processes => {
                    println!("Processes:");
                    for process in state.process_list.processes.iter() {
                        println!(
                            "#{:-10}  {}",
                            process.definition.process_id, process.definition.folder_name
                        );
                    }
                }
                Display::Model => {
                    println!("{}", state.model.name)
                }
            },
            Commands::Run(Run { path, commands }) => {
                let mut a = ControlFlow::Continue(());
                if let Some(mut new_run_history) = path.map(|p| RunHistory::load(p)).transpose()? {
                    *global_settings = new_run_history.global_settings.clone();
                    *default_runtime_settings = new_run_history.default_runtime_settings.clone();
                    let res =
                        new_run_history.run(self, state, global_settings, default_runtime_settings);

                    run_history.merge(new_run_history);
                    a = res?;
                };

                if let Some(c) = commands {
                    for c in c.split(';') {
                        if c.trim().is_empty() {
                            continue;
                        }
                        info!("Running command from --commands: {}", c);
                        let cmd = CommandHistory::from_raw_string(c)?;
                        let res = self.run_command(
                            cmd.command,
                            state,
                            run_history,
                            global_settings,
                            default_runtime_settings,
                        )?;
                        a = res;

                        if let ControlFlow::Break(()) = res {
                            return Ok(res);
                        }
                    }
                }
                return Ok(a);
            }
            Commands::Reset(Reset::Processes {}) => {
                let n_processes = state.process_list.processes.len();
                state.process_list = ProcessList::default();
                status_info!(
                    "All {} processes have been cleared from the current state.",
                    n_processes
                );
            }
            Commands::Batch {
                process_file: _process_file,
                batch_input_file: _batch_input_file,
                name: _name,
                output_name: _output_name,
            } => {
                todo!("Batch command not implemented yet");
            }
        }
        Ok(ControlFlow::Continue(()))
    }

    pub fn run(mut self, raw: String) -> Result<()> {
        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        };
        initialise()?;

        let mut state = match State::load(
            self.state_folder.clone(),
            self.model_file.clone(),
            self.trace_logs_filename.clone(),
        ) {
            Ok(state) => state,
            Err(e) => {
                info!(
                    "No valid state folder found ({e}) :{} model_file:{}. Creating new default state.",
                    self.state_folder.display(),
                    self.model_file
                        .as_ref()
                        .map(|p| p.display().to_string())
                        .unwrap_or("None".to_string())
                );
                State::new(self.state_folder.clone(), self.trace_logs_filename.clone())
            }
        };
        let mut run_history = self.get_run_history()?;
        let mut global_settings = run_history.global_settings.clone();
        let mut default_runtime_settings = run_history.default_runtime_settings.clone();

        // let mut state = State::load(&cli.state_file,);

        if let Some(a) = self.command.take() {
            let _ = self.run_command(
                a.clone(),
                &mut state,
                &mut run_history,
                &mut global_settings,
                &mut default_runtime_settings,
            )?;

            run_history.push_with_raw(a, Some(raw));
        } else {
            print_banner();
            let prompt = DefaultPrompt {
                left_prompt: DefaultPromptSegment::Basic(format!(
                    "{} | γloop ",
                    self.state_folder.display()
                )),
                ..DefaultPrompt::default()
            };

            // 2. Build the REPL – clap‑repl takes ownership and configures rustyline.
            let mut repl = ClapEditor::<Cli>::builder().with_prompt(Box::new(prompt));

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
                    ReadCommandOutput::Command(mut cli, raw_input) => {
                        if let Some(c) = cli.command.take() {
                            match cli.run_command(
                                c.clone(),
                                &mut state,
                                &mut run_history,
                                &mut global_settings,
                                &mut default_runtime_settings,
                            ) {
                                Err(e) => {
                                    eprintln!("{e:?}");
                                    self.override_settings(cli);
                                }
                                Ok(ControlFlow::Break(())) => {
                                    run_history.push_with_raw(c, Some(raw_input.clone()));
                                    self.override_settings(cli);
                                    break;
                                }
                                _ => {
                                    run_history.push_with_raw(c, Some(raw_input));
                                    self.override_settings(cli);
                                }
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
            self.save(
                &mut state,
                &run_history,
                &default_runtime_settings,
                &global_settings,
                None,
                false,
            )?;
        }
        Ok(())
    }

    // pub fn initialize(&self) {}

    fn get_run_history(&self) -> Result<RunHistory> {
        let default_path = self.state_folder.join("run.toml");
        let path = self.run_history.as_ref().unwrap_or(&default_path);
        if !path.exists() {
            status_info!("Loading default state at {}", path.display());
            Ok(RunHistory::default())
        } else {
            match RunHistory::load(path) {
                Ok(r) => Ok(r),
                Err(e) => {
                    if path != &default_path {
                        Err(eyre!(
                            "Could not load run history from {}: {}",
                            path.display(),
                            e
                        ))
                    } else {
                        info!(
                            "No valid run history found at {}, creating new default run history",
                            path.display()
                        );
                        Ok(RunHistory::default())
                    }
                }
            }
        }
    }

    pub fn save(
        &self,
        state: &mut State,
        run_history: &RunHistory,
        default_runtime_settings: &RuntimeSettings,
        global_settings: &GlobalSettings,
        root_folder: Option<PathBuf>,
        strict: bool,
    ) -> Result<()> {
        // let root_folder = root_folder.join("gammaloop_state");

        // check if the export root exists, if not create it, if it does return error
        let mut selected_root_folder = root_folder.unwrap_or(self.state_folder.clone());
        if !selected_root_folder.exists() {
            fs::create_dir_all(&selected_root_folder)?;
        } else {
            if strict {
                return Err(eyre!(
                    "Export root already exists, please choose a different path or remove the existing directory",
                ));
            }

            if !self.override_state {
                while selected_root_folder.clone().exists() {
                    eprint!(
                        "Gammaloop export root {} already exists. Specify '{}' for overwriting, '{}' for not saving, or '{}' to specify where to save current state to:\n > ",
                        selected_root_folder.display().to_string().green(),
                        "o".red().bold(),
                        "n".blue().bold(),
                        "<NEW_PATH>".green().bold()
                    );
                    let mut user_input = String::new();
                    std::io::stdin()
                        .read_line(&mut user_input)
                        .expect("Could not read user-specified gammaloop state export destination");
                    //user_input = user_input.trim().into();
                    match user_input.trim() {
                        "o" => {
                            status_info!(
                                "Overwriting existing gammaloop state at {}",
                                selected_root_folder.display().to_string().green()
                            );
                            break;
                        }
                        "n" => {
                            return Ok(());
                        }
                        new_path => {
                            selected_root_folder = new_path.into();
                            continue;
                        }
                    }
                }
            }
        }

        state.save(&selected_root_folder, true, false)?;
        if self.no_skip_default {
            SHOWDEFAULTS.store(true, std::sync::atomic::Ordering::Relaxed);
        }

        set_serialize_commands_as_strings(self.try_strings);
        run_history.save_toml(&selected_root_folder, true, false)?;
        set_serialize_commands_as_strings(false);

        default_runtime_settings.to_file(
            &selected_root_folder.join("default_runtime_settings.toml"),
            true,
        )?;
        global_settings.to_file(&selected_root_folder.join("global_settings.toml"), true)?;

        if self.no_skip_default {
            SHOWDEFAULTS.store(false, std::sync::atomic::Ordering::Relaxed);
        }

        Ok(())
    }
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Import {
    Model(import_model::ImportModel),
    Amplitude {
        // #[arg(short = 'p')]
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: PathBuf,

        #[arg(short = 'i')]
        process_id: Option<usize>,

        #[arg(short = 'n')]
        integrand_name: Option<String>,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Save {
    Dot {
        #[arg(value_hint = clap::ValueHint::FilePath)]
        path: Option<PathBuf>,
    },
    State {
        /// Path to save the state to, by default is the current state folder
        #[arg(short = 'p', long, value_hint = clap::ValueHint::FilePath)]
        path: Option<PathBuf>,

        /// Do not skip fields that are default
        #[arg(short = 'd', long)]
        no_skip_default: bool,
    },
    /// regenerate the schema files
    Schema {},
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Reset {
    Processes {},
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Display {
    Model,
    Processes,
    Integrands { process_id: Option<usize> },
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Run {
    /// Path to a run file to execute
    #[arg(value_hint = clap::ValueHint::FilePath)]
    path: Option<PathBuf>,

    #[arg(short = 'c', long)]
    commands: Option<String>,
}

pub enum Log {
    Level(LevelFilter),
    Format(LogFormat),
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Commands {
    #[clap(subcommand)]
    Display(Display),
    #[clap(subcommand)]
    Set(Set),
    #[clap(subcommand)]
    Import(Import),
    #[clap(subcommand)]
    Save(Save),

    Run(Run),

    #[clap(subcommand)]
    Reset(Reset),

    Integrate(Integrate),

    Generate(generate::Generate),

    /// Quit gammaloop
    Quit {
        // #[arg(short = 'o', long, default_value_t = false)]
        // override_state: bool,
    },
    /// Inspect a single phase‑space point / momentum configuration
    // #[clap(subcommand)]
    Inspect(inspect::Inspect),

    /// Benchmark raw integrand evaluation speed
    Bench {
        /// Number of random samples to evaluate
        #[arg(short = 's', long, value_name = "SAMPLES")]
        samples: usize,
        /// The process id to inspect
        #[arg(short = 'i', long = "process-id", value_name = "ID")]
        process_id: usize,

        /// The name of the process to inspect
        #[arg(short = 'n', long = "name", value_name = "NAME")]
        process_name: String,
        /// Number of cores to parallelize over
        #[arg(short = 'c', long)]
        n_cores: usize,
    },

    /// HPC batch evaluation branch
    Batch {
        #[arg(value_name = "PROCESS_FILE", value_hint = clap::ValueHint::FilePath)]
        process_file: PathBuf,
        #[arg(value_name = "BATCH_INPUT_FILE", value_hint = clap::ValueHint::FilePath)]
        batch_input_file: PathBuf,
        #[arg(short = 'n', long, value_name = "NAME")]
        name: String,
        #[arg(value_name = "NAME")]
        output_name: String,
    },
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

impl Commands {}

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
