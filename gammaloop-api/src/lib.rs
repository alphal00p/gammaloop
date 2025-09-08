use ::tracing::{debug, info, level_filters::LevelFilter};
use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand};

use clap_repl::{
    reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
    ClapEditor, ReadCommandOutput,
};
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
use state::{current_log_spec, RunHistory, State};
use std::{ffi::OsString, fs, path::PathBuf};
use std::{fs::File, ops::ControlFlow};
use symbolica::activate_oem_license;
// use tracing::LogLevel;

pub mod generate;
pub mod import_model;
pub mod inspect;
pub mod integrate;
pub mod state;
pub mod tracing;

#[derive(Parser, Debug)]
#[command(name = "gammaLoop", version, about)]
#[command(next_line_help = true)]
pub struct Cli {
    /// Path to the a run file as history, and as settings, by default is: `./gammaloop_state/run.yaml`
    #[arg(short = 'r', long)]
    run_history: Option<PathBuf>,

    /// Path to the state folder
    #[arg(short = 's', long, default_value = "./gammaloop_state")]
    pub state_folder: PathBuf,

    /// Path to the model file
    #[arg(short = 'm', long)]
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

    // /// Debug level
    // #[arg(short = 'd', long, value_enum, default_value_t = LogLevel::Info)]
    // debug_level: LogLevel,
    /// Optional sub‑command
    #[command(subcommand)]
    pub command: Option<Commands>,
}

pub struct Parsed {
    pub argv: Vec<OsString>, // raw argv from the OS (including program name)
    pub matches: clap::ArgMatches, // only what the user actually supplied
    pub cli: Cli,            // your typed struct, built from matches
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
        let cli =
            <Cli as FromArgMatches>::from_arg_matches(&matches).map_err(|e| e.format(&mut cmd))?;

        Ok(Parsed { argv, matches, cli })
    }

    fn override_settings(&mut self, other: Cli) {
        self.state_folder = other.state_folder;
        self.model_file = other.model_file;
        self.override_state = other.override_state;
        self.no_save_state = other.no_save_state;
    }

    fn run_command(
        &mut self,
        command: Commands,
        state: &mut State,
        run_history: &mut RunHistory,
        global_settings: &mut GlobalSettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<ControlFlow<()>, Report> {
        if let Some(level) = self.level {
            let spec = level.to_env_spec();
            state.set_log_spec(spec)?;
        }
        if self.debug {
            state.set_log_spec(LogLevel::Debug.to_env_spec())?;
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
                Import::Amplitude { path } => {
                    state.import_amplitude(path, None)?;
                }
                Import::Model(im) => {
                    im.run(state)?;
                }
            },
            Commands::Save(s) => match s {
                Save::Dot { path } => {
                    state.export_dots(path.unwrap_or(self.state_folder.clone()))?;
                }
                Save::State {
                    path,
                    no_skip_default,
                } => {
                    self.save(
                        state,
                        run_history,
                        &default_runtime_settings,
                        &global_settings,
                        path,
                        false,
                        no_skip_default,
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
                    let process = &state.process_list.processes[process_id];
                    println!("Integrands for process {}:", process_id);
                    for integrand in process.get_integrand_names() {
                        println!("  {}", integrand);
                    }
                }
                Display::Processes => {
                    println!("Processes:");
                    for (index, process) in state.process_list.processes.iter().enumerate() {
                        println!("  {}", process.definition.folder_name(&state.model, index));
                    }
                }
                Display::Model => {
                    println!("{}", state.model.name)
                }
            },
            Commands::Run(Run { path }) => {
                let mut new_run_history = RunHistory::new(path)?;
                let res =
                    new_run_history.run(self, state, global_settings, default_runtime_settings);

                run_history.merge(new_run_history);
                return res;
            }
            _ => {}
        }
        Ok(ControlFlow::Continue(()))
    }

    pub fn run(mut self) -> Result<()> {
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

            run_history.push(a);
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
                    ReadCommandOutput::Command(mut cli) => {
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
                                    run_history.push(c);
                                    self.override_settings(cli);
                                    break;
                                }
                                _ => {
                                    run_history.push(c);
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
                false,
            )?;
        }
        Ok(())
    }

    // pub fn initialize(&self) {}

    fn get_run_history(&self) -> Result<RunHistory> {
        let default_path = self.state_folder.join("run.toml");
        let path = self.run_history.as_ref().unwrap_or(&default_path);
        RunHistory::new(path)
    }

    pub fn save(
        &self,
        state: &mut State,
        run_history: &RunHistory,
        default_runtime_settings: &RuntimeSettings,
        global_settings: &GlobalSettings,
        root_folder: Option<PathBuf>,
        strict: bool,
        no_skip_default: bool,
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
        run_history.save_toml(&selected_root_folder, true, false)?;
        if no_skip_default {
            SHOWDEFAULTS.store(true, std::sync::atomic::Ordering::Relaxed);
        }

        default_runtime_settings.to_file(
            &selected_root_folder.join("default_runtime_settings.toml"),
            true,
        )?;
        global_settings.to_file(&selected_root_folder.join("global_settings.toml"), true)?;

        if no_skip_default {
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
        path: PathBuf,
        // format: String,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Save {
    Dot {
        path: Option<PathBuf>,
    },
    State {
        /// Path to save the state to, by default is the current state folder
        #[arg(short = 'p', long)]
        path: Option<PathBuf>,

        /// Do not skip fields that are default
        #[arg(short = 'd', long)]
        no_skip_default: bool,
    },
    /// regenerate the schema files
    Schema {},
}

pub mod set;

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Display {
    Model,
    Processes,
    Integrands { process_id: usize },
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Run {
    /// Path to a run file to execute
    path: PathBuf,
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
        #[arg(value_name = "PROCESS_FILE")]
        process_file: PathBuf,
        #[arg(value_name = "BATCH_INPUT_FILE")]
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
    let spec = current_log_spec();
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
    let folder = get_schema_folder()?;

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
