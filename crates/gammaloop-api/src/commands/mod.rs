//! Typed command definitions and command-execution results.
//!
//! Commands share the CLI's stateful behavior even when executed through
//! [`crate::session::CliSession`]. Most communicate through state/settings mutations, logs, or
//! explicitly requested files. [`CommandOutput`] is reserved for results that the embedded Rust
//! caller should consume directly.

use std::{ops::ControlFlow, path::PathBuf, str::FromStr};

use clap::{builder::ArgExt, Arg, Subcommand};
use color_eyre::Report;
use gammalooprs::settings::RuntimeSettings;
use save::SaveState;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    state::{CommandHistory, RunHistory, State},
    CLISettings,
};
use symbolica::atom::Atom;
pub mod approach;
pub use approach::Approach;
pub mod bench;
pub use bench::Bench;
pub mod commands_block;
pub use commands_block::StartCommandsBlock;
pub mod display;
pub use display::Display;
pub mod duplicate;
pub use duplicate::Duplicate;
pub mod generate;
pub use generate::Generate;
pub mod import;
pub use import::Import;
pub mod inspect;
pub use inspect::Inspect;
pub mod integrate;
pub use integrate::{Integrate, IntegrationOutput};
pub mod remove;
pub use remove::Remove;
pub mod save;
pub use save::Save;
pub mod select;
pub use select::Select;
pub mod set;
pub mod shell;
pub use set::Set;
pub use shell::Shell;
pub mod threedreps;
pub use threedreps::ThreeDRep;
pub mod run;
pub use run::Run;
pub mod evaluate;
pub mod evaluate_samples;
pub use evaluate::Evaluate;
pub mod renormalize;
pub use renormalize::Renormalize;
pub mod profile;
pub use profile::Profile;
pub(crate) mod process_settings;

/// Export metadata for Clap settings without public introspection APIs.
#[doc(hidden)]
#[derive(Clone, Debug, Default)]
pub struct CliArgumentMetadata {
    pub requires: Vec<&'static str>,
    pub default_missing_values: Vec<&'static str>,
}

impl ArgExt for CliArgumentMetadata {}

pub(crate) trait CliArgumentMetadataExt {
    fn cli_requires(self, id: &'static str) -> Self;
    fn cli_default_missing_value(self, value: &'static str) -> Self;
}

impl CliArgumentMetadataExt for Arg {
    fn cli_requires(self, id: &'static str) -> Self {
        let mut metadata = self
            .get::<CliArgumentMetadata>()
            .cloned()
            .unwrap_or_default();
        metadata.requires.push(id);
        self.requires(id).add(metadata)
    }

    fn cli_default_missing_value(self, value: &'static str) -> Self {
        let mut metadata = self
            .get::<CliArgumentMetadata>()
            .cloned()
            .unwrap_or_default();
        metadata.default_missing_values.push(value);
        self.default_missing_value(value).add(metadata)
    }
}

/// Structured value returned by a command, when that command has an embedded result.
///
/// Most commands return [`None`](CommandOutput::None); this does not mean that they had no effect.
/// They may have changed the in-memory state or settings, emitted diagnostics, or written an
/// explicitly requested file.
#[derive(Debug, Clone, Default)]
pub enum CommandOutput {
    /// No direct Rust value; inspect the mutated state, settings, logs, or requested output file.
    #[default]
    None,
    /// Symbolic result produced by an [`Evaluate`](Commands::Evaluate) command.
    Evaluate(Atom),
    /// Estimate, uncertainty, and integration metadata produced by an
    /// [`Integrate`](Commands::Integrate) command.
    Integrate(IntegrationOutput),
}

/// Control-flow request and optional structured value produced by one command.
///
/// Callers must inspect `flow` even when they do not need `output`. A break carries the
/// [`SaveState`] policy requested by save/quit handling; it does not persist the state merely by
/// being dropped.
#[derive(Debug, Clone)]
pub struct CommandExecution {
    /// Continue the session or return the requested save/quit policy to the embedding caller.
    pub flow: ControlFlow<SaveState>,
    /// Embedded result for evaluation/integration, otherwise [`CommandOutput::None`].
    pub output: CommandOutput,
}

impl CommandExecution {
    pub fn continue_with(output: CommandOutput) -> Self {
        Self {
            flow: ControlFlow::Continue(()),
            output,
        }
    }

    pub fn continue_without_output() -> Self {
        Self::continue_with(CommandOutput::None)
    }

    pub fn break_with(save_state: SaveState) -> Self {
        Self {
            flow: ControlFlow::Break(save_state),
            output: CommandOutput::None,
        }
    }
}

#[allow(clippy::large_enum_variant)]
#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Commands {
    /// Inspect models, processes, integrands, settings, and named session data.
    #[clap(subcommand)]
    Display(Display),
    /// Copy generated data to a new process or integrand name.
    #[clap(subcommand)]
    Duplicate(Duplicate),
    /// Change global, runtime, model, or process settings.
    #[clap(subcommand)]
    Set(Set),
    /// Import models or previously generated graphs into the active state.
    #[clap(subcommand)]
    Import(Import),
    /// Export graph data, standalone evaluators, schemas, or the current state.
    #[clap(subcommand)]
    Save(Save),

    /// Execute a stored command block or an inline command list.
    Run(Run),
    /// Begin recording subsequent commands under a reusable block name.
    #[command(name = "start_commands_block")]
    StartCommandsBlock(StartCommandsBlock),
    /// Stop recording commands and store the active command block.
    #[command(name = "finish_commands_block")]
    FinishCommandsBlock,

    /// Remove generated processes from the active state.
    #[command(name = "remove")]
    #[clap(subcommand)]
    Remove(Remove),

    /// Integrate one or more generated process integrands.
    Integrate(Integrate),

    Generate(Generate),

    /// Select graph groups and write the selection to a new process or integrand.
    Select(Select),

    /// End the interactive GammaLoop session without executing further commands.
    Quit(SaveState),
    /// Inspect integrand contributions at one phase-space point or momentum configuration.
    // #[clap(subcommand)]
    Inspect(Inspect),

    /// Sample an integrand while approaching a phase-space point along selected scaling axes.
    Approach(Approach),

    /// Evaluate a vacuum amplitude analytically or numerically with Vakint.
    Evaluate(Evaluate),

    /// Compute and export ultraviolet renormalization contributions for an amplitude.
    Renormalize(Renormalize),

    /// Measure evaluation throughput at a fixed point over a target duration.
    Bench(Bench),
    /// Profile ultraviolet or infrared scaling limits of a generated integrand.
    #[clap(subcommand)]
    Profile(Profile),

    /// Validate a graph or build its diagnostic three-dimensional energy representation.
    #[command(name = "3Drep")]
    #[clap(subcommand)]
    ThreeDRep(ThreeDRep),

    /// Evaluate a process-definition file against an HPC batch input and write named results.
    Batch {
        /// Process-definition file used to construct the batch evaluator.
        #[arg(value_name = "PROCESS_FILE", value_hint = clap::ValueHint::FilePath)]
        process_file: PathBuf,
        /// Batch input file containing the evaluation points or jobs to execute.
        #[arg(value_name = "BATCH_INPUT_FILE", value_hint = clap::ValueHint::FilePath)]
        batch_input_file: PathBuf,
        /// Logical batch name used in generated output metadata.
        #[arg(short = 'n', long, value_name = "NAME")]
        name: String,
        /// Base name assigned to the batch result output.
        #[arg(value_name = "NAME")]
        output_name: String,
    },
    /// Run an operating-system shell command from the interactive GammaLoop session.
    #[command(name = "!")]
    Shell(Shell),

    #[doc(hidden)]
    #[command(skip)]
    CommandTemplate(String),
}

impl FromStr for Commands {
    type Err = Report;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        Ok(CommandHistory::from_raw_string(s)?.command)
    }
}

impl Commands {
    pub fn run(
        self,
        state: &mut State,
        run_history: &mut RunHistory,
        global_cli_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<CommandExecution, Report> {
        match self {
            Commands::Profile(mut p) => {
                if let profile::ProfileResult::UltraViolet(analysis) =
                    p.run(state, global_cli_settings)?
                {
                    let verdict = analysis.pass_fail(gammalooprs::uv::profile::UV_PROFILE_MAX_DOD);
                    if verdict.failed > 0 {
                        let Profile::UltraViolet(options) = &mut p else {
                            unreachable!()
                        };
                        let (process_id, integrand_name) = state.find_integrand_ref(
                            options.process.as_ref(),
                            options.integrand_name.as_ref(),
                        )?;
                        options.process = Some(crate::state::ProcessRef::Id(process_id));
                        options.integrand_name = Some(integrand_name.clone());
                        options.seed = Some(options.seed.unwrap_or(42));
                        if !options.uv_ray_directions.is_empty() && options.uv_ray_norms.is_empty()
                        {
                            options.uv_ray_norms.push(
                                state
                                    .process_list
                                    .get_integrand_mut(process_id, &integrand_name)?
                                    .get_settings()
                                    .kinematics
                                    .e_cm,
                            );
                        }
                        let scope = if analysis.stopped_early {
                            "Stopped after the first failing limit; only completed limits are reported."
                        } else {
                            "All selected limits were profiled."
                        };
                        // Reuse the command serialization accepted by run cards.
                        // Resolve defaults that otherwise depend on session selection.
                        let reproduction =
                            toml::to_string_pretty(&std::collections::BTreeMap::from([(
                                "commands",
                                vec![Commands::Profile(p)],
                            )]))?;
                        return Err(eyre::eyre!(
                            "{verdict}\n{scope}\nReproduce in the same generated state with this run-card fragment:\n{reproduction}"
                        ));
                    }
                }
            }
            Commands::ThreeDRep(command) => {
                command.run(state, global_cli_settings)?;
            }
            Commands::Quit(s) => {
                return Ok(CommandExecution::break_with(s));
            }
            Commands::Inspect(inspect) => {
                if let Some(path) = &inspect.json_output {
                    global_cli_settings.ensure_write_target_outside_active_state(
                        path,
                        "write inspect JSON output",
                    )?;
                }
                let _ = inspect.run(state)?;
            }
            Commands::Approach(approach) => {
                let _ = approach.run(state, global_cli_settings)?;
            }
            Commands::Bench(bench) => {
                if let Some(path) = &bench.json_output {
                    global_cli_settings.ensure_write_target_outside_active_state(
                        path,
                        "write benchmark JSON output",
                    )?;
                }
                bench.run(state)?;
            }
            Commands::Import(s) => s.run(state, global_cli_settings)?,
            Commands::Save(s) => s.run(
                state,
                run_history,
                default_runtime_settings,
                global_cli_settings,
            )?,
            Commands::Set(s) => s.run(state, global_cli_settings, default_runtime_settings)?,
            Commands::Generate(g) => {
                let would_compile_into_active_state = global_cli_settings.session.read_only_state
                    && global_cli_settings.global.generation.evaluator.compile
                    && matches!(
                        g.mode.as_ref(),
                        None | Some(generate::GenerateCmd::Existing(_))
                    );
                if would_compile_into_active_state {
                    return Err(Report::msg(format!(
                        "Cannot compile generated integrands into '{}' because this session was started with --read-only-state. Disable `global.generation.evaluator.compile`, restart without --read-only-state, or save the state elsewhere first.",
                        global_cli_settings.state.folder.display()
                    )));
                }

                g.run(
                    state,
                    &global_cli_settings.state.folder,
                    global_cli_settings.override_state,
                    &global_cli_settings.global,
                    default_runtime_settings,
                )?
            }
            Commands::Select(s) => {
                s.run(state, global_cli_settings)?;
            }
            Commands::Integrate(g) => {
                return Ok(CommandExecution::continue_with(CommandOutput::Integrate(
                    g.run(state, global_cli_settings)?,
                )));
            }
            Commands::Evaluate(g) => {
                return Ok(CommandExecution::continue_with(CommandOutput::Evaluate(
                    g.run(state, global_cli_settings, default_runtime_settings)?,
                )));
            }
            Commands::Renormalize(r) => {
                _ = r.run(state, global_cli_settings)?;
            }
            Commands::Display(l) => {
                l.run(
                    state,
                    global_cli_settings,
                    default_runtime_settings,
                    run_history,
                )?;
            }
            Commands::Duplicate(command) => {
                command.run(state)?;
            }
            Commands::Run(r) => {
                return r.run(
                    state,
                    global_cli_settings,
                    default_runtime_settings,
                    run_history,
                );
            }
            Commands::StartCommandsBlock(_) | Commands::FinishCommandsBlock => {
                return Err(Report::msg(
                    "Command block definition commands must be handled by the CLI session",
                ));
            }
            Commands::Remove(r) => r.run(state)?,
            Commands::Batch {
                process_file: _process_file,
                batch_input_file: _batch_input_file,
                name: _name,
                output_name: _output_name,
            } => {
                todo!("Batch command not implemented yet");
            }
            Commands::Shell(s) => {
                s.run()?;
            }
            Commands::CommandTemplate(_) => {
                return Err(Report::msg(
                    "Command templates can only be executed through `run` with `-D/--define` variables",
                ));
            }
        }
        Ok(CommandExecution::continue_without_output())
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::{Approach, Bench, Commands, Inspect, Renormalize};
    use crate::{
        commands::generate::{Generate, GenerateCmd, ProcessArgs},
        state::{ProcessRef, RunHistory, State},
        CLISettings,
    };
    use gammalooprs::settings::RuntimeSettings;

    #[test]
    fn select_command_parses_graph_names() {
        let command: Commands =
            "select -p #0 -i default --with-only-graph-names GL02 GL03 --with-graph-names GL04 GL05 --without-graph-names GL06"
                .parse()
                .unwrap();
        match command {
            Commands::Select(select) => {
                assert_eq!(select.process, Some(ProcessRef::Id(0)));
                assert_eq!(select.integrand_name.as_deref(), Some("default"));
                assert_eq!(
                    select.with_graph_names,
                    vec!["GL04".to_string(), "GL05".to_string()]
                );
                assert_eq!(
                    select.with_only_graph_names,
                    vec!["GL02".to_string(), "GL03".to_string()]
                );
                assert_eq!(select.without_graph_names, vec!["GL06".to_string()]);
            }
            other => panic!("expected select command, got {other:?}"),
        }
    }

    #[test]
    fn select_command_requires_a_with_only_graph_name() {
        let error = "select -p #0 -i default --with-only-graph-names"
            .parse::<Commands>()
            .unwrap_err();
        assert!(error.to_string().contains("--with-only-graph-names"));
    }

    #[test]
    fn select_command_parses_filter_options() {
        let command: Commands = "select -p #0 -i default --amplitude-graphs --with-raised-propagator-signatures '[2]' --without-massive-raised-propagator-signatures '[]' --with-cycle-signatures '[(3,21)]' --without-cycle-signatures '[(ghost)]' --with-vertices '[V_6,V_9]' --without-vertices '[V_36]' --with-particles '[t,b]' '[g]' --without-particles '(e+,e-)'"
            .parse()
            .unwrap();
        match command {
            Commands::Select(select) => {
                assert!(select.amplitude_graphs);
                assert_eq!(
                    select.with_raised_propagator_signatures,
                    vec!["[2]".to_string()]
                );
                assert_eq!(
                    select.without_massive_raised_propagator_signatures,
                    vec!["[]".to_string()]
                );
                assert_eq!(select.with_cycle_signatures, vec!["[(3,21)]".to_string()]);
                assert_eq!(
                    select.without_cycle_signatures,
                    vec!["[(ghost)]".to_string()]
                );
                assert_eq!(select.with_vertices, vec!["[V_6,V_9]".to_string()]);
                assert_eq!(select.without_vertices, vec!["[V_36]".to_string()]);
                assert_eq!(
                    select.with_particles,
                    vec!["[t,b]".to_string(), "[g]".to_string()]
                );
                assert_eq!(select.without_particles, vec!["(e+,e-)".to_string()]);
            }
            other => panic!("expected select command, got {other:?}"),
        }
    }

    #[test]
    fn select_command_parses_output_targets() {
        let command: Commands = "select -p #0 -i default --with-graph-names GL04 --output_process selected_proc --output_integrand selected_itg --clear-existing-processes"
            .parse()
            .unwrap();
        match command {
            Commands::Select(select) => {
                assert_eq!(select.process, Some(ProcessRef::Id(0)));
                assert_eq!(select.integrand_name.as_deref(), Some("default"));
                assert_eq!(select.output_process.as_deref(), Some("selected_proc"));
                assert_eq!(select.output_integrand.as_deref(), Some("selected_itg"));
                assert!(select.clear_existing_processes);
            }
            other => panic!("expected select command, got {other:?}"),
        }
    }

    #[test]
    fn approach_command_parses_axes_spacing_and_output() {
        let command: Commands = "approach -p #0 -i default -x 0.5 0.25 --approach-axis 1.0,0.0 --approach-axis 0.0,1.0 --n-points 3 --skip-midpoint --logarithmic --min-abs-t 1e-4 --n-cores 2 --output-results approach.json"
            .parse()
            .unwrap();
        match command {
            Commands::Approach(approach) => {
                assert_eq!(approach.process, Some(ProcessRef::Id(0)));
                assert_eq!(approach.integrand_name.as_deref(), Some("default"));
                assert_eq!(approach.point, vec![0.5, 0.25]);
                assert_eq!(
                    approach.approach_axes,
                    vec!["1.0,0.0".to_string(), "0.0,1.0".to_string()]
                );
                assert_eq!(approach.n_points, 3);
                assert!(approach.skip_midpoint);
                assert!(approach.logarithmic);
                assert!(!approach.linear);
                assert_eq!(approach.min_abs_t, 1.0e-4);
                assert_eq!(approach.n_cores, Some(2));
                assert_eq!(
                    approach.output_results.as_deref(),
                    Some("approach.json".as_ref())
                );
            }
            other => panic!("expected approach command, got {other:?}"),
        }
    }

    #[test]
    fn approach_command_parses_momentum_space_selectors() {
        let command: Commands = "approach -p #0 -i default -x 1.0,-7.0e-2,0.0,-0.2,0.3,0.4 --approach-axis=-1.0e-3,0.0,0.0,0.0,0.0,0.0 --n-points 1 --linear --n-cores 1 --momentum-space --graph-id 2 --orientation-id 1"
            .parse()
            .unwrap();
        match command {
            Commands::Approach(approach) => {
                assert_eq!(approach.point, vec![1.0, -7.0e-2, 0.0, -0.2, 0.3, 0.4]);
                assert_eq!(approach.approach_axes, vec!["-1.0e-3,0.0,0.0,0.0,0.0,0.0"]);
                assert!(approach.momentum_space);
                assert!(approach.linear);
                assert_eq!(approach.graph_id, Some(2));
                assert_eq!(approach.orientation_id, Some(1));
                assert_eq!(approach.n_cores, Some(1));
            }
            other => panic!("expected approach command, got {other:?}"),
        }
    }

    #[test]
    fn generate_rejects_compilation_into_active_state_in_read_only_mode() {
        let mut state = State::new_test();
        let mut run_history = RunHistory::default();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();
        cli_settings.session.read_only_state = true;
        cli_settings.global.generation.evaluator.compile = true;

        let err = Commands::Generate(Generate {
            keep_sources: false,
            mode: Some(GenerateCmd::Existing(ProcessArgs {
                process: None,
                integrand_name: None,
            })),
        })
        .run(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut runtime_settings,
        )
        .unwrap_err();

        assert!(format!("{err:?}").contains("--read-only-state"));
    }

    #[test]
    fn inspect_rejects_json_output_inside_read_only_state_before_integrand_lookup() {
        let mut state = State::new_test();
        let mut run_history = RunHistory::default();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();
        cli_settings.state.folder = PathBuf::from("/tmp/read_only_state");
        cli_settings.session.read_only_state = true;

        let err = Commands::Inspect(Inspect {
            point: vec![0.1, 0.2],
            json_output: Some(PathBuf::from("/tmp/read_only_state/inspect.json")),
            ..Inspect::default()
        })
        .run(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut runtime_settings,
        )
        .unwrap_err();

        let rendered = format!("{err:?}");
        assert!(
            rendered.contains("Cannot write inspect JSON output"),
            "{rendered}"
        );
        assert!(rendered.contains("--read-only-state"), "{rendered}");
    }

    #[test]
    fn bench_rejects_json_output_inside_read_only_state_before_integrand_lookup() {
        let mut state = State::new_test();
        let mut run_history = RunHistory::default();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();
        cli_settings.state.folder = PathBuf::from("/tmp/read_only_state");
        cli_settings.session.read_only_state = true;

        let err = Commands::Bench(Bench {
            point: vec![0.1, 0.2],
            duration: "1ms".to_string(),
            json_output: Some(PathBuf::from("/tmp/read_only_state/bench.json")),
            ..Bench::default()
        })
        .run(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut runtime_settings,
        )
        .unwrap_err();

        let rendered = format!("{err:?}");
        assert!(
            rendered.contains("Cannot write benchmark JSON output"),
            "{rendered}"
        );
        assert!(rendered.contains("--read-only-state"), "{rendered}");
    }

    #[test]
    fn approach_rejects_output_inside_read_only_state_before_integrand_lookup() {
        let mut state = State::new_test();
        let mut run_history = RunHistory::default();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();
        cli_settings.state.folder = PathBuf::from("/tmp/read_only_state");
        cli_settings.session.read_only_state = true;

        let err = Commands::Approach(Approach {
            point: vec![0.1, 0.2],
            output_results: Some(PathBuf::from("/tmp/read_only_state/approach.json")),
            ..Approach::default()
        })
        .run(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut runtime_settings,
        )
        .unwrap_err();

        let rendered = format!("{err:?}");
        assert!(
            rendered.contains("Cannot write approach results"),
            "{rendered}"
        );
        assert!(rendered.contains("--read-only-state"), "{rendered}");
    }

    #[test]
    fn renormalize_rejects_output_inside_read_only_state_before_integrand_lookup() {
        let mut state = State::new_test();
        let mut run_history = RunHistory::default();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();
        cli_settings.state.folder = PathBuf::from("/tmp/read_only_state");
        cli_settings.session.read_only_state = true;

        let err = Commands::Renormalize(Renormalize {
            result_path: Some(PathBuf::from("/tmp/read_only_state/renormalization")),
            ..Renormalize::default()
        })
        .run(
            &mut state,
            &mut run_history,
            &mut cli_settings,
            &mut runtime_settings,
        )
        .unwrap_err();

        let rendered = format!("{err:?}");
        assert!(
            rendered.contains("Cannot write renormalization results"),
            "{rendered}"
        );
        assert!(rendered.contains("--read-only-state"), "{rendered}");
    }
}
