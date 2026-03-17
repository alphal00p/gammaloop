use std::{
    fs,
    io::{self, IsTerminal, Write},
    path::PathBuf,
};

use clap::Args;
use crossterm::{
    cursor::{MoveDown, MoveToColumn, MoveUp},
    queue,
    terminal::{self, Clear, ClearType},
};
use gammalooprs::utils::serde_utils::SmartSerde;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use unicode_width::UnicodeWidthChar;

use color_eyre::Result;
use colored::Colorize;
use gammalooprs::{
    integrate::{
        IntegrationState, IntegrationStatusKind, emit_integration_status_via_tracing,
        havana_integrate, print_integral_result,
    },
    settings::{RuntimeSettings, runtime::IntegrationResult},
    utils::F,
};
use tracing::{info, warn};

use crate::{
    CLISettings,
    completion::CompletionArgExt,
    state::{ProcessRef, State},
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, unsendable, name = "IntegrationSettings")
)]
#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Integrate {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// The integrand name to integrate
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// The path to store results in
    #[arg(short = 's', long, value_hint = clap::ValueHint::FilePath)]
    pub result_path: Option<PathBuf>,

    /// Number of cores to parallelize over
    #[arg(short = 'c', long)]
    pub n_cores: Option<usize>,

    /// The path to run the integrationg within
    #[arg(short = 'w', long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    /// Specify the target integration result to compare against
    #[arg(short = 't', num_args = 2, long,
          value_delimiter = ',',          // allow --vals 1,-2,3
          allow_negative_numbers = true,  // treat -2 as a value, not a flag
          allow_hyphen_values=true
    )]
    pub target: Option<Vec<f64>>,

    /// Whether to restart the integration from scratch, or continue from a previous run if possible
    #[arg(short = 'r', long)]
    pub restart: bool,

    /// Hide max-weight information during intermediate iteration updates
    #[arg(long = "no-show-max-weights")]
    pub no_show_max_weights: bool,

    /// Hide integration statistics during intermediate iteration updates
    #[arg(long = "no-show-integration-statistics")]
    pub no_show_integration_statistics: bool,

    /// Disable live streaming of intermediate integration updates in interactive terminals
    #[arg(long = "no-stream")]
    pub no_stream: bool,
}

struct StreamRenderer {
    stderr: io::Stderr,
    rendered_line_count: usize,
}

impl StreamRenderer {
    fn new() -> Self {
        Self {
            stderr: io::stderr(),
            rendered_line_count: 0,
        }
    }

    fn render(&mut self, block: &str) -> Result<()> {
        let prepared = prepare_stream_block(block, stream_terminal_width());
        if self.rendered_line_count > 0 {
            self.clear_rendered_block()?;
        }

        queue!(self.stderr, MoveToColumn(0))?;
        write!(self.stderr, "{prepared}")?;
        self.stderr.flush()?;
        self.rendered_line_count = prepared.lines().count().max(1);
        Ok(())
    }

    fn clear(&mut self) -> Result<()> {
        if self.rendered_line_count == 0 {
            return Ok(());
        }

        self.clear_rendered_block()?;
        self.rendered_line_count = 0;
        self.stderr.flush()?;
        Ok(())
    }

    fn clear_rendered_block(&mut self) -> Result<()> {
        queue!(self.stderr, MoveToColumn(0))?;
        if self.rendered_line_count > 1 {
            queue!(
                self.stderr,
                MoveUp((self.rendered_line_count.saturating_sub(1)) as u16)
            )?;
        }

        for line_index in 0..self.rendered_line_count {
            queue!(self.stderr, MoveToColumn(0), Clear(ClearType::CurrentLine))?;
            if line_index + 1 < self.rendered_line_count {
                queue!(self.stderr, MoveDown(1))?;
            }
        }
        queue!(self.stderr, MoveToColumn(0))?;
        if self.rendered_line_count > 1 {
            queue!(
                self.stderr,
                MoveUp((self.rendered_line_count.saturating_sub(1)) as u16)
            )?;
        }
        Ok(())
    }
}

impl Drop for StreamRenderer {
    fn drop(&mut self) {
        let _ = self.clear();
    }
}

fn stream_terminal_width() -> usize {
    terminal::size()
        .map(|(width, _)| width.saturating_sub(1).max(1) as usize)
        .unwrap_or(120)
}

fn prepare_stream_block(block: &str, max_width: usize) -> String {
    block
        .lines()
        .map(|line| truncate_ansi_line(line, max_width))
        .collect::<Vec<_>>()
        .join("\n")
}

fn truncate_ansi_line(line: &str, max_width: usize) -> String {
    if max_width == 0 {
        return String::new();
    }

    let mut truncated = String::new();
    let mut chars = line.chars().peekable();
    let mut visible_width = 0usize;
    let mut saw_escape = false;
    let mut was_truncated = false;

    while let Some(ch) = chars.next() {
        if ch == '\u{1b}' && chars.peek() == Some(&'[') {
            saw_escape = true;
            truncated.push(ch);
            truncated.push(chars.next().unwrap());
            for code in chars.by_ref() {
                truncated.push(code);
                if ('@'..='~').contains(&code) {
                    break;
                }
            }
            continue;
        }

        let char_width = UnicodeWidthChar::width(ch).unwrap_or(0);
        if visible_width + char_width > max_width {
            was_truncated = true;
            break;
        }

        truncated.push(ch);
        visible_width += char_width;
    }

    if was_truncated && saw_escape && !truncated.ends_with("\u{1b}[0m") {
        truncated.push_str("\u{1b}[0m");
    }

    truncated
}

impl Integrate {
    fn default_workspace_path(&self, global_cli_settings: &CLISettings) -> PathBuf {
        if global_cli_settings.session.read_only_state {
            let workspace_name = global_cli_settings
                .state
                .name
                .as_deref()
                .map(str::trim)
                .filter(|name| !name.is_empty())
                .map(|name| format!("{name}_integration_workspace"))
                .unwrap_or_else(|| "integration_workspace".to_string());
            PathBuf::from(".").join(workspace_name)
        } else {
            global_cli_settings
                .state
                .folder
                .join("integration_workspace")
        }
    }

    pub fn run(
        &self,
        state: &mut State,
        global_cli_settings: &CLISettings,
    ) -> Result<IntegrationResult> {
        let default_workspace_path = self.default_workspace_path(global_cli_settings);
        let workspace_path = if let Some(p) = self.workspace_path.clone() {
            p
        } else {
            default_workspace_path.clone()
        };

        let result_path = if let Some(p) = self.result_path.clone() {
            p
        } else {
            default_workspace_path.join("integration_result.json")
        };

        let target = self.target.clone().map(|t| Complex::new(F(t[0]), F(t[1])));

        let (process_id, integrand_name) =
            state.find_integrand_ref(self.process.as_ref(), self.integrand_name.as_ref())?;

        if self.restart && workspace_path.exists() {
            fs::remove_dir_all(&workspace_path)?;
        }

        let gloop_integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;

        gloop_integrand.warm_up(&state.model)?;

        info!("Gammaloop now integrates {}", integrand_name.green().bold());

        let path_to_state = workspace_path.join("integration_state");

        let integration_state = match fs::read(path_to_state) {
            Ok(state_bytes) => {
                info!(
                    "{}",
                    "Found integration state, result of previous integration:".yellow()
                );
                info!("");

                let state: IntegrationState = bincode::decode_from_slice::<IntegrationState, _>(
                    &state_bytes,
                    bincode::config::standard(),
                )
                .expect("Could not deserialize state")
                .0;

                let path_to_workspace_settings = workspace_path.join("settings.yaml");

                let workspace_settings: RuntimeSettings =
                    RuntimeSettings::from_file(path_to_workspace_settings, "workspace settings")?;
                // force the settings to be the same as the ones used in the previous integration
                if *gloop_integrand.get_mut_settings() != workspace_settings.clone() {
                    warn!("settings have changed with respect to workspace, reverting changes");
                    *gloop_integrand.get_mut_settings() = workspace_settings.clone();
                }

                print_integral_result(
                    &state.integral.re,
                    1,
                    state.iter,
                    "re",
                    target.map(|c| c.re),
                );

                print_integral_result(
                    &state.integral.im,
                    2,
                    state.iter,
                    "im",
                    target.map(|c| c.im),
                );
                info!("");
                warn!(
                    "Any changes to the settings will be ignored, integrate with the {} option for changes to take effect",
                    "--restart".blue()
                );
                info!("{}", "Resuming integration".yellow());

                Some(state)
            }

            Err(_) => {
                info!("No integration state found, starting new integration");
                None
            }
        };

        if !workspace_path.exists() {
            fs::create_dir_all(&workspace_path)?;
            info!(
                "Created workspace directory at {}",
                workspace_path.display()
            );
        }
        let settings = gloop_integrand.get_settings().clone();

        let n_cores = self
            .n_cores
            .unwrap_or(global_cli_settings.global.n_cores.integrate);

        let enable_streaming = if self.no_stream {
            false
        } else if io::stderr().is_terminal() {
            true
        } else {
            info!(
                "Streaming integration updates disabled because stderr is not a TTY; falling back to append-only logging."
            );
            false
        };
        let mut stream_renderer = enable_streaming.then(StreamRenderer::new);

        let result = havana_integrate(
            &settings,
            &state.model,
            |set| gloop_integrand.user_data_generator(n_cores, set),
            target,
            integration_state,
            Some(workspace_path.clone()),
            !self.no_show_max_weights,
            !self.no_show_integration_statistics,
            move |kind, status_block| {
                if let Some(renderer) = stream_renderer.as_mut() {
                    match kind {
                        IntegrationStatusKind::Iteration => renderer.render(&status_block)?,
                        IntegrationStatusKind::Final => {
                            renderer.clear()?;
                            emit_integration_status_via_tracing(kind, &status_block)?;
                        }
                    }
                } else {
                    emit_integration_status_via_tracing(kind, &status_block)?;
                }

                Ok(())
            },
        )?;

        fs::write(&result_path, serde_json::to_string(&result)?)?;

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::Integrate;
    use crate::{CLISettings, SessionSettings, StateSettings};
    use std::path::PathBuf;

    #[test]
    fn read_only_state_uses_cwd_workspace_default() {
        let integrate = Integrate::default();
        let mut settings = CLISettings::default();
        settings.state = StateSettings {
            folder: PathBuf::from("/tmp/saved_state"),
            name: Some("gg_hhh_1l".to_string()),
        };
        settings.session = SessionSettings {
            read_only_state: true,
            ..SessionSettings::default()
        };

        assert_eq!(
            integrate.default_workspace_path(&settings),
            PathBuf::from("./gg_hhh_1l_integration_workspace")
        );

        settings.state.name = None;
        assert_eq!(
            integrate.default_workspace_path(&settings),
            PathBuf::from("./integration_workspace")
        );
    }

    #[test]
    fn truncate_ansi_line_limits_visible_width() {
        let line = "\u{1b}[32mhello\u{1b}[0m world";
        let truncated = super::truncate_ansi_line(&line, 7);

        assert!(truncated.contains('\u{1b}'));
        assert!(truncated.ends_with("\u{1b}[0m"));
        assert!(!truncated.contains("world"));
    }

    #[test]
    fn prepare_stream_block_truncates_each_line_independently() {
        let block = "123456789\nabcdefghi";
        let prepared = super::prepare_stream_block(block, 5);

        assert_eq!(prepared, "12345\nabcde");
    }
}
