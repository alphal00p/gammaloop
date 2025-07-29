//! **clap‑repl**‑based REPL for *gammaLoop*
//!
//! Replaces the hand‑rolled rustyline helper with the off‑the‑shelf
//! [`clap-repl`](https://crates.io/crates/clap-repl) crate.  This gives you
//! completions, history, and command dispatch with far less code – and no more
//! borrow‑checker gymnastics.
//!
//! ---------------------------------------------------------------------------
//! ### `Cargo.toml` additions
//! ```toml
//! clap-repl   = "0.4"         # or the latest 0.x version
//! dirs        = "5"           # cross‑platform home‑dir helper (optional)
//! ```
//! (You can drop `rustyline` and `shell-words` from *your* dependency list – the
//! `clap-repl` crate re‑exports everything it needs.)
//!
//! ---------------------------------------------------------------------------
//! Assumptions:
//! * `_gammaloop::cli_functions` exposes
//!   ```rust
//!   pub use crate::cli_functions::{Cli, run_cli};
//!   ```
//!   where `Cli` is your `#[derive(Parser)]` root and `run_cli(parsed)` contains
//!   the heavy dispatch logic.
//!
//! * You **optionally** collect global flags (cores, config, …) into `pre_args`
//!   in `main.rs` and forward them here.  They’ll be inserted **before** the
//!   user’s input so every command inherits them.
//!
//! ---------------------------------------------------------------------------
//! There’s essentially no state: `clap-repl` handles history loading/saving and
//! completion for you.

use std::ops::ControlFlow;

use crate::{cli::Cli, Settings};
use clap::{CommandFactory, Parser};
use clap_repl::{
    reedline::{DefaultPrompt, DefaultPromptSegment, FileBackedHistory},
    ClapEditor, ReadCommandOutput,
};
use color_eyre::Report;
use console::style;
use dirs::home_dir;
use log::warn;

use super::State;

/// Launch the REPL.
/// * `prompt`   – prompt prefix (e.g. "gammaLoop").
/// * `pre_args` – global flags to prepend to every command (may be empty).
pub fn start(mut settings: Settings) -> Result<(), Report> {
    // 1. Create the base clap::Command from your `Cli` type.

    let prompt = DefaultPrompt {
        left_prompt: DefaultPromptSegment::Basic("γloop".to_owned()),
        ..DefaultPrompt::default()
    };
    let mut state = State::default();
    warn!("Starting REPL");

    // 2. Build the REPL – clap‑repl takes ownership and configures rustyline.
    let mut repl = ClapEditor::<Cli>::builder().with_prompt(Box::new(prompt));

    if let Some(home) = home_dir() {
        repl = repl.with_editor_hook(move |reed| {
            // Do custom things with `Reedline` instance here
            reed.with_history(Box::new(
                FileBackedHistory::with_file(10000, home.join(".gammaLoop_history")).unwrap(),
            ))
        })
    }
    let mut r = repl.build();

    loop {
        match r.read_command() {
            ReadCommandOutput::Command(c) => match c.run_with_settings(&mut settings, &mut state) {
                Err(e) => eprintln!("{e}"),
                Ok(ControlFlow::Break(())) => {
                    break;
                }
                _ => {}
            },
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

    Ok(())
}
