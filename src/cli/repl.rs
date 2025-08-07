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

    Ok(())
}
