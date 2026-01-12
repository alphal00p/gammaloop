use clap::Args;
use color_eyre::Result;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct Shell {
    // Capture everything after '!' as raw tokens
    #[arg(
        trailing_var_arg = true,
        allow_hyphen_values = true,
        value_name = "CMD",
        num_args = 1..,
    )]
    cmd: Vec<std::ffi::OsString>,
}

impl Shell {
    pub fn run(self) -> Result<()> {
        use std::process::Command;

        // Pick a shell. macOS typically has /bin/zsh; fallback to /bin/sh
        let shell = std::env::var("SHELL").unwrap_or_else(|_| "/bin/sh".to_string());

        // Turn args back into a safe shell string: 'abc' -> 'abc', foo bar -> 'foo' 'bar'
        fn sh_quote(arg: &std::ffi::OsString) -> String {
            // Convert OsStr lossily and single-quote. Replace ' with '\'' (close-quote + escaped quote + reopen)
            let s = arg.to_string_lossy();
            if s.is_empty() {
                "''".to_string()
            } else {
                format!("'{}'", s.replace('\'', "'\"'\"'"))
            }
        }

        let cmd_string = self.cmd.iter().map(sh_quote).collect::<Vec<_>>().join(" ");

        let status = Command::new(shell)
            .arg("-lc") // login + run command; drop -l if you don’t want login semantics
            .arg(cmd_string)
            .status()?;

        if !status.success() {
            return Err(eyre::eyre!(
                "Shell command exited with status {:?}",
                status.code()
            ));
        }
        Ok(())
    }
}
