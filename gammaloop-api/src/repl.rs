use std::{ffi::OsString, marker::PhantomData, path::Path};

use clap::Parser;
use console::style;

use reedline::{Prompt, Reedline, Signal, Span};
use shlex::Shlex;

mod builder {
    use std::marker::PhantomData;

    use clap::Parser;
    use nu_ansi_term::{Color, Style};
    use reedline::{
        default_emacs_keybindings, DefaultHinter, DefaultPrompt, EditMode, Emacs, IdeMenu,
        KeyModifiers, MenuBuilder, Prompt, Reedline, ReedlineEvent, ReedlineMenu,
    };

    use crate::repl::{ClapEditor, ReedCompleter};

    pub struct ClapEditorBuilder<C: Parser + Send + Sync + 'static> {
        prompt: Box<dyn Prompt>,
        edit_mode: Box<dyn EditMode>,
        hook: Box<dyn FnOnce(Reedline) -> Reedline>,
        c_phantom: PhantomData<C>,
    }

    impl<C: Parser + Send + Sync + 'static> ClapEditorBuilder<C> {
        pub(crate) fn new() -> Self {
            Self {
                prompt: Box::<DefaultPrompt>::default(),
                edit_mode: {
                    let mut keybindings = default_emacs_keybindings();
                    keybindings.add_binding(
                        KeyModifiers::NONE,
                        reedline::KeyCode::Tab,
                        ReedlineEvent::UntilFound(vec![
                            ReedlineEvent::Menu("completion_menu".to_string()),
                            ReedlineEvent::MenuNext,
                        ]),
                    );
                    Box::new(Emacs::new(keybindings))
                },
                hook: Box::new(|e| e),
                c_phantom: PhantomData,
            }
        }

        pub fn with_prompt(mut self, prompt: Box<dyn Prompt>) -> Self {
            self.prompt = prompt;
            self
        }

        pub fn with_edit_mode(mut self, edit_mode: Box<dyn EditMode>) -> Self {
            self.edit_mode = edit_mode;
            self
        }

        pub fn with_editor_hook(
            mut self,
            hook: impl FnOnce(Reedline) -> Reedline + 'static,
        ) -> Self {
            self.hook = Box::new(hook);
            self
        }

        pub fn build(self) -> ClapEditor<C> {
            let completion_menu = Box::new(
                IdeMenu::default()
                    .with_default_border()
                    .with_name("completion_menu"),
            );

            let rl = Reedline::create()
                .with_completer(Box::new(ReedCompleter::<C> {
                    c_phantom: PhantomData,
                }))
                .with_menu(ReedlineMenu::EngineCompleter(completion_menu))
                .with_hinter(Box::new(
                    DefaultHinter::default().with_style(Style::new().italic().fg(Color::DarkGray)),
                ))
                .with_edit_mode(self.edit_mode);
            let rl = (self.hook)(rl);
            ClapEditor {
                rl,
                prompt: self.prompt,
                c_phantom: PhantomData,
            }
        }
    }
}

pub use builder::ClapEditorBuilder;

pub struct ClapEditor<C: Parser + Send + Sync + 'static> {
    rl: Reedline,
    prompt: Box<dyn Prompt>,
    c_phantom: PhantomData<C>,
}

struct ReedCompleter<C: Parser + Send + Sync + 'static> {
    c_phantom: PhantomData<C>,
}

impl<C: Parser + Send + Sync + 'static> reedline::Completer for ReedCompleter<C> {
    fn complete(&mut self, line: &str, pos: usize) -> Vec<reedline::Suggestion> {
        let args = Shlex::new(line);
        let mut args: Vec<OsString> = std::iter::once(OsString::new())
            .chain(args.map(OsString::from))
            .collect();

        if line.ends_with(' ') {
            args.push(OsString::new());
        }

        let mut suggestions = Vec::new();
        let root_cmd = C::command();
        let mut cmd = &root_cmd;

        // Navigate through completed subcommands to find current context
        let mut cmd_args = &args[1..]; // Skip the program name
        let mut remaining_args = cmd_args; // Track remaining arguments for current command context

        // Process all completed arguments (not including the one we're currently typing)
        while cmd_args.len() > 1 {
            let arg = &cmd_args[0];
            let arg_str = arg.to_string_lossy();

            // Check if this completed argument is a subcommand
            if let Some(subcmd) = cmd.get_subcommands().find(|sc| sc.get_name() == arg_str) {
                cmd = subcmd;
                cmd_args = &cmd_args[1..];
                remaining_args = cmd_args; // Update remaining args for the new command context
            } else {
                // This argument is not a subcommand, so we stop navigating
                break;
            }
        }

        // Handle the special case where we just completed a subcommand and added a space
        if line.ends_with(' ') && cmd_args.len() == 1 {
            let arg = &cmd_args[0];
            let arg_str = arg.to_string_lossy();

            // If the last completed argument is a subcommand, navigate into it
            if let Some(subcmd) = cmd.get_subcommands().find(|sc| sc.get_name() == arg_str) {
                cmd = subcmd;
                cmd_args = &[]; // Now we're starting fresh in the new subcommand
                remaining_args = &[]; // No remaining args in the new context
            }
        }

        // Determine what we're currently completing
        let (current_arg, start_pos) = if cmd_args.is_empty() {
            // We're starting a new argument
            (&OsString::new(), pos)
        } else {
            // We're completing the last argument
            let current_arg = &cmd_args[cmd_args.len() - 1];
            let start_pos = pos.saturating_sub(current_arg.len());
            (current_arg, start_pos)
        };

        let current_arg_str = current_arg.to_string_lossy();

        // If current argument starts with '-', suggest flags
        if current_arg_str.starts_with('-') {
            for arg in cmd.get_arguments() {
                // Short flags
                if let Some(short) = arg.get_short() {
                    let flag = format!("-{}", short);
                    if flag.starts_with(current_arg_str.as_ref()) {
                        suggestions.push(reedline::Suggestion {
                            value: flag,
                            description: arg.get_help().map(|s| s.to_string()),
                            style: None,
                            extra: None,
                            span: Span::new(start_pos, pos),
                            append_whitespace: true,
                        });
                    }
                }

                // Long flags
                if let Some(long) = arg.get_long() {
                    let flag = format!("--{}", long);
                    if flag.starts_with(current_arg_str.as_ref()) {
                        suggestions.push(reedline::Suggestion {
                            value: flag,
                            description: arg.get_help().map(|s| s.to_string()),
                            style: None,
                            extra: None,
                            span: Span::new(start_pos, pos),
                            append_whitespace: true,
                        });
                    }
                }
            }
        } else {
            // Check if we should complete paths first
            if should_complete_paths(cmd, remaining_args) {
                if let Some(path_suggestions) = complete_path(&current_arg_str, start_pos, pos) {
                    suggestions.extend(path_suggestions);
                }
            } else {
                // Only suggest subcommands if we're not completing a flag value
                for subcommand in cmd.get_subcommands() {
                    let name = subcommand.get_name();
                    if name.starts_with(current_arg_str.as_ref()) {
                        suggestions.push(reedline::Suggestion {
                            value: name.to_string(),
                            description: subcommand.get_about().map(|s| s.to_string()),
                            style: None,
                            extra: None,
                            span: Span::new(start_pos, pos),
                            append_whitespace: true,
                        });
                    }
                }
            }
        }

        suggestions
    }
}

/// Check if we should complete paths based on the current argument position
fn should_complete_paths(cmd: &clap::Command, cmd_args: &[std::ffi::OsString]) -> bool {
    // If we're currently completing after a flag that expects a path value
    if cmd_args.len() >= 2 {
        // Check if the second-to-last argument was a flag that expects a path
        let flag_arg = &cmd_args[cmd_args.len() - 2];
        let flag_arg_str = flag_arg.to_string_lossy();

        if flag_arg_str.starts_with('-') {
            // We're completing the value after a flag
            let flag_expects_path = cmd
                .get_arguments()
                .find(|a| {
                    // Check long flag
                    let long_matches = if let Some(long) = a.get_long() {
                        let long_flag = format!("--{}", long);
                        flag_arg_str == long_flag
                    } else {
                        false
                    };

                    // Check short flag
                    let short_matches = if let Some(short) = a.get_short() {
                        let short_flag = format!("-{}", short);
                        flag_arg_str == short_flag
                    } else {
                        false
                    };

                    long_matches || short_matches
                })
                .map(|arg| is_path_argument(arg))
                .unwrap_or(false);
            if flag_expects_path {
                return true;
            }
        }
    }

    // Count non-flag positional arguments to determine which positional argument we're completing
    let mut positional_count: usize = 0;
    for arg in cmd_args {
        let arg_str = arg.to_string_lossy();
        if !arg_str.starts_with('-') {
            positional_count += 1;
        }
    }

    // Get the positional argument definition for this position (0-based)
    let positional_args: Vec<_> = cmd
        .get_arguments()
        .filter(|arg| !arg.get_long().is_some() && !arg.get_short().is_some())
        .collect();

    if let Some(positional_arg) = positional_args.get(positional_count.saturating_sub(1)) {
        is_path_argument(positional_arg)
    } else {
        false
    }
}

/// Determine if an argument expects a path based on its type and hints
fn is_path_argument(arg: &clap::Arg) -> bool {
    // Check value hint for path-related hints
    match arg.get_value_hint() {
        clap::ValueHint::FilePath
        | clap::ValueHint::DirPath
        | clap::ValueHint::AnyPath
        | clap::ValueHint::ExecutablePath => true,
        _ => {
            // Check if the argument name suggests it's a path
            let name = arg.get_id().as_str().to_lowercase();
            name.contains("path")
                || name.contains("file")
                || name.contains("dir")
                || name.contains("directory")
                || name == "input"
                || name == "output"
        }
    }
}

/// Complete file and directory paths
fn complete_path(
    partial_path: &str,
    start_pos: usize,
    pos: usize,
) -> Option<Vec<reedline::Suggestion>> {
    let mut suggestions = Vec::new();

    // Determine the directory to search and the prefix to match
    let path = Path::new(partial_path);
    let (search_dir, file_prefix) = if partial_path.ends_with('/') || partial_path.ends_with('\\') {
        // Complete from the specified directory
        (path.to_path_buf(), String::new())
    } else if let Some(parent) = path.parent() {
        // Complete from parent directory with filename prefix
        if parent.as_os_str().is_empty() {
            (
                std::path::PathBuf::from("."),
                path.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string(),
            )
        } else {
            (
                parent.to_path_buf(),
                path.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string(),
            )
        }
    } else {
        // Complete from current directory
        (std::path::PathBuf::from("."), partial_path.to_string())
    };

    // Try to read the directory
    match std::fs::read_dir(&search_dir) {
        Ok(entries) => {
            for entry in entries.flatten() {
                let file_name = entry.file_name();
                let file_name_str = file_name.to_string_lossy();

                // Skip hidden files unless we're explicitly typing them
                if file_name_str.starts_with('.') && !file_prefix.starts_with('.') {
                    continue;
                }

                // Check if this entry matches our prefix
                if file_name_str.starts_with(&file_prefix) {
                    let is_dir = entry.file_type().map(|ft| ft.is_dir()).unwrap_or(false);
                    let mut full_path = if search_dir.to_string_lossy() == "." {
                        file_name_str.to_string()
                    } else {
                        search_dir.join(&file_name).to_string_lossy().to_string()
                    };

                    // Add trailing slash for directories
                    if is_dir {
                        full_path.push('/');
                    }

                    let description = if is_dir {
                        Some("Directory".to_string())
                    } else {
                        // Try to get file size for description
                        entry
                            .metadata()
                            .ok()
                            .map(|meta| format!("File ({} bytes)", meta.len()))
                    };

                    suggestions.push(reedline::Suggestion {
                        value: full_path,
                        description,
                        style: None,
                        extra: None,
                        span: Span::new(start_pos, pos),
                        append_whitespace: !is_dir, // Don't add space after directories
                    });
                }
            }
        }
        Err(_e) => {
            // Directory doesn't exist or can't be read - this is normal
        }
    }

    if suggestions.is_empty() {
        None
    } else {
        // Sort suggestions: directories first, then files, both alphabetically
        suggestions.sort_by(|a, b| {
            let a_is_dir = a.value.ends_with('/');
            let b_is_dir = b.value.ends_with('/');

            match (a_is_dir, b_is_dir) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                _ => a.value.cmp(&b.value),
            }
        });

        Some(suggestions)
    }
}

/// Result of reading a command from the REPL.
///
pub enum ReadCommandOutput<C> {
    /// Input parsed successfully. Contains the parsed command and the raw input string.
    Command(C, String),

    /// Input was empty.
    EmptyLine,

    /// Clap parse error happened. You should print the error manually.
    ClapError(clap::error::Error),

    /// Input was not lexically valid, for example it had odd number of `"`
    ShlexError,

    /// Reedline failed to work with stdio.
    ReedlineError(std::io::Error),

    /// User pressed ctrl+C
    CtrlC,

    /// User pressed ctrl+D
    CtrlD,
}

impl<C> ReadCommandOutput<C> {
    /// Extract just the parsed command, discarding the raw input.
    /// Returns `None` if this is not a `Command` variant.
    pub fn command(self) -> Option<C> {
        match self {
            ReadCommandOutput::Command(cmd, _) => Some(cmd),
            _ => None,
        }
    }

    /// Extract both the parsed command and raw input.
    /// Returns `None` if this is not a `Command` variant.
    pub fn command_with_raw(self) -> Option<(C, String)> {
        match self {
            ReadCommandOutput::Command(cmd, raw) => Some((cmd, raw)),
            _ => None,
        }
    }

    /// Get the raw input string if this is a `Command` variant.
    /// Returns `None` if this is not a `Command` variant.
    pub fn raw_input(&self) -> Option<&str> {
        match self {
            ReadCommandOutput::Command(_, raw) => Some(raw),
            _ => None,
        }
    }
}

impl<C: Parser + Send + Sync + 'static> ClapEditor<C> {
    pub fn builder() -> ClapEditorBuilder<C> {
        ClapEditorBuilder::<C>::new()
    }

    pub fn get_editor(&mut self) -> &mut Reedline {
        &mut self.rl
    }

    pub fn set_prompt(&mut self, prompt: Box<dyn Prompt>) {
        self.prompt = prompt;
    }

    pub fn read_command(&mut self) -> ReadCommandOutput<C> {
        let line = match self.rl.read_line(&*self.prompt) {
            Ok(Signal::Success(buffer)) => buffer,
            Ok(Signal::CtrlC) => return ReadCommandOutput::CtrlC,
            Ok(Signal::CtrlD) => return ReadCommandOutput::CtrlD,
            Err(e) => return ReadCommandOutput::ReedlineError(e),
        };
        if line.trim().is_empty() {
            return ReadCommandOutput::EmptyLine;
        }

        // _ = self.rl.add_history_entry(line.as_str());

        match shlex::split(&line) {
            Some(split) => {
                match C::try_parse_from(std::iter::once("").chain(split.iter().map(String::as_str)))
                {
                    Ok(c) => ReadCommandOutput::Command(c, line),
                    Err(e) => ReadCommandOutput::ClapError(e),
                }
            }
            None => ReadCommandOutput::ShlexError,
        }
    }

    pub fn repl(mut self, mut handler: impl FnMut(C, String)) {
        loop {
            match self.read_command() {
                ReadCommandOutput::Command(c, raw) => handler(c, raw),
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

    /// Convenience method for backward compatibility - runs REPL with handler that only receives parsed command
    #[allow(unused_mut)]
    pub fn repl_simple(mut self, handler: impl FnMut(C)) {
        let mut handler = handler;
        self.repl(|cmd, _raw| handler(cmd))
    }

    pub async fn repl_async<F, Fut>(mut self, mut handler: F)
    where
        F: FnMut(C, String) -> Fut,
        Fut: std::future::Future<Output = ()>,
    {
        loop {
            match self.read_command() {
                ReadCommandOutput::Command(c, raw) => handler(c, raw).await,
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

    /// Convenience method for backward compatibility - runs async REPL with handler that only receives parsed command
    #[allow(unused_mut)]
    pub async fn repl_simple_async<F, Fut>(mut self, handler: F)
    where
        F: FnMut(C) -> Fut,
        Fut: std::future::Future<Output = ()>,
    {
        let mut handler = handler;
        self.repl_async(|cmd, _raw| handler(cmd)).await
    }
}
