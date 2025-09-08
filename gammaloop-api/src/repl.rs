use std::{ffi::OsString, marker::PhantomData};

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

        let arg_index = args.len() - 1;
        let current_arg = if arg_index > 0 && arg_index < args.len() {
            &args[arg_index]
        } else {
            &OsString::new()
        };

        // Use basic completion for subcommands and flags
        let mut suggestions = Vec::new();
        let cmd = C::command();

        // Simple completion: suggest subcommands if no args yet
        if args.len() <= 1 || line.trim().is_empty() {
            for subcommand in cmd.get_subcommands() {
                let name = subcommand.get_name();
                if name.starts_with(current_arg.to_string_lossy().as_ref()) {
                    let start_pos = pos.saturating_sub(current_arg.len());
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

        suggestions
    }
}

/// Result of reading a command from the REPL.
///
/// # Example
///
/// ```rust
/// use clap::Parser;
/// use gammaloop_api::repl::ClapEditor;
///
/// #[derive(Parser)]
/// struct MyCli {
///     #[arg(short, long)]
///     verbose: bool,
///     command: String,
/// }
///
/// let editor = ClapEditor::<MyCli>::builder().build();
///
/// // Option 1: Handle both parsed command and raw input
/// editor.repl(|parsed_cmd, raw_input| {
///     println!("User typed: '{}'", raw_input);
///     println!("Parsed verbose flag: {}", parsed_cmd.verbose);
///     println!("Command: {}", parsed_cmd.command);
/// });
///
/// // Option 2: Handle only parsed command (backward compatibility)
/// editor.repl_simple(|parsed_cmd| {
///     println!("Command: {}", parsed_cmd.command);
/// });
///
/// // Option 3: Manual handling with ReadCommandOutput
/// loop {
///     match editor.read_command() {
///         ReadCommandOutput::Command(cmd, raw) => {
///             // Do something with both cmd and raw
///         },
///         ReadCommandOutput::EmptyLine => continue,
///         ReadCommandOutput::CtrlD => break,
///         // ... handle other cases
///         _ => {}
///     }
/// }
/// ```
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
