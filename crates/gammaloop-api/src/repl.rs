use std::{
    collections::{BTreeMap, HashSet},
    marker::PhantomData,
    path::Path,
    sync::OnceLock,
    sync::{Arc, RwLock},
};

use clap::Parser;
use console::style;
use gammalooprs::model::ParameterType;
use gammalooprs::settings::RuntimeSettings;
use serde_json::Value as JsonValue;

use reedline::{Prompt, Reedline, Signal, Span, Suggestion, ValidationResult};

use crate::{
    command_parser::split_command_line,
    commands::import::model::{builtin_json_model_names, builtin_json_model_restriction_names},
    commands::process_settings::{
        observable_completion_root, observable_schema, quantity_completion_root_for_kind,
        quantity_kind_names, quantity_schema, selector_completion_root_for_kind,
        selector_kind_names, selector_schema, NamedProcessSettingKind,
        ProcessSettingsCompletionEntry,
    },
    completion::{arg_value_completion, ArgValueCompletion, SelectorKind},
    session::CliSession,
    settings_tree::{
        schema_at_path, schema_enum_values, serialize_schema, serialize_settings_with_defaults,
        value_at_path,
    },
    CLISettings,
};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CompletionState {
    commands_block_names: Vec<String>,
    process_entries: Vec<ProcessCompletionEntry>,
    process_settings_entries: Vec<ProcessSettingsCompletionEntry>,
    ir_profile_entries: Vec<IrProfileCompletionEntry>,
    model_parameter_entries: Vec<ModelParameterCompletionEntry>,
    model_particle_names: Vec<String>,
    model_coupling_names: Vec<String>,
    model_vertex_names: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelParameterCompletionEntry {
    pub name: String,
    pub parameter_type: ParameterType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcessKind {
    Amplitude,
    CrossSection,
}

impl ProcessKind {
    fn label(self) -> &'static str {
        match self {
            ProcessKind::Amplitude => "amplitudes",
            ProcessKind::CrossSection => "cross sections",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcessCompletionEntry {
    pub id: usize,
    pub name: String,
    pub kind: ProcessKind,
    pub integrand_names: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IrProfileCompletionEntry {
    pub process_name: String,
    pub integrand_name: String,
    pub graph_names: Vec<String>,
    pub graph_limit_entries: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct SharedCompletionState(Arc<RwLock<CompletionState>>);

impl SharedCompletionState {
    pub fn new() -> Self {
        Self(Arc::new(RwLock::new(CompletionState::default())))
    }

    fn snapshot(&self) -> CompletionState {
        self.0.read().map(|state| state.clone()).unwrap_or_default()
    }

    fn write(&self, update: impl FnOnce(&mut CompletionState)) {
        if let Ok(mut state) = self.0.write() {
            update(&mut state);
        }
    }

    pub fn update_from_session(&self, session: &CliSession<'_>) {
        self.write(|state| {
            state.commands_block_names = session.current_commands_block_names();
            state.process_entries = session.current_process_entries();
            state.process_settings_entries = session.current_process_settings_entries();
            state.ir_profile_entries = session.current_ir_profile_entries();
            state.model_parameter_entries = session.current_model_parameter_entries();
            state.model_particle_names = session.current_model_particle_names();
            state.model_coupling_names = session.current_model_coupling_names();
            state.model_vertex_names = session.current_model_vertex_names();
        });
    }
}

mod builder {
    use std::marker::PhantomData;

    use clap::Parser;
    use nu_ansi_term::{Color, Style};
    use reedline::{
        default_emacs_keybindings, DefaultHinter, DefaultPrompt, EditMode, Emacs, IdeMenu,
        KeyModifiers, MenuBuilder, Prompt, Reedline, ReedlineEvent, ReedlineMenu,
    };

    use crate::repl::{ClapEditor, ReedCompleter, SharedCompletionState, ShellLikeValidator};

    pub struct ClapEditorBuilder<C: Parser + Send + Sync + 'static> {
        prompt: Box<dyn Prompt>,
        edit_mode: Box<dyn EditMode>,
        hook: Box<dyn FnOnce(Reedline) -> Reedline>,
        completion_state: SharedCompletionState,
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
                completion_state: SharedCompletionState::new(),
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

        pub fn with_completion_state(mut self, completion_state: SharedCompletionState) -> Self {
            self.completion_state = completion_state;
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
                    completion_state: self.completion_state,
                    c_phantom: PhantomData,
                }))
                .with_menu(ReedlineMenu::EngineCompleter(completion_menu))
                .with_hinter(Box::new(
                    DefaultHinter::default().with_style(Style::new().italic().fg(Color::DarkGray)),
                ))
                .with_validator(Box::new(ShellLikeValidator))
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QuoteStyle {
    None,
    Single,
    Double,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CompletionToken {
    raw: String,
    cooked: String,
    start: usize,
    end: usize,
    quote_style: QuoteStyle,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PathCompletionRequest {
    partial_path: String,
    span_start: usize,
    quote_style: QuoteStyle,
}

#[derive(Debug, Clone)]
struct FlagValueContext<'a> {
    request: PathCompletionRequest,
    arg: &'a clap::Arg,
    consumed_values: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SettingsCatalogKind {
    Global,
    Runtime,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProcessSettingsMutationKind {
    Add,
    Update,
    Remove,
    Display,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ProcessSettingsCommandKind {
    mutation: ProcessSettingsMutationKind,
    setting: NamedProcessSettingKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SettingsValueKind {
    Boolean,
    Integer,
    Number,
    String,
    Array,
    Object,
    Null,
}

impl SettingsValueKind {
    fn description(self) -> &'static str {
        match self {
            SettingsValueKind::Boolean => "expects a boolean",
            SettingsValueKind::Integer => "expects an integer",
            SettingsValueKind::Number => "expects a number",
            SettingsValueKind::String => "expects a string",
            SettingsValueKind::Array => "expects an array",
            SettingsValueKind::Object => "expects an object",
            SettingsValueKind::Null => "expects null",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SettingsPathKind {
    Container,
    Leaf(SettingsValueKind),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ValueCompletionRequest {
    partial_value: String,
    span_start: usize,
    quote_style: QuoteStyle,
}

#[derive(Debug, Clone)]
struct CommandContext<'a> {
    root_cmd: &'a clap::Command,
    cmd: &'a clap::Command,
    matched_path: Vec<&'a str>,
    all_completed_tokens: &'a [CompletionToken],
    completed_tokens: &'a [CompletionToken],
    current_token: &'a CompletionToken,
}

#[derive(Debug, Clone)]
struct ReedCompleter<C: Parser + Send + Sync + 'static> {
    completion_state: SharedCompletionState,
    c_phantom: PhantomData<C>,
}

struct ShellLikeValidator;

impl reedline::Validator for ShellLikeValidator {
    fn validate(&self, line: &str) -> ValidationResult {
        if split_command_line(line).is_ok() {
            ValidationResult::Complete
        } else {
            ValidationResult::Incomplete
        }
    }
}

impl<C: Parser + Send + Sync + 'static> reedline::Completer for ReedCompleter<C> {
    fn complete(&mut self, line: &str, pos: usize) -> Vec<reedline::Suggestion> {
        let completion_state = self.completion_state.snapshot();
        collect_completions::<C>(line, pos, &completion_state)
    }
}

/// Determine if an argument expects a path based on its type and hints
fn is_path_argument(arg: &clap::Arg) -> bool {
    matches!(
        arg.get_value_hint(),
        clap::ValueHint::FilePath
            | clap::ValueHint::DirPath
            | clap::ValueHint::AnyPath
            | clap::ValueHint::ExecutablePath
    )
}

fn global_settings_root() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_settings_with_defaults(&CLISettings::default(), "CLI settings completion")
            .expect("default CLI settings must serialize for completion")
    })
}

fn global_settings_schema() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_schema::<CLISettings>("CLI settings completion")
            .expect("CLI settings schema must serialize for completion")
    })
}

fn runtime_settings_root() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_settings_with_defaults(&RuntimeSettings::default(), "runtime settings completion")
            .expect("default runtime settings must serialize for completion")
    })
}

fn runtime_settings_schema() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_schema::<RuntimeSettings>("runtime settings completion")
            .expect("runtime settings schema must serialize for completion")
    })
}

fn matches_command_path(context: &CommandContext<'_>, path: &[&str]) -> bool {
    context.matched_path.as_slice() == path
}

fn looks_like_filesystem_path(fragment: &str) -> bool {
    fragment.starts_with("./")
        || fragment.starts_with("../")
        || fragment.starts_with('/')
        || fragment.starts_with('~')
        || fragment.contains('/')
        || fragment.contains('\\')
}

fn looks_like_option_prefix(fragment: &str) -> bool {
    if matches!(fragment, "-" | "--") {
        return true;
    }

    if let Some(long) = fragment.strip_prefix("--") {
        return long.chars().next().is_some_and(|c| c.is_ascii_alphabetic());
    }

    if let Some(short) = fragment.strip_prefix('-') {
        let mut chars = short.chars();
        return matches!(
            (chars.next(), chars.next()),
            (Some(first), None) if first.is_ascii_alphabetic()
        );
    }

    false
}

fn collect_completions<C: Parser + Send + Sync + 'static>(
    line: &str,
    pos: usize,
    completion_state: &CompletionState,
) -> Vec<Suggestion> {
    let root_cmd = C::command();
    let tokens = tokenize_completion_input(line, pos);
    let Some(context) = resolve_command_context(&root_cmd, &tokens) else {
        return Vec::new();
    };

    let mut suggestions = Vec::new();
    let mut seen = HashSet::new();
    let used_arg_ids = collect_used_arg_ids(context.cmd, context.completed_tokens);

    if let Some(flag_value_context) = flag_value_completion_request(
        context.cmd,
        context.completed_tokens,
        context.current_token,
        pos,
    ) {
        if !flag_value_context.consumed_values.is_empty()
            && looks_like_option_prefix(context.current_token.cooked.as_str())
        {
            add_flag_suggestions(
                context.cmd,
                context.current_token,
                pos,
                &mut suggestions,
                &mut seen,
                &used_arg_ids,
            );
            if !suggestions.is_empty() {
                return suggestions;
            }
        }

        let value_request = &flag_value_context.request;
        let arg = flag_value_context.arg;
        if is_path_argument(arg) {
            add_path_suggestions(&mut suggestions, &mut seen, value_request, pos);
        } else if let Some(completion) = arg_value_completion(arg) {
            match completion {
                ArgValueCompletion::ProcessSelector(kind) => {
                    add_process_suggestions(
                        completion_state,
                        value_request,
                        kind,
                        pos,
                        &mut suggestions,
                        &mut seen,
                    );
                }
                ArgValueCompletion::IntegrandSelector(kind) => {
                    add_integrand_suggestions(
                        completion_state,
                        context.all_completed_tokens,
                        context.root_cmd,
                        value_request,
                        kind,
                        pos,
                        &mut suggestions,
                        &mut seen,
                    );
                }
                ArgValueCompletion::SelectedIntegrandTarget => {
                    add_selected_integrand_target_suggestions(
                        completion_state,
                        context.all_completed_tokens,
                        context.root_cmd,
                        value_request,
                        pos,
                        &mut suggestions,
                        &mut seen,
                    );
                }
                ArgValueCompletion::Disabled => {}
            }
        } else if is_ir_profile_select_argument(&context, arg) {
            add_ir_profile_select_suggestions(
                completion_state,
                context.all_completed_tokens,
                context.root_cmd,
                value_request,
                pos,
                &mut suggestions,
                &mut seen,
            );
        } else if is_generate_vertex_interactions_argument(arg) {
            add_generate_vertex_suggestions(
                completion_state,
                value_request,
                &flag_value_context.consumed_values,
                pos,
                &mut suggestions,
                &mut seen,
            );
        } else if add_possible_value_suggestions(
            arg,
            value_request,
            pos,
            &mut suggestions,
            &mut seen,
        ) {
        }
        return suggestions;
    }

    if context.current_token.cooked.starts_with('-') {
        add_flag_suggestions(
            context.cmd,
            context.current_token,
            pos,
            &mut suggestions,
            &mut seen,
            &used_arg_ids,
        );
        return suggestions;
    }

    if context.current_token.cooked.is_empty() && process_settings_command_kind(&context).is_none()
    {
        add_flag_suggestions(
            context.cmd,
            context.current_token,
            pos,
            &mut suggestions,
            &mut seen,
            &used_arg_ids,
        );
    }

    if let Some(path_request) = positional_path_completion_request(
        context.cmd,
        context.completed_tokens,
        context.current_token,
    ) {
        add_value_suggestions(&context, pos, completion_state, &mut suggestions, &mut seen);
        add_path_suggestions(&mut suggestions, &mut seen, &path_request, pos);
    } else {
        add_value_suggestions(&context, pos, completion_state, &mut suggestions, &mut seen);
        add_subcommand_suggestions(
            context.cmd,
            context.current_token,
            pos,
            &mut suggestions,
            &mut seen,
        );
    }

    suggestions
}

fn resolve_command_context<'a>(
    root_cmd: &'a clap::Command,
    tokens: &'a [CompletionToken],
) -> Option<CommandContext<'a>> {
    let current_token = tokens.last()?;
    let completed_tokens = &tokens[..tokens.len().saturating_sub(1)];
    let mut cmd = root_cmd;
    let mut index = 0;
    let mut remaining_start = 0;
    let mut matched_path = Vec::new();

    while index < completed_tokens.len() {
        let token = completed_tokens[index].cooked.as_str();
        if let Some(subcmd) = find_subcommand(cmd, token) {
            matched_path.push(subcmd.get_name());
            cmd = subcmd;
            index += 1;
            remaining_start = index;
            continue;
        }

        if let Some(arg) = find_flag(cmd, token) {
            index += 1;
            if arg_takes_value(arg) && inline_flag_value_start(token).is_none() {
                index = index.saturating_add(1);
            }
            continue;
        }
        break;
    }

    Some(CommandContext {
        root_cmd,
        cmd,
        matched_path,
        all_completed_tokens: completed_tokens,
        completed_tokens: &completed_tokens[remaining_start..],
        current_token,
    })
}

fn find_subcommand<'a>(cmd: &'a clap::Command, token: &str) -> Option<&'a clap::Command> {
    cmd.get_subcommands()
        .find(|subcmd| subcommand_matches(subcmd, token))
}

fn subcommand_matches(subcommand: &clap::Command, token: &str) -> bool {
    subcommand.get_name() == token || subcommand.get_all_aliases().any(|alias| alias == token)
}

fn find_flag<'a>(cmd: &'a clap::Command, token: &str) -> Option<&'a clap::Arg> {
    if let Some(long) = token.strip_prefix("--") {
        let long = long.split_once('=').map(|(name, _)| name).unwrap_or(long);
        return cmd.get_arguments().find(|arg| arg.get_long() == Some(long));
    }

    if token.starts_with('-') && !token.starts_with("--") {
        let mut chars = token.chars();
        let _dash = chars.next();
        let short = chars.next()?;
        if chars.next().is_none() {
            return cmd
                .get_arguments()
                .find(|arg| arg.get_short() == Some(short));
        }
    }

    None
}

fn arg_takes_value(arg: &clap::Arg) -> bool {
    arg.get_action().takes_values()
}

fn arg_id(arg: &clap::Arg) -> String {
    arg.get_id().to_string()
}

fn arg_max_values(arg: &clap::Arg) -> usize {
    arg.get_num_args()
        .map(|range| range.max_values())
        .unwrap_or(1)
}

fn arg_accepts_more_values(arg: &clap::Arg, consumed_values: usize) -> bool {
    consumed_values < arg_max_values(arg)
}

fn inline_flag_value_start(token: &str) -> Option<usize> {
    token
        .strip_prefix("--")
        .and_then(|long| long.find('=').map(|offset| offset + 3))
}

fn flag_value_completion_request<'a>(
    cmd: &'a clap::Command,
    completed_tokens: &[CompletionToken],
    current_token: &CompletionToken,
    _pos: usize,
) -> Option<FlagValueContext<'a>> {
    if let Some(flag_value_context) = inline_flag_value_completion_request(cmd, current_token) {
        return Some(flag_value_context);
    }

    let mut active_arg = None;
    let mut consumed_values = Vec::new();
    let mut consumed_count = 0usize;

    for token in completed_tokens {
        if let Some(arg) = find_flag(cmd, token.cooked.as_str()) {
            if !arg_takes_value(arg) {
                active_arg = None;
                consumed_values.clear();
                consumed_count = 0;
                continue;
            }

            consumed_values.clear();
            consumed_count = 0;
            if inline_flag_value_start(token.cooked.as_str()).is_some() {
                consumed_count = 1;
                if let Some((_, value)) = token.cooked.split_once('=') {
                    consumed_values.push(value.to_string());
                }
            }
            active_arg = arg_accepts_more_values(arg, consumed_count).then_some(arg);
            continue;
        }

        if token.cooked.starts_with('-') {
            active_arg = None;
            consumed_values.clear();
            consumed_count = 0;
            continue;
        }

        let Some(arg) = active_arg else {
            continue;
        };
        consumed_values.push(token.cooked.clone());
        consumed_count += 1;
        if !arg_accepts_more_values(arg, consumed_count) {
            active_arg = None;
            consumed_values.clear();
            consumed_count = 0;
        }
    }

    Some(FlagValueContext {
        request: PathCompletionRequest {
            partial_path: current_token.cooked.clone(),
            span_start: current_token.start,
            quote_style: current_token.quote_style,
        },
        arg: active_arg?,
        consumed_values,
    })
}

fn inline_flag_value_completion_request<'a>(
    cmd: &'a clap::Command,
    current_token: &CompletionToken,
) -> Option<FlagValueContext<'a>> {
    let separator = inline_flag_value_start(current_token.cooked.as_str())?;
    let arg = find_flag(cmd, current_token.cooked.as_str())?;
    if !arg_takes_value(arg) {
        return None;
    }

    let raw_separator = current_token.raw.find('=')? + 1;
    let fragment = parse_shell_fragment(
        &current_token.raw[raw_separator..],
        current_token.start + raw_separator,
    );
    let cooked_value = current_token.cooked[separator..].to_string();

    Some(FlagValueContext {
        request: PathCompletionRequest {
            partial_path: cooked_value,
            span_start: fragment.start,
            quote_style: fragment.quote_style,
        },
        arg,
        consumed_values: Vec::new(),
    })
}

fn collect_used_arg_ids(
    cmd: &clap::Command,
    completed_tokens: &[CompletionToken],
) -> HashSet<String> {
    let mut used_arg_ids = HashSet::new();
    for token in completed_tokens {
        if let Some(arg) = find_flag(cmd, token.cooked.as_str()) {
            if arg_allows_multiple_occurrences(arg) {
                continue;
            }
            used_arg_ids.insert(arg_id(arg));
        }
    }
    used_arg_ids
}

fn arg_allows_multiple_occurrences(arg: &clap::Arg) -> bool {
    matches!(
        arg.get_action(),
        clap::ArgAction::Append | clap::ArgAction::Count
    )
}

fn positional_path_completion_request(
    cmd: &clap::Command,
    completed_tokens: &[CompletionToken],
    current_token: &CompletionToken,
) -> Option<PathCompletionRequest> {
    let mut positional_index = 0usize;
    let mut index = 0usize;

    while index < completed_tokens.len() {
        let token = completed_tokens[index].cooked.as_str();
        if let Some(arg) = find_flag(cmd, token) {
            index += 1;
            if arg_takes_value(arg) && inline_flag_value_start(token).is_none() {
                index = index.saturating_add(1);
            }
            continue;
        }

        if !token.starts_with('-') {
            positional_index += 1;
        }
        index += 1;
    }

    let mut positional_args: Vec<_> = cmd
        .get_arguments()
        .filter(|arg| arg.is_positional())
        .collect();
    positional_args.sort_by_key(|arg| arg.get_index().unwrap_or(usize::MAX));

    let positional_arg = positional_args.get(positional_index).copied().or_else(|| {
        positional_args.last().copied().filter(|arg| {
            arg.get_num_args()
                .map(|range| range.max_values() > 1)
                .unwrap_or(false)
        })
    })?;

    if !is_path_argument(positional_arg) {
        return None;
    }

    Some(PathCompletionRequest {
        partial_path: current_token.cooked.clone(),
        span_start: current_token.start,
        quote_style: current_token.quote_style,
    })
}

fn add_flag_suggestions(
    cmd: &clap::Command,
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
    used_arg_ids: &HashSet<String>,
) {
    let prefix = current_token.cooked.as_str();

    for arg in cmd.get_arguments() {
        if used_arg_ids.contains(&arg_id(arg)) {
            continue;
        }
        if let Some(short) = arg.get_short() {
            let flag = format!("-{}", short);
            if (prefix.is_empty() || flag.starts_with(prefix)) && seen.insert(flag.clone()) {
                suggestions.push(reedline::Suggestion {
                    value: flag,
                    description: arg.get_help().map(|help| help.to_string()),
                    style: None,
                    extra: None,
                    span: Span::new(current_token.start, pos),
                    append_whitespace: true,
                });
            }
        }

        if let Some(long) = arg.get_long() {
            let flag = format!("--{}", long);
            if (prefix.is_empty() || flag.starts_with(prefix)) && seen.insert(flag.clone()) {
                suggestions.push(reedline::Suggestion {
                    value: flag,
                    description: arg.get_help().map(|help| help.to_string()),
                    style: None,
                    extra: None,
                    span: Span::new(current_token.start, pos),
                    append_whitespace: true,
                });
            }
        }
    }
}

fn is_generate_vertex_interactions_argument(arg: &clap::Arg) -> bool {
    matches!(
        arg.get_long(),
        Some("veto-vertex-interactions" | "allowed-vertex-interactions")
    )
}

fn add_possible_value_suggestions(
    arg: &clap::Arg,
    request: &PathCompletionRequest,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) -> bool {
    let prefix = request.partial_path.as_str();
    let mut added = false;

    for possible_value in arg.get_possible_values() {
        if possible_value.is_hide_set() {
            continue;
        }

        let description = possible_value.get_help().map(|help| help.to_string());
        for candidate in possible_value.get_name_and_aliases() {
            if !prefix.is_empty() && !candidate.starts_with(prefix) {
                continue;
            }
            let rendered = render_value_completion(candidate, request.quote_style);
            if !seen.insert(rendered.clone()) {
                continue;
            }
            suggestions.push(reedline::Suggestion {
                value: rendered.clone(),
                description: description.clone(),
                style: None,
                extra: None,
                span: Span::new(request.span_start, pos),
                append_whitespace: value_completion_appends_whitespace(&rendered),
            });
            added = true;
        }
    }

    added
}

fn add_generate_vertex_suggestions(
    completion_state: &CompletionState,
    request: &PathCompletionRequest,
    consumed_values: &[String],
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = request.partial_path.as_str();
    let consumed_values = consumed_values
        .iter()
        .map(String::as_str)
        .collect::<HashSet<_>>();

    for vertex_name in &completion_state.model_vertex_names {
        if consumed_values.contains(vertex_name.as_str()) {
            continue;
        }
        if !prefix.is_empty() && !vertex_name.starts_with(prefix) {
            continue;
        }
        let rendered = render_value_completion(vertex_name, request.quote_style);
        if !seen.insert(rendered.clone()) {
            continue;
        }
        suggestions.push(reedline::Suggestion {
            value: rendered.clone(),
            description: Some("Model vertex interaction".to_string()),
            style: None,
            extra: None,
            span: Span::new(request.span_start, pos),
            append_whitespace: value_completion_appends_whitespace(&rendered),
        });
    }
}

fn add_subcommand_suggestions(
    cmd: &clap::Command,
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = current_token.cooked.as_str();
    for subcommand in cmd.get_subcommands() {
        add_subcommand_name_suggestion(
            subcommand.get_name(),
            subcommand,
            prefix,
            current_token.start,
            pos,
            suggestions,
            seen,
        );

        for alias in subcommand.get_all_aliases() {
            add_subcommand_name_suggestion(
                alias,
                subcommand,
                prefix,
                current_token.start,
                pos,
                suggestions,
                seen,
            );
        }
    }
}

fn add_subcommand_name_suggestion(
    name: &str,
    subcommand: &clap::Command,
    prefix: &str,
    span_start: usize,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    if !prefix.is_empty() && !name.starts_with(prefix) {
        return;
    }

    if !seen.insert(name.to_string()) {
        return;
    }

    suggestions.push(reedline::Suggestion {
        value: name.to_string(),
        description: subcommand.get_about().map(|about| about.to_string()),
        style: None,
        extra: None,
        span: Span::new(span_start, pos),
        append_whitespace: true,
    });
}

fn add_value_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    if context.cmd.get_name() == "run"
        || matches_command_path(context, &["display", "command_block"])
    {
        add_run_block_suggestions(context, pos, completion_state, suggestions, seen);
    }

    if matches_command_path(context, &["import", "model"]) {
        add_builtin_model_suggestions(context, pos, suggestions, seen);
    }

    if matches_command_path(context, &["set", "model"]) {
        add_model_parameter_suggestions(context, pos, completion_state, suggestions, seen);
    }

    if matches_generate_command_path(context) {
        add_generate_spec_suggestions(context, pos, completion_state, suggestions, seen);
    }

    if add_process_settings_suggestions(context, pos, completion_state, suggestions, seen) {
        return;
    }

    if let Some(catalog_kind) = settings_catalog_kind(context) {
        add_settings_kv_suggestions(context, pos, catalog_kind, suggestions, seen);
    }
}

fn add_run_block_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = context.current_token.cooked.as_str();
    for name in &completion_state.commands_block_names {
        let rendered = render_value_completion(name, context.current_token.quote_style);
        if (prefix.is_empty() || name.starts_with(prefix)) && seen.insert(rendered.clone()) {
            suggestions.push(Suggestion {
                value: rendered,
                description: Some("Active command block".to_string()),
                style: None,
                extra: None,
                span: Span::new(context.current_token.start, pos),
                append_whitespace: true,
            });
        }
    }
}

fn add_builtin_model_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = context.current_token.cooked.as_str();
    if looks_like_filesystem_path(prefix) {
        return;
    }

    if let Some((model_name, restriction_prefix)) = builtin_model_restriction_request(prefix) {
        add_builtin_model_restriction_suggestions(
            model_name,
            restriction_prefix,
            context.current_token.quote_style,
            context.current_token.start,
            pos,
            suggestions,
            seen,
        );
        return;
    }

    for name in builtin_json_model_names() {
        let rendered = render_value_completion(name, context.current_token.quote_style);
        if (prefix.is_empty() || name.starts_with(prefix)) && seen.insert(rendered.clone()) {
            suggestions.push(Suggestion {
                value: rendered,
                description: Some("Built-in JSON model".to_string()),
                style: None,
                extra: None,
                span: Span::new(context.current_token.start, pos),
                append_whitespace: true,
            });
        }
    }
}

fn builtin_model_restriction_request(prefix: &str) -> Option<(&str, &str)> {
    let (model_name, restriction_prefix) = prefix.split_once('-')?;
    (!model_name.is_empty()).then_some((model_name, restriction_prefix))
}

fn add_builtin_model_restriction_suggestions(
    model_name: &str,
    restriction_prefix: &str,
    quote_style: QuoteStyle,
    span_start: usize,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    if !builtin_json_model_names()
        .iter()
        .any(|candidate| candidate == model_name)
    {
        return;
    }

    let mut candidates = builtin_json_model_restriction_names(model_name)
        .unwrap_or(&[])
        .to_vec();
    candidates.push("full".to_string());
    candidates.sort();
    candidates.dedup();

    for restriction in candidates {
        if !restriction_prefix.is_empty() && !restriction.starts_with(restriction_prefix) {
            continue;
        }
        let value = format!("{model_name}-{restriction}");
        let rendered = render_value_completion(&value, quote_style);
        if !seen.insert(rendered.clone()) {
            continue;
        }
        suggestions.push(Suggestion {
            value: rendered,
            description: Some("Built-in JSON model restriction".to_string()),
            style: None,
            extra: None,
            span: Span::new(span_start, pos),
            append_whitespace: true,
        });
    }
}

fn add_model_parameter_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    if context.completed_tokens.last().is_some_and(|token| {
        matches!(
            token.cooked.as_str(),
            "-p" | "--process" | "-i" | "--integrand-name"
        )
    }) {
        return;
    }

    let current = context.current_token.cooked.as_str();
    if current.contains('=') {
        add_model_parameter_value_hint(context, completion_state, pos, suggestions, seen);
        return;
    }

    if current == "defaults" {
        return;
    }

    let mut assigned = context
        .completed_tokens
        .iter()
        .filter_map(|token| {
            token
                .cooked
                .split_once('=')
                .map(|(key, _)| key.trim().to_string())
        })
        .collect::<HashSet<_>>();
    if let Some((key, _)) = current.split_once('=') {
        assigned.insert(key.trim().to_string());
    }

    if !current.contains('=') && ("defaults".starts_with(current) || current.is_empty()) {
        let rendered = render_value_completion("defaults", context.current_token.quote_style);
        if seen.insert(rendered.clone()) {
            suggestions.push(Suggestion {
                value: rendered,
                description: Some("Reset targeted model settings to defaults".to_string()),
                style: None,
                extra: None,
                span: Span::new(context.current_token.start, pos),
                append_whitespace: true,
            });
        }
    }

    for entry in &completion_state.model_parameter_entries {
        let name = &entry.name;
        if assigned.contains(name) || !(current.is_empty() || name.starts_with(current)) {
            continue;
        }
        let rendered =
            render_value_completion(&format!("{name}="), context.current_token.quote_style);
        if !seen.insert(rendered.clone()) {
            continue;
        }
        suggestions.push(Suggestion {
            value: rendered,
            description: Some("External model parameter".to_string()),
            style: None,
            extra: None,
            span: Span::new(context.current_token.start, pos),
            append_whitespace: false,
        });
    }
}

fn add_model_parameter_value_hint(
    context: &CommandContext<'_>,
    completion_state: &CompletionState,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let parameter_type = context
        .current_token
        .cooked
        .split_once('=')
        .and_then(|(name, _)| {
            completion_state
                .model_parameter_entries
                .iter()
                .find(|entry| entry.name == name.trim())
                .map(|entry| entry.parameter_type.clone())
        });
    let hint = crate::model_parameters::model_value_format_hint(parameter_type);
    let rendered = render_value_completion(
        context.current_token.cooked.as_str(),
        context.current_token.quote_style,
    );
    let marker = format!("{}::{hint}", rendered);
    if !seen.insert(marker) {
        return;
    }

    suggestions.push(Suggestion {
        value: rendered,
        description: Some(hint.to_string()),
        style: None,
        extra: None,
        span: Span::new(context.current_token.start, pos),
        append_whitespace: false,
    });
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenerateMode {
    Amp,
    Xs,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenerateSpecRegion {
    Initial,
    Final,
    Veto,
    Only,
    ProcessOptions,
    Perturbative,
}

#[derive(Debug, Clone)]
struct GenerateCompletionRequest {
    mode: GenerateMode,
    region: GenerateSpecRegion,
    spec_tokens: Vec<String>,
}

fn matches_generate_command_path(context: &CommandContext<'_>) -> bool {
    matches_command_path(context, &["generate", "amp"])
        || matches_command_path(context, &["generate", "xs"])
}

fn generate_mode_for_context(context: &CommandContext<'_>) -> Option<GenerateMode> {
    if matches_command_path(context, &["generate", "amp"]) {
        Some(GenerateMode::Amp)
    } else if matches_command_path(context, &["generate", "xs"]) {
        Some(GenerateMode::Xs)
    } else {
        None
    }
}

fn add_generate_spec_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let Some(request) = generate_completion_request(context) else {
        return;
    };

    match request.region {
        GenerateSpecRegion::Initial => {
            add_generate_particle_suggestions(
                &completion_state.model_particle_names,
                context.current_token,
                pos,
                "Model particle",
                suggestions,
                seen,
            );
            add_generate_literal_suggestion(
                ">",
                "Initial/final state separator",
                context.current_token,
                pos,
                suggestions,
                seen,
            );
            add_generate_literal_suggestion(
                "to",
                "Initial/final state separator",
                context.current_token,
                pos,
                suggestions,
                seen,
            );
        }
        GenerateSpecRegion::Final => {
            add_generate_particle_suggestions(
                &completion_state.model_particle_names,
                context.current_token,
                pos,
                "Model particle",
                suggestions,
                seen,
            );
            add_generate_process_option_starters(context.current_token, pos, suggestions, seen);
            add_generate_coupling_suggestions(
                request.mode,
                &completion_state.model_coupling_names,
                context.current_token,
                pos,
                suggestions,
                seen,
            );
        }
        GenerateSpecRegion::Veto | GenerateSpecRegion::Only => {
            add_generate_particle_suggestions(
                &completion_state.model_particle_names,
                context.current_token,
                pos,
                "Model particle",
                suggestions,
                seen,
            );
            add_generate_coupling_suggestions(
                request.mode,
                &completion_state.model_coupling_names,
                context.current_token,
                pos,
                suggestions,
                seen,
            );
        }
        GenerateSpecRegion::ProcessOptions => {
            add_generate_process_option_starters(context.current_token, pos, suggestions, seen);
            add_generate_coupling_suggestions(
                request.mode,
                &completion_state.model_coupling_names,
                context.current_token,
                pos,
                suggestions,
                seen,
            );
            add_generate_perturbative_block_starters(context.current_token, pos, suggestions, seen);
        }
        GenerateSpecRegion::Perturbative => {
            add_generate_perturbative_suggestions(
                &request.spec_tokens,
                &completion_state.model_coupling_names,
                context.current_token,
                pos,
                suggestions,
                seen,
            );
        }
    }
}

fn generate_completion_request(context: &CommandContext<'_>) -> Option<GenerateCompletionRequest> {
    let mode = generate_mode_for_context(context)?;
    if context.current_token.cooked.starts_with('-') {
        return None;
    }

    let mut spec_tokens = Vec::new();
    for token in context.completed_tokens {
        if token.cooked.starts_with('-') {
            return None;
        }
        spec_tokens.push(token.cooked.clone());
    }

    Some(GenerateCompletionRequest {
        mode,
        region: generate_spec_region(&spec_tokens, context.current_token.cooked.as_str()),
        spec_tokens,
    })
}

fn generate_spec_region(spec_tokens: &[String], current_prefix: &str) -> GenerateSpecRegion {
    if has_unclosed_perturbative_block(spec_tokens, current_prefix) {
        return GenerateSpecRegion::Perturbative;
    }

    let Some(arrow_index) = find_generate_arrow(spec_tokens) else {
        return GenerateSpecRegion::Initial;
    };

    let mut region = GenerateSpecRegion::Final;
    for token in spec_tokens.iter().skip(arrow_index + 1) {
        if matches!(token.as_str(), "/" | "veto") {
            region = GenerateSpecRegion::Veto;
            continue;
        }
        if matches!(token.as_str(), "|" | "only") {
            region = GenerateSpecRegion::Only;
            continue;
        }
        if token.starts_with('[') || is_generate_coupling_token(token) {
            region = GenerateSpecRegion::ProcessOptions;
        }
    }

    region
}

fn has_unclosed_perturbative_block(spec_tokens: &[String], current_prefix: &str) -> bool {
    let mut depth = 0i32;
    for token in spec_tokens
        .iter()
        .map(String::as_str)
        .chain([current_prefix])
    {
        for ch in token.chars() {
            match ch {
                '[' => depth += 1,
                ']' => depth -= 1,
                _ => {}
            }
        }
    }
    depth > 0 || current_prefix.starts_with('[')
}

fn find_generate_arrow(spec_tokens: &[String]) -> Option<usize> {
    spec_tokens
        .iter()
        .position(|token| token == ">" || token.eq_ignore_ascii_case("to"))
}

fn is_generate_coupling_token(token: &str) -> bool {
    generate_coupling_name_prefix(token).is_some()
        && (token.contains("==") || token.contains(">=") || token.contains("<="))
}

fn generate_coupling_name_prefix(token: &str) -> Option<&str> {
    let mut chars = token.char_indices();
    let (_, first) = chars.next()?;
    if !(first.is_ascii_alphabetic() || first == '_') {
        return None;
    }

    let mut end = first.len_utf8();
    for (index, ch) in chars {
        if ch.is_ascii_alphanumeric() || ch == '_' {
            end = index + ch.len_utf8();
            continue;
        }
        break;
    }
    Some(&token[..end])
}

fn add_generate_particle_suggestions(
    particle_names: &[String],
    current_token: &CompletionToken,
    pos: usize,
    description: &str,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = current_token.cooked.as_str();
    for particle_name in particle_names {
        if !prefix.is_empty() && !particle_name.starts_with(prefix) {
            continue;
        }
        add_generate_literal_suggestion(
            particle_name,
            description,
            current_token,
            pos,
            suggestions,
            seen,
        );
    }
}

fn add_generate_process_option_starters(
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    for (value, description) in [
        ("/", "Start vetoed-particle list"),
        ("veto", "Start vetoed-particle list"),
        ("|", "Start allowed-particle list"),
        ("only", "Start allowed-particle list"),
    ] {
        add_generate_literal_suggestion(value, description, current_token, pos, suggestions, seen);
    }
}

fn add_generate_coupling_suggestions(
    mode: GenerateMode,
    coupling_names: &[String],
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = current_token.cooked.as_str();
    for coupling_name in coupling_names {
        for candidate in [
            format!("{coupling_name}=="),
            format!("{coupling_name}>="),
            format!("{coupling_name}<="),
        ] {
            if !prefix.is_empty() && !candidate.starts_with(prefix) {
                continue;
            }
            add_generate_literal_suggestion(
                &candidate,
                "Coupling-order constraint",
                current_token,
                pos,
                suggestions,
                seen,
            );
        }

        if mode == GenerateMode::Xs {
            let candidate = format!("{coupling_name}^");
            if prefix.is_empty() || candidate.starts_with(prefix) {
                add_generate_literal_suggestion(
                    &candidate,
                    "Cross-section coupling-order power",
                    current_token,
                    pos,
                    suggestions,
                    seen,
                );
            }
        }
    }
}

fn add_generate_perturbative_block_starters(
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    for (value, description) in [
        ("[", "Start perturbative block"),
        ("[{1}]", "One-loop amplitude/summed graphs"),
        ("[{{1}}]", "One-loop forward graphs"),
    ] {
        add_generate_literal_suggestion(value, description, current_token, pos, suggestions, seen);
    }
}

fn add_generate_perturbative_suggestions(
    spec_tokens: &[String],
    coupling_names: &[String],
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let needs_opening_bracket = !spec_tokens.iter().any(|token| token.contains('['))
        && current_token.cooked.starts_with('[')
        && !current_token.cooked.contains(']');
    for (value, description) in [
        ("{1}", "One-loop amplitude/summed graphs"),
        ("{{1}}", "One-loop forward graphs"),
    ] {
        let candidate = if needs_opening_bracket {
            format!("[{value}")
        } else {
            value.to_string()
        };
        add_generate_literal_suggestion(
            &candidate,
            description,
            current_token,
            pos,
            suggestions,
            seen,
        );
    }

    for coupling_name in coupling_names {
        for candidate in [coupling_name.clone(), format!("{coupling_name}=")] {
            let candidate = if needs_opening_bracket {
                format!("[{candidate}")
            } else {
                candidate
            };
            add_generate_literal_suggestion(
                &candidate,
                "Perturbative coupling order",
                current_token,
                pos,
                suggestions,
                seen,
            );
        }
    }
}

fn add_generate_literal_suggestion(
    value: &str,
    description: &str,
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = current_token.cooked.as_str();
    if !prefix.is_empty() && !value.starts_with(prefix) {
        return;
    }

    let rendered = render_value_completion(value, current_token.quote_style);
    if !seen.insert(rendered.clone()) {
        return;
    }

    suggestions.push(Suggestion {
        value: rendered.clone(),
        description: Some(description.to_string()),
        style: None,
        extra: None,
        span: Span::new(current_token.start, pos),
        append_whitespace: value_completion_appends_whitespace(&rendered),
    });
}

fn settings_catalog_kind(context: &CommandContext<'_>) -> Option<SettingsCatalogKind> {
    if matches_command_path(context, &["set", "global", "kv"]) {
        Some(SettingsCatalogKind::Global)
    } else if matches_command_path(context, &["set", "default-runtime", "kv"])
        || matches_command_path(context, &["set", "process", "kv"])
    {
        Some(SettingsCatalogKind::Runtime)
    } else {
        None
    }
}

fn process_settings_command_kind(
    context: &CommandContext<'_>,
) -> Option<ProcessSettingsCommandKind> {
    let cases = [
        (
            &["set", "process", "add", "quantity"][..],
            ProcessSettingsMutationKind::Add,
            NamedProcessSettingKind::Quantity,
        ),
        (
            &["set", "process", "add", "observable"][..],
            ProcessSettingsMutationKind::Add,
            NamedProcessSettingKind::Observable,
        ),
        (
            &["set", "process", "add", "selector"][..],
            ProcessSettingsMutationKind::Add,
            NamedProcessSettingKind::Selector,
        ),
        (
            &["set", "process", "update", "quantity"][..],
            ProcessSettingsMutationKind::Update,
            NamedProcessSettingKind::Quantity,
        ),
        (
            &["set", "process", "update", "observable"][..],
            ProcessSettingsMutationKind::Update,
            NamedProcessSettingKind::Observable,
        ),
        (
            &["set", "process", "update", "selector"][..],
            ProcessSettingsMutationKind::Update,
            NamedProcessSettingKind::Selector,
        ),
        (
            &["set", "process", "remove", "quantity"][..],
            ProcessSettingsMutationKind::Remove,
            NamedProcessSettingKind::Quantity,
        ),
        (
            &["set", "process", "remove", "observable"][..],
            ProcessSettingsMutationKind::Remove,
            NamedProcessSettingKind::Observable,
        ),
        (
            &["set", "process", "remove", "selector"][..],
            ProcessSettingsMutationKind::Remove,
            NamedProcessSettingKind::Selector,
        ),
        (
            &["display", "quantities"][..],
            ProcessSettingsMutationKind::Display,
            NamedProcessSettingKind::Quantity,
        ),
        (
            &["display", "observables"][..],
            ProcessSettingsMutationKind::Display,
            NamedProcessSettingKind::Observable,
        ),
        (
            &["display", "selectors"][..],
            ProcessSettingsMutationKind::Display,
            NamedProcessSettingKind::Selector,
        ),
    ];

    cases.iter().find_map(|(path, mutation, setting)| {
        matches_command_path(context, path).then_some(ProcessSettingsCommandKind {
            mutation: *mutation,
            setting: *setting,
        })
    })
}

fn add_settings_kv_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    catalog_kind: SettingsCatalogKind,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let (key_request, value_request) = split_kv_completion_request(context.current_token);
    let root = match catalog_kind {
        SettingsCatalogKind::Global => global_settings_root(),
        SettingsCatalogKind::Runtime => runtime_settings_root(),
    };
    let schema = match catalog_kind {
        SettingsCatalogKind::Global => global_settings_schema(),
        SettingsCatalogKind::Runtime => runtime_settings_schema(),
    };

    if let Some(value_request) = value_request {
        add_settings_value_suggestions(
            root,
            schema,
            &key_request,
            &value_request,
            pos,
            suggestions,
            seen,
        );
        return;
    }

    add_settings_path_suggestions(root, &key_request, pos, suggestions, seen);
}

fn add_process_settings_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) -> bool {
    let Some(command_kind) = process_settings_command_kind(context) else {
        return false;
    };

    match command_kind.mutation {
        ProcessSettingsMutationKind::Add => add_process_settings_add_suggestions(
            context,
            pos,
            completion_state,
            command_kind.setting,
            suggestions,
            seen,
        ),
        ProcessSettingsMutationKind::Update => add_process_settings_update_suggestions(
            context,
            pos,
            completion_state,
            command_kind.setting,
            suggestions,
            seen,
        ),
        ProcessSettingsMutationKind::Remove | ProcessSettingsMutationKind::Display => {
            add_process_settings_name_suggestions_for_existing_entries(
                context,
                pos,
                completion_state,
                command_kind.setting,
                suggestions,
                seen,
            )
        }
    }

    true
}

fn add_process_settings_add_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    setting_kind: NamedProcessSettingKind,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    if context.completed_tokens.is_empty() {
        let request = ValueCompletionRequest {
            partial_value: context.current_token.cooked.clone(),
            span_start: context.current_token.start,
            quote_style: context.current_token.quote_style,
        };
        add_type_hint_suggestion(SettingsValueKind::String, &request, pos, suggestions, seen);
        return;
    }

    match setting_kind {
        NamedProcessSettingKind::Quantity | NamedProcessSettingKind::Selector
            if context.completed_tokens.len() == 1 =>
        {
            let kinds = match setting_kind {
                NamedProcessSettingKind::Quantity => quantity_kind_names(),
                NamedProcessSettingKind::Selector => selector_kind_names(),
                NamedProcessSettingKind::Observable => Vec::new(),
            };
            add_literal_value_suggestions(
                kinds.iter().map(String::as_str),
                &format!("Available {} kind", setting_kind.singular()),
                context.current_token,
                pos,
                suggestions,
                seen,
            );
            return;
        }
        _ => {}
    }

    let scopes = selected_process_settings_entries(
        completion_state,
        context.root_cmd,
        context.all_completed_tokens,
    );
    let quantity_candidates = common_quantity_names(&scopes);
    let (root, schema, accepts_quantity_reference) = match setting_kind {
        NamedProcessSettingKind::Quantity if context.completed_tokens.len() >= 2 => (
            quantity_completion_root_for_kind(context.completed_tokens[1].cooked.as_str()),
            Some(quantity_schema()),
            false,
        ),
        NamedProcessSettingKind::Observable if !context.completed_tokens.is_empty() => (
            Some(observable_completion_root()),
            Some(observable_schema()),
            true,
        ),
        NamedProcessSettingKind::Selector if context.completed_tokens.len() >= 2 => (
            selector_completion_root_for_kind(context.completed_tokens[1].cooked.as_str()),
            Some(selector_schema()),
            true,
        ),
        _ => (None, None, false),
    };

    if let (Some(root), Some(schema)) = (root, schema) {
        add_process_settings_kv_suggestions(
            root,
            schema,
            accepts_quantity_reference.then_some(quantity_candidates.as_slice()),
            context.current_token,
            pos,
            suggestions,
            seen,
        );
    }
}

fn add_process_settings_update_suggestions(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    setting_kind: NamedProcessSettingKind,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    if context.completed_tokens.is_empty() {
        add_process_settings_name_suggestions_for_existing_entries(
            context,
            pos,
            completion_state,
            setting_kind,
            suggestions,
            seen,
        );
        return;
    }

    let scopes = selected_process_settings_entries(
        completion_state,
        context.root_cmd,
        context.all_completed_tokens,
    );
    let selected_name = context.completed_tokens[0].cooked.as_str();
    let roots = common_named_setting_roots(&scopes, setting_kind, selected_name);
    let Some(root) = intersect_settings_roots(&roots) else {
        return;
    };

    let schema = match setting_kind {
        NamedProcessSettingKind::Quantity => quantity_schema(),
        NamedProcessSettingKind::Observable => observable_schema(),
        NamedProcessSettingKind::Selector => selector_schema(),
    };
    let quantity_candidates = common_quantity_names(&scopes);
    let accepts_quantity_reference = matches!(
        setting_kind,
        NamedProcessSettingKind::Observable | NamedProcessSettingKind::Selector
    );

    add_process_settings_kv_suggestions(
        &root,
        schema,
        accepts_quantity_reference.then_some(quantity_candidates.as_slice()),
        context.current_token,
        pos,
        suggestions,
        seen,
    );
}

fn add_process_settings_name_suggestions_for_existing_entries(
    context: &CommandContext<'_>,
    pos: usize,
    completion_state: &CompletionState,
    setting_kind: NamedProcessSettingKind,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let scopes = selected_process_settings_entries(
        completion_state,
        context.root_cmd,
        context.all_completed_tokens,
    );
    let names = common_named_setting_names(&scopes, setting_kind);
    add_literal_value_suggestions(
        names.iter().map(String::as_str),
        &format!(
            "{} present in {} selected integrand(s)",
            setting_kind.singular(),
            scopes.len()
        ),
        context.current_token,
        pos,
        suggestions,
        seen,
    );
}

fn add_process_settings_kv_suggestions(
    root: &JsonValue,
    schema: &JsonValue,
    quantity_name_candidates: Option<&[String]>,
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let (key_request, value_request) = split_kv_completion_request(current_token);
    if let Some(value_request) = value_request {
        if key_request == "quantity" {
            if let Some(candidates) = quantity_name_candidates {
                let before = suggestions.len();
                add_path_request_literal_suggestions(
                    candidates.iter().map(String::as_str),
                    "Existing quantity name",
                    &value_request,
                    pos,
                    suggestions,
                    seen,
                );
                if suggestions.len() > before {
                    return;
                }
            }
        }
        add_settings_value_suggestions(
            root,
            schema,
            &key_request,
            &value_request,
            pos,
            suggestions,
            seen,
        );
        return;
    }

    add_settings_path_suggestions(root, &key_request, pos, suggestions, seen);
}

fn add_literal_value_suggestions<'a>(
    candidates: impl IntoIterator<Item = &'a str>,
    description: &str,
    current_token: &CompletionToken,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = current_token.cooked.as_str();
    for candidate in candidates {
        if !prefix.is_empty() && !candidate.starts_with(prefix) {
            continue;
        }
        let rendered = render_value_completion(candidate, current_token.quote_style);
        if !seen.insert(rendered.clone()) {
            continue;
        }
        suggestions.push(Suggestion {
            value: rendered.clone(),
            description: Some(description.to_string()),
            style: None,
            extra: None,
            span: Span::new(current_token.start, pos),
            append_whitespace: value_completion_appends_whitespace(&rendered),
        });
    }
}

fn add_path_request_literal_suggestions<'a>(
    candidates: impl IntoIterator<Item = &'a str>,
    description: &str,
    request: &ValueCompletionRequest,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    for candidate in candidates {
        if !request.partial_value.is_empty() && !candidate.starts_with(&request.partial_value) {
            continue;
        }
        let rendered = render_value_completion(candidate, request.quote_style);
        if !seen.insert(rendered.clone()) {
            continue;
        }
        suggestions.push(Suggestion {
            value: rendered,
            description: Some(description.to_string()),
            style: None,
            extra: None,
            span: Span::new(request.span_start, pos),
            append_whitespace: true,
        });
    }
}

fn selected_process_settings_entries<'a>(
    completion_state: &'a CompletionState,
    root_cmd: &clap::Command,
    completed_tokens: &[CompletionToken],
) -> Vec<&'a ProcessSettingsCompletionEntry> {
    let selected_process =
        find_selected_process_entry(completion_state, root_cmd, completed_tokens).or_else(|| {
            (completion_state.process_entries.len() == 1)
                .then(|| &completion_state.process_entries[0])
        });
    let Some(selected_process) = selected_process else {
        return Vec::new();
    };

    let integrand_filter = find_last_flag_value(root_cmd, completed_tokens, "integrand-name");
    completion_state
        .process_settings_entries
        .iter()
        .filter(|entry| entry.process_id == selected_process.id)
        .filter(|entry| {
            integrand_filter
                .as_deref()
                .is_none_or(|integrand_name| entry.integrand_name == integrand_name)
        })
        .collect()
}

fn named_settings_map(
    entry: &ProcessSettingsCompletionEntry,
    setting_kind: NamedProcessSettingKind,
) -> &BTreeMap<String, JsonValue> {
    match setting_kind {
        NamedProcessSettingKind::Quantity => &entry.quantities,
        NamedProcessSettingKind::Observable => &entry.observables,
        NamedProcessSettingKind::Selector => &entry.selectors,
    }
}

fn common_named_setting_names(
    scopes: &[&ProcessSettingsCompletionEntry],
    setting_kind: NamedProcessSettingKind,
) -> Vec<String> {
    let Some(first) = scopes.first() else {
        return Vec::new();
    };

    named_settings_map(first, setting_kind)
        .keys()
        .filter(|name| {
            scopes
                .iter()
                .all(|scope| named_settings_map(scope, setting_kind).contains_key(*name))
        })
        .cloned()
        .collect()
}

fn common_quantity_names(scopes: &[&ProcessSettingsCompletionEntry]) -> Vec<String> {
    common_named_setting_names(scopes, NamedProcessSettingKind::Quantity)
}

fn common_named_setting_roots<'a>(
    scopes: &[&'a ProcessSettingsCompletionEntry],
    setting_kind: NamedProcessSettingKind,
    name: &str,
) -> Vec<&'a JsonValue> {
    let roots = scopes
        .iter()
        .filter_map(|scope| named_settings_map(scope, setting_kind).get(name))
        .collect::<Vec<_>>();
    if roots.len() == scopes.len() {
        roots
    } else {
        Default::default()
    }
}

fn intersect_settings_roots(roots: &[&JsonValue]) -> Option<JsonValue> {
    let mut iter = roots.iter();
    let mut intersection = (*iter.next()?).clone();

    for root in iter {
        intersection = intersect_json_values(&intersection, root)?;
    }

    Some(intersection)
}

fn intersect_json_values(left: &JsonValue, right: &JsonValue) -> Option<JsonValue> {
    match (left, right) {
        (JsonValue::Object(left_map), JsonValue::Object(right_map)) => {
            let mut merged = serde_json::Map::new();
            for (key, left_value) in left_map {
                let Some(right_value) = right_map.get(key) else {
                    continue;
                };
                if let Some(value) = intersect_json_values(left_value, right_value) {
                    merged.insert(key.clone(), value);
                }
            }
            Some(JsonValue::Object(merged))
        }
        (JsonValue::Array(left_items), JsonValue::Array(right_items))
            if left_items.len() == right_items.len() =>
        {
            Some(JsonValue::Array(left_items.clone()))
        }
        (left_value, right_value) if json_values_have_matching_kind(left_value, right_value) => {
            Some(left_value.clone())
        }
        _ => None,
    }
}

fn json_values_have_matching_kind(left: &JsonValue, right: &JsonValue) -> bool {
    matches!(
        (left, right),
        (JsonValue::Null, JsonValue::Null)
            | (JsonValue::Bool(_), JsonValue::Bool(_))
            | (JsonValue::Number(_), JsonValue::Number(_))
            | (JsonValue::String(_), JsonValue::String(_))
            | (JsonValue::Array(_), JsonValue::Array(_))
            | (JsonValue::Object(_), JsonValue::Object(_))
    )
}

fn add_process_suggestions(
    completion_state: &CompletionState,
    request: &PathCompletionRequest,
    selector_kind: SelectorKind,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = request.partial_path.as_str();
    let process_entries = completion_state
        .process_entries
        .iter()
        .filter(|entry| matches_selector_kind(entry.kind, selector_kind));

    if let Some(name_prefix) = prefix.strip_prefix("name:") {
        for entry in process_entries {
            if !entry.name.starts_with(name_prefix) {
                continue;
            }
            let value =
                render_value_completion(&format!("name:{}", entry.name), request.quote_style);
            if !seen.insert(value.clone()) {
                continue;
            }
            suggestions.push(reedline::Suggestion {
                value,
                description: Some(format!("Process #{} ({})", entry.id, entry.kind.label())),
                style: None,
                extra: None,
                span: Span::new(request.span_start, pos),
                append_whitespace: true,
            });
        }
        return;
    }

    if let Some(id_prefix) = prefix.strip_prefix('#') {
        for entry in process_entries {
            let id = entry.id.to_string();
            if !id.starts_with(id_prefix) {
                continue;
            }
            let value = render_value_completion(&format!("#{}", id), request.quote_style);
            if !seen.insert(value.clone()) {
                continue;
            }
            suggestions.push(reedline::Suggestion {
                value,
                description: Some(format!(
                    "Process {} ({}, {})",
                    entry.name,
                    entry.kind.label(),
                    entry.id
                )),
                style: None,
                extra: None,
                span: Span::new(request.span_start, pos),
                append_whitespace: true,
            });
        }
        return;
    }

    for entry in process_entries {
        if !entry.name.starts_with(prefix) {
            continue;
        }
        let value = render_value_completion(&entry.name, request.quote_style);
        if !seen.insert(value.clone()) {
            continue;
        }
        suggestions.push(reedline::Suggestion {
            value,
            description: Some(format!("Process #{} ({})", entry.id, entry.kind.label())),
            style: None,
            extra: None,
            span: Span::new(request.span_start, pos),
            append_whitespace: true,
        });
    }
}

fn add_integrand_suggestions(
    completion_state: &CompletionState,
    completed_tokens: &[CompletionToken],
    root_cmd: &clap::Command,
    request: &PathCompletionRequest,
    selector_kind: SelectorKind,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = request.partial_path.as_str();
    let process_filter = find_selected_process_entry(completion_state, root_cmd, completed_tokens);
    let mut matching_processes = Vec::new();

    for entry in &completion_state.process_entries {
        if !matches_selector_kind(entry.kind, selector_kind) {
            continue;
        }
        if let Some(selected_process) = process_filter {
            if entry.id != selected_process.id {
                continue;
            }
        }
        matching_processes.push(entry);
    }

    for entry in matching_processes {
        for integrand_name in &entry.integrand_names {
            if !integrand_name.starts_with(prefix) {
                continue;
            }
            let rendered = render_value_completion(integrand_name, request.quote_style);
            if !seen.insert(rendered.clone()) {
                continue;
            }

            let description = Some(if let Some(selected_process) = process_filter {
                format!(
                    "{} integrand in process {}",
                    selected_process.kind.label(),
                    selected_process.name
                )
            } else {
                let processes = completion_state
                    .process_entries
                    .iter()
                    .filter(|candidate| {
                        matches_selector_kind(candidate.kind, selector_kind)
                            && candidate
                                .integrand_names
                                .iter()
                                .any(|name| name == integrand_name)
                    })
                    .map(|candidate| candidate.name.as_str())
                    .collect::<Vec<_>>();
                if processes.len() == 1 {
                    format!("Process {}", processes[0])
                } else {
                    format!("Processes: {}", processes.join(", "))
                }
            });

            suggestions.push(reedline::Suggestion {
                value: rendered,
                description,
                style: None,
                extra: None,
                span: Span::new(request.span_start, pos),
                append_whitespace: true,
            });
        }
    }
}

fn add_selected_integrand_target_suggestions(
    completion_state: &CompletionState,
    completed_tokens: &[CompletionToken],
    root_cmd: &clap::Command,
    request: &PathCompletionRequest,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    let prefix = request.partial_path.as_str();
    if prefix.contains('=') {
        return;
    }

    let selected_keys =
        selected_integrand_keys_for_target_completion(completion_state, completed_tokens, root_cmd);
    let specified_target_keys = specified_target_keys(completed_tokens, root_cmd);

    for slot_key in selected_keys {
        if specified_target_keys.contains(&slot_key) {
            continue;
        }
        let suggestion = format!("{slot_key}=");
        if !suggestion.starts_with(prefix) {
            continue;
        }
        let rendered = render_value_completion(&suggestion, request.quote_style);
        if !seen.insert(rendered.clone()) {
            continue;
        }
        suggestions.push(reedline::Suggestion {
            value: rendered,
            description: Some("Target for selected integrand".to_string()),
            style: None,
            extra: None,
            span: Span::new(request.span_start, pos),
            append_whitespace: false,
        });
    }
}

fn is_ir_profile_select_argument(context: &CommandContext<'_>, arg: &clap::Arg) -> bool {
    matches_command_path(context, &["profile", "infra-red"]) && arg.get_long() == Some("select")
}

fn add_ir_profile_select_suggestions(
    completion_state: &CompletionState,
    completed_tokens: &[CompletionToken],
    root_cmd: &clap::Command,
    request: &PathCompletionRequest,
    pos: usize,
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
) {
    let process_filter = find_last_flag_value(root_cmd, completed_tokens, "process");
    let integrand_filter = find_last_flag_value(root_cmd, completed_tokens, "integrand-name");
    let prefix = request.partial_path.as_str();

    for entry in &completion_state.ir_profile_entries {
        if let Some(process_value) = process_filter.as_deref() {
            let Some(process_entry) =
                resolve_process_reference(&completion_state.process_entries, process_value)
            else {
                continue;
            };
            if entry.process_name != process_entry.name {
                continue;
            }
        }

        if let Some(integrand_name) = integrand_filter.as_deref() {
            if entry.integrand_name != integrand_name {
                continue;
            }
        }

        for graph_name in &entry.graph_names {
            if !prefix.is_empty() && !graph_name.starts_with(prefix) {
                continue;
            }
            let rendered = render_value_completion(graph_name, request.quote_style);
            if !seen.insert(rendered.clone()) {
                continue;
            }
            suggestions.push(reedline::Suggestion {
                value: rendered,
                description: Some(format!(
                    "IR profile graph in {} / {}",
                    entry.process_name, entry.integrand_name
                )),
                style: None,
                extra: None,
                span: Span::new(request.span_start, pos),
                append_whitespace: true,
            });
        }

        for value in &entry.graph_limit_entries {
            if !prefix.is_empty() && !value.starts_with(prefix) {
                continue;
            }
            let rendered = render_value_completion(value, request.quote_style);
            if !seen.insert(rendered.clone()) {
                continue;
            }
            suggestions.push(reedline::Suggestion {
                value: rendered,
                description: Some(format!(
                    "IR profile selection in {} / {}",
                    entry.process_name, entry.integrand_name
                )),
                style: None,
                extra: None,
                span: Span::new(request.span_start, pos),
                append_whitespace: true,
            });
        }
    }
}

fn split_kv_completion_request(
    current_token: &CompletionToken,
) -> (String, Option<ValueCompletionRequest>) {
    let Some(separator) = current_token.raw.find('=') else {
        return (current_token.cooked.clone(), None);
    };

    let key = current_token.raw[..separator].trim().to_string();
    let fragment = parse_shell_fragment(
        &current_token.raw[separator + 1..],
        current_token.start + separator + 1,
    );

    (
        key,
        Some(ValueCompletionRequest {
            partial_value: fragment.cooked,
            span_start: fragment.start,
            quote_style: fragment.quote_style,
        }),
    )
}

fn add_settings_path_suggestions(
    root: &JsonValue,
    key_request: &str,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let (container_path, child_prefix) = settings_container_request(key_request);
    let Ok(container) = value_at_path(root, &container_path) else {
        return;
    };
    let JsonValue::Object(map) = container else {
        return;
    };

    let mut entries = map.iter().collect::<Vec<_>>();
    entries.sort_by(|(left_key, left_value), (right_key, right_value)| {
        let left_is_container = matches!(left_value, JsonValue::Object(_));
        let right_is_container = matches!(right_value, JsonValue::Object(_));
        right_is_container
            .cmp(&left_is_container)
            .then_with(|| left_key.cmp(right_key))
    });

    for (key, value) in entries {
        if !child_prefix.is_empty() && !key.starts_with(&child_prefix) {
            continue;
        }

        let path_kind = settings_path_kind(value);
        let replacement = render_settings_path_completion(&container_path, key, path_kind);
        if !seen.insert(replacement.clone()) {
            continue;
        }

        suggestions.push(Suggestion {
            value: replacement.clone(),
            description: Some(settings_path_description(path_kind).to_string()),
            style: None,
            extra: None,
            span: Span::new(pos.saturating_sub(key_request.len()), pos),
            append_whitespace: value_completion_appends_whitespace(&replacement),
        });
    }
}

fn add_settings_value_suggestions(
    root: &JsonValue,
    schema: &JsonValue,
    key_request: &str,
    value_request: &ValueCompletionRequest,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let Ok(value) = value_at_path(root, key_request) else {
        return;
    };

    if let Some(schema_node) = schema_at_path(schema, key_request) {
        let enum_values = schema_enum_values(schema, schema_node);
        if !enum_values.is_empty() {
            for candidate in enum_values {
                if !candidate.starts_with(&value_request.partial_value) {
                    continue;
                }
                let rendered = render_value_completion(&candidate, value_request.quote_style);
                if !seen.insert(rendered.clone()) {
                    continue;
                }
                suggestions.push(Suggestion {
                    value: rendered,
                    description: Some("expects one of the allowed enum variants".to_string()),
                    style: None,
                    extra: None,
                    span: Span::new(value_request.span_start, pos),
                    append_whitespace: true,
                });
            }
            return;
        }
    }

    let kind = settings_value_kind(value);

    if kind == SettingsValueKind::Boolean {
        for candidate in ["false", "true"] {
            if !candidate.starts_with(&value_request.partial_value) {
                continue;
            }
            let rendered = render_value_completion(candidate, value_request.quote_style);
            if !seen.insert(rendered.clone()) {
                continue;
            }
            suggestions.push(Suggestion {
                value: rendered,
                description: Some(kind.description().to_string()),
                style: None,
                extra: None,
                span: Span::new(value_request.span_start, pos),
                append_whitespace: true,
            });
        }
        return;
    }

    add_type_hint_suggestion(kind, value_request, pos, suggestions, seen);
}

fn add_type_hint_suggestion(
    kind: SettingsValueKind,
    value_request: &ValueCompletionRequest,
    pos: usize,
    suggestions: &mut Vec<Suggestion>,
    seen: &mut HashSet<String>,
) {
    let rendered = render_value_completion(&value_request.partial_value, value_request.quote_style);
    let marker = format!("{}::{}", rendered, kind.description());
    if !seen.insert(marker) {
        return;
    }

    suggestions.push(Suggestion {
        value: rendered,
        description: Some(kind.description().to_string()),
        style: None,
        extra: None,
        span: Span::new(value_request.span_start, pos),
        append_whitespace: false,
    });
}

fn settings_container_request(key_request: &str) -> (String, String) {
    let trimmed = key_request.trim();
    if trimmed.is_empty() {
        return (String::new(), String::new());
    }

    if let Some(container) = trimmed.strip_suffix('.') {
        return (container.to_string(), String::new());
    }

    if let Some((container, child_prefix)) = trimmed.rsplit_once('.') {
        return (container.to_string(), child_prefix.to_string());
    }

    (String::new(), trimmed.to_string())
}

fn render_settings_path_completion(
    container_path: &str,
    key: &str,
    path_kind: SettingsPathKind,
) -> String {
    let mut value = if container_path.is_empty() {
        key.to_string()
    } else {
        format!("{container_path}.{key}")
    };
    match path_kind {
        SettingsPathKind::Container => value.push('.'),
        SettingsPathKind::Leaf(_) => value.push('='),
    }
    value
}

fn settings_path_kind(value: &JsonValue) -> SettingsPathKind {
    match value {
        JsonValue::Object(_) => SettingsPathKind::Container,
        _ => SettingsPathKind::Leaf(settings_value_kind(value)),
    }
}

fn settings_path_description(path_kind: SettingsPathKind) -> &'static str {
    match path_kind {
        SettingsPathKind::Container => "Settings container",
        SettingsPathKind::Leaf(kind) => kind.description(),
    }
}

fn settings_value_kind(value: &JsonValue) -> SettingsValueKind {
    match value {
        JsonValue::Bool(_) => SettingsValueKind::Boolean,
        JsonValue::Number(number) if number.is_i64() || number.is_u64() => {
            SettingsValueKind::Integer
        }
        JsonValue::Number(_) => SettingsValueKind::Number,
        JsonValue::String(_) => SettingsValueKind::String,
        JsonValue::Array(_) => SettingsValueKind::Array,
        JsonValue::Object(_) => SettingsValueKind::Object,
        JsonValue::Null => SettingsValueKind::Null,
    }
}

fn matches_selector_kind(kind: ProcessKind, selector_kind: SelectorKind) -> bool {
    match selector_kind {
        SelectorKind::Any => true,
        SelectorKind::Amplitude => kind == ProcessKind::Amplitude,
        SelectorKind::CrossSection => kind == ProcessKind::CrossSection,
    }
}

fn find_selected_process_entry<'a>(
    completion_state: &'a CompletionState,
    root_cmd: &clap::Command,
    completed_tokens: &[CompletionToken],
) -> Option<&'a ProcessCompletionEntry> {
    let process_value = find_last_flag_value(root_cmd, completed_tokens, "process")?;
    resolve_process_reference(&completion_state.process_entries, &process_value)
}

fn find_last_flag_value(
    root_cmd: &clap::Command,
    completed_tokens: &[CompletionToken],
    long_name: &str,
) -> Option<String> {
    let mut result = None;
    let mut cmd = root_cmd;
    let mut index = 0usize;

    while index < completed_tokens.len() {
        let token = completed_tokens[index].cooked.as_str();
        if let Some(subcmd) = find_subcommand(cmd, token) {
            cmd = subcmd;
            index += 1;
            continue;
        }

        if let Some(arg) = find_flag(cmd, token) {
            let is_target = arg.get_long() == Some(long_name);
            let inline_start = inline_flag_value_start(token);
            if is_target {
                if let Some(start) = inline_start {
                    result = Some(token[start..].to_string());
                } else if arg_takes_value(arg) {
                    result = completed_tokens
                        .get(index + 1)
                        .map(|next| next.cooked.clone());
                }
            }

            index += 1;
            if arg_takes_value(arg) && inline_start.is_none() {
                index = index.saturating_add(1);
            }
            continue;
        }
        index += 1;
    }

    result
}

fn selected_integrand_keys_for_target_completion(
    completion_state: &CompletionState,
    completed_tokens: &[CompletionToken],
    root_cmd: &clap::Command,
) -> Vec<String> {
    let mut keys = Vec::new();
    let mut pending_process: Option<String> = None;
    let mut cmd = root_cmd;
    let mut index = 0usize;

    while index < completed_tokens.len() {
        let token = completed_tokens[index].cooked.as_str();
        if let Some(subcmd) = find_subcommand(cmd, token) {
            cmd = subcmd;
            index += 1;
            continue;
        }

        let Some(arg) = find_flag(cmd, token) else {
            index += 1;
            continue;
        };

        let inline_start = inline_flag_value_start(token);
        let value = if let Some(start) = inline_start {
            Some(token[start..].to_string())
        } else if arg_takes_value(arg) {
            completed_tokens
                .get(index + 1)
                .map(|next| next.cooked.clone())
        } else {
            None
        };

        match arg.get_long() {
            Some("process") => pending_process = value,
            Some("integrand-name") => {
                let Some(integrand_name) = value else {
                    index += 1;
                    if arg_takes_value(arg) && inline_start.is_none() {
                        index = index.saturating_add(1);
                    }
                    continue;
                };

                let Some(process_value) = pending_process.as_deref() else {
                    index += 1;
                    if arg_takes_value(arg) && inline_start.is_none() {
                        index = index.saturating_add(1);
                    }
                    continue;
                };

                let Some(process_entry) =
                    resolve_process_reference(&completion_state.process_entries, process_value)
                else {
                    index += 1;
                    if arg_takes_value(arg) && inline_start.is_none() {
                        index = index.saturating_add(1);
                    }
                    continue;
                };

                if process_entry
                    .integrand_names
                    .iter()
                    .any(|candidate| candidate == &integrand_name)
                {
                    keys.push(format!("{}@{}", process_entry.name, integrand_name));
                }
            }
            _ => {}
        }

        index += 1;
        if arg_takes_value(arg) && inline_start.is_none() {
            index = index.saturating_add(1);
        }
    }

    keys
}

fn specified_target_keys(
    completed_tokens: &[CompletionToken],
    root_cmd: &clap::Command,
) -> HashSet<String> {
    let mut specified = HashSet::new();
    let mut cmd = root_cmd;
    let mut index = 0usize;

    while index < completed_tokens.len() {
        let token = completed_tokens[index].cooked.as_str();
        if let Some(subcmd) = find_subcommand(cmd, token) {
            cmd = subcmd;
            index += 1;
            continue;
        }

        let Some(arg) = find_flag(cmd, token) else {
            index += 1;
            continue;
        };

        let inline_start = inline_flag_value_start(token);
        let value = if let Some(start) = inline_start {
            Some(token[start..].to_string())
        } else if arg_takes_value(arg) {
            completed_tokens
                .get(index + 1)
                .map(|next| next.cooked.clone())
        } else {
            None
        };

        if arg.get_long() == Some("target") {
            if let Some(value) = value {
                if let Some((key, _)) = value.split_once('=') {
                    specified.insert(key.to_string());
                }
            }
        }

        index += 1;
        if arg_takes_value(arg) && inline_start.is_none() {
            index = index.saturating_add(1);
        }
    }

    specified
}

fn resolve_process_reference<'a>(
    process_entries: &'a [ProcessCompletionEntry],
    value: &str,
) -> Option<&'a ProcessCompletionEntry> {
    if let Some(id) = value
        .strip_prefix('#')
        .and_then(|rest| rest.parse::<usize>().ok())
    {
        return process_entries.iter().find(|entry| entry.id == id);
    }

    if let Some(name) = value.strip_prefix("name:") {
        return process_entries.iter().find(|entry| entry.name == name);
    }

    let name_match = process_entries.iter().find(|entry| entry.name == value);
    let id_match = value
        .parse::<usize>()
        .ok()
        .and_then(|id| process_entries.iter().find(|entry| entry.id == id));

    match (id_match, name_match) {
        (Some(_), Some(_)) => None,
        (Some(entry), None) | (None, Some(entry)) => Some(entry),
        (None, None) => None,
    }
}

fn tokenize_completion_input(line: &str, pos: usize) -> Vec<CompletionToken> {
    let prefix = &line[..pos];
    let mut tokens = tokenize_shell_prefix(prefix, 0);

    if prefix.is_empty() || prefix_ends_with_separator_whitespace(prefix) {
        tokens.push(CompletionToken {
            raw: String::new(),
            cooked: String::new(),
            start: pos,
            end: pos,
            quote_style: QuoteStyle::None,
        });
    }

    tokens
}

fn prefix_ends_with_separator_whitespace(prefix: &str) -> bool {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum ParseMode {
        Unquoted,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut mode = ParseMode::Unquoted;
    let mut chars = prefix.char_indices().peekable();
    let mut ended_with_separator = false;

    while let Some((_, ch)) = chars.next() {
        match mode {
            ParseMode::Unquoted => match ch {
                '\\' => {
                    ended_with_separator = false;
                    let _ = chars.next();
                }
                '\'' => {
                    ended_with_separator = false;
                    mode = ParseMode::SingleQuoted;
                }
                '"' => {
                    ended_with_separator = false;
                    mode = ParseMode::DoubleQuoted;
                }
                _ if ch.is_whitespace() => {
                    ended_with_separator = true;
                }
                _ => {
                    ended_with_separator = false;
                }
            },
            ParseMode::SingleQuoted => {
                ended_with_separator = false;
                if ch == '\'' {
                    mode = ParseMode::Unquoted;
                }
            }
            ParseMode::DoubleQuoted => match ch {
                '\\' => {
                    ended_with_separator = false;
                    let _ = chars.next();
                }
                '"' => {
                    ended_with_separator = false;
                    mode = ParseMode::Unquoted;
                }
                _ => {
                    ended_with_separator = false;
                }
            },
        }
    }

    ended_with_separator
}

fn tokenize_shell_prefix(prefix: &str, base_offset: usize) -> Vec<CompletionToken> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum ParseMode {
        Unquoted,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut tokens = Vec::new();
    let mut mode = ParseMode::Unquoted;
    let mut raw = String::new();
    let mut cooked = String::new();
    let mut start = None;
    let mut quote_style = QuoteStyle::None;
    let mut chars = prefix.char_indices().peekable();

    let finalize = |tokens: &mut Vec<CompletionToken>,
                    raw: &mut String,
                    cooked: &mut String,
                    start: &mut Option<usize>,
                    quote_style: &mut QuoteStyle,
                    end: usize| {
        if let Some(start) = start.take() {
            tokens.push(CompletionToken {
                raw: std::mem::take(raw),
                cooked: std::mem::take(cooked),
                start,
                end,
                quote_style: *quote_style,
            });
            *quote_style = QuoteStyle::None;
        }
    };

    while let Some((offset, ch)) = chars.next() {
        let absolute_offset = base_offset + offset;
        if mode == ParseMode::Unquoted && ch.is_whitespace() {
            finalize(
                &mut tokens,
                &mut raw,
                &mut cooked,
                &mut start,
                &mut quote_style,
                absolute_offset,
            );
            continue;
        }

        if start.is_none() {
            start = Some(absolute_offset);
        }

        match mode {
            ParseMode::Unquoted => match ch {
                '\\' => {
                    raw.push(ch);
                    if let Some((_, next)) = chars.next() {
                        raw.push(next);
                        cooked.push(next);
                    }
                }
                '\'' => {
                    if raw.is_empty() {
                        quote_style = QuoteStyle::Single;
                    }
                    raw.push(ch);
                    mode = ParseMode::SingleQuoted;
                }
                '"' => {
                    if raw.is_empty() {
                        quote_style = QuoteStyle::Double;
                    }
                    raw.push(ch);
                    mode = ParseMode::DoubleQuoted;
                }
                _ => {
                    raw.push(ch);
                    cooked.push(ch);
                }
            },
            ParseMode::SingleQuoted => {
                raw.push(ch);
                if ch == '\'' {
                    mode = ParseMode::Unquoted;
                } else {
                    cooked.push(ch);
                }
            }
            ParseMode::DoubleQuoted => match ch {
                '\\' => {
                    raw.push(ch);
                    if let Some((_, next)) = chars.next() {
                        raw.push(next);
                        cooked.push(next);
                    }
                }
                '"' => {
                    raw.push(ch);
                    mode = ParseMode::Unquoted;
                }
                _ => {
                    raw.push(ch);
                    cooked.push(ch);
                }
            },
        }
    }

    finalize(
        &mut tokens,
        &mut raw,
        &mut cooked,
        &mut start,
        &mut quote_style,
        base_offset + prefix.len(),
    );

    tokens
}

fn parse_shell_fragment(fragment: &str, base_offset: usize) -> CompletionToken {
    tokenize_shell_prefix(fragment, base_offset)
        .into_iter()
        .next()
        .unwrap_or(CompletionToken {
            raw: String::new(),
            cooked: String::new(),
            start: base_offset,
            end: base_offset,
            quote_style: QuoteStyle::None,
        })
}

fn add_path_suggestions(
    suggestions: &mut Vec<reedline::Suggestion>,
    seen: &mut HashSet<String>,
    request: &PathCompletionRequest,
    pos: usize,
) {
    let Some(path_candidates) = complete_path(&request.partial_path) else {
        return;
    };

    for candidate in path_candidates {
        let value = render_path_completion(&candidate.path, request.quote_style, candidate.is_dir);
        if !seen.insert(value.clone()) {
            continue;
        }

        suggestions.push(reedline::Suggestion {
            value,
            description: candidate.description,
            style: None,
            extra: None,
            span: Span::new(request.span_start, pos),
            append_whitespace: !candidate.is_dir,
        });
    }
}

#[derive(Debug, Clone)]
struct PathCompletionCandidate {
    path: String,
    is_dir: bool,
    description: Option<String>,
}

/// Complete file and directory paths
fn complete_path(partial_path: &str) -> Option<Vec<PathCompletionCandidate>> {
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

                    suggestions.push(PathCompletionCandidate {
                        path: full_path,
                        is_dir,
                        description,
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
            let a_is_dir = a.is_dir;
            let b_is_dir = b.is_dir;

            match (a_is_dir, b_is_dir) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                _ => a.path.cmp(&b.path),
            }
        });

        Some(suggestions)
    }
}

fn render_path_completion(path: &str, quote_style: QuoteStyle, is_dir: bool) -> String {
    match quote_style {
        QuoteStyle::None => shell_escape_unquoted(path),
        QuoteStyle::Double => {
            let mut rendered = "\"".to_string();
            for ch in path.chars() {
                if matches!(ch, '\\' | '"' | '$' | '`') {
                    rendered.push('\\');
                }
                rendered.push(ch);
            }
            if !is_dir {
                rendered.push('"');
            }
            rendered
        }
        QuoteStyle::Single => {
            let escaped = path.replace('\'', "'\\''");
            let mut rendered = format!("'{}", escaped);
            if !is_dir {
                rendered.push('\'');
            }
            rendered
        }
    }
}

fn render_value_completion(value: &str, quote_style: QuoteStyle) -> String {
    match quote_style {
        QuoteStyle::None => shell_escape_unquoted(value),
        QuoteStyle::Double => {
            let mut rendered = "\"".to_string();
            for ch in value.chars() {
                if matches!(ch, '\\' | '"' | '$' | '`') {
                    rendered.push('\\');
                }
                rendered.push(ch);
            }
            rendered.push('"');
            rendered
        }
        QuoteStyle::Single => {
            let escaped = value.replace('\'', "'\\''");
            format!("'{}'", escaped)
        }
    }
}

fn value_completion_appends_whitespace(value: &str) -> bool {
    !(value.ends_with('=') || value.ends_with('.'))
}

fn shell_escape_unquoted(path: &str) -> String {
    let mut escaped = String::with_capacity(path.len());
    for ch in path.chars() {
        if matches!(
            ch,
            ' ' | '\t'
                | '\n'
                | '\\'
                | '\''
                | '"'
                | ';'
                | '&'
                | '|'
                | '<'
                | '>'
                | '('
                | ')'
                | '['
                | ']'
                | '{'
                | '}'
                | '$'
                | '`'
                | '!'
        ) {
            escaped.push('\\');
        }
        escaped.push(ch);
    }
    escaped
}

#[cfg(test)]
mod tests {
    use std::{collections::BTreeMap, fs};

    use clap::CommandFactory;
    use gammalooprs::model::ParameterType;
    use serde_json::json;
    use tempfile::tempdir;

    use crate::{
        commands::process_settings::ProcessSettingsCompletionEntry,
        completion::{arg_value_completion, ArgValueCompletion},
        Repl,
    };

    use super::{
        collect_completions, CompletionState, IrProfileCompletionEntry,
        ModelParameterCompletionEntry, ProcessCompletionEntry, ProcessKind,
    };

    fn sample_process_entries() -> Vec<ProcessCompletionEntry> {
        vec![
            ProcessCompletionEntry {
                id: 0,
                name: "triangle".to_string(),
                kind: ProcessKind::Amplitude,
                integrand_names: vec!["LO".to_string(), "NLO".to_string()],
            },
            ProcessCompletionEntry {
                id: 1,
                name: "epem_a_tth".to_string(),
                kind: ProcessKind::Amplitude,
                integrand_names: vec!["LO".to_string(), "virtual".to_string()],
            },
            ProcessCompletionEntry {
                id: 2,
                name: "epem_xs".to_string(),
                kind: ProcessKind::CrossSection,
                integrand_names: vec!["LO".to_string(), "subtracted".to_string()],
            },
        ]
    }

    fn sample_ir_profile_entries() -> Vec<IrProfileCompletionEntry> {
        vec![IrProfileCompletionEntry {
            process_name: "epem_xs".to_string(),
            integrand_name: "subtracted".to_string(),
            graph_names: vec!["GL1".to_string(), "GL2".to_string()],
            graph_limit_entries: vec![
                "GL1 C[1,2]".to_string(),
                "GL1 S(3)".to_string(),
                "GL2 C[4,5]".to_string(),
            ],
        }]
    }

    fn sample_process_settings_entries() -> Vec<ProcessSettingsCompletionEntry> {
        vec![
            ProcessSettingsCompletionEntry {
                process_id: 0,
                process_name: "triangle".to_string(),
                integrand_name: "LO".to_string(),
                quantities: BTreeMap::from([(
                    "top_pt".to_string(),
                    json!({
                        "type": "particle_scalar",
                        "pdgs": [6, -6],
                        "quantity": "PT"
                    }),
                )]),
                observables: BTreeMap::from([(
                    "top_pt_hist".to_string(),
                    json!({
                        "quantity": "top_pt",
                        "entry_selection": "all",
                        "entry_index": 0,
                        "value_transform": "identity",
                        "phase": "real",
                        "misbinning_max_normalized_distance": null,
                        "x_min": 0.0,
                        "x_max": 500.0,
                        "n_bins": 50,
                        "log_x_axis": false,
                        "log_y_axis": true
                    }),
                )]),
                selectors: BTreeMap::from([(
                    "top_cut".to_string(),
                    json!({
                        "quantity": "top_pt",
                        "entry_selection": "all",
                        "entry_index": 0,
                        "selector": "value_range",
                        "min": 10.0,
                        "max": 500.0,
                        "reduction": "any_in_range"
                    }),
                )]),
            },
            ProcessSettingsCompletionEntry {
                process_id: 0,
                process_name: "triangle".to_string(),
                integrand_name: "NLO".to_string(),
                quantities: BTreeMap::from([(
                    "top_pt".to_string(),
                    json!({
                        "type": "particle_scalar",
                        "pdgs": [6, -6],
                        "quantity": "PT"
                    }),
                )]),
                observables: BTreeMap::from([(
                    "top_pt_hist".to_string(),
                    json!({
                        "quantity": "top_pt",
                        "entry_selection": "all",
                        "entry_index": 0,
                        "value_transform": "identity",
                        "phase": "real",
                        "misbinning_max_normalized_distance": null,
                        "x_min": 0.0,
                        "x_max": 500.0,
                        "n_bins": 50,
                        "log_x_axis": false,
                        "log_y_axis": true
                    }),
                )]),
                selectors: BTreeMap::from([(
                    "top_cut".to_string(),
                    json!({
                        "quantity": "top_pt",
                        "entry_selection": "all",
                        "entry_index": 0,
                        "selector": "value_range",
                        "min": 10.0,
                        "max": 500.0,
                        "reduction": "any_in_range"
                    }),
                )]),
            },
            ProcessSettingsCompletionEntry {
                process_id: 2,
                process_name: "epem_xs".to_string(),
                integrand_name: "subtracted".to_string(),
                quantities: BTreeMap::from([(
                    "mll".to_string(),
                    json!({
                        "type": "particle_scalar",
                        "pdgs": [11, -11],
                        "quantity": "E"
                    }),
                )]),
                observables: BTreeMap::new(),
                selectors: BTreeMap::new(),
            },
        ]
    }

    fn generate_completion_state() -> CompletionState {
        CompletionState {
            process_entries: sample_process_entries(),
            model_particle_names: vec![
                "g".to_string(),
                "h".to_string(),
                "e+".to_string(),
                "e-".to_string(),
            ],
            model_coupling_names: vec!["QCD".to_string(), "QED".to_string()],
            model_vertex_names: vec!["V_6".to_string(), "V_9".to_string(), "V_36".to_string()],
            process_settings_entries: sample_process_settings_entries(),
            ..CompletionState::default()
        }
    }

    fn completion_values(line: &str, completion_state: &CompletionState) -> Vec<String> {
        collect_completions::<Repl>(line, line.len(), completion_state)
            .into_iter()
            .map(|suggestion| suggestion.value)
            .collect()
    }

    fn completion_suggestions(
        line: &str,
        completion_state: &CompletionState,
    ) -> Vec<reedline::Suggestion> {
        collect_completions::<Repl>(line, line.len(), completion_state)
    }

    #[test]
    fn completion_offers_leaf_flags_on_empty_argument() {
        let values = completion_values("display model ", &CompletionState::default());

        assert!(values.contains(&"-a".to_string()));
        assert!(values.contains(&"--show-all".to_string()));
    }

    #[test]
    fn completion_uses_subcommand_aliases() {
        let values = completion_values("display settings d", &CompletionState::default());

        assert!(values.contains(&"defaults".to_string()));
    }

    #[test]
    fn completion_offers_run_block_names_and_commands_flag() {
        let completion_state = CompletionState {
            commands_block_names: vec!["alpha".to_string(), "beta".to_string()],
            ..CompletionState::default()
        };

        let values = completion_values("run ", &completion_state);

        assert!(values.contains(&"-c".to_string()));
        assert!(values.contains(&"--commands".to_string()));
        assert!(values.contains(&"alpha".to_string()));
        assert!(values.contains(&"beta".to_string()));
    }

    #[test]
    fn completion_does_not_offer_run_blocks_while_typing_commands_value() {
        let completion_state = CompletionState {
            commands_block_names: vec!["alpha".to_string()],
            ..CompletionState::default()
        };

        let values = completion_values("run -c ", &completion_state);

        assert!(values.is_empty());
    }

    #[test]
    fn completion_uses_path_hints_instead_of_name_heuristics() {
        let values = completion_values("display settings global ", &CompletionState::default());

        assert!(values.is_empty());
    }

    #[test]
    fn completion_completes_quoted_paths_safely() {
        let temp = tempdir().unwrap();
        let file_path = temp.path().join("my file.txt");
        fs::write(&file_path, "test").unwrap();

        let line = format!("import graphs \"{}", temp.path().join("my f").display());
        let values = completion_values(&line, &CompletionState::default());

        let expected = format!("\"{}", file_path.display()) + "\"";
        assert!(values.contains(&expected));
    }

    #[test]
    fn completion_completes_path_values_for_long_flags_with_equals() {
        let temp = tempdir().unwrap();
        let file_path = temp.path().join("profile.json");
        fs::write(&file_path, "test").unwrap();

        let line = format!(
            "integrate --workspace-path={}",
            temp.path().join("pro").display()
        );
        let values = completion_values(&line, &CompletionState::default());

        let expected = format!("{}", file_path.display());
        assert!(values.contains(&expected));
    }

    #[test]
    fn completion_offers_builtin_json_models_for_import_model() {
        let values = completion_values("import model s", &CompletionState::default());

        assert!(values.contains(&"sm".to_string()));
        assert!(values.contains(&"scalars".to_string()));
    }

    #[test]
    fn completion_offers_external_model_parameters_for_set_model() {
        let completion_state = CompletionState {
            model_parameter_entries: vec![
                ModelParameterCompletionEntry {
                    name: "alpha".to_string(),
                    parameter_type: ParameterType::Real,
                },
                ModelParameterCompletionEntry {
                    name: "beta".to_string(),
                    parameter_type: ParameterType::Imaginary,
                },
            ],
            ..CompletionState::default()
        };

        let values = completion_values("set model a", &completion_state);

        assert_eq!(values, vec!["alpha=".to_string()]);
    }

    #[test]
    fn completion_skips_model_parameter_suggestions_for_target_flags() {
        let completion_state = CompletionState {
            model_parameter_entries: vec![ModelParameterCompletionEntry {
                name: "alpha".to_string(),
                parameter_type: ParameterType::Real,
            }],
            ..CompletionState::default()
        };

        let values = completion_values("set model -p pro", &completion_state);

        assert!(!values.contains(&"alpha=".to_string()));
    }

    #[test]
    fn completion_offers_defaults_for_set_model() {
        let values = completion_values("set model d", &CompletionState::default());

        assert_eq!(values, vec!["defaults".to_string()]);
    }

    #[test]
    fn completion_skips_already_assigned_model_parameters() {
        let completion_state = CompletionState {
            model_parameter_entries: vec![
                ModelParameterCompletionEntry {
                    name: "alpha".to_string(),
                    parameter_type: ParameterType::Real,
                },
                ModelParameterCompletionEntry {
                    name: "beta".to_string(),
                    parameter_type: ParameterType::Real,
                },
                ModelParameterCompletionEntry {
                    name: "gamma".to_string(),
                    parameter_type: ParameterType::Imaginary,
                },
            ],
            ..CompletionState::default()
        };

        let values = completion_values("set model alpha=1 b", &completion_state);

        assert_eq!(values, vec!["beta=".to_string()]);
        assert!(!values.contains(&"alpha=".to_string()));
    }

    #[test]
    fn completion_shows_real_format_hint_for_set_model_values() {
        let suggestions = completion_suggestions(
            "set model alpha=",
            &CompletionState {
                model_parameter_entries: vec![ModelParameterCompletionEntry {
                    name: "alpha".to_string(),
                    parameter_type: ParameterType::Real,
                }],
                ..CompletionState::default()
            },
        );

        assert!(suggestions.iter().any(|suggestion| {
            suggestion.description.as_deref()
                == Some(crate::model_parameters::MODEL_REAL_VALUE_FORMAT_HINT)
                && suggestion.value == "alpha="
        }));
    }

    #[test]
    fn completion_shows_complex_format_hint_for_complex_set_model_values() {
        let suggestions = completion_suggestions(
            "set model beta=",
            &CompletionState {
                model_parameter_entries: vec![ModelParameterCompletionEntry {
                    name: "beta".to_string(),
                    parameter_type: ParameterType::Imaginary,
                }],
                ..CompletionState::default()
            },
        );

        assert!(suggestions.iter().any(|suggestion| {
            suggestion.description.as_deref()
                == Some(crate::model_parameters::MODEL_COMPLEX_VALUE_FORMAT_HINT)
                && suggestion.value == "beta="
        }));
    }

    #[test]
    fn completion_offers_process_names_for_process_selectors() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("integrate -p e", &completion_state);

        assert!(values.contains(&"epem_a_tth".to_string()));
        assert!(values.contains(&"epem_xs".to_string()));
        assert!(!values.contains(&"triangle".to_string()));
    }

    #[test]
    fn completion_offers_process_ids_when_process_prefix_requests_them() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("integrate --process #", &completion_state);

        assert!(values.contains(&"#0".to_string()));
        assert!(values.contains(&"#1".to_string()));
        assert!(values.contains(&"#2".to_string()));
    }

    #[test]
    fn completion_offers_processes_for_set_process() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("set process -p e", &completion_state);

        assert!(values.contains(&"epem_a_tth".to_string()));
        assert!(values.contains(&"epem_xs".to_string()));
        assert!(!values.contains(&"triangle".to_string()));
    }

    #[test]
    fn completion_offers_integrands_for_set_process() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("set process -p #0 -i ", &completion_state);

        assert!(values.contains(&"LO".to_string()));
        assert!(values.contains(&"NLO".to_string()));
        assert!(!values.contains(&"virtual".to_string()));
    }

    #[test]
    fn completion_keeps_repeatable_integrate_selectors_available() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("integrate -p triangle -i LO -", &completion_state);

        assert!(values.contains(&"-p".to_string()));
        assert!(values.contains(&"-i".to_string()));
    }

    #[test]
    fn completion_offers_integrands_for_reset_processes() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values(
            "reset processes --process epem_a_tth --integrand-name ",
            &completion_state,
        );

        assert!(values.contains(&"LO".to_string()));
        assert!(values.contains(&"virtual".to_string()));
        assert!(!values.contains(&"subtracted".to_string()));
    }

    #[test]
    fn completion_offers_processes_for_reset_processes() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("reset processes -p e", &completion_state);

        assert!(values.contains(&"epem_a_tth".to_string()));
        assert!(values.contains(&"epem_xs".to_string()));
        assert!(!values.contains(&"triangle".to_string()));
    }

    #[test]
    fn completion_filters_integrands_by_selected_process() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values(
            "integrate --process epem_a_tth --integrand-name ",
            &completion_state,
        );

        assert!(values.contains(&"LO".to_string()));
        assert!(values.contains(&"virtual".to_string()));
        assert!(!values.contains(&"subtracted".to_string()));
    }

    #[test]
    fn completion_filters_integrands_by_selector_kind() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let uv_values = completion_values("profile ultra-violet -i ", &completion_state);
        assert!(uv_values.contains(&"LO".to_string()));
        assert!(uv_values.contains(&"virtual".to_string()));
        assert!(!uv_values.contains(&"subtracted".to_string()));

        let ir_values = completion_values("profile infra-red -i ", &completion_state);
        assert!(ir_values.contains(&"LO".to_string()));
        assert!(ir_values.contains(&"subtracted".to_string()));
        assert!(!ir_values.contains(&"virtual".to_string()));
    }

    #[test]
    fn completion_filters_profile_processes_by_selector_kind() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let uv_values = completion_values("profile ultra-violet -p e", &completion_state);
        assert!(uv_values.contains(&"epem_a_tth".to_string()));
        assert!(!uv_values.contains(&"epem_xs".to_string()));

        let ir_values = completion_values("profile infra-red -p e", &completion_state);
        assert!(ir_values.contains(&"epem_xs".to_string()));
        assert!(!ir_values.contains(&"epem_a_tth".to_string()));
    }

    #[test]
    fn completion_filters_processes_by_selector_kind() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("evaluate --process ", &completion_state);

        assert!(values.contains(&"triangle".to_string()));
        assert!(values.contains(&"epem_a_tth".to_string()));
        assert!(!values.contains(&"epem_xs".to_string()));
    }

    #[test]
    fn completion_supports_inline_process_value_when_filtering_integrands() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values(
            "integrate --process=name:epem_xs --integrand-name=s",
            &completion_state,
        );

        assert_eq!(values, vec!["subtracted".to_string()]);
    }

    #[test]
    fn completion_offers_selected_integrand_target_keys() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values(
            "integrate -p triangle -i LO -p epem_xs -i subtracted --target ",
            &completion_state,
        );

        assert!(values.contains(&"triangle@LO=".to_string()));
        assert!(values.contains(&"epem_xs@subtracted=".to_string()));
    }

    #[test]
    fn completion_hides_already_targeted_slot_keys() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values(
            "integrate -p triangle -i LO -p epem_xs -i subtracted --target triangle@LO=0,1 --target ",
            &completion_state,
        );

        assert!(!values.contains(&"triangle@LO=".to_string()));
        assert!(values.contains(&"epem_xs@subtracted=".to_string()));
    }

    #[test]
    fn completion_does_not_offer_existing_integrands_for_free_form_names() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ..CompletionState::default()
        };

        let values = completion_values("generate amp --integrand-name ", &completion_state);

        assert!(values.is_empty());
    }

    #[test]
    fn completion_offers_global_settings_paths() {
        let values = completion_values(
            "set global kv global.generation.eva",
            &CompletionState::default(),
        );

        assert!(values.contains(&"global.generation.evaluator.".to_string()));
    }

    #[test]
    fn completion_offers_global_log_directive_paths() {
        let values = completion_values("set global kv global.", &CompletionState::default());

        assert!(values.contains(&"global.display_directive=".to_string()));
        assert!(values.contains(&"global.logfile_directive=".to_string()));
    }

    #[test]
    fn completion_does_not_offer_removed_dummy_global_setting() {
        let values = completion_values("set global kv ", &CompletionState::default());

        assert!(!values.contains(&"dummy=".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_runtime_settings_paths_for_set_process_after_flags() {
        let values = completion_values(
            "set process -p triangle -i LO kv integ",
            &CompletionState::default(),
        );

        assert!(values.contains(&"integrator.".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_observable_update_paths_after_name() {
        let values = completion_values(
            "set process -p triangle -i LO update observable top_pt_hist x_",
            &generate_completion_state(),
        );

        assert!(values.contains(&"x_max=".to_string()), "{values:?}");
        assert!(values.contains(&"x_min=".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_observable_update_enum_values() {
        let values = completion_values(
            "set process -p triangle -i LO update observable top_pt_hist phase=",
            &generate_completion_state(),
        );

        assert!(values.contains(&"real".to_string()), "{values:?}");
        assert!(values.contains(&"imag".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_selector_update_paths_after_name() {
        let values = completion_values(
            "set process -p triangle -i LO update selector top_cut entry_",
            &generate_completion_state(),
        );

        assert!(values.contains(&"entry_index=".to_string()), "{values:?}");
        assert!(
            values.contains(&"entry_selection=".to_string()),
            "{values:?}"
        );
    }

    #[test]
    fn completion_offers_quantity_update_paths_after_name() {
        let values = completion_values(
            "set process -p triangle -i LO update quantity top_pt p",
            &generate_completion_state(),
        );

        assert!(values.contains(&"pdgs=".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_existing_quantity_names_for_process_update() {
        let values = completion_values(
            "set process -p triangle update quantity ",
            &generate_completion_state(),
        );

        assert_eq!(values, vec!["top_pt".to_string()]);
    }

    #[test]
    fn completion_offers_quantity_kinds_for_process_add() {
        let values = completion_values(
            "set process -p triangle add quantity my_quantity ",
            &generate_completion_state(),
        );

        assert!(
            values.contains(&"particle_scalar".to_string()),
            "{values:?}"
        );
        assert!(values.contains(&"jet_pt".to_string()), "{values:?}");
        assert!(values.contains(&"cross_section".to_string()), "{values:?}");
    }

    #[test]
    fn completion_shows_type_hint_for_process_add_name() {
        let suggestions = completion_suggestions(
            "set process -p triangle add quantity ",
            &generate_completion_state(),
        );

        assert!(
            suggestions
                .iter()
                .any(|suggestion| suggestion.description.as_deref() == Some("expects a string")),
            "{suggestions:?}"
        );
    }

    #[test]
    fn completion_offers_selector_kinds_for_process_add() {
        let values = completion_values(
            "set process -p triangle add selector my_cut ",
            &generate_completion_state(),
        );

        assert!(values.contains(&"value_range".to_string()), "{values:?}");
        assert!(values.contains(&"count_range".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_quantity_names_for_observable_quantity_references() {
        let values = completion_values(
            "set process -p triangle -i LO add observable my_hist quantity=",
            &generate_completion_state(),
        );

        assert_eq!(values, vec!["top_pt".to_string()]);
    }

    #[test]
    fn completion_offers_existing_selector_names_for_process_remove() {
        let values = completion_values(
            "set process -p triangle remove selector ",
            &generate_completion_state(),
        );

        assert_eq!(values, vec!["top_cut".to_string()]);
    }

    #[test]
    fn completion_offers_display_quantity_names() {
        let values = completion_values(
            "display quantities -p triangle ",
            &generate_completion_state(),
        );

        assert_eq!(values, vec!["top_pt".to_string()]);
    }

    #[test]
    fn completion_offers_process_and_integrand_selectors_for_display_quantities() {
        let completion_state = generate_completion_state();

        let process_values = completion_values("display quantities -p t", &completion_state);
        assert_eq!(process_values, vec!["triangle".to_string()]);

        let integrand_values =
            completion_values("display quantities -p triangle -i ", &completion_state);
        assert!(integrand_values.contains(&"LO".to_string()));
        assert!(integrand_values.contains(&"NLO".to_string()));
    }

    #[test]
    fn completion_offers_command_block_names_for_display_command_block() {
        let completion_state = CompletionState {
            commands_block_names: vec!["alpha".to_string(), "beta".to_string()],
            ..CompletionState::default()
        };

        let values = completion_values("display command_block ", &completion_state);

        assert!(values.contains(&"alpha".to_string()));
        assert!(values.contains(&"beta".to_string()));
    }

    #[test]
    fn completion_offers_boolean_values_for_settings_leaf() {
        let values = completion_values(
            "set global kv global.generation.evaluator.iterative_orientation_optimization=",
            &CompletionState::default(),
        );

        assert!(values.contains(&"false".to_string()));
        assert!(values.contains(&"true".to_string()));
    }

    #[test]
    fn completion_offers_enum_values_for_settings_leaf() {
        let values = completion_values(
            "set global kv global.log_style.log_format=",
            &CompletionState::default(),
        );

        assert!(values.contains(&"Long".to_string()), "{values:?}");
        assert!(values.contains(&"Short".to_string()), "{values:?}");
        assert!(values.contains(&"Min".to_string()), "{values:?}");
        assert!(values.contains(&"None".to_string()), "{values:?}");
    }

    #[test]
    fn completion_shows_type_hint_for_non_boolean_settings_leaf() {
        let suggestions = completion_suggestions(
            "set default-runtime kv integrator.n_start=",
            &CompletionState::default(),
        );

        assert!(
            suggestions
                .iter()
                .any(|suggestion| suggestion.description.as_deref() == Some("expects an integer")),
            "{suggestions:?}"
        );
    }

    #[test]
    fn completion_does_not_repeat_used_options_by_short_or_long_form() {
        let values =
            completion_values("integrate --process triangle ", &CompletionState::default());

        assert!(!values.contains(&"-p".to_string()));
        assert!(!values.contains(&"--process".to_string()));
    }

    #[test]
    fn completion_does_not_append_whitespace_for_settings_assignment_keys() {
        let suggestions = completion_suggestions(
            "set global kv global.generation.evaluator.iterative_orientation_optimization",
            &CompletionState::default(),
        );

        let suggestion = suggestions
            .iter()
            .find(|suggestion| {
                suggestion.value
                    == "global.generation.evaluator.iterative_orientation_optimization="
            })
            .unwrap();
        assert!(!suggestion.append_whitespace);
    }

    #[test]
    fn completion_prefers_options_first_when_new_flags_can_start() {
        let suggestions =
            completion_suggestions("generate amp g g > ", &generate_completion_state());

        assert!(suggestions.first().unwrap().value.starts_with('-'));
    }

    #[test]
    fn completion_offers_generate_particles_and_arrow() {
        let values = completion_values("generate amp g ", &generate_completion_state());

        assert!(values.contains(&"\\>".to_string()));
        assert!(values.contains(&"to".to_string()));
        assert!(values.contains(&"g".to_string()));
        assert!(values.contains(&"h".to_string()));
    }

    #[test]
    fn completion_offers_generate_process_option_starters_and_couplings() {
        let values = completion_values("generate amp g g > h ", &generate_completion_state());

        assert!(values.contains(&"/".to_string()));
        assert!(values.contains(&"\\|".to_string()));
        assert!(values.contains(&"QCD==".to_string()));
        assert!(values.contains(&"QED\\>=".to_string()));
    }

    #[test]
    fn completion_offers_generate_xs_powered_coupling_skeletons() {
        let values = completion_values("generate xs g g > h Q", &generate_completion_state());

        assert!(values.contains(&"QCD^".to_string()));
        assert!(values.contains(&"QED^".to_string()));
    }

    #[test]
    fn completion_offers_generate_perturbative_block_suggestions() {
        let values = completion_values("generate amp g g > h [", &generate_completion_state());

        assert!(values.contains(&"\\[\\{1\\}".to_string()));
        assert!(values.contains(&"\\[\\{\\{1\\}\\}".to_string()));
        assert!(values.contains(&"\\[QED".to_string()));
        assert!(values.contains(&"\\[QED=".to_string()));
    }

    #[test]
    fn completion_offers_generate_value_enum_variants() {
        let values = completion_values(
            "generate amp g g > h --numerator-grouping ",
            &generate_completion_state(),
        );

        assert!(values.contains(&"no_grouping".to_string()));
        assert!(values.contains(&"only_detect_zeroes".to_string()));
        assert!(values.contains(&"group_identical_graphs_up_to_sign".to_string()));
    }

    #[test]
    fn completion_offers_generate_vertex_names_for_list_flags() {
        let values = completion_values(
            "generate amp g g > h --veto-vertex-interactions V_",
            &generate_completion_state(),
        );

        assert!(values.contains(&"V_6".to_string()));
        assert!(values.contains(&"V_9".to_string()));
        assert!(values.contains(&"V_36".to_string()));
    }

    #[test]
    fn completion_keeps_generate_vertex_list_active_and_filters_existing_values() {
        let values = completion_values(
            "generate amp g g > h --allowed-vertex-interactions V_6 ",
            &generate_completion_state(),
        );

        assert!(!values.contains(&"V_6".to_string()));
        assert!(values.contains(&"V_9".to_string()));
        assert!(values.contains(&"V_36".to_string()));
    }

    #[test]
    fn completion_switches_back_to_flags_after_variadic_generate_values() {
        let values = completion_values(
            "generate amp g g > h --allowed-vertex-interactions V_6 V_9 -",
            &generate_completion_state(),
        );

        assert!(
            values.contains(&"--only-diagrams".to_string()),
            "{values:?}"
        );
        assert!(
            values.contains(&"--symmetrize-final-states".to_string()),
            "{values:?}"
        );
    }

    #[test]
    fn completion_offers_builtin_model_restriction_suffixes() {
        let values = completion_values("import model sm-d", &CompletionState::default());

        assert!(values.contains(&"sm-default".to_string()), "{values:?}");
    }

    #[test]
    fn completion_offers_ir_profile_select_graphs_and_limits() {
        let completion_state = CompletionState {
            process_entries: sample_process_entries(),
            ir_profile_entries: sample_ir_profile_entries(),
            ..CompletionState::default()
        };

        let graph_values = completion_values(
            "profile infra-red -p epem_xs -i subtracted --select GL",
            &completion_state,
        );
        assert!(
            graph_values.contains(&"GL1".to_string()),
            "{graph_values:?}"
        );
        assert!(
            graph_values.contains(&"GL2".to_string()),
            "{graph_values:?}"
        );

        let limit_values = completion_values(
            "profile infra-red -p epem_xs -i subtracted --select GL1\\ ",
            &completion_state,
        );
        assert!(
            limit_values.contains(&"GL1\\ C\\[1,2\\]".to_string()),
            "{limit_values:?}"
        );
        assert!(
            limit_values.contains(&"GL1\\ S\\(3\\)".to_string()),
            "{limit_values:?}"
        );
    }

    fn visit_args(
        cmd: &clap::Command,
        command_path: &mut Vec<String>,
        visit: &mut impl FnMut(&[String], &clap::Arg),
    ) {
        for arg in cmd.get_arguments() {
            visit(command_path, arg);
        }

        for subcommand in cmd.get_subcommands() {
            command_path.push(subcommand.get_name().to_string());
            visit_args(subcommand, command_path, visit);
            command_path.pop();
        }
    }

    #[test]
    fn completion_metadata_covers_all_process_and_integrand_name_args() {
        let command = Repl::command();
        let mut command_path = vec![command.get_name().to_string()];
        let mut missing = Vec::new();
        let mut invalid_process_metadata = Vec::new();

        visit_args(
            &command,
            &mut command_path,
            &mut |path, arg| match arg.get_long() {
                Some("process") => match arg_value_completion(arg) {
                    Some(ArgValueCompletion::ProcessSelector(_)) => {}
                    other => invalid_process_metadata.push(format!(
                        "{} --process => {:?}",
                        path.join(" "),
                        other
                    )),
                },
                Some("integrand-name") => {
                    if arg_value_completion(arg).is_none() {
                        missing.push(format!("{} --integrand-name", path.join(" ")));
                    }
                }
                _ => {}
            },
        );

        assert!(
            missing.is_empty(),
            "missing completion metadata: {missing:?}"
        );
        assert!(
            invalid_process_metadata.is_empty(),
            "invalid process completion metadata: {invalid_process_metadata:?}"
        );
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

        match split_command_line(&line) {
            Ok(split) => {
                match C::try_parse_from(std::iter::once("").chain(split.iter().map(String::as_str)))
                {
                    Ok(c) => ReadCommandOutput::Command(c, line),
                    Err(e) => ReadCommandOutput::ClapError(e),
                }
            }
            Err(_) => ReadCommandOutput::ShlexError,
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
