use std::{
    fs::{self, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
    sync::{Arc, LazyLock, Mutex, OnceLock},
    time::Duration,
};

use chrono::{Datelike, Local, SecondsFormat, Timelike};
use colored::{ColoredString, Colorize};
use eyre::Context;
use gammaloop_tracing_filter::GammaLogFilter;
use gammalooprs::utils::tracing::{LogFormat, LogStyle};
use indicatif::ProgressState;
use itertools::Itertools;
use tabled::{builder::Builder, settings::Style};
use tracing::{field::Visit, level_filters::LevelFilter, Event, Subscriber};
use tracing_indicatif::{filter::IndicatifFilter, style::ProgressStyle, IndicatifLayer};
use tracing_subscriber::field::RecordFields;
use tracing_subscriber::Layer;
use tracing_subscriber::{
    filter::Filtered,
    fmt::FormattedFields,
    fmt::{
        self,
        format::Writer,
        writer::{BoxMakeWriter, MakeWriter},
        FmtContext, FormatEvent, FormatFields,
    },
    layer::SubscriberExt,
    registry::LookupSpan,
    reload,
    util::SubscriberInitExt,
};

use color_eyre::Result;

fn file_filter_from(user_spec: &str) -> Result<GammaLogFilter> {
    GammaLogFilter::parse_with_default(user_spec, Some(LevelFilter::WARN)).map_err(Into::into)
}

fn display_filter_from(user_spec: &str) -> Result<GammaLogFilter> {
    GammaLogFilter::parse_with_default(user_spec, Some(LevelFilter::OFF))
        .context(format!("Trying to parse filter spec: {user_spec}"))
}

const ENV_FILE_LOG_FILTER: &str = "GL_LOGFILE_FILTER";
const ENV_NO_GL_HARD_WARNINGS: &str = "GL_NO_HARD_WARNINGS";
const ENV_DISPLAY_LOG_FILTER: &str = "GL_DISPLAY_FILTER";
const ENV_ALL_LOG_FILTER: &str = "GL_ALL_LOG_FILTER";

#[derive(Debug, Clone)]
struct StderrLogSpecState {
    base_spec: String,
    override_spec: Option<String>,
}

impl Default for StderrLogSpecState {
    fn default() -> Self {
        let base_spec = if let Ok(all) = std::env::var(ENV_ALL_LOG_FILTER) {
            all
        } else {
            std::env::var(ENV_DISPLAY_LOG_FILTER)
                .ok()
                .unwrap_or_default()
        };
        Self {
            base_spec,
            override_spec: None,
        }
    }
}

impl StderrLogSpecState {
    fn effective_spec(&self) -> &str {
        self.override_spec
            .as_deref()
            .unwrap_or(self.base_spec.as_str())
    }

    fn effective_spec_string(&self) -> String {
        self.effective_spec().to_string()
    }
}

#[derive(Debug, Clone, Default)]
struct FileLogSpecState {
    base_spec: String,
    override_spec: Option<String>,
    hard_disabled_reason: Option<String>,
}

impl FileLogSpecState {
    fn effective_spec(&self) -> &str {
        self.override_spec
            .as_deref()
            .unwrap_or(self.base_spec.as_str())
    }

    fn effective_spec_string(&self) -> String {
        self.effective_spec().to_string()
    }

    fn hard_disabled(&self) -> bool {
        self.hard_disabled_reason.is_some()
    }
}

// Statics to hold the current log specifications
static FILE_LOG_SPEC: LazyLock<Mutex<FileLogSpecState>> = LazyLock::new(|| {
    let base_spec = if let Ok(all) = std::env::var(ENV_ALL_LOG_FILTER) {
        all
    } else {
        std::env::var(ENV_FILE_LOG_FILTER).ok().unwrap_or_default()
    };
    Mutex::new(FileLogSpecState {
        base_spec,
        ..FileLogSpecState::default()
    })
});
static STDERR_LOG_SPEC: LazyLock<Mutex<StderrLogSpecState>> =
    LazyLock::new(|| Mutex::new(StderrLogSpecState::default()));

type ReloadFilterFn = Box<dyn Fn(GammaLogFilter) -> Result<()> + Send + Sync>;

struct FilterHandles {
    file_reload: Option<ReloadFilterFn>,
    stderr_reload: ReloadFilterFn,
}

static FILTER_HANDLES: OnceLock<FilterHandles> = OnceLock::new();

fn detached_filter_handles(file_reload_enabled: bool) -> FilterHandles {
    FilterHandles {
        file_reload: file_reload_enabled.then(|| Box::new(|_| Ok(())) as ReloadFilterFn),
        stderr_reload: Box::new(|_| Ok(())),
    }
}

fn log_filter_env_override(specific_env: &str) -> Option<String> {
    if let Ok(all) = std::env::var(ENV_ALL_LOG_FILTER) {
        Some(all)
    } else {
        std::env::var(specific_env).ok()
    }
}

struct LazyJsonWriterState {
    dir: PathBuf,
    filename: String,
    file: Mutex<Option<fs::File>>,
}

#[derive(Clone)]
struct LazyJsonMakeWriter {
    state: Arc<LazyJsonWriterState>,
}

struct LazyJsonWriter {
    state: Arc<LazyJsonWriterState>,
    buffer: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum LogSink {
    #[default]
    Display,
    File,
}

const DISPLAY_ONLY_FIELD_PREFIX: &str = "display.";
const FILE_ONLY_FIELD_PREFIX: &str = "file.";

fn route_field_name(name: &str, sink: LogSink) -> Option<String> {
    if let Some(stripped) = name.strip_prefix(DISPLAY_ONLY_FIELD_PREFIX) {
        return (sink == LogSink::Display).then(|| stripped.to_string());
    }
    if let Some(stripped) = name.strip_prefix(FILE_ONLY_FIELD_PREFIX) {
        return (sink == LogSink::File).then(|| stripped.to_string());
    }
    Some(name.to_string())
}

#[derive(Default)]
struct RoutedFieldsVisitor {
    sink: LogSink,
    message: Option<String>,
    fields: serde_json::Map<String, serde_json::Value>,
}

impl RoutedFieldsVisitor {
    fn for_sink(sink: LogSink) -> Self {
        Self {
            sink,
            ..Self::default()
        }
    }

    fn record_json_value(&mut self, field: &tracing::field::Field, value: serde_json::Value) {
        let Some(name) = route_field_name(field.name(), self.sink) else {
            return;
        };
        if name == "message" {
            self.message = match value {
                serde_json::Value::String(string) => Some(string),
                other => Some(other.to_string()),
            };
        } else {
            self.fields.insert(name, value);
        }
    }
}

impl Visit for RoutedFieldsVisitor {
    fn record_bool(&mut self, field: &tracing::field::Field, value: bool) {
        self.record_json_value(field, serde_json::Value::Bool(value));
    }

    fn record_i64(&mut self, field: &tracing::field::Field, value: i64) {
        self.record_json_value(field, value.into());
    }

    fn record_u64(&mut self, field: &tracing::field::Field, value: u64) {
        self.record_json_value(field, value.into());
    }

    fn record_i128(&mut self, field: &tracing::field::Field, value: i128) {
        self.record_json_value(field, serde_json::Value::String(value.to_string()));
    }

    fn record_u128(&mut self, field: &tracing::field::Field, value: u128) {
        self.record_json_value(field, serde_json::Value::String(value.to_string()));
    }

    fn record_f64(&mut self, field: &tracing::field::Field, value: f64) {
        let value = serde_json::Number::from_f64(value)
            .map(serde_json::Value::Number)
            .unwrap_or_else(|| serde_json::Value::String(value.to_string()));
        self.record_json_value(field, value);
    }

    fn record_str(&mut self, field: &tracing::field::Field, value: &str) {
        self.record_json_value(field, serde_json::Value::String(value.to_string()));
    }

    fn record_error(
        &mut self,
        field: &tracing::field::Field,
        value: &(dyn std::error::Error + 'static),
    ) {
        self.record_json_value(field, serde_json::Value::String(value.to_string()));
    }

    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        self.record_json_value(field, serde_json::Value::String(format!("{value:?}")));
    }
}

fn render_display_message(message: &str) -> String {
    message.to_string()
}

const DISPLAY_TAG_FIELDS: &[&str] = &[
    "generation",
    "integration",
    "profile",
    "persistence",
    "uv",
    "ir",
    "subtraction",
    "sampling",
    "stability",
    "observables",
    "selectors",
    "cache",
    "graph",
    "group",
    "orientation",
    "channel",
    "cut",
    "event",
    "sample",
    "iteration",
    "term",
    "solver",
    "compile",
    "inspect",
    "summary",
    "dump",
];

#[derive(Debug, Clone, Eq, PartialEq)]
struct DisplayTag {
    name: String,
    enabled: bool,
}

fn extract_display_tags(
    fields: &mut serde_json::Map<String, serde_json::Value>,
) -> Vec<DisplayTag> {
    let mut tags = fields
        .iter()
        .filter_map(|(name, value)| match value {
            serde_json::Value::Bool(enabled) => Some(DisplayTag {
                name: name.clone(),
                enabled: *enabled,
            }),
            _ => None,
        })
        .collect_vec();

    for tag in &tags {
        fields.remove(&tag.name);
    }

    tags.sort_by_key(|tag| {
        DISPLAY_TAG_FIELDS
            .iter()
            .position(|known| known == &tag.name.as_str())
            .map_or_else(
                || (usize::MAX, tag.name.clone()),
                |index| (index, String::new()),
            )
    });
    tags
}

fn render_display_tags(tags: &[DisplayTag]) -> String {
    if tags.is_empty() {
        String::new()
    } else {
        let mut rendered_tags = tags.iter().map(|tag| {
            if tag.enabled {
                format!("#{}", tag.name)
            } else {
                format!("#!{}", tag.name)
            }
        });
        format!("[{{{}}}]", rendered_tags.join(", "))
            .bright_black()
            .to_string()
    }
}

fn render_display_field_table(
    style: &LogStyle,
    fields: &serde_json::Map<String, serde_json::Value>,
) -> Option<String> {
    if !style.include_fields || fields.is_empty() {
        return None;
    }

    let mut builder = Builder::new();
    builder.push_record([
        "field".bold().blue().to_string(),
        "value".bold().blue().to_string(),
    ]);

    for (key, value) in fields {
        builder.push_record([key.clone(), display_field_value(value)]);
    }

    let table = builder.build().with(Style::blank()).to_string();
    Some(indent_block(&table, "    "))
}

fn display_field_value(value: &serde_json::Value) -> String {
    match value {
        serde_json::Value::String(string) => string.clone(),
        other => other.to_string(),
    }
}

fn indent_block(block: &str, prefix: &str) -> String {
    block
        .lines()
        .map(|line| format!("{prefix}{line}"))
        .join("\n")
}

impl LazyJsonMakeWriter {
    fn new(dir: PathBuf, filename: String) -> Self {
        Self {
            state: Arc::new(LazyJsonWriterState {
                dir,
                filename,
                file: Mutex::new(None),
            }),
        }
    }
}

impl<'a> MakeWriter<'a> for LazyJsonMakeWriter {
    type Writer = LazyJsonWriter;

    fn make_writer(&'a self) -> Self::Writer {
        LazyJsonWriter {
            state: Arc::clone(&self.state),
            buffer: Vec::new(),
        }
    }
}

impl io::Write for LazyJsonWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.buffer.extend_from_slice(buf);

        while let Some(newline_index) = self.buffer.iter().position(|byte| *byte == b'\n') {
            let line = self.buffer.drain(..=newline_index).collect::<Vec<_>>();
            self.write_sanitized_line(&line)?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        if !self.buffer.is_empty() {
            let line = std::mem::take(&mut self.buffer);
            self.write_sanitized_line(&line)?;
        }
        if let Some(file) = self.state.file.lock().unwrap().as_mut() {
            file.flush()?;
        }
        Ok(())
    }
}

impl LazyJsonWriter {
    fn write_sanitized_line(&self, line: &[u8]) -> io::Result<()> {
        let mut file = self.state.file.lock().unwrap();
        if file.is_none() {
            fs::create_dir_all(&self.state.dir)?;
            let opened = OpenOptions::new()
                .create(true)
                .append(true)
                .open(self.state.dir.join(&self.state.filename))?;
            *file = Some(opened);
        }

        let file = file.as_mut().unwrap();
        let line = sanitize_json_log_line(line);
        file.write_all(&line)
    }
}

fn sanitize_json_log_line(line: &[u8]) -> Vec<u8> {
    let mut trailing_newline = false;
    let line = if let Some(stripped) = line.strip_suffix(b"\n") {
        trailing_newline = true;
        stripped
    } else {
        line
    };

    let sanitized = match serde_json::from_slice::<serde_json::Value>(line) {
        Ok(mut value) => {
            strip_ansi_from_json_value(&mut value);
            serde_json::to_vec(&value).unwrap_or_else(|_| line.to_vec())
        }
        Err(_) => line.to_vec(),
    };

    if trailing_newline {
        let mut with_newline = sanitized;
        with_newline.push(b'\n');
        with_newline
    } else {
        sanitized
    }
}

fn strip_ansi_from_json_value(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::String(string) => {
            *string = strip_ansi_escape_codes(string);
        }
        serde_json::Value::Array(items) => {
            for item in items {
                strip_ansi_from_json_value(item);
            }
        }
        serde_json::Value::Object(map) => {
            for value in map.values_mut() {
                strip_ansi_from_json_value(value);
            }
        }
        serde_json::Value::Null | serde_json::Value::Bool(_) | serde_json::Value::Number(_) => {}
    }
}

fn strip_ansi_escape_codes(line: &str) -> String {
    let mut stripped = String::with_capacity(line.len());
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        if ch == '\u{1b}' && chars.peek() == Some(&'[') {
            chars.next();
            for code in chars.by_ref() {
                if ('@'..='~').contains(&code) {
                    break;
                }
            }
            continue;
        }

        stripped.push(ch);
    }

    stripped
}

fn directive_is_effectively_off(spec: &str) -> bool {
    GammaLogFilter::is_effectively_off(spec)
}

fn file_log_disabled_error(reason: &str) -> color_eyre::Report {
    eyre::eyre!(
        "Cannot enable logfile logging because this session was started with the logfile logger disabled ({reason})"
    )
}

pub fn configure_file_log_boot_mode(disabled: bool, reason: Option<&str>) -> Result<()> {
    let mut state = FILE_LOG_SPEC.lock().unwrap();
    state.hard_disabled_reason = if disabled {
        Some(reason.unwrap_or("CLI boot option").to_string())
    } else {
        None
    };
    if state.hard_disabled() && !directive_is_effectively_off(state.effective_spec()) {
        return Err(file_log_disabled_error(
            state
                .hard_disabled_reason
                .as_deref()
                .unwrap_or("CLI boot option"),
        ));
    }
    Ok(())
}

pub fn file_log_boot_disabled_reason() -> Option<String> {
    FILE_LOG_SPEC.lock().unwrap().hard_disabled_reason.clone()
}

/// Set the file log filter at runtime while preserving required targets.
pub fn set_file_log_filter(user_spec: impl AsRef<str>) -> Result<()> {
    let user_spec = if let Some(file) = log_filter_env_override(ENV_FILE_LOG_FILTER) {
        if std::env::var(ENV_NO_GL_HARD_WARNINGS).is_err() {
            println!(
                "WARNING, file log filter is set to {file}, will override settings {}",
                user_spec.as_ref()
            );
        }
        file
    } else {
        user_spec.as_ref().to_string()
    };
    let effective_spec = {
        let mut state = FILE_LOG_SPEC.lock().unwrap();
        let old_base = state.base_spec.clone();
        state.base_spec = user_spec;
        let effective = state.effective_spec_string();
        if let Some(reason) = state.hard_disabled_reason.clone() {
            if !directive_is_effectively_off(&effective) {
                state.base_spec = old_base;
                return Err(file_log_disabled_error(&reason));
            }
        }
        effective
    };

    let Some(handles) = FILTER_HANDLES.get() else {
        return Ok(());
    };
    let Some(reload) = handles.file_reload.as_ref() else {
        return Ok(());
    };
    let full_spec = file_filter_from(&effective_spec)?;
    reload(full_spec)?;
    Ok(())
}

pub fn set_file_log_filter_override(user_spec: Option<String>) -> Result<()> {
    let effective_spec = {
        let mut state = FILE_LOG_SPEC.lock().unwrap();
        let old_override = state.override_spec.clone();
        state.override_spec = user_spec;
        let effective = state.effective_spec_string();
        if let Some(reason) = state.hard_disabled_reason.clone() {
            if !directive_is_effectively_off(&effective) {
                state.override_spec = old_override;
                return Err(file_log_disabled_error(&reason));
            }
        }
        effective
    };

    let Some(handles) = FILTER_HANDLES.get() else {
        return Ok(());
    };
    let Some(reload) = handles.file_reload.as_ref() else {
        return Ok(());
    };
    let full_spec = file_filter_from(&effective_spec)?;
    reload(full_spec)?;
    Ok(())
}

pub fn clear_file_log_filter_override_on_settings_change() -> Result<()> {
    set_file_log_filter_override(None)
}

/// Set the stderr log filter at runtime.
pub fn set_stderr_log_filter(user_spec: impl AsRef<str>) -> Result<()> {
    let user_spec = if let Some(display) = log_filter_env_override(ENV_DISPLAY_LOG_FILTER) {
        if std::env::var(ENV_NO_GL_HARD_WARNINGS).is_err() {
            println!(
                "WARNING, display log filter is set to {display}, will override settings {}",
                user_spec.as_ref()
            );
        }
        display
    } else {
        user_spec.as_ref().to_string()
    };

    let effective_spec = {
        let mut state = STDERR_LOG_SPEC.lock().unwrap();
        state.base_spec = user_spec;
        state.effective_spec_string()
    };

    let Some(handles) = FILTER_HANDLES.get() else {
        return Ok(());
    };
    let full_spec = display_filter_from(&effective_spec)?;
    // println!("Modifying display filter to: {}", full_spec);
    (handles.stderr_reload)(full_spec)?;
    Ok(())
}

pub fn set_stderr_log_filter_override(user_spec: Option<String>) -> Result<()> {
    let effective_spec = {
        let mut state = STDERR_LOG_SPEC.lock().unwrap();
        state.override_spec = user_spec;
        state.effective_spec_string()
    };

    let Some(handles) = FILTER_HANDLES.get() else {
        return Ok(());
    };
    let full_spec = display_filter_from(&effective_spec)?;
    (handles.stderr_reload)(full_spec)?;
    Ok(())
}

/// Get the current file log filter specification.
pub fn get_file_log_filter() -> String {
    FILE_LOG_SPEC.lock().unwrap().effective_spec_string()
}

/// Get the current stderr log filter specification.
pub fn get_stderr_log_filter() -> String {
    STDERR_LOG_SPEC.lock().unwrap().effective_spec_string()
}

pub fn get_stderr_log_filter_label() -> String {
    let spec = get_stderr_log_filter();
    collapse_scoped_gamma_level(&spec).unwrap_or(spec)
}

fn collapse_scoped_gamma_level(spec: &str) -> Option<String> {
    let (left, right) = spec.split_once(',')?;
    let left = left.trim().strip_prefix("gammaloop_api=")?;
    let right = right.trim().strip_prefix("gammalooprs=")?;
    (left == right).then(|| left.to_string())
}

fn elapsed_subsec(state: &ProgressState, writer: &mut dyn std::fmt::Write) {
    let seconds = state.elapsed().as_secs();
    let sub_seconds = (state.elapsed().as_millis() % 1000) / 100;
    let _ = writer.write_str(&format!("{}.{}s", seconds, sub_seconds));
}

pub(crate) fn init_tracing(dir: impl AsRef<Path>, log_file_name: Option<String>) {
    symbolica::GLOBAL_SETTINGS
        .initialize_tracing
        .store(false, std::sync::atomic::Ordering::Relaxed);
    FILTER_HANDLES.get_or_init(|| {
        let file_state = FILE_LOG_SPEC.lock().unwrap().clone();
        let file_spec = file_state.effective_spec_string();
        // Keep startup filter semantics aligned with runtime reloads.
        let file_filter = file_filter_from(&file_spec).unwrap_or_else(|err| {
            eprintln!(
                "Invalid initial file log filter '{file_spec}': {err}. Falling back to warn."
            );
            file_filter_from("").expect("default file log filter must parse")
        });
        // println!(
        //     "File filter: {file_filter},:{}",
        //     FILE_LOG_SPEC.lock().unwrap().as_str()
        // );
        let stderr_spec = STDERR_LOG_SPEC.lock().unwrap().effective_spec_string();
        // Keep startup filter semantics aligned with runtime reloads.
        let stderr_filter = display_filter_from(&stderr_spec).unwrap_or_else(|err| {
            eprintln!(
                "Invalid initial display log filter '{stderr_spec}': {err}. Falling back to off."
            );
            display_filter_from("").expect("default display log filter must parse")
        });
        // println!(
        //     "Stderr filter: {stderr_filter},:{}",
        //     STDERR_LOG_SPEC.lock().unwrap().as_str()
        // );

        let (file_filter_layer, file_handle) = reload::Layer::new(file_filter.clone());
        let (stderr_filter_layer, stderr_handle) = reload::Layer::new(stderr_filter.clone());

        // e.g. "2025-08-26T21-05-33.123"
        let ts = Local::now()
            .to_rfc3339_opts(SecondsFormat::Millis, true)
            .replace([':', '+'], "-"); // keep it filename-friendly

        let filename = if let Some(log_name) = log_file_name {
            format!("gammalog-{log_name}.jsonl")
        } else {
            format!("gammalog-{ts}.jsonl")
        };
        let file_writer = BoxMakeWriter::new(LazyJsonMakeWriter::new(
            dir.as_ref().to_path_buf(),
            filename,
        ));

        let json = fmt::layer()
            .event_format(FileJsonFmt)
            .with_writer(file_writer);

        let indicatif_layer = IndicatifLayer::new()
            .with_span_field_formatter(IndicatifPbMsgFields{})
            .with_max_progress_bars(1024, None) .with_progress_style( ProgressStyle::with_template(
                      "{color_start}{span_child_prefix}{span_fields} -- {wide_msg} {elapsed_subsec}{color_end}",
                  )
                  .unwrap()
                  .with_key(
                      "elapsed_subsec",
                      elapsed_subsec,
                  )
                  .with_key(
                      "color_start",
                      |state: &ProgressState, writer: &mut dyn std::fmt::Write| {
                          let elapsed = state.elapsed();

                          if elapsed > Duration::from_secs(60) {
                              // Red
                              let _ = write!(writer, "\x1b[{}m", 1 + 30);
                          } else if elapsed > Duration::from_secs(4) {
                              // Yellow
                              let _ = write!(writer, "\x1b[{}m", 3 + 30);
                          }
                      },
                  )
                  .with_key(
                      "color_end",
                      |state: &ProgressState, writer: &mut dyn std::fmt::Write| {
                          if state.elapsed() > Duration::from_secs(4) {
                              let _ = write!(writer, "\x1b[0m");
                          }
                      },
                  ),
              )
              .with_span_child_prefix_symbol("↳ ")
              .with_span_child_prefix_indent(" ");

        // Pretty status layer for stderr - only show status events
        let status_layer = fmt::layer()
            .with_target(false)
            .event_format(StatusFmt)
            .with_writer(indicatif_layer.get_stderr_writer());

        let subscriber = tracing_subscriber::registry()
            .with(Filtered::new(status_layer, stderr_filter_layer))
            .with(indicatif_layer.with_filter(IndicatifFilter::new(false)));

        let init_result = if file_state.hard_disabled() {
            subscriber.try_init()
        } else {
            subscriber
                .with(Filtered::new(json, file_filter_layer))
                .try_init()
        };
        if init_result.is_err() {
            return detached_filter_handles(!file_state.hard_disabled());
        }

        let stderr_handle = stderr_handle.clone();
        let file_reload = (!file_state.hard_disabled()).then_some(Box::new(move |filter| {
            file_handle
                .modify(|f| *f = filter)
                .map_err(color_eyre::Report::from)
        }) as ReloadFilterFn);
        let stderr_reload = Box::new(move |filter| {
            stderr_handle
                .modify(|f| *f = filter)
                .map_err(color_eyre::Report::from)
        }) as ReloadFilterFn;

        FilterHandles {
            file_reload,
            stderr_reload,
        }
    });
}

#[derive(Debug, Clone, Default)]
struct LogStyleState {
    base_style: LogStyle,
    format_override: Option<LogFormat>,
}

impl LogStyleState {
    fn effective_style(&self) -> LogStyle {
        let mut style = self.base_style.clone();
        if let Some(log_format) = self.format_override {
            style.log_format = log_format;
        }
        style
    }
}

static LOG_STYLE_STATE: LazyLock<Mutex<LogStyleState>> =
    LazyLock::new(|| Mutex::new(LogStyleState::default()));

pub fn set_log_style(style: LogStyle) {
    LOG_STYLE_STATE.lock().unwrap().base_style = style;
}

pub fn set_log_format_override(log_format: Option<LogFormat>) {
    LOG_STYLE_STATE.lock().unwrap().format_override = log_format;
}

#[derive(Clone, Default)]
struct IndicatifPbMsgFields;

impl<'writer> FormatFields<'writer> for IndicatifPbMsgFields {
    fn format_fields<R: RecordFields>(
        &self,
        mut writer: Writer<'writer>,
        fields: R,
    ) -> std::fmt::Result {
        struct Visitor {
            msg: Option<String>,
        }

        impl tracing::field::Visit for Visitor {
            fn record_debug(&mut self, f: &tracing::field::Field, v: &dyn std::fmt::Debug) {
                if f.name() == "indicatif.pb_msg" {
                    self.msg = Some(format!("{v:?}"));
                }
            }

            fn record_str(&mut self, f: &tracing::field::Field, v: &str) {
                if f.name() == "indicatif.pb_msg" {
                    self.msg = Some(v.to_string());
                }
            }
        }

        let mut v = Visitor { msg: None };
        fields.record(&mut v);
        if let Some(msg) = v.msg {
            write!(writer, "{msg}")?;
        }
        Ok(())
    }
}

/// Pretty formatter for the *status* layer; mimics your fern formatting.
#[derive(Clone, Default)]
struct StatusFmt;
impl<S, N> FormatEvent<S, N> for StatusFmt
where
    S: Subscriber + for<'a> LookupSpan<'a>,
    N: for<'writer> FormatFields<'writer> + 'static,
{
    fn format_event(
        &self,
        _ctx: &FmtContext<'_, S, N>,
        mut w: Writer<'_>,
        event: &Event<'_>,
    ) -> std::fmt::Result {
        let now = Local::now();
        let meta = event.metadata();
        let style = LOG_STYLE_STATE.lock().unwrap().effective_style();

        let mut v = RoutedFieldsVisitor::for_sink(LogSink::Display);
        event.record(&mut v);
        let msg = v.message.as_deref().unwrap_or("");
        let mut display_fields = v.fields;
        let display_tags = extract_display_tags(&mut display_fields);
        let rendered = render_display_message(msg);
        let rendered_field_table = render_display_field_table(&style, &display_fields);

        let ts_long = format!(
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:03}",
            now.year(),
            now.month(),
            now.day(),
            now.hour(),
            now.minute(),
            now.second(),
            now.timestamp_subsec_millis()
        );
        let ts_short_ms = format!(
            "{:02}:{:02}:{:02}.{:03}",
            now.hour(),
            now.minute(),
            now.second(),
            now.timestamp_subsec_millis()
        );
        let ts_short = format!("{:02}:{:02}:{:02}", now.hour(), now.minute(), now.second());

        let long_source = if style.full_line_source {
            format_full_source(meta)
        } else {
            format_target(meta.module_path().unwrap_or("").to_string(), *meta.level())
        };
        let full_source = if style.full_line_source {
            format_full_source(meta)
        } else {
            format_target_full(meta.module_path().unwrap_or("").to_string())
        };
        let long_source = format!("{}{}", long_source, render_display_tags(&display_tags));
        let full_source = format!("{}{}", full_source, render_display_tags(&display_tags));

        match style.log_format {
            LogFormat::Long => {
                let timestamp = if style.short_timestamp {
                    &ts_short_ms
                } else {
                    &ts_long
                };
                write!(
                    w,
                    "[{}] @{} {}: {}",
                    timestamp,
                    long_source,
                    format_level(*meta.level()),
                    rendered
                )?;
            }
            LogFormat::Full => {
                let timestamp = if style.short_timestamp {
                    &ts_short_ms
                } else {
                    &ts_long
                };
                write!(
                    w,
                    "[{}] @{} {}: {}",
                    timestamp,
                    full_source,
                    format_level(*meta.level()),
                    rendered
                )?;
            }
            LogFormat::Short => {
                write!(
                    w,
                    "[{}] {}: {}",
                    ts_short,
                    format_level(*meta.level()),
                    rendered
                )?;
            }
            LogFormat::Min => {
                write!(w, "{}: {}", format_level(*meta.level()), rendered)?;
            }
            LogFormat::None => {
                write!(w, "{}", rendered)?;
            }
        }
        if let Some(field_table) = rendered_field_table {
            write!(w, "\n{field_table}")?;
        }
        writeln!(w)
    }
}

#[derive(Clone, Default)]
struct FileJsonFmt;

impl<S, N> FormatEvent<S, N> for FileJsonFmt
where
    S: Subscriber + for<'a> LookupSpan<'a>,
    N: for<'writer> FormatFields<'writer> + 'static,
{
    fn format_event(
        &self,
        ctx: &FmtContext<'_, S, N>,
        mut w: Writer<'_>,
        event: &Event<'_>,
    ) -> std::fmt::Result {
        let now = Local::now().to_rfc3339_opts(SecondsFormat::Millis, true);
        let meta = event.metadata();

        let mut visitor = RoutedFieldsVisitor::for_sink(LogSink::File);
        event.record(&mut visitor);

        let mut payload = serde_json::Map::new();
        payload.insert("timestamp".to_string(), serde_json::Value::String(now));
        payload.insert(
            "level".to_string(),
            serde_json::Value::String(meta.level().to_string()),
        );
        payload.insert(
            "target".to_string(),
            serde_json::Value::String(meta.target().to_string()),
        );

        if let Some(module_path) = meta.module_path() {
            payload.insert(
                "module_path".to_string(),
                serde_json::Value::String(module_path.to_string()),
            );
        }
        if let Some(file) = meta.file() {
            payload.insert(
                "file".to_string(),
                serde_json::Value::String(file.to_string()),
            );
        }
        if let Some(line) = meta.line() {
            payload.insert("line".to_string(), line.into());
        }
        if let Some(message) = visitor.message {
            payload.insert("message".to_string(), serde_json::Value::String(message));
        }

        for (key, value) in visitor.fields {
            payload.insert(key, value);
        }

        if let Some(span) = ctx.lookup_current() {
            let mut current_span = serde_json::Map::new();
            current_span.insert(
                "name".to_string(),
                serde_json::Value::String(span.name().to_string()),
            );
            if let Some(fields) = span.extensions().get::<FormattedFields<N>>() {
                let formatted = fields.to_string();
                if !formatted.is_empty() {
                    current_span.insert("fields".to_string(), serde_json::Value::String(formatted));
                }
            }
            payload.insert("span".to_string(), serde_json::Value::Object(current_span));
        }

        let line = serde_json::to_string(&payload).map_err(|_| std::fmt::Error)?;
        w.write_str(&line)?;
        w.write_char('\n')
    }
}

pub(crate) fn format_level(level: tracing::Level) -> ColoredString {
    match level {
        tracing::Level::ERROR => format!("{:<8}", "ERROR").red(),
        tracing::Level::WARN => format!("{:<8}", "WARNING").yellow(),
        tracing::Level::INFO => format!("{:<8}", "INFO").into(),
        tracing::Level::DEBUG => format!("{:<8}", "DEBUG").bright_black(),
        tracing::Level::TRACE => format!("{:<8}", "TRACE").into(),
    }
}

pub(crate) fn format_target(target: String, level: tracing::Level) -> ColoredString {
    // println!("Target: {}", target);
    let split_targets = target.split("::").collect::<Vec<_>>();
    //[-2..].iter().join("::");
    let start = split_targets.len().saturating_sub(2);
    let mut shortened_path = split_targets[start..].join("::");
    if level < tracing::Level::DEBUG && shortened_path.len() > 20 {
        shortened_path = format!("{}...", shortened_path.chars().take(17).collect::<String>());
    }
    format!("{:<20}", shortened_path).bright_blue()
}

pub(crate) fn format_target_full(target: String) -> ColoredString {
    target.bright_blue()
}

pub(crate) fn format_full_source(meta: &tracing::Metadata<'_>) -> ColoredString {
    let module = meta.module_path().unwrap_or("");
    let file = meta.file().unwrap_or("");
    let line = meta.line();

    let mut parts = Vec::new();
    if !module.is_empty() {
        parts.push(module.to_string());
    }
    if !file.is_empty() {
        let mut file_part = file.to_string();
        if let Some(line) = line {
            file_part.push(':');
            file_part.push_str(&line.to_string());
        }
        parts.push(file_part);
    } else if let Some(line) = line {
        parts.push(line.to_string());
    }

    let source = if parts.is_empty() {
        meta.target().to_string()
    } else {
        parts.join(" ")
    };

    source.bright_blue()
}

#[cfg(test)]
mod tests {
    use tracing_subscriber::filter::LevelFilter;

    use crate::tracing::strip_ansi_escape_codes;

    use super::{
        collapse_scoped_gamma_level, display_filter_from, extract_display_tags, file_filter_from,
        format_target, format_target_full, get_stderr_log_filter, get_stderr_log_filter_label,
        indent_block, render_display_field_table, render_display_message, render_display_tags,
        route_field_name, set_stderr_log_filter, set_stderr_log_filter_override, DisplayTag,
        LogSink, StderrLogSpecState, STDERR_LOG_SPEC,
    };

    struct StderrStateRestore(StderrLogSpecState);

    impl Drop for StderrStateRestore {
        fn drop(&mut self) {
            *STDERR_LOG_SPEC.lock().unwrap() = self.0.clone();
        }
    }

    fn preserve_stderr_filter_state() -> StderrStateRestore {
        StderrStateRestore(STDERR_LOG_SPEC.lock().unwrap().clone())
    }

    #[test]
    fn stderr_filter_override_supersedes_base_and_can_be_cleared() {
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        let _restore = preserve_stderr_filter_state();

        set_stderr_log_filter("warn").unwrap();
        set_stderr_log_filter_override(Some("debug".to_string())).unwrap();
        assert_eq!(get_stderr_log_filter(), "debug");

        set_stderr_log_filter("error").unwrap();
        assert_eq!(get_stderr_log_filter(), "debug");

        set_stderr_log_filter_override(None).unwrap();
        assert_eq!(get_stderr_log_filter(), "error");
    }

    #[test]
    fn scoped_gamma_directive_collapses_to_plain_level_label() {
        assert_eq!(
            collapse_scoped_gamma_level("gammaloop_api=debug,gammalooprs=debug"),
            Some("debug".to_string())
        );
        assert_eq!(
            collapse_scoped_gamma_level("gammaloop_api=debug,gammalooprs=info"),
            None
        );
        assert_eq!(
            collapse_scoped_gamma_level("symbolica=debug,gammalooprs=debug"),
            None
        );
    }

    #[test]
    fn stderr_filter_label_prefers_collapsed_scoped_gamma_directive() {
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());
        let _restore = preserve_stderr_filter_state();

        set_stderr_log_filter("warn").unwrap();
        set_stderr_log_filter_override(Some("gammaloop_api=trace,gammalooprs=trace".to_string()))
            .unwrap();
        assert_eq!(get_stderr_log_filter_label(), "trace");

        set_stderr_log_filter_override(None).unwrap();
        assert_eq!(get_stderr_log_filter_label(), "warn");
    }

    #[test]
    fn display_filter_empty_spec_defaults_to_off() {
        assert_eq!(
            display_filter_from("").unwrap().max_level_hint(),
            Some(LevelFilter::OFF)
        );
    }

    #[test]
    fn file_filter_empty_spec_defaults_to_warn() {
        assert_eq!(
            file_filter_from("").unwrap().max_level_hint(),
            Some(LevelFilter::WARN)
        );
    }

    #[test]
    fn display_filter_examples_with_multi_field_groups_parse() {
        for spec in [
            "gammalooprs::uv::forest[{generation,uv}]=debug",
            "gammalooprs::uv::forest[{generation,uv,orientation,dump}]=debug",
            "gammalooprs::uv::forest[{generation,uv,term}]=debug",
            "gammalooprs::integrands::process::amplitude[{integration,subtraction}]=debug",
            "gammalooprs::integrands::process::cross_section[{integration,sample,inspect}]=debug",
            "gammalooprs::integrands::process::cross_section[{integration,cut,solver}]=debug",
        ] {
            assert!(
                display_filter_from(spec).is_ok(),
                "display filter example should parse: {spec}"
            );
        }
    }

    #[test]
    fn file_filter_examples_with_multi_field_groups_parse() {
        let spec = "gammalooprs::uv::forest[{generation,uv,orientation,dump}]=debug";
        assert!(file_filter_from(spec).is_ok());
    }

    #[test]
    fn full_target_format_preserves_long_module_paths() {
        let target = "gammaloop_api::commands::very_long_module_name".to_string();
        assert_eq!(
            strip_ansi_escape_codes(&format_target_full(target.clone()).to_string()),
            target
        );
        assert!(format_target(target, tracing::Level::INFO)
            .to_string()
            .contains("..."));
    }

    #[test]
    fn route_field_name_hides_display_only_fields_from_file_sink() {
        assert_eq!(
            route_field_name("display.progress", LogSink::Display),
            Some("progress".to_string())
        );
        assert_eq!(route_field_name("display.progress", LogSink::File), None);
        assert_eq!(
            route_field_name("plain", LogSink::Display),
            Some("plain".to_string())
        );
        assert_eq!(
            route_field_name("plain", LogSink::File),
            Some("plain".to_string())
        );
    }

    #[test]
    fn route_field_name_hides_file_only_fields_from_display_sink() {
        assert_eq!(route_field_name("file.json", LogSink::Display), None);
        assert_eq!(
            route_field_name("file.json", LogSink::File),
            Some("json".to_string())
        );
    }

    #[test]
    fn display_message_keeps_only_message_text() {
        assert_eq!(render_display_message("hello"), "hello");
    }

    #[test]
    fn display_formatter_extracts_all_boolean_fields_as_tags() {
        let mut fields = serde_json::Map::from_iter([
            ("generation".to_string(), serde_json::Value::Bool(true)),
            ("uv".to_string(), serde_json::Value::Bool(true)),
            ("inspect".to_string(), serde_json::Value::Bool(false)),
            ("custom".to_string(), serde_json::Value::Bool(true)),
            (
                "progress".to_string(),
                serde_json::Value::String("step-1".to_string()),
            ),
        ]);

        let tags = extract_display_tags(&mut fields);
        assert_eq!(
            tags,
            vec![
                DisplayTag {
                    name: "generation".to_string(),
                    enabled: true,
                },
                DisplayTag {
                    name: "uv".to_string(),
                    enabled: true,
                },
                DisplayTag {
                    name: "inspect".to_string(),
                    enabled: false,
                },
                DisplayTag {
                    name: "custom".to_string(),
                    enabled: true,
                },
            ]
        );
        assert_eq!(fields.get("generation"), None);
        assert_eq!(fields.get("uv"), None);
        assert_eq!(fields.get("inspect"), None);
        assert_eq!(fields.get("custom"), None);
        assert_eq!(
            fields.get("progress"),
            Some(&serde_json::Value::String("step-1".to_string()))
        );
    }

    #[test]
    fn display_tags_render_next_to_source() {
        assert_eq!(
            strip_ansi_escape_codes(&render_display_tags(&[
                DisplayTag {
                    name: "generation".to_string(),
                    enabled: true,
                },
                DisplayTag {
                    name: "inspect".to_string(),
                    enabled: false,
                },
                DisplayTag {
                    name: "custom".to_string(),
                    enabled: true,
                },
            ])),
            "[{#generation, #!inspect, #custom}]"
        );
    }

    #[test]
    fn full_source_with_display_tags_is_copy_pasteable_as_directive() {
        let directive_head = format!(
            "{}{}",
            strip_ansi_escape_codes(
                &format_target_full("gammalooprs::uv::forest".to_string()).to_string()
            ),
            strip_ansi_escape_codes(&render_display_tags(&[
                DisplayTag {
                    name: "generation".to_string(),
                    enabled: true,
                },
                DisplayTag {
                    name: "uv".to_string(),
                    enabled: true,
                },
                DisplayTag {
                    name: "inspect".to_string(),
                    enabled: false,
                },
            ])),
        );

        assert_eq!(
            directive_head,
            "gammalooprs::uv::forest[{#generation, #uv, #!inspect}]"
        );
        assert!(display_filter_from(&format!("{directive_head}=debug")).is_ok());
    }

    #[test]
    fn display_field_table_renders_only_non_boolean_fields() {
        let style = gammalooprs::utils::tracing::LogStyle {
            include_fields: true,
            ..Default::default()
        };
        let mut fields = serde_json::Map::from_iter([
            ("inspect".to_string(), serde_json::Value::Bool(false)),
            (
                "progress".to_string(),
                serde_json::Value::String("step-1".to_string()),
            ),
        ]);

        let _tags = extract_display_tags(&mut fields);
        let rendered = render_display_field_table(&style, &fields).unwrap();
        let rendered = strip_ansi_escape_codes(&rendered);
        assert!(rendered.contains("field"));
        assert!(rendered.contains("value"));
        assert!(!rendered.contains("inspect"));
        assert!(!rendered.contains("false"));
        assert!(rendered.contains("progress"));
        assert!(rendered.contains("step-1"));
        assert!(rendered.lines().all(|line| line.starts_with("    ")));
    }

    #[test]
    fn indent_block_prefixes_each_line() {
        assert_eq!(indent_block("a\nb", "  "), "  a\n  b");
    }
}
