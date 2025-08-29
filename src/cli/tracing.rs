use crate::cli::state::{format_level, format_target};

use super::state::LOG_SPEC;
use super::LogFormat;
use bincode_trait_derive::{Decode, Encode};
use chrono::{Datelike, Local, SecondsFormat, Timelike};
use clap::ValueEnum;
use momtrop::log::Logger;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::sync::{LazyLock, Mutex, OnceLock};
use tracing::level_filters::LevelFilter;
use tracing_appender::non_blocking::{NonBlockingBuilder, WorkerGuard};
use tracing_appender::rolling::{RollingFileAppender, Rotation};
use tracing_subscriber::filter::filter_fn;
use tracing_subscriber::fmt::format::Writer;
use tracing_subscriber::fmt::{FmtContext, FormatEvent, FormatFields};
use tracing_subscriber::registry::LookupSpan;
use tracing_subscriber::{fmt, prelude::*, registry::Registry, reload, EnvFilter};

#[repr(usize)]
#[derive(
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Debug,
    Hash,
    ValueEnum,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    JsonSchema,
)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
pub enum LogLevel {
    /// A level lower than all log levels.
    Off,
    /// Corresponds to the `Error` log level.
    Error,
    /// Corresponds to the `Warn` log level.
    Warn,
    /// Corresponds to the `Info` log level.
    Info,
    /// Corresponds to the `Debug` log level.
    Debug,
    /// Corresponds to the `Trace` log level.
    Trace,
}

impl Default for LogLevel {
    fn default() -> Self {
        LogLevel::Info
    }
}

impl From<LogLevel> for LevelFilter {
    fn from(value: LogLevel) -> Self {
        match value {
            LogLevel::Debug => LevelFilter::DEBUG,
            LogLevel::Trace => LevelFilter::TRACE,
            LogLevel::Warn => LevelFilter::WARN,
            LogLevel::Error => LevelFilter::ERROR,
            LogLevel::Off => LevelFilter::OFF,
            LogLevel::Info => LevelFilter::INFO,
        }
    }
}

impl LogLevel {
    pub fn to_env_spec(self) -> &'static str {
        match self {
            LogLevel::Off => "off",
            LogLevel::Error => "error",
            LogLevel::Warn => "warn",
            LogLevel::Info => "info",
            LogLevel::Debug => "debug",
            LogLevel::Trace => "trace",
        }
    }
}

// Global one-time slots
pub(super) static FILTER_HANDLE: OnceLock<reload::Handle<EnvFilter, Registry>> = OnceLock::new();
pub(super) static LOG_GUARD: OnceLock<WorkerGuard> = OnceLock::new();

pub(crate) fn init_tracing(
    default_spec: impl AsRef<str>,
    dir: impl AsRef<Path>,
) -> reload::Handle<EnvFilter, Registry> {
    FILTER_HANDLE
        .get_or_init(|| {
            // 2) reloadable filter
            let spec =
                std::env::var("RUST_LOG").unwrap_or_else(|_| default_spec.as_ref().to_string());
            LOG_SPEC.set(Mutex::new(spec.clone())).ok();
            let (filter_layer, handle) = reload::Layer::new(EnvFilter::new(spec));

            let _ = std::fs::create_dir_all(dir.as_ref());

            // e.g. "2025-08-26T21-05-33.123"
            let ts = Local::now()
                .to_rfc3339_opts(SecondsFormat::Millis, true)
                .replace(':', "-")
                .replace('+', "-"); // keep it filename-friendly

            let filename = format!("gammalog-{ts}.jsonl");
            // One file per state (per process), no rotation
            let file = RollingFileAppender::new(Rotation::NEVER, dir.as_ref(), &filename);
            let (nb, guard) = NonBlockingBuilder::default()
                .buffered_lines_limit(200_000)
                .lossy(false)
                .finish(file);
            LOG_GUARD.set(guard).ok();

            let json = fmt::layer()
                .json()
                .flatten_event(true)
                .with_current_span(true)
                .with_span_list(true)
                .with_writer(nb);

            // 4) pretty status to stderr, opt-in via `target="status"`
            let status_layer = fmt::layer()
                .with_target(false)
                .event_format(StatusFmt::default())
                .with_writer(std::io::stderr)
                .with_filter(filter_fn(|m| m.target() == "status"));

            tracing_subscriber::registry()
                .with(filter_layer)
                .with(json)
                .with(status_layer)
                .init();

            handle
        })
        .clone()
}

pub fn init_test_tracing() -> reload::Handle<EnvFilter, Registry> {
    FILTER_HANDLE
        .get_or_init(|| {
            // 2) reloadable filter
            let spec = std::env::var("RUST_LOG").unwrap_or_else(|_| "debug".to_string());

            let (filter_layer, handle) = reload::Layer::new(EnvFilter::new(spec));

            // 4) pretty status to stderr, opt-in via `target="status"`
            let status_layer = fmt::layer()
                .with_target(false)
                .pretty()
                .with_writer(std::io::stderr);

            tracing_subscriber::registry()
                .with(filter_layer)
                .with(status_layer)
                .init();

            handle
        })
        .clone()
}

pub fn init_bench_tracing() -> reload::Handle<EnvFilter, Registry> {
    FILTER_HANDLE
        .get_or_init(|| {
            let (filter_layer, handle) = reload::Layer::new(EnvFilter::new("warn"));

            // 4) pretty status to stderr, opt-in via `target="status"`
            let status_layer = fmt::layer()
                .with_target(false)
                .pretty()
                .with_writer(std::io::stderr);

            tracing_subscriber::registry()
                .with(filter_layer)
                .with(status_layer)
                .init();

            handle
        })
        .clone()
}

use tracing::{event, Event, Level, Subscriber};

pub enum Target {
    Lib,
    Status,
}

pub struct TracingLogger {
    level: Level,
    sink: Target, // chooses a literal target at compile time
}

impl TracingLogger {
    pub fn to_lib(level: Level) -> Self {
        Self {
            level,
            sink: Target::Lib,
        }
    }
    pub fn to_status(level: Level) -> Self {
        Self {
            level,
            sink: Target::Status,
        }
    }
}

impl Logger for TracingLogger {
    fn write<T: Serialize>(&self, msg: &str, data: &T) {
        let json = match serde_json::to_string(data) {
            Ok(s) => s,
            Err(e) => format!(r#"{{"serde_error":"{}"}}"#, e),
        };

        // helper to pick the literal target
        macro_rules! log_to {
            ($lvl:ident, $($rest:tt)*) => {
                match self.sink {
                    Target::Lib    => event!(Level::$lvl,        $($rest)*),
                    Target::Status => event!(target: "status", Level::$lvl, $($rest)*),
                }
            }
        }

        match self.level {
            Level::TRACE => log_to!(TRACE, data = %json, "{msg}"),
            Level::DEBUG => log_to!(DEBUG, data = %json, "{msg}"),
            Level::INFO => log_to!(INFO,  data = %json, "{msg}"),
            Level::WARN => log_to!(WARN,  data = %json, "{msg}"),
            Level::ERROR => log_to!(ERROR, data = %json, "{msg}"),
        }
    }
}

use serde_json::Value;
use std::borrow::Cow;

/// Anything you want to show in status + attach as structured JSON.
pub trait StatusRenderable {
    /// Human-friendly, possibly multi-line (e.g. a table).
    fn status_pretty(&self) -> Cow<'_, str>;

    /// Machine-readable payload (stable shape for analysis).
    fn status_json(&self) -> Value;
}

// If you already implement Display + Serialize on a type,
// you get a decent default for free.
impl<T> StatusRenderable for T
where
    T: Serialize + std::fmt::Display,
{
    fn status_pretty(&self) -> Cow<'_, str> {
        Cow::Owned(self.to_string())
    }
    fn status_json(&self) -> Value {
        serde_json::to_value(self).unwrap_or(Value::Null)
    }
}
// io/macros.rs
#[macro_export]
macro_rules! status_event {
    // With structured payload: message expr + data expr
    ($level:ident, $fmt:expr ; data = $data:expr $(,)?) => {{
        let __d = &$data;
        let __json_s = serde_json::to_string(
            &$crate::cli::tracing::StatusRenderable::status_json(__d)
        ).unwrap_or_else(|e| format!(r#"{{"serde_error":"{e}"}}"#));

        // Always emit a structured status event with JSON payload
        ::tracing::event!(target:"status_data",::tracing::Level::$level, data=%__json_s, $fmt);

        // And, if interactive, emit a pretty rendering through the status target
        if $crate::cli::tracing::stderr_is_tty() {
            let __pretty = $crate::cli::tracing::StatusRenderable::status_pretty(__d);
            ::tracing::event!(target:"status", ::tracing::Level::$level, "{}", __pretty);
        }
    }};

    // Plain status event (no payload)
    ($level:ident, $($t:tt)*) => {{
        ::tracing::event!(target:"status", ::tracing::Level::$level, $($t)*);
    }};
}

#[macro_export]
macro_rules! status_info  { ($($x:tt)*) => { $crate::status_event!(INFO,  $($x)*); }; }
#[macro_export]
macro_rules! status_warn  { ($($x:tt)*) => { $crate::status_event!(WARN,  $($x)*); }; }
#[macro_export]
macro_rules! status_debug { ($($x:tt)*) => { $crate::status_event!(DEBUG, $($x)*); }; }

use std::io::IsTerminal;
pub fn stderr_is_tty() -> bool {
    std::io::stderr().is_terminal()
}

/// Collect the event's formatted "message" field.
struct MessageVisitor {
    message: Option<String>,
}
impl tracing::field::Visit for MessageVisitor {
    fn record_debug(&mut self, f: &tracing::field::Field, v: &dyn std::fmt::Debug) {
        if f.name() == "message" {
            self.message = Some(format!("{v:?}"));
        }
    }
}

pub static LOG_FORMAT: LazyLock<Mutex<LogFormat>> = LazyLock::new(|| Mutex::new(LogFormat::Long));

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

        let mut v = MessageVisitor { message: None };
        event.record(&mut v);
        let msg = v.message.as_deref().unwrap_or("");

        match *LOG_FORMAT.lock().unwrap() {
            LogFormat::Long => {
                write!(
                    w,
                    "[{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:03}] @{} {}: {}",
                    now.year(),
                    now.month(),
                    now.day(),
                    now.hour(),
                    now.minute(),
                    now.second(),
                    now.timestamp_subsec_millis(),
                    format_target(meta.module_path().unwrap_or("").to_string(), *meta.level()),
                    format_level(*meta.level()),
                    msg
                )?;
            }
            LogFormat::Short => {
                write!(
                    w,
                    "[{:02}:{:02}:{:02}] {}: {}",
                    now.hour(),
                    now.minute(),
                    now.second(),
                    format_level(*meta.level()),
                    msg
                )?;
            }
            LogFormat::Min => {
                write!(w, "{}: {}", format_level(*meta.level()), msg)?;
            }
            LogFormat::None => {
                write!(w, "{}", msg)?;
            }
        }
        writeln!(w)
    }
}
