use bincode_trait_derive::{Decode, Encode};
use clap::ValueEnum;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;
use tracing::level_filters::LevelFilter;
use tracing_appender::non_blocking::WorkerGuard;
use tracing_subscriber::filter::filter_fn;
use tracing_subscriber::{fmt, prelude::*, registry::Registry, reload, EnvFilter};

#[repr(usize)]
#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, ValueEnum, Serialize, Deserialize,
)]
pub enum LogFormat {
    Long,
    Short,
    Min,
    None,
}

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
            LogLevel::Off => "gammalooprs=off",
            LogLevel::Error => "gammalooprs=error",
            LogLevel::Warn => "gammalooprs=warn",
            LogLevel::Info => "gammalooprs=info",
            LogLevel::Debug => "gammalooprs=debug",
            LogLevel::Trace => "gammalooprs=trace",
        }
    }
}

// Global one-time slots
pub static FILTER_HANDLE: OnceLock<reload::Handle<EnvFilter, Registry>> = OnceLock::new();
pub static LOG_GUARD: OnceLock<WorkerGuard> = OnceLock::new();

pub fn init_test_tracing() -> reload::Handle<EnvFilter, Registry> {
    FILTER_HANDLE
        .get_or_init(|| {
            // 2) reloadable filter - only allow gammaloop crate logs
            let _spec = std::env::var("RUST_LOG").unwrap_or_else(|_| "debug".to_string());

            // Use EnvFilter for reload compatibility
            let env_filter = EnvFilter::new("gammalooprs=debug,status=trace,status_data=trace");
            let (filter_layer, handle) = reload::Layer::new(env_filter);

            // 4) pretty status to stderr, opt-in via `target="status"`
            // Add strict filtering to status layer
            let test_status_filter = filter_fn(|metadata| {
                let target = metadata.target();
                target.starts_with("gammalooprs") || target == "status" || target == "status_data"
            });

            let status_layer = fmt::layer()
                .with_target(false)
                .pretty()
                .with_writer(std::io::stderr)
                .with_filter(test_status_filter);

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
            // Use EnvFilter for reload compatibility
            let env_filter = EnvFilter::new("_gammaloop=warn,status=trace,status_data=trace");
            let (filter_layer, handle) = reload::Layer::new(env_filter);

            // 4) pretty status to stderr, opt-in via `target="status"`
            // Add strict filtering to status layer
            let bench_status_filter = filter_fn(|metadata| {
                let target = metadata.target();
                target.starts_with("_gammaloop") || target == "status" || target == "status_data"
            });

            let status_layer = fmt::layer()
                .with_target(false)
                .pretty()
                .with_writer(std::io::stderr)
                .with_filter(bench_status_filter);

            tracing_subscriber::registry()
                .with(filter_layer)
                .with(status_layer)
                .init();

            handle
        })
        .clone()
}

use tracing::{event, Level};

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
        // Ask tracing whether these events would be recorded.
        let __data_enabled =
            ::tracing::event_enabled!(target: "status_data", ::tracing::Level::$level);
        let __pretty_enabled =
            $crate::utils::tracing::stderr_is_tty()
            && ::tracing::event_enabled!(target: "status", ::tracing::Level::$level);

        // Fast path: nothing to do.
        if __data_enabled || __pretty_enabled {
            // Evaluate `$data` exactly once if anything is enabled.
            let __d_val = $data;

            if __data_enabled {
                // Build JSON only when the structured event is enabled.
                let __json_s = serde_json::to_string(
                    &$crate::utils::tracing::StatusRenderable::status_json(&__d_val)
                ).unwrap_or_else(|e| format!(r#"{{"serde_error":"{e}"}}"#));

                ::tracing::event!(
                    target: "status_data",
                    ::tracing::Level::$level,
                    data = %__json_s,
                    $fmt
                );
            }

            if __pretty_enabled {
                // Pretty render only when an interactive TTY AND the event is enabled.
                let __pretty =
                    $crate::utils::tracing::StatusRenderable::status_pretty(&__d_val);
                ::tracing::event!(
                    target: "status",
                    ::tracing::Level::$level,
                    "{}", __pretty
                );
            }
        }
    }};

    ($level:ident, $($t:tt)*) => {{
        ::tracing::event!(target:"status", ::tracing::Level::$level, $($t)*);
    }};
}
#[macro_export]
macro_rules! status_info  { ($($x:tt)*) => { $crate::status_event!(INFO,  $($x)*); }; }
#[macro_export]
macro_rules! status_warn  { ($($x:tt)*) => { $crate::status_event!(WARN,  $($x)*); }; }
#[macro_export]
macro_rules! status_error  { ($($x:tt)*) => { $crate::status_event!(ERROR,  $($x)*); }; }
#[macro_export]
macro_rules! status_debug { ($($x:tt)*) => { $crate::status_event!(DEBUG, $($x)*); }; }
#[macro_export]
macro_rules! status_trace { ($($x:tt)*) => { $crate::status_event!(TRACE, $($x)*); }; }

use std::io::IsTerminal;
pub fn stderr_is_tty() -> bool {
    std::io::stderr().is_terminal()
}
