use bincode_trait_derive::{Decode, Encode};
use clap::ValueEnum;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;
use tracing::level_filters::LevelFilter;
use tracing_appender::non_blocking::WorkerGuard;
use tracing_subscriber::filter::filter_fn;
use tracing_subscriber::{EnvFilter, fmt, prelude::*, registry::Registry, reload};
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
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
    Encode,
    Decode,
    JsonSchema,
    ValueEnum,
    Serialize,
    Deserialize,
    Default,
)]
pub enum LogFormat {
    #[default]
    Long,
    Short,
    Min,
    None,
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, JsonSchema, PartialEq, Default)]
#[serde(default, deny_unknown_fields)]
pub struct LogStyle {
    #[serde(skip_serializing_if = "crate::utils::serde_utils::IsDefault::is_default")]
    pub log_format: LogFormat,
    #[serde(skip_serializing_if = "crate::utils::serde_utils::is_false")]
    pub short_timestamp: bool,
    #[serde(skip_serializing_if = "crate::utils::serde_utils::is_false")]
    pub full_line_source: bool,
    #[serde(skip_serializing_if = "crate::utils::serde_utils::is_false")]
    pub include_fields: bool,
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
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Default)]
pub enum LogLevel {
    /// A level lower than all log levels.
    Off,
    /// Corresponds to the `Error` log level.
    Error,
    /// Corresponds to the `Warn` log level.
    Warn,
    /// Corresponds to the `Info` log level.
    #[default]
    Info,
    /// Corresponds to the `Debug` log level.
    Debug,
    /// Corresponds to the `Trace` log level.
    Trace,
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

    pub fn to_cli_display_directive_spec(self) -> &'static str {
        match self {
            LogLevel::Off => "gammaloop_api=off,gammalooprs=off",
            LogLevel::Error => "gammaloop_api=error,gammalooprs=error",
            LogLevel::Warn => "gammaloop_api=warn,gammalooprs=warn",
            LogLevel::Info => "gammaloop_api=info,gammalooprs=info",
            LogLevel::Debug => "gammaloop_api=debug,gammalooprs=debug",
            LogLevel::Trace => "gammaloop_api=trace,gammalooprs=trace",
        }
    }

    pub fn to_cli_logfile_directive_spec(self) -> &'static str {
        self.to_cli_display_directive_spec()
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

            _ = tracing_subscriber::registry()
                .with(filter_layer)
                .with(status_layer)
                .try_init();

            handle
        })
        .clone()
}

pub fn init_bench_tracing() -> reload::Handle<EnvFilter, Registry> {
    FILTER_HANDLE
        .get_or_init(|| {
            // Use EnvFilter for reload compatibility
            let env_filter = EnvFilter::new("gammaloop_api=warn,status=trace,status_data=trace");
            let (filter_layer, handle) = reload::Layer::new(env_filter);

            // 4) pretty status to stderr, opt-in via `target="status"`
            // Add strict filtering to status layer
            let bench_status_filter = filter_fn(|metadata| {
                let target = metadata.target();
                target.starts_with("gammaloop_api") || target == "status" || target == "status_data"
            });

            let status_layer = fmt::layer()
                .with_target(false)
                .pretty()
                .with_writer(std::io::stderr)
                .with_filter(bench_status_filter);

            _ = tracing_subscriber::registry()
                .with(filter_layer)
                .with(status_layer)
                .try_init();

            handle
        })
        .clone()
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

use std::io::IsTerminal;
pub fn stderr_is_tty() -> bool {
    std::io::stderr().is_terminal()
}
