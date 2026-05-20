use bincode_trait_derive::{Decode, Encode};
use clap::ValueEnum;
use gammaloop_tracing_filter::GammaLogFilter;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{path::PathBuf, sync::OnceLock};
use tracing::level_filters::LevelFilter;
use tracing_appender::non_blocking::WorkerGuard;
use tracing_subscriber::{fmt, prelude::*};
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
    Full,
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

// Global one-time slots.
static TEST_TRACING_INITIALISED: OnceLock<()> = OnceLock::new();
pub static LOG_GUARD: OnceLock<WorkerGuard> = OnceLock::new();

const ENV_FILE_LOG_FILTER: &str = "GL_LOGFILE_FILTER";
const ENV_DISPLAY_LOG_FILTER: &str = "GL_DISPLAY_FILTER";
const ENV_ALL_LOG_FILTER: &str = "GL_ALL_LOG_FILTER";
const ENV_TEST_LOG_DIR: &str = "GL_TEST_LOG_DIR";

pub fn init_test_tracing() {
    TEST_TRACING_INITIALISED.get_or_init(|| {
        init_test_tracing_with_defaults("info", "off");
    });
}

pub fn init_bench_tracing() {
    TEST_TRACING_INITIALISED.get_or_init(|| {
        init_test_tracing_with_defaults("warn", "off");
    });
}

fn init_test_tracing_with_defaults(display_default: &str, file_default: &str) {
    let display_spec = log_filter_env_override(ENV_DISPLAY_LOG_FILTER)
        .unwrap_or_else(|| display_default.to_string());
    let file_spec =
        log_filter_env_override(ENV_FILE_LOG_FILTER).unwrap_or_else(|| file_default.to_string());

    let (_display_spec, display_filter) =
        parse_log_filter_or_default(&display_spec, display_default, "display");
    let (file_spec, file_filter) = parse_log_filter_or_default(&file_spec, file_default, "file");

    let display_layer = fmt::layer()
        .pretty()
        .with_writer(std::io::stderr)
        .with_filter(display_filter);

    let subscriber = tracing_subscriber::registry().with(display_layer);
    if GammaLogFilter::is_effectively_off(&file_spec) {
        _ = subscriber.try_init();
    } else {
        let file_appender = tracing_appender::rolling::never(
            test_log_dir(),
            format!("gammaloop-test-{}.jsonl", std::process::id()),
        );
        let (file_writer, file_guard) = tracing_appender::non_blocking(file_appender);
        let _ = LOG_GUARD.set(file_guard);

        let file_layer = fmt::layer()
            .json()
            .with_writer(file_writer)
            .with_filter(file_filter);

        _ = subscriber.with(file_layer).try_init();
    }
}

fn log_filter_env_override(specific_env: &str) -> Option<String> {
    if let Ok(all) = std::env::var(ENV_ALL_LOG_FILTER) {
        Some(all)
    } else {
        std::env::var(specific_env).ok()
    }
}

fn parse_log_filter_or_default(
    spec: &str,
    default_spec: &str,
    sink: &str,
) -> (String, GammaLogFilter) {
    let spec = if spec.trim().is_empty() {
        default_spec
    } else {
        spec
    };

    match GammaLogFilter::parse(spec) {
        Ok(filter) => (spec.to_string(), filter),
        Err(err) => {
            eprintln!(
                "Invalid test {sink} log filter '{spec}': {err}. Falling back to {default_spec}."
            );
            (
                default_spec.to_string(),
                GammaLogFilter::parse(default_spec).expect("test log default filter must parse"),
            )
        }
    }
}

fn test_log_dir() -> PathBuf {
    std::env::var(ENV_TEST_LOG_DIR)
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../..")
                .join("target/test-logs")
        })
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
