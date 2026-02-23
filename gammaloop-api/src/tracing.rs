use std::{
    path::Path,
    sync::{LazyLock, Mutex, OnceLock},
    time::Duration,
};

use chrono::{Datelike, Local, SecondsFormat, Timelike};
use colored::{ColoredString, Colorize};
use eyre::Context;
use gammalooprs::utils::tracing::{LogFormat, LogStyle, LOG_GUARD};
use indicatif::ProgressState;
use tracing::{level_filters::LevelFilter, Event, Subscriber};
use tracing_appender::{
    non_blocking::NonBlockingBuilder,
    rolling::{RollingFileAppender, Rotation},
};
use tracing_indicatif::{filter::IndicatifFilter, style::ProgressStyle, IndicatifLayer};
use tracing_subscriber::field::RecordFields;
use tracing_subscriber::{
    filter::Filtered,
    fmt::{self, format::Writer, FmtContext, FormatEvent, FormatFields},
    layer::SubscriberExt,
    registry::LookupSpan,
    reload,
    util::SubscriberInitExt,
    EnvFilter, Registry,
};
use tracing_subscriber::{fmt::format::FmtSpan, Layer};

use color_eyre::Result;

fn file_filter_from(user_spec: &str) -> Result<EnvFilter> {
    // Start from a strict global default…
    let mut filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::WARN.into()) // global floor
        .parse_lossy(""); // no user rules yet

    // for req in ["gammalooprs=debug", "_gammaloop=info", "symbolica=off"] {
    //     let Ok(d) = req.parse() else {
    //         continue;
    //     };
    //     filter = filter.add_directive(d);
    // }
    for a in user_spec.split(",") {
        if a.trim().is_empty() {
            continue;
        }
        filter = filter.add_directive(a.trim().parse()?);
    }
    Ok(filter)
}

fn display_filter_from(user_spec: &str) -> Result<EnvFilter> {
    // Default: only show `status` target, everything else OFF.
    // Users can override/expand with their spec.
    let mut f = EnvFilter::builder()
        .with_default_directive(LevelFilter::OFF.into())
        .parse_lossy("");

    for a in user_spec.split(",") {
        if a.trim().is_empty() {
            continue;
        }
        f = f.add_directive(
            a.trim()
                .parse()
                .context(format!("Trying to get directive from :{}", a))?,
        );
    }

    Ok(f)
}

const ENV_FILE_LOG_FILTER: &str = "GL_LOGFILE_FILTER";
const ENV_NO_GL_HARD_WARNINGS: &str = "GL_NO_HARD_WARNINGS";
const ENV_DISPLAY_LOG_FILTER: &str = "GL_DISPLAY_FILTER";
const ENV_ALL_LOG_FILTER: &str = "GL_ALL_LOG_FILTER";

// Statics to hold the current log specifications
static FILE_LOG_SPEC: LazyLock<Mutex<String>> = LazyLock::new(|| {
    if let Ok(all) = std::env::var(ENV_ALL_LOG_FILTER) {
        // println!("All:{all}");
        return Mutex::new(all);
    }
    let directive = std::env::var(ENV_FILE_LOG_FILTER)
        .ok()
        .unwrap_or("".to_string());
    Mutex::new(directive)
});
static STDERR_LOG_SPEC: LazyLock<Mutex<String>> = LazyLock::new(|| {
    if let Ok(all) = std::env::var(ENV_ALL_LOG_FILTER) {
        return Mutex::new(all);
    }
    let directive = std::env::var(ENV_DISPLAY_LOG_FILTER)
        .ok()
        .unwrap_or("".to_string());
    Mutex::new(directive)
});

struct FilterHandles {
    file_handle: reload::Handle<EnvFilter, Registry>,
    stderr_handle: reload::Handle<
        EnvFilter,
        tracing_subscriber::layer::Layered<
            Filtered<
                fmt::Layer<
                    Registry,
                    fmt::format::JsonFields,
                    fmt::format::Format<fmt::format::Json>,
                    tracing_appender::non_blocking::NonBlocking,
                >,
                reload::Layer<EnvFilter, Registry>,
                Registry,
            >,
            Registry,
        >,
    >,
}

static FILTER_HANDLES: OnceLock<FilterHandles> = OnceLock::new();

/// Set the file log filter at runtime while preserving required targets.
pub fn set_file_log_filter(user_spec: impl AsRef<str>) -> Result<()> {
    let Some(handles) = FILTER_HANDLES.get() else {
        return Ok(());
    };
    let user_spec = if let Ok(file) = std::env::var(ENV_FILE_LOG_FILTER) {
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
    *FILE_LOG_SPEC.lock().unwrap() = user_spec.to_string();

    let full_spec = file_filter_from(&user_spec)?;
    // println!(
    //     "Modifying file filter to: {} with userspec {}",
    //     full_spec, user_spec
    // );
    handles.file_handle.modify(|f| *f = full_spec)?;
    Ok(())
}

/// Set the stderr log filter at runtime.
pub fn set_stderr_log_filter(user_spec: impl AsRef<str>) -> Result<()> {
    let Some(handles) = FILTER_HANDLES.get() else {
        return Ok(());
    };
    let user_spec = if let Ok(display) = std::env::var(ENV_DISPLAY_LOG_FILTER) {
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
    let full_spec = display_filter_from(&user_spec)?;
    *STDERR_LOG_SPEC.lock().unwrap() = user_spec;
    // println!("Modifying display filter to: {}", full_spec);
    handles.stderr_handle.modify(|f| *f = full_spec)?;
    Ok(())
}

/// Get the current file log filter specification.
pub fn get_file_log_filter() -> String {
    FILE_LOG_SPEC.lock().unwrap().clone()
}

/// Get the current stderr log filter specification.
pub fn get_stderr_log_filter() -> String {
    STDERR_LOG_SPEC.lock().unwrap().clone()
}

fn elapsed_subsec(state: &ProgressState, writer: &mut dyn std::fmt::Write) {
    let seconds = state.elapsed().as_secs();
    let sub_seconds = (state.elapsed().as_millis() % 1000) / 100;
    let _ = writer.write_str(&format!("{}.{}s", seconds, sub_seconds));
}

pub(crate) fn init_tracing(
    dir: impl AsRef<Path>,
    log_file_name: Option<String>,
) -> reload::Handle<EnvFilter, Registry> {
    symbolica::GLOBAL_SETTINGS
        .initialize_tracing
        .store(false, std::sync::atomic::Ordering::Relaxed);
    // println!("Init tracing");
    let handles = FILTER_HANDLES.get_or_init(|| {
        let file_filter = EnvFilter::new(FILE_LOG_SPEC.lock().unwrap().as_str());
        // println!(
        //     "File filter: {file_filter},:{}",
        //     FILE_LOG_SPEC.lock().unwrap().as_str()
        // );
        let stderr_filter = EnvFilter::new(STDERR_LOG_SPEC.lock().unwrap().as_str());
        // println!(
        //     "Stderr filter: {stderr_filter},:{}",
        //     STDERR_LOG_SPEC.lock().unwrap().as_str()
        // );

        let (file_filter_layer, file_handle) = reload::Layer::new(file_filter.clone());
        let (stderr_filter_layer, stderr_handle) = reload::Layer::new(stderr_filter.clone());

        let _ = std::fs::create_dir_all(dir.as_ref());

        // e.g. "2025-08-26T21-05-33.123"
        let ts = Local::now()
            .to_rfc3339_opts(SecondsFormat::Millis, true)
            .replace([':', '+'], "-"); // keep it filename-friendly

        let filename = if let Some(log_name) = log_file_name {
            format!("gammalog-{log_name}.jsonl")
        } else {
            format!("gammalog-{ts}.jsonl")
        };

        // One file per state (per process), no rotation
        let file = RollingFileAppender::new(Rotation::NEVER, dir.as_ref(), &filename);
        let (nb, guard) = NonBlockingBuilder::default()
            .buffered_lines_limit(200_000)
            .lossy(false)
            .finish(file);
        LOG_GUARD.set(guard).ok();

        // JSON layer for file output with its own filter
        let json = fmt::layer()
            .json()
            .flatten_event(true)
            .with_current_span(true)
            .with_writer(nb);

        let indicatif_layer = IndicatifLayer::new()
            .with_span_field_formatter(IndicatifPbMsgFields::default())
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
                              let _ =write!(writer, "\x1b[0m");
                          }
                      },
                  ),
              ).with_span_child_prefix_symbol("↳ ").with_span_child_prefix_indent(" ");
;

        // Pretty status layer for stderr - only show status events
        let status_layer = fmt::layer()
            .with_target(false)
            .event_format(StatusFmt)
            .with_writer(indicatif_layer.get_stderr_writer());

        tracing_subscriber::registry()
            .with(Filtered::new(json, file_filter_layer))
            .with(Filtered::new(status_layer, stderr_filter_layer))
            .with(indicatif_layer.with_filter(IndicatifFilter::new(false)))
            .init();

        FilterHandles {
            file_handle,
            stderr_handle,
        }
    });

    // Return the file handle for backward compatibility
    handles.file_handle.clone()
}

pub static LOG_STYLE: LazyLock<Mutex<LogStyle>> = LazyLock::new(|| Mutex::new(LogStyle::default()));

pub fn set_log_style(style: LogStyle) {
    *LOG_STYLE.lock().unwrap() = style;
}

/// Collect the event's formatted "message" field.
#[derive(Default)]
struct MessageVisitor {
    message: Option<String>,
    fields: Vec<(String, String)>,
}
impl tracing::field::Visit for MessageVisitor {
    fn record_debug(&mut self, f: &tracing::field::Field, v: &dyn std::fmt::Debug) {
        self.record_value(f, format!("{v:?}"));
    }

    fn record_str(&mut self, f: &tracing::field::Field, v: &str) {
        self.record_value(f, v.to_string());
    }

    fn record_bool(&mut self, f: &tracing::field::Field, v: bool) {
        self.record_value(f, v.to_string());
    }

    fn record_i64(&mut self, f: &tracing::field::Field, v: i64) {
        self.record_value(f, v.to_string());
    }

    fn record_u64(&mut self, f: &tracing::field::Field, v: u64) {
        self.record_value(f, v.to_string());
    }

    fn record_f64(&mut self, f: &tracing::field::Field, v: f64) {
        self.record_value(f, v.to_string());
    }

    fn record_error(&mut self, f: &tracing::field::Field, v: &(dyn std::error::Error + 'static)) {
        self.record_value(f, v.to_string());
    }
}
impl MessageVisitor {
    fn record_value(&mut self, f: &tracing::field::Field, value: String) {
        if f.name() == "message" {
            self.message = Some(value);
        } else {
            self.fields.push((f.name().to_string(), value));
        }
    }
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
        let style = LOG_STYLE.lock().unwrap().clone();

        let mut v = MessageVisitor::default();
        event.record(&mut v);
        let msg = v.message.as_deref().unwrap_or("");
        let field_suffix = if style.include_fields && !v.fields.is_empty() {
            let mut s = String::new();
            if !msg.is_empty() {
                s.push(' ');
            }
            s.push('{');
            for (i, (k, v)) in v.fields.iter().enumerate() {
                if i > 0 {
                    s.push_str(", ");
                }
                s.push_str(k);
                s.push('=');
                s.push_str(v);
            }
            s.push('}');
            s
        } else {
            String::new()
        };
        let rendered = format!("{msg}{field_suffix}");

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
        let short_source = if style.full_line_source {
            Some(format_full_source(meta))
        } else {
            None
        };

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
            LogFormat::Short => {
                if let Some(source) = short_source {
                    write!(
                        w,
                        "[{}] @{} {}: {}",
                        ts_short,
                        source,
                        format_level(*meta.level()),
                        rendered
                    )?;
                } else {
                    write!(
                        w,
                        "[{}] {}: {}",
                        ts_short,
                        format_level(*meta.level()),
                        rendered
                    )?;
                }
            }
            LogFormat::Min => {
                if let Some(source) = short_source {
                    write!(
                        w,
                        "@{} {}: {}",
                        source,
                        format_level(*meta.level()),
                        rendered
                    )?;
                } else {
                    write!(w, "{}: {}", format_level(*meta.level()), rendered)?;
                }
            }
            LogFormat::None => {
                if let Some(source) = short_source {
                    write!(w, "@{} {}", source, rendered)?;
                } else {
                    write!(w, "{}", rendered)?;
                }
            }
        }
        writeln!(w)
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
