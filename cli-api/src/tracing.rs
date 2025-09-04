use std::{
    path::Path,
    sync::{LazyLock, Mutex},
};

use chrono::{Datelike, Local, SecondsFormat, Timelike};
use colored::{ColoredString, Colorize};
use gammalooprs::utils::tracing::{LogFormat, FILTER_HANDLE, LOG_GUARD};
use tracing::{Event, Subscriber};
use tracing_appender::{
    non_blocking::NonBlockingBuilder,
    rolling::{RollingFileAppender, Rotation},
};
use tracing_subscriber::{
    fmt::{self, format::Writer, FmtContext, FormatEvent, FormatFields},
    layer::SubscriberExt,
    registry::LookupSpan,
    reload,
    util::SubscriberInitExt,
    EnvFilter, Registry,
};

use crate::state::LOG_SPEC;

pub(crate) fn init_tracing(
    default_spec: impl AsRef<str>,
    dir: impl AsRef<Path>,
) -> reload::Handle<EnvFilter, Registry> {
    FILTER_HANDLE
        .get_or_init(|| {
            // 2) reloadable filter - only allow gammaloop and status targets
            let spec =
                std::env::var("RUST_LOG").unwrap_or_else(|_| default_spec.as_ref().to_string());
            LOG_SPEC.set(Mutex::new(spec.clone())).ok();

            // Use EnvFilter for reload compatibility, but add per-layer filtering
            let env_filter = EnvFilter::new(&format!(
                "_gammaloop={},status=trace,status_data=trace",
                spec
            ));
            let (filter_layer, handle) = reload::Layer::new(env_filter);

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

            // Add strict filtering to JSON layer
            // let json_filter = filter_fn(|metadata| {
            //     let target = metadata.target();
            //     target.starts_with("_gammaloop") || target == "status" || target == "status_data"
            // });

            let json = fmt::layer()
                .json()
                .flatten_event(true)
                .with_current_span(true)
                .with_span_list(true)
                .with_writer(nb);
            // .with_filter(json_filter);

            // 4) pretty status to stderr, opt-in via `target="status"`
            let status_layer = fmt::layer()
                .with_target(false)
                .event_format(StatusFmt::default())
                .with_writer(std::io::stderr);
            // .with_filter(filter_fn(|m| m.target() == "status"));

            tracing_subscriber::registry()
                .with(filter_layer)
                .with(json)
                .with(status_layer)
                .init();

            handle
        })
        .clone()
}

pub static LOG_FORMAT: LazyLock<Mutex<LogFormat>> = LazyLock::new(|| Mutex::new(LogFormat::Long));

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
