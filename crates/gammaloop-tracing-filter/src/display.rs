use std::sync::Arc;

use chrono::{Datelike, Local, Timelike};
use colored::{ColoredString, Colorize};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use tabled::{builder::Builder, settings::Style};
use tracing::{Event, Subscriber, field::Visit};
use tracing_subscriber::{
    fmt::{FmtContext, FormatEvent, FormatFields, format::Writer},
    registry::LookupSpan,
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[repr(usize)]
#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Serialize, Deserialize, Default,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
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
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default, deny_unknown_fields)]
pub struct LogStyle {
    #[serde(skip_serializing_if = "is_default")]
    pub log_format: LogFormat,
    #[serde(skip_serializing_if = "is_false")]
    pub short_timestamp: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub full_line_source: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub include_fields: bool,
}

fn is_default<T: Default + PartialEq>(value: &T) -> bool {
    value == &T::default()
}

fn is_false(value: &bool) -> bool {
    !*value
}

#[derive(Clone)]
pub struct GammaDisplayFormat {
    style: Arc<dyn Fn() -> LogStyle + Send + Sync>,
}

impl Default for GammaDisplayFormat {
    fn default() -> Self {
        Self::new(LogStyle::default())
    }
}

impl GammaDisplayFormat {
    pub fn new(style: LogStyle) -> Self {
        Self {
            style: Arc::new(move || style.clone()),
        }
    }

    pub fn full_with_fields() -> Self {
        Self::new(LogStyle {
            log_format: LogFormat::Full,
            include_fields: true,
            ..LogStyle::default()
        })
    }

    pub fn dynamic(style: impl Fn() -> LogStyle + Send + Sync + 'static) -> Self {
        Self {
            style: Arc::new(style),
        }
    }
}

impl<S, N> FormatEvent<S, N> for GammaDisplayFormat
where
    S: Subscriber + for<'a> LookupSpan<'a>,
    N: for<'writer> FormatFields<'writer> + 'static,
{
    fn format_event(
        &self,
        _ctx: &FmtContext<'_, S, N>,
        mut writer: Writer<'_>,
        event: &Event<'_>,
    ) -> std::fmt::Result {
        let now = Local::now();
        let meta = event.metadata();
        let style = (self.style)();

        let mut visitor = RoutedFieldsVisitor::for_sink(LogSink::Display);
        event.record(&mut visitor);
        let msg = visitor.message.as_deref().unwrap_or("");
        let mut display_fields = visitor.fields;
        let display_tags = extract_display_tags(&mut display_fields);
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

        let target = meta.target().to_string();
        let long_source = format!(
            "{}{}",
            format_target(target.clone(), *meta.level()),
            render_display_tags(&display_tags)
        );
        let full_source = format!(
            "{}{}",
            format_target_full(target),
            render_display_tags(&display_tags)
        );
        let line_source = render_line_source_line(&style, meta);

        match style.log_format {
            LogFormat::Long => {
                let timestamp = if style.short_timestamp {
                    &ts_short_ms
                } else {
                    &ts_long
                };
                write!(
                    writer,
                    "[{}] @{} {}: {}",
                    timestamp,
                    long_source,
                    format_level(*meta.level()),
                    msg
                )?;
            }
            LogFormat::Full => {
                let timestamp = if style.short_timestamp {
                    &ts_short_ms
                } else {
                    &ts_long
                };
                write!(
                    writer,
                    "[{}] @{} {}: {}",
                    timestamp,
                    full_source,
                    format_level(*meta.level()),
                    msg
                )?;
            }
            LogFormat::Short => {
                write!(
                    writer,
                    "[{}] {}: {}",
                    ts_short,
                    format_level(*meta.level()),
                    msg
                )?;
            }
            LogFormat::Min => {
                write!(writer, "{}: {}", format_level(*meta.level()), msg)?;
            }
            LogFormat::None => {
                write!(writer, "{msg}")?;
            }
        }

        if let Some(field_table) = rendered_field_table {
            write!(writer, "\n{field_table}")?;
        }
        if let Some(line_source) = line_source {
            write!(writer, "\n{line_source}")?;
        }
        writeln!(writer)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LogSink {
    #[default]
    Display,
    File,
}

const DISPLAY_ONLY_FIELD_PREFIX: &str = "display.";
const FILE_ONLY_FIELD_PREFIX: &str = "file.";

pub fn route_field_name(name: &str, sink: LogSink) -> Option<String> {
    if let Some(stripped) = name.strip_prefix(DISPLAY_ONLY_FIELD_PREFIX) {
        return (sink == LogSink::Display).then(|| stripped.to_string());
    }
    if let Some(stripped) = name.strip_prefix(FILE_ONLY_FIELD_PREFIX) {
        return (sink == LogSink::File).then(|| stripped.to_string());
    }
    Some(name.to_string())
}

#[derive(Default)]
pub struct RoutedFieldsVisitor {
    sink: LogSink,
    pub message: Option<String>,
    pub fields: serde_json::Map<String, serde_json::Value>,
}

impl RoutedFieldsVisitor {
    pub fn for_sink(sink: LogSink) -> Self {
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

fn format_level(level: tracing::Level) -> ColoredString {
    match level {
        tracing::Level::ERROR => format!("{:<8}", "ERROR").red(),
        tracing::Level::WARN => format!("{:<8}", "WARNING").yellow(),
        tracing::Level::INFO => format!("{:<8}", "INFO").into(),
        tracing::Level::DEBUG => format!("{:<8}", "DEBUG").bright_black(),
        tracing::Level::TRACE => format!("{:<8}", "TRACE").into(),
    }
}

fn format_target(target: String, level: tracing::Level) -> ColoredString {
    let split_targets = target.split("::").collect::<Vec<_>>();
    let start = split_targets.len().saturating_sub(2);
    let mut shortened_path = split_targets[start..].join("::");
    if level < tracing::Level::DEBUG && shortened_path.len() > 20 {
        shortened_path = format!("{}...", shortened_path.chars().take(17).collect::<String>());
    }
    format!("{shortened_path:<20}").bright_blue()
}

fn format_target_full(target: String) -> ColoredString {
    target.bright_blue()
}

fn render_line_source_line(style: &LogStyle, meta: &tracing::Metadata<'_>) -> Option<String> {
    if !style.full_line_source {
        return None;
    }
    render_line_source_line_from_parts(meta.file(), meta.line())
}

fn render_line_source_line_from_parts(file: Option<&str>, line: Option<u32>) -> Option<String> {
    format_line_source(file, line)
        .map(|source| format!("    {}", format!("at {source}").bright_black()))
}

fn format_line_source(file: Option<&str>, line: Option<u32>) -> Option<String> {
    let file = file.filter(|file| !file.is_empty())?;
    Some(match line {
        Some(line) => format!("{file}:{line}"),
        None => file.to_string(),
    })
}

pub fn strip_ansi_escape_codes(line: &str) -> String {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_target_format_preserves_long_module_paths() {
        let target = "gammaloop_api::commands::very_long_module_name".to_string();
        assert_eq!(
            strip_ansi_escape_codes(&format_target_full(target.clone()).to_string()),
            target
        );
        assert!(
            format_target(target, tracing::Level::INFO)
                .to_string()
                .contains("...")
        );
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
    fn full_target_with_display_tags_is_copy_pasteable_as_directive() {
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
        assert!(crate::GammaLogFilter::parse(&format!("{directive_head}=debug")).is_ok());
    }

    #[test]
    fn line_source_renders_as_trailing_line() {
        let line = strip_ansi_escape_codes(
            &render_line_source_line_from_parts(Some("src/lib.rs"), Some(42)).unwrap(),
        );

        assert_eq!(line, "    at src/lib.rs:42");
    }

    #[test]
    fn display_field_table_renders_only_non_boolean_fields() {
        let style = LogStyle {
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
