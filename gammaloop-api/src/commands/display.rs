use clap::Subcommand;
use colored::Colorize;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value as JsonValue;
use tabled::{builder::Builder, settings::Style};
use tracing::info;

use color_eyre::Result;
use eyre::{Context, eyre};
use gammalooprs::processes::{Amplitude, CrossSection, Process, ProcessCollection};
use gammalooprs::settings::RuntimeSettings;
use gammalooprs::utils::serde_utils::SHOWDEFAULTS;
use std::sync::atomic::Ordering;

use crate::{
    CLISettings,
    commands::generate::ProcessArgs,
    state::{ProcessRef, State},
};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Display {
    Model {
        #[arg(short = 'c', long = "show-couplings", default_value_t = false)]
        show_couplings: bool,
        #[arg(short = 'v', long = "show-vertices", default_value_t = false)]
        show_vertices: bool,
        #[arg(short = 'r', long = "show-parameters", default_value_t = false)]
        show_parameters: bool,
        #[arg(short = 'p', long = "show-particles", default_value_t = false)]
        show_particles: bool,
        #[arg(short = 'a', long = "show-all", default_value_t = false)]
        show_all: bool,
    },
    Processes,
    Integrands {
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(short = 'p', long = "process", value_name = "PROCESS")]
        process: Option<ProcessRef>,
    },
    Settings {
        #[command(subcommand)]
        target: DisplaySettingsTarget,
    },
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum DisplaySettingsTarget {
    Global {
        /// Optional dotted path to a specific setting
        #[arg(value_name = "PATH")]
        path: Option<String>,
    },
    #[command(alias = "defaults")]
    DefaultRuntime {
        /// Optional dotted path to a specific setting
        #[arg(value_name = "PATH")]
        path: Option<String>,
    },
    Process {
        #[command(flatten)]
        process: ProcessArgs,
        /// Optional dotted path to a specific setting
        #[arg(value_name = "PATH")]
        path: Option<String>,
    },
}

impl Display {
    pub fn run(
        &self,
        state: &State,
        global_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        match self {
            Display::Integrands { process } => {
                let process_id = process
                    .as_ref()
                    .map(|process_ref| state.resolve_process_ref(Some(process_ref)))
                    .transpose()?;
                render_integrands_table(state, process_id)?;
            }
            Display::Processes => {
                render_processes_table(state)?;
            }
            Display::Model {
                show_couplings,
                show_vertices,
                show_parameters,
                show_particles,
                show_all,
            } => {
                info!(
                    "\n{}",
                    state.model.get_description(
                        *show_particles || *show_all,
                        *show_parameters || *show_all,
                        *show_vertices || *show_all,
                        *show_couplings || *show_all,
                    )
                )
            }
            Display::Settings { target } => {
                target.run(state, global_settings, default_runtime_settings)?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct IntegrandMetrics {
    name: String,
    graphs: usize,
    graph_groups: usize,
    record_size_bytes: usize,
}

fn render_processes_table(state: &State) -> Result<()> {
    if state.process_list.processes.is_empty() {
        info!("{}", "No processes generated yet.".yellow());
        return Ok(());
    }

    let mut builder = Builder::new();
    builder.push_record([
        "process #".bold().blue().to_string(),
        "process".bold().blue().to_string(),
        "# integrands".bold().blue().to_string(),
        "integrand names".bold().blue().to_string(),
        "graphs / integrand".bold().blue().to_string(),
        "graph groups / integrand".bold().blue().to_string(),
        ".bin size / integrand".bold().blue().to_string(),
    ]);

    for process in &state.process_list.processes {
        let metrics = collect_integrand_metrics(process)?;
        builder.push_record([
            format!("#{}", process.definition.process_id)
                .blue()
                .to_string(),
            process
                .definition
                .folder_name
                .clone()
                .green()
                .bold()
                .to_string(),
            metrics.len().to_string().yellow().to_string(),
            format_integrand_names(&metrics),
            format_metrics_column(&metrics, |metric| metric.graphs.to_string()),
            format_metrics_column(&metrics, |metric| metric.graph_groups.to_string()),
            format_metrics_column(&metrics, |metric| format_bytes(metric.record_size_bytes)),
        ]);
    }

    let mut table = builder.build();
    table.with(Style::rounded());
    info!("{}", "Processes".bold().blue());
    info!("\n{table}");
    Ok(())
}

fn render_integrands_table(state: &State, process_filter: Option<usize>) -> Result<()> {
    if state.process_list.processes.is_empty() {
        info!("{}", "No processes generated yet.".yellow());
        return Ok(());
    }

    let mut builder = Builder::new();
    builder.push_record([
        "process #".bold().blue().to_string(),
        "process".bold().blue().to_string(),
        "integrand".bold().blue().to_string(),
        "graphs".bold().blue().to_string(),
        "graph groups".bold().blue().to_string(),
        ".bin size".bold().blue().to_string(),
    ]);

    let mut num_rows = 0usize;
    for (process_index, process) in state.process_list.processes.iter().enumerate() {
        if process_filter.is_some_and(|selected| selected != process_index) {
            continue;
        }
        let metrics = collect_integrand_metrics(process)?;
        for metric in metrics {
            num_rows += 1;
            builder.push_record([
                format!("#{}", process.definition.process_id)
                    .blue()
                    .to_string(),
                process
                    .definition
                    .folder_name
                    .clone()
                    .green()
                    .bold()
                    .to_string(),
                metric.name.cyan().to_string(),
                metric.graphs.to_string().yellow().to_string(),
                metric.graph_groups.to_string().yellow().to_string(),
                format_bytes(metric.record_size_bytes).magenta().to_string(),
            ]);
        }
    }

    if num_rows == 0 {
        if let Some(process_id) = process_filter {
            let process = &state.process_list.processes[process_id];
            info!(
                "{}",
                format!(
                    "No integrands found for process #{} ({})",
                    process.definition.process_id, process.definition.folder_name
                )
                .yellow()
            );
        } else {
            info!("{}", "No integrands found.".yellow());
        }
        return Ok(());
    }

    let mut table = builder.build();
    table.with(Style::rounded());

    if let Some(process_id) = process_filter {
        let process = &state.process_list.processes[process_id];
        info!(
            "{}",
            format!(
                "Integrands for process #{} ({})",
                process.definition.process_id, process.definition.folder_name
            )
            .bold()
            .blue()
        );
    } else {
        info!("{}", "Integrands".bold().blue());
    }
    info!("\n{table}");
    Ok(())
}

fn collect_integrand_metrics(process: &Process) -> Result<Vec<IntegrandMetrics>> {
    match &process.collection {
        ProcessCollection::Amplitudes(amplitudes) => amplitudes
            .values()
            .map(integrand_metrics_from_amplitude)
            .collect(),
        ProcessCollection::CrossSections(cross_sections) => cross_sections
            .values()
            .map(integrand_metrics_from_cross_section)
            .collect(),
    }
}

fn integrand_metrics_from_amplitude(amplitude: &Amplitude) -> Result<IntegrandMetrics> {
    Ok(IntegrandMetrics {
        name: amplitude.name.clone(),
        graphs: amplitude.graphs.len(),
        graph_groups: amplitude.graph_group_structure.len(),
        record_size_bytes: integrand_record_size_from_amplitude(amplitude)?,
    })
}

fn integrand_metrics_from_cross_section(cross_section: &CrossSection) -> Result<IntegrandMetrics> {
    Ok(IntegrandMetrics {
        name: cross_section.name.clone(),
        graphs: cross_section.supergraphs.len(),
        graph_groups: cross_section.graph_group_structure.len(),
        record_size_bytes: integrand_record_size_from_cross_section(cross_section)?,
    })
}

fn integrand_record_size_from_amplitude(amplitude: &Amplitude) -> Result<usize> {
    let mut record = amplitude.clone();
    record.integrand = None;
    let encoded =
        bincode::encode_to_vec(&record, bincode::config::standard()).with_context(|| {
            format!(
                "While serializing amplitude '{}' for display",
                amplitude.name
            )
        })?;
    Ok(encoded.len())
}

fn integrand_record_size_from_cross_section(cross_section: &CrossSection) -> Result<usize> {
    let mut record = cross_section.clone();
    record.integrand = None;
    let encoded =
        bincode::encode_to_vec(&record, bincode::config::standard()).with_context(|| {
            format!(
                "While serializing cross-section '{}' for display",
                cross_section.name
            )
        })?;
    Ok(encoded.len())
}

fn format_integrand_names(metrics: &[IntegrandMetrics]) -> String {
    if metrics.is_empty() {
        return "(none)".dimmed().to_string();
    }
    metrics
        .iter()
        .map(|metric| metric.name.clone().cyan().to_string())
        .collect::<Vec<_>>()
        .join("\n")
}

fn format_metrics_column(
    metrics: &[IntegrandMetrics],
    format_value: impl Fn(&IntegrandMetrics) -> String,
) -> String {
    if metrics.is_empty() {
        return "(none)".dimmed().to_string();
    }
    metrics
        .iter()
        .map(|metric| format!("{}: {}", metric.name.cyan(), format_value(metric)))
        .collect::<Vec<_>>()
        .join("\n")
}

fn format_bytes(size: usize) -> String {
    const KIB: usize = 1024;
    const MIB: usize = KIB * 1024;
    const GIB: usize = MIB * 1024;

    if size < KIB {
        format!("{size} B")
    } else if size < MIB {
        format!("{:.2} KiB", size as f64 / KIB as f64)
    } else if size < GIB {
        format!("{:.2} MiB", size as f64 / MIB as f64)
    } else {
        format!("{:.2} GiB", size as f64 / GIB as f64)
    }
}

impl DisplaySettingsTarget {
    fn run(
        &self,
        state: &State,
        global_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        match self {
            DisplaySettingsTarget::Global { path } => {
                let root = serialize_settings_with_defaults(
                    global_settings,
                    "global settings for display",
                )?;
                render_settings_table("global settings", &root, path.as_deref())?;
            }
            DisplaySettingsTarget::DefaultRuntime { path } => {
                let root = serialize_settings_with_defaults(
                    default_runtime_settings,
                    "default runtime settings for display",
                )?;
                render_settings_table("default runtime settings", &root, path.as_deref())?;
            }
            DisplaySettingsTarget::Process { process, path } => {
                let process_id = state.resolve_process_ref(process.process.as_ref())?;
                let process_ref = &state.process_list.processes[process_id];
                if let Some(name) = &process.integrand_name {
                    let integrand = state.process_list.get_integrand(process_id, name)?;
                    let settings = serialize_settings_with_defaults(
                        integrand.get_settings(),
                        "process runtime settings for display",
                    )?;
                    render_settings_table(
                        &format!(
                            "process settings for #{} ({}) integrand '{}'",
                            process_ref.definition.process_id,
                            process_ref.definition.folder_name,
                            name
                        ),
                        &settings,
                        path.as_deref(),
                    )?;
                } else {
                    for integrand_name in process_ref.get_integrand_names() {
                        let integrand = state
                            .process_list
                            .get_integrand(process_id, integrand_name)?;
                        let settings = serialize_settings_with_defaults(
                            integrand.get_settings(),
                            "process runtime settings for display",
                        )?;
                        render_settings_table(
                            &format!(
                                "process settings for #{} ({}) integrand '{}'",
                                process_ref.definition.process_id,
                                process_ref.definition.folder_name,
                                integrand_name
                            ),
                            &settings,
                            path.as_deref(),
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}

fn render_settings_table(title: &str, root: &JsonValue, key_path: Option<&str>) -> Result<()> {
    let key_path = key_path.map(str::trim).filter(|path| !path.is_empty());
    let selected = if let Some(path) = key_path {
        value_at_path(root, path)?
    } else {
        root
    };
    let mut builder = Builder::new();
    builder.push_record(["key", "value"]);

    match selected {
        JsonValue::Object(map) => {
            let mut keys: Vec<_> = map.keys().cloned().collect();
            keys.sort();
            for key in keys {
                let value = map
                    .get(&key)
                    .expect("sorted keys were derived from the same map");
                builder.push_record([key, format_json_value(value)]);
            }
        }
        JsonValue::Array(array) => {
            if array.is_empty() {
                builder.push_record(["(empty)".to_string(), "[]".to_string()]);
            } else {
                for (index, value) in array.iter().enumerate() {
                    builder.push_record([index.to_string(), format_json_value(value)]);
                }
            }
        }
        primitive => {
            builder.push_record(["value".to_string(), format_json_value(primitive)]);
        }
    }

    let mut table = builder.build();
    table.with(Style::rounded());
    if let Some(path) = key_path {
        info!("{title} (path: {path}):");
    } else {
        info!("{title}:");
    }
    info!("\n{table}");
    Ok(())
}

fn serialize_settings_with_defaults<T: Serialize>(
    settings: &T,
    context: &str,
) -> Result<JsonValue> {
    let _show_defaults_guard = ShowDefaultsGuard::new(true);
    serde_json::to_value(settings).context(format!("While serializing {context}"))
}

struct ShowDefaultsGuard {
    previous: bool,
}

impl ShowDefaultsGuard {
    fn new(show_defaults: bool) -> Self {
        let previous = SHOWDEFAULTS.swap(show_defaults, Ordering::Relaxed);
        Self { previous }
    }
}

impl Drop for ShowDefaultsGuard {
    fn drop(&mut self) {
        SHOWDEFAULTS.store(self.previous, Ordering::Relaxed);
    }
}

fn value_at_path<'a>(root: &'a JsonValue, path: &str) -> Result<&'a JsonValue> {
    let mut current = root;
    for segment in path.split('.') {
        if segment.is_empty() {
            continue;
        }
        current = match current {
            JsonValue::Object(map) => map.get(segment).ok_or_else(|| {
                let mut available: Vec<_> = map.keys().cloned().collect();
                available.sort();
                eyre!(
                    "Key segment '{}' not found while resolving '{}'. Available keys: {}",
                    segment,
                    path,
                    available.join(", ")
                )
            })?,
            JsonValue::Array(items) => {
                let index = segment.parse::<usize>().map_err(|_| {
                    eyre!(
                        "Array segment '{}' is not a valid index while resolving '{}'",
                        segment,
                        path
                    )
                })?;
                items.get(index).ok_or_else(|| {
                    eyre!(
                        "Array index {} out of bounds while resolving '{}', len={}",
                        index,
                        path,
                        items.len()
                    )
                })?
            }
            _ => {
                return Err(eyre!(
                    "Cannot descend into segment '{}' while resolving '{}': current value is {}",
                    segment,
                    path,
                    json_type_name(current)
                ));
            }
        };
    }
    Ok(current)
}

fn json_type_name(value: &JsonValue) -> &'static str {
    match value {
        JsonValue::Null => "null",
        JsonValue::Bool(_) => "boolean",
        JsonValue::Number(_) => "number",
        JsonValue::String(_) => "string",
        JsonValue::Array(_) => "array",
        JsonValue::Object(_) => "object",
    }
}

fn format_json_value(value: &JsonValue) -> String {
    match value {
        JsonValue::Null => "null".to_string(),
        JsonValue::Bool(_) | JsonValue::Number(_) => value.to_string(),
        JsonValue::String(s) => s.clone(),
        _ => serde_json::to_string_pretty(value).unwrap_or_else(|_| value.to_string()),
    }
}

#[cfg(test)]
mod test {
    use clap::Parser;
    use serde_json::json;

    use crate::{CLISettings, Repl, commands::Commands, state::ProcessRef};
    use gammalooprs::settings::RuntimeSettings;

    use super::{
        Display, DisplaySettingsTarget, format_bytes, serialize_settings_with_defaults,
        value_at_path,
    };

    #[test]
    fn parse_display_settings_process_with_path() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "settings",
            "process",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
            "integrator.n_max",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Settings { target }) => match target {
                DisplaySettingsTarget::Process { process, path } => {
                    assert_eq!(
                        process.process,
                        Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                    );
                    assert_eq!(process.integrand_name, Some("LO".to_string()));
                    assert_eq!(path.as_deref(), Some("integrator.n_max"));
                }
                other => panic!("Expected process target, got {other:?}"),
            },
            other => panic!("Expected display settings command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_settings_process_without_path() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "settings",
            "process",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Settings { target }) => match target {
                DisplaySettingsTarget::Process { path, .. } => {
                    assert_eq!(path, None);
                }
                other => panic!("Expected process target, got {other:?}"),
            },
            other => panic!("Expected display settings command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_settings_global_without_query() {
        let repl = Repl::try_parse_from(["gammaloop", "display", "settings", "global"]).unwrap();

        match repl.command {
            Commands::Display(Display::Settings { target }) => match target {
                DisplaySettingsTarget::Global { path } => {
                    assert_eq!(path, None);
                }
                other => panic!("Expected global target, got {other:?}"),
            },
            other => panic!("Expected display settings command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_settings_defaults_alias() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "settings",
            "defaults",
            "integrator.n_start",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Settings { target }) => match target {
                DisplaySettingsTarget::DefaultRuntime { path } => {
                    assert_eq!(path.as_deref(), Some("integrator.n_start"));
                }
                other => panic!("Expected default runtime target, got {other:?}"),
            },
            other => panic!("Expected display settings command, got {other:?}"),
        }
    }

    #[test]
    fn value_lookup_supports_nested_paths() {
        let value = json!({
            "integrator": {
                "n_max": 42
            }
        });
        let found = value_at_path(&value, "integrator.n_max").unwrap();
        assert_eq!(found, &json!(42));
    }

    #[test]
    fn serialized_settings_include_default_values() {
        let global = CLISettings::default();
        let global_json = serialize_settings_with_defaults(&global, "global settings").unwrap();
        assert!(value_at_path(&global_json, "global.n_cores.generate").is_ok());

        let runtime = RuntimeSettings::default();
        let runtime_json =
            serialize_settings_with_defaults(&runtime, "default runtime settings").unwrap();
        assert!(value_at_path(&runtime_json, "integrator.n_start").is_ok());
    }

    #[test]
    fn parse_display_integrands_without_process() {
        let repl = Repl::try_parse_from(["gammaloop", "display", "integrands"]).unwrap();

        match repl.command {
            Commands::Display(Display::Integrands { process }) => {
                assert_eq!(process, None);
            }
            other => panic!("Expected display integrands command, got {other:?}"),
        }
    }

    #[test]
    fn format_bytes_uses_binary_units() {
        assert_eq!(format_bytes(999), "999 B");
        assert_eq!(format_bytes(1024), "1.00 KiB");
        assert_eq!(format_bytes(1024 * 1024), "1.00 MiB");
    }
}
