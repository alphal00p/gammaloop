use clap::{Args, Subcommand, ValueEnum};
use colored::Colorize;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value as JsonValue;
use std::{
    ffi::OsStr,
    fs,
    path::{Path, PathBuf},
};
use tabled::{
    builder::Builder,
    settings::{style::HorizontalLine, themes::Theme, Style},
    Tabled,
};
use tracing::info;
use walkdir::WalkDir;

use color_eyre::Result;
use eyre::{eyre, Context};
use gammalooprs::{
    processes::{Amplitude, CrossSection, Process, ProcessCollection},
    settings::RuntimeSettings,
};

use crate::{
    commands::generate::{render_generation_summary, ProcessArgs},
    commands::process_settings::{
        observable_kind, quantity_kind, selector_kind, serialize_runtime_named_settings,
        summarize_observable, summarize_quantity, summarize_selector, NamedProcessSettingKind,
    },
    completion::CompletionArgExt,
    integrand_info::{IntegrandGraphGroupInfo, IntegrandInfo, IntegrandKind},
    session::display_command,
    settings_tree::{serialize_settings_with_defaults, value_at_path},
    state::{CommandsBlock, IntegrandGenerationSummary, ProcessRef, RunHistory, State},
    CLISettings,
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
        #[arg(long = "show-particles", default_value_t = false)]
        show_particles: bool,
        #[arg(short = 'a', long = "show-all", default_value_t = false)]
        show_all: bool,
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            requires = "integrand_name",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: Option<ProcessRef>,
        /// Integrand name inside the selected process
        #[arg(
            short = 'i',
            long = "integrand-name",
            value_name = "NAME",
            requires = "process",
            completion_integrand_selector(crate::completion::SelectorKind::Any)
        )]
        integrand_name: Option<String>,
    },
    Processes,
    #[command(name = "integrand")]
    Integrands {
        /// Process reference: #<id>, name:<name>, or <id>/<name>
        #[arg(
            short = 'p',
            long = "process",
            value_name = "PROCESS",
            completion_process_selector(crate::completion::SelectorKind::Any)
        )]
        process: Option<ProcessRef>,
        /// Integrand name inside the selected process
        #[arg(
            short = 'i',
            long = "integrand-name",
            value_name = "NAME",
            requires = "process",
            completion_integrand_selector(crate::completion::SelectorKind::Any)
        )]
        integrand_name: Option<String>,
        /// Restrict the detailed view to the selected master graph names
        #[arg(
            short = 'g',
            long = "graph",
            value_name = "MASTER_GRAPH",
            num_args = 1..,
            requires = "integrand_name",
            completion_selected_master_graph()
        )]
        graphs: Vec<String>,
        /// Restrict the detailed view to the selected category tables
        #[arg(
            long = "category",
            value_name = "CATEGORY",
            num_args = 1..,
            requires = "integrand_name",
            completion_selected_integrand_category()
        )]
        categories: Vec<IntegrandDisplayCategory>,
    },
    Quantities {
        #[command(flatten)]
        target: DisplayProcessNamedSettingsArgs,
    },
    Observables {
        #[command(flatten)]
        target: DisplayProcessNamedSettingsArgs,
    },
    Selectors {
        #[command(flatten)]
        target: DisplayProcessNamedSettingsArgs,
    },
    #[command(name = "command_block")]
    CommandBlock {
        #[arg(value_name = "NAME")]
        name: Option<String>,
    },
    Settings {
        #[command(subcommand)]
        target: DisplaySettingsTarget,
    },
}

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct DisplayProcessNamedSettingsArgs {
    #[command(flatten)]
    process: ProcessArgs,
    #[arg(value_name = "NAME")]
    name: Option<String>,
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
        run_history: &RunHistory,
    ) -> Result<()> {
        match self {
            Display::Integrands {
                process,
                integrand_name,
                graphs,
                categories,
            } => {
                let process_id = process
                    .as_ref()
                    .map(|process_ref| state.resolve_process_ref(Some(process_ref)))
                    .transpose()?;
                if let Some(integrand_name) = integrand_name.as_deref() {
                    render_integrand_detail(
                        state,
                        &global_settings.state.folder,
                        process_id.expect("clap requires --process when --integrand-name is set"),
                        integrand_name,
                        graphs,
                        categories,
                    )?;
                } else {
                    render_integrands_table(state, &global_settings.state.folder, process_id)?;
                }
            }
            Display::Processes => {
                render_processes_table(state, &global_settings.state.folder)?;
            }
            Display::Quantities { target } => {
                render_named_process_settings(state, target, NamedProcessSettingKind::Quantity)?;
            }
            Display::Observables { target } => {
                render_named_process_settings(state, target, NamedProcessSettingKind::Observable)?;
            }
            Display::Selectors { target } => {
                render_named_process_settings(state, target, NamedProcessSettingKind::Selector)?;
            }
            Display::CommandBlock { name } => {
                render_command_blocks(run_history, name.as_deref())?;
            }
            Display::Model {
                show_couplings,
                show_vertices,
                show_parameters,
                show_particles,
                show_all,
                process,
                integrand_name,
            } => {
                let model = if let (Some(process), Some(integrand_name)) =
                    (process.as_ref(), integrand_name.as_ref())
                {
                    let process_id = state.resolve_process_ref(Some(process))?;
                    state.resolve_model_for_integrand(process_id, integrand_name)?
                } else {
                    state.model.clone()
                };
                info!(
                    "\n{}",
                    model.get_description(
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
    bin_disk_size_bytes: Option<u64>,
    so_disk_size_bytes: Option<u64>,
}

#[derive(Clone, Debug, Default)]
struct IntegrandArtifactSizes {
    bin_disk_size_bytes: Option<u64>,
    so_disk_size_bytes: Option<u64>,
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, ValueEnum, JsonSchema, PartialEq, Eq, Hash,
)]
#[serde(rename_all = "snake_case")]
pub enum IntegrandDisplayCategory {
    #[value(name = "generation")]
    Generation,
    #[value(name = "orientation")]
    Orientation,
    #[value(name = "loop_momentum_basis")]
    LoopMomentumBasis,
    #[value(name = "cuts")]
    Cuts,
}

impl IntegrandDisplayCategory {
    fn title(self) -> &'static str {
        match self {
            Self::Generation => "Generation",
            Self::Orientation => "Orientations",
            Self::LoopMomentumBasis => "Loop momentum bases",
            Self::Cuts => "Cuts",
        }
    }

    fn value_header(self) -> &'static str {
        match self {
            Self::Generation => "generation",
            Self::Orientation => "orientation",
            Self::LoopMomentumBasis => "loop momentum basis",
            Self::Cuts => "cut",
        }
    }
}

#[derive(Tabled)]
struct DetailSummaryRow {
    field: String,
    value: String,
}

fn orientation_column_widths(edge_ids: &[usize], fallback_count: usize) -> Vec<usize> {
    if edge_ids.is_empty() {
        vec![1; fallback_count]
    } else {
        edge_ids
            .iter()
            .map(|edge_id| edge_id.to_string().len().max(1))
            .collect()
    }
}

fn render_orientation_edge_ids(edge_ids: &[usize]) -> String {
    if edge_ids.is_empty() {
        return String::new();
    }

    let widths = orientation_column_widths(edge_ids, edge_ids.len());
    edge_ids
        .iter()
        .zip(widths.iter())
        .map(|(edge_id, width)| {
            format!("{edge_id:<width$}", width = *width)
                .magenta()
                .to_string()
        })
        .join(" ")
}

fn render_orientation_signature(signature: &[i8], edge_ids: &[usize]) -> String {
    let widths = orientation_column_widths(edge_ids, signature.len());
    signature
        .iter()
        .zip(widths.iter())
        .map(|(sign, width)| {
            let padded = format!(
                "{:<width$}",
                match sign {
                    1 => "+",
                    -1 => "-",
                    0 => "0",
                    _ => "?",
                },
                width = *width
            );
            match sign {
                1 => padded.green().bold().to_string(),
                -1 => padded.red().bold().to_string(),
                0 => padded.dimmed().to_string(),
                _ => padded.to_string(),
            }
        })
        .join(" ")
}

fn render_edge_ids_with_highlight(edge_ids: &[usize], highlight: bool) -> String {
    if edge_ids.is_empty() {
        return String::new();
    }

    let rendered = format!(
        "({})",
        edge_ids.iter().map(|edge_id| edge_id.to_string()).join(",")
    );
    if highlight {
        rendered.green().bold().to_string()
    } else {
        rendered.magenta().to_string()
    }
}

fn render_edge_ids(edge_ids: &[usize]) -> String {
    render_edge_ids_with_highlight(edge_ids, false)
}

fn render_threshold_esurface_ids(esurface_ids: &[usize]) -> String {
    if esurface_ids.is_empty() {
        return "none".dimmed().to_string();
    }
    format!(
        "[{}]",
        esurface_ids
            .iter()
            .map(|esurface_id| esurface_id.to_string())
            .join(",")
    )
    .cyan()
    .to_string()
}

fn render_threshold_esurface(edge_ids: &[usize]) -> String {
    render_edge_ids(edge_ids)
}

fn render_cut(edge_ids: &[usize], raising_power: usize) -> String {
    let edge_ids = render_edge_ids(edge_ids);
    if edge_ids.is_empty() || raising_power <= 1 {
        edge_ids
    } else {
        format!("{edge_ids}^{raising_power}").magenta().to_string()
    }
}

#[derive(Clone)]
struct CategoryRow {
    primary_id: String,
    secondary_id: String,
    tertiary_id: String,
    value: String,
}

fn render_category_id(id: usize, highlight: bool) -> String {
    let rendered = format!("#{id}");
    if highlight {
        rendered.green().bold().to_string()
    } else {
        rendered.yellow().to_string()
    }
}

fn render_optional_category_id(id: Option<usize>, highlight: bool) -> String {
    id.map(|id| render_category_id(id, highlight))
        .unwrap_or_default()
}

fn master_graph(
    group: &IntegrandGraphGroupInfo,
) -> Result<&crate::integrand_info::IntegrandGraphInfo> {
    group
        .graphs
        .iter()
        .find(|graph| graph.is_master)
        .ok_or_else(|| eyre!("Graph group #{} has no master graph", group.group_id))
}

fn available_integrand_categories(detail: &IntegrandInfo) -> Vec<IntegrandDisplayCategory> {
    match detail.kind {
        IntegrandKind::Amplitude => vec![
            IntegrandDisplayCategory::Generation,
            IntegrandDisplayCategory::Orientation,
            IntegrandDisplayCategory::LoopMomentumBasis,
            IntegrandDisplayCategory::Cuts,
        ],
        IntegrandKind::CrossSection => vec![
            IntegrandDisplayCategory::Generation,
            IntegrandDisplayCategory::Orientation,
            IntegrandDisplayCategory::LoopMomentumBasis,
            IntegrandDisplayCategory::Cuts,
        ],
    }
}

fn selected_integrand_categories(
    detail: &IntegrandInfo,
    requested: &[IntegrandDisplayCategory],
) -> Result<Vec<IntegrandDisplayCategory>> {
    let available = available_integrand_categories(detail);
    if requested.is_empty() {
        return Ok(available);
    }

    let mut selected = Vec::new();
    for category in requested.iter().copied() {
        if !available.contains(&category) {
            return Err(eyre!(
                "Category '{}' is not available for {} integrands.",
                category.to_possible_value().expect("value enum").get_name(),
                detail.kind
            ));
        }
        if !selected.contains(&category) {
            selected.push(category);
        }
    }
    Ok(selected)
}

fn filtered_graph_groups<'a>(
    detail: &'a IntegrandInfo,
    requested_graphs: &[String],
) -> Result<Vec<&'a IntegrandGraphGroupInfo>> {
    if requested_graphs.is_empty() {
        return Ok(detail.graph_groups.iter().collect());
    }

    let available_master_names = detail
        .graph_groups
        .iter()
        .map(master_graph)
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .map(|graph| graph.name.clone())
        .collect::<Vec<_>>();
    let missing = requested_graphs
        .iter()
        .filter(|graph_name| {
            !available_master_names
                .iter()
                .any(|name| name == *graph_name)
        })
        .cloned()
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        return Err(eyre!(
            "Unknown master graph filter(s): {}. Available master graphs: {}",
            missing.join(", "),
            available_master_names.join(", ")
        ));
    }

    Ok(detail
        .graph_groups
        .iter()
        .filter(|group| {
            master_graph(group)
                .map(|graph| requested_graphs.iter().any(|name| name == &graph.name))
                .unwrap_or(false)
        })
        .collect())
}

fn render_integrand_detail(
    state: &State,
    state_folder: &Path,
    process_id: usize,
    integrand_name: &str,
    requested_graphs: &[String],
    requested_categories: &[IntegrandDisplayCategory],
) -> Result<()> {
    let detail = state.get_integrand_info(
        Some(&ProcessRef::Id(process_id)),
        Some(&integrand_name.to_string()),
    )?;
    let process = &state.process_list.processes[process_id];
    let artifact_sizes =
        collect_integrand_artifact_sizes(state_folder, process, &detail.integrand_name)?;
    let generation_summary = state.generation_summary(process_id, &detail.integrand_name);
    render_integrand_detail_from_info(
        &detail,
        &artifact_sizes,
        generation_summary,
        requested_graphs,
        requested_categories,
    )
}

fn category_rows(
    group: &IntegrandGraphGroupInfo,
    category: IntegrandDisplayCategory,
) -> Vec<CategoryRow> {
    match category {
        IntegrandDisplayCategory::Generation => Vec::new(),
        IntegrandDisplayCategory::Orientation => group
            .orientations
            .iter()
            .map(|orientation| CategoryRow {
                primary_id: render_category_id(orientation.orientation_id, false),
                secondary_id: String::new(),
                tertiary_id: String::new(),
                value: render_orientation_signature(
                    &orientation.signature,
                    &group.orientation_edge_ids,
                ),
            })
            .collect(),
        IntegrandDisplayCategory::LoopMomentumBasis => group
            .loop_momentum_bases
            .iter()
            .map(|basis| CategoryRow {
                primary_id: render_category_id(basis.basis_id, basis.matches_generation_basis),
                secondary_id: render_optional_category_id(
                    basis.channel_id,
                    basis.matches_generation_basis,
                ),
                tertiary_id: String::new(),
                value: render_edge_ids_with_highlight(
                    &basis.edge_ids,
                    basis.matches_generation_basis,
                ),
            })
            .collect(),
        IntegrandDisplayCategory::Cuts => group
            .cuts
            .iter()
            .map(|cut| CategoryRow {
                primary_id: render_category_id(cut.cut_id, false),
                secondary_id: render_threshold_esurface_ids(&cut.left_threshold_esurface_ids),
                tertiary_id: render_threshold_esurface_ids(&cut.right_threshold_esurface_ids),
                value: render_cut(&cut.edge_ids, cut.raising_power),
            })
            .collect(),
    }
}

fn category_group_header(
    group: &IntegrandGraphGroupInfo,
    category: IntegrandDisplayCategory,
) -> String {
    match category {
        IntegrandDisplayCategory::Orientation => {
            render_orientation_edge_ids(&group.orientation_edge_ids)
        }
        IntegrandDisplayCategory::Generation
        | IntegrandDisplayCategory::LoopMomentumBasis
        | IntegrandDisplayCategory::Cuts => String::new(),
    }
}

fn render_integrand_thresholds_table(
    groups: &[&IntegrandGraphGroupInfo],
) -> Result<Option<String>> {
    let mut builder = Builder::new();
    builder.push_record([
        "Group #".bold().blue().to_string(),
        "Graph".bold().blue().to_string(),
        "Esurface".bold().blue().to_string(),
        "edges".bold().blue().to_string(),
    ]);

    let mut inserted_blocks = 0usize;
    let mut separator_rows = vec![1usize];
    let mut next_row = 1usize;

    for group in groups {
        if group.threshold_esurfaces.is_empty() {
            continue;
        }

        let master = master_graph(group)?;
        inserted_blocks += 1;
        builder.push_record([
            format!("#{}", group.group_id).blue().to_string(),
            format!("#{} : {}", master.graph_id, master.name)
                .green()
                .bold()
                .to_string(),
            String::new(),
            String::new(),
        ]);
        next_row += 1;
        separator_rows.push(next_row);

        for threshold in &group.threshold_esurfaces {
            builder.push_record([
                String::new(),
                String::new(),
                format!("#{}", threshold.esurface_id).yellow().to_string(),
                render_threshold_esurface(&threshold.edge_ids),
            ]);
            next_row += 1;
        }
        separator_rows.push(next_row);
    }

    if inserted_blocks == 0 {
        return Ok(None);
    }

    let mut table = builder.build();
    let mut style = Theme::from_style(Style::rounded().remove_horizontals());
    let mut unique_separator_rows = separator_rows.into_iter().unique().collect_vec();
    if let Some(last_row) = unique_separator_rows.pop() {
        for row in unique_separator_rows {
            style.insert_horizontal_line(
                row,
                HorizontalLine::new('─')
                    .intersection('┼')
                    .left('├')
                    .right('┤'),
            );
        }
        style.insert_horizontal_line(
            last_row,
            HorizontalLine::new('─')
                .intersection('┴')
                .left('╰')
                .right('╯'),
        );
    }
    table.with(style);
    Ok(Some(table.to_string()))
}

fn render_integrand_category_table(
    groups: &[&IntegrandGraphGroupInfo],
    category: IntegrandDisplayCategory,
) -> Result<Option<String>> {
    let mut builder = Builder::new();
    match category {
        IntegrandDisplayCategory::LoopMomentumBasis => {
            builder.push_record([
                "Group #".bold().blue().to_string(),
                "Graph".bold().blue().to_string(),
                "basis ID".bold().blue().to_string(),
                "channel ID".bold().blue().to_string(),
                category.value_header().bold().blue().to_string(),
            ]);
        }
        IntegrandDisplayCategory::Cuts => {
            builder.push_record([
                "Group #".bold().blue().to_string(),
                "Graph".bold().blue().to_string(),
                "ID".bold().blue().to_string(),
                category.value_header().bold().blue().to_string(),
                "left thresholds".bold().blue().to_string(),
                "right thresholds".bold().blue().to_string(),
            ]);
        }
        _ => {
            builder.push_record([
                "Group #".bold().blue().to_string(),
                "Graph".bold().blue().to_string(),
                "ID".bold().blue().to_string(),
                category.value_header().bold().blue().to_string(),
            ]);
        }
    }

    let mut inserted_blocks = 0usize;
    let mut separator_rows = vec![1usize];
    let mut next_row = 1usize;

    for group in groups {
        let rows = category_rows(group, category);
        if rows.is_empty() {
            continue;
        }

        let master = master_graph(group)?;
        inserted_blocks += 1;
        match category {
            IntegrandDisplayCategory::LoopMomentumBasis => {
                builder.push_record([
                    format!("#{}", group.group_id).blue().to_string(),
                    format!("#{} : {}", master.graph_id, master.name)
                        .green()
                        .bold()
                        .to_string(),
                    String::new(),
                    String::new(),
                    category_group_header(group, category),
                ]);
            }
            IntegrandDisplayCategory::Cuts => {
                builder.push_record([
                    format!("#{}", group.group_id).blue().to_string(),
                    format!("#{} : {}", master.graph_id, master.name)
                        .green()
                        .bold()
                        .to_string(),
                    String::new(),
                    category_group_header(group, category),
                    String::new(),
                    String::new(),
                ]);
            }
            _ => {
                builder.push_record([
                    format!("#{}", group.group_id).blue().to_string(),
                    format!("#{} : {}", master.graph_id, master.name)
                        .green()
                        .bold()
                        .to_string(),
                    String::new(),
                    category_group_header(group, category),
                ]);
            }
        }
        next_row += 1;
        separator_rows.push(next_row);

        for row in rows {
            match category {
                IntegrandDisplayCategory::LoopMomentumBasis => {
                    builder.push_record([
                        String::new(),
                        String::new(),
                        row.primary_id,
                        row.secondary_id,
                        row.value,
                    ]);
                }
                IntegrandDisplayCategory::Cuts => {
                    builder.push_record([
                        String::new(),
                        String::new(),
                        row.primary_id,
                        row.value,
                        row.secondary_id,
                        row.tertiary_id,
                    ]);
                }
                _ => {
                    builder.push_record([String::new(), String::new(), row.primary_id, row.value]);
                }
            }
            next_row += 1;
        }
        separator_rows.push(next_row);
    }

    if inserted_blocks == 0 {
        return Ok(None);
    }

    let mut table = builder.build();
    let mut style = Theme::from_style(Style::rounded().remove_horizontals());
    let mut unique_separator_rows = separator_rows.into_iter().unique().collect_vec();
    if let Some(last_row) = unique_separator_rows.pop() {
        for row in unique_separator_rows {
            style.insert_horizontal_line(
                row,
                HorizontalLine::new('─')
                    .intersection('┼')
                    .left('├')
                    .right('┤'),
            );
        }
        style.insert_horizontal_line(
            last_row,
            HorizontalLine::new('─')
                .intersection('┴')
                .left('╰')
                .right('╯'),
        );
    }
    table.with(style);
    Ok(Some(table.to_string()))
}

fn render_generation_category(summary: Option<&IntegrandGenerationSummary>) -> String {
    summary
        .and_then(|summary| {
            render_generation_summary(&summary.reports, summary.peak_ram_bytes, None, None)
        })
        .unwrap_or_else(|| {
            "No generation summary has been recorded for this integrand."
                .yellow()
                .to_string()
        })
}

fn render_integrand_detail_from_info(
    detail: &IntegrandInfo,
    artifact_sizes: &IntegrandArtifactSizes,
    generation_summary: Option<&IntegrandGenerationSummary>,
    requested_graphs: &[String],
    requested_categories: &[IntegrandDisplayCategory],
) -> Result<()> {
    info!(
        "{}",
        format!(
            "Integrand '{}' for process #{} ({})",
            detail.integrand_name, detail.process_id, detail.process_name
        )
        .bold()
        .blue()
    );

    let summary_rows = vec![
        DetailSummaryRow {
            field: "kind".to_string(),
            value: detail.kind.to_string().yellow().to_string(),
        },
        DetailSummaryRow {
            field: "compile enabled".to_string(),
            value: detail
                .generation_compilation
                .compile_enabled()
                .to_string()
                .yellow()
                .to_string(),
        },
        DetailSummaryRow {
            field: "generation backend".to_string(),
            value: detail
                .generation_compilation
                .active_backend_name()
                .yellow()
                .to_string(),
        },
        DetailSummaryRow {
            field: "generation compile options".to_string(),
            value: detail
                .generation_compilation
                .external_options()
                .map(ToString::to_string)
                .unwrap_or_else(|| "(none)".to_string())
                .magenta()
                .to_string(),
        },
        DetailSummaryRow {
            field: "active f64 backend".to_string(),
            value: detail.active_f64_backend.to_string().yellow().to_string(),
        },
        DetailSummaryRow {
            field: "graphs".to_string(),
            value: detail.graph_count.to_string().yellow().to_string(),
        },
        DetailSummaryRow {
            field: "graph groups".to_string(),
            value: detail.graph_group_count.to_string().yellow().to_string(),
        },
        DetailSummaryRow {
            field: ".bin size".to_string(),
            value: format_artifact_size(artifact_sizes.bin_disk_size_bytes),
        },
        DetailSummaryRow {
            field: ".so size".to_string(),
            value: format_artifact_size(artifact_sizes.so_disk_size_bytes),
        },
    ];
    let mut summary_table = tabled::Table::new(summary_rows);
    summary_table.with(Style::rounded());
    info!(
        "
{summary_table}"
    );

    let groups = filtered_graph_groups(detail, requested_graphs)?;
    let categories = selected_integrand_categories(detail, requested_categories)?;
    for category in categories {
        if category == IntegrandDisplayCategory::Generation {
            let table = render_generation_category(generation_summary);
            info!(
                "
{}
{table}",
                category.title().bold().blue()
            );
            continue;
        }

        if let Some(table) = render_integrand_category_table(&groups, category)? {
            info!(
                "
{}
{table}",
                category.title().bold().blue()
            );
        }

        if category == IntegrandDisplayCategory::Cuts {
            if let Some(table) = render_integrand_thresholds_table(&groups)? {
                info!(
                    "
{}
{table}",
                    "Threshold esurfaces".bold().blue()
                );
            }
        }
    }

    Ok(())
}

fn render_processes_table(state: &State, state_folder: &Path) -> Result<()> {
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
        ".so size / integrand".bold().blue().to_string(),
    ]);

    for process in &state.process_list.processes {
        let metrics = collect_integrand_metrics(state_folder, process)?;
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
            format_metrics_column(&metrics, |metric| {
                format_artifact_size(metric.bin_disk_size_bytes)
            }),
            format_metrics_column(&metrics, |metric| {
                format_artifact_size(metric.so_disk_size_bytes)
            }),
        ]);
    }

    let mut table = builder.build();
    table.with(Style::rounded());
    info!("{}", "Processes".bold().blue());
    info!("\n{table}");
    Ok(())
}

fn render_integrands_table(
    state: &State,
    state_folder: &Path,
    process_filter: Option<usize>,
) -> Result<()> {
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
        ".so size".bold().blue().to_string(),
    ]);

    let mut num_rows = 0usize;
    for (process_index, process) in state.process_list.processes.iter().enumerate() {
        if process_filter.is_some_and(|selected| selected != process_index) {
            continue;
        }
        let metrics = collect_integrand_metrics(state_folder, process)?;
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
                format_artifact_size(metric.bin_disk_size_bytes),
                format_artifact_size(metric.so_disk_size_bytes),
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

fn collect_integrand_metrics(
    state_folder: &Path,
    process: &Process,
) -> Result<Vec<IntegrandMetrics>> {
    match &process.collection {
        ProcessCollection::Amplitudes(amplitudes) => amplitudes
            .values()
            .map(|amplitude| integrand_metrics_from_amplitude(state_folder, process, amplitude))
            .collect(),
        ProcessCollection::CrossSections(cross_sections) => cross_sections
            .values()
            .map(|cross_section| {
                integrand_metrics_from_cross_section(state_folder, process, cross_section)
            })
            .collect(),
    }
}

fn integrand_metrics_from_amplitude(
    state_folder: &Path,
    process: &Process,
    amplitude: &Amplitude,
) -> Result<IntegrandMetrics> {
    let artifact_sizes = collect_integrand_artifact_sizes(state_folder, process, &amplitude.name)?;
    Ok(IntegrandMetrics {
        name: amplitude.name.clone(),
        graphs: amplitude.graphs.len(),
        graph_groups: amplitude.graph_group_structure.len(),
        bin_disk_size_bytes: artifact_sizes.bin_disk_size_bytes,
        so_disk_size_bytes: artifact_sizes.so_disk_size_bytes,
    })
}

fn integrand_metrics_from_cross_section(
    state_folder: &Path,
    process: &Process,
    cross_section: &CrossSection,
) -> Result<IntegrandMetrics> {
    let artifact_sizes =
        collect_integrand_artifact_sizes(state_folder, process, &cross_section.name)?;
    Ok(IntegrandMetrics {
        name: cross_section.name.clone(),
        graphs: cross_section.supergraphs.len(),
        graph_groups: cross_section.graph_group_structure.len(),
        bin_disk_size_bytes: artifact_sizes.bin_disk_size_bytes,
        so_disk_size_bytes: artifact_sizes.so_disk_size_bytes,
    })
}

fn collect_integrand_artifact_sizes(
    state_folder: &Path,
    process: &Process,
    integrand_name: &str,
) -> Result<IntegrandArtifactSizes> {
    let artifact_root = integrand_artifact_root(state_folder, process, integrand_name);
    Ok(IntegrandArtifactSizes {
        bin_disk_size_bytes: disk_usage_for_extension(&artifact_root, "bin")?,
        so_disk_size_bytes: disk_usage_for_extension(&artifact_root, "so")?,
    })
}

fn integrand_artifact_root(
    state_folder: &Path,
    process: &Process,
    integrand_name: &str,
) -> PathBuf {
    let collection_dir = match &process.collection {
        ProcessCollection::Amplitudes(_) => "amplitudes",
        ProcessCollection::CrossSections(_) => "cross_sections",
    };
    state_folder
        .join("processes")
        .join(collection_dir)
        .join(&process.definition.folder_name)
        .join(integrand_name)
}

fn disk_usage_for_extension(root: &Path, extension: &str) -> Result<Option<u64>> {
    if !root.try_exists().with_context(|| {
        format!(
            "While checking saved integrand artifacts under {}",
            root.display()
        )
    })? {
        return Ok(None);
    }

    let mut total = 0u64;
    for entry in WalkDir::new(root) {
        let entry = entry.with_context(|| {
            format!(
                "While walking saved integrand artifacts under {}",
                root.display()
            )
        })?;
        if !entry.file_type().is_file() || entry.path().extension() != Some(OsStr::new(extension)) {
            continue;
        }

        let metadata = entry
            .metadata()
            .with_context(|| format!("While reading metadata for {}", entry.path().display()))?;
        total = total
            .checked_add(file_disk_usage_bytes(&metadata))
            .ok_or_else(|| {
                eyre!(
                    "Artifact size overflow while summing .{} files under {}",
                    extension,
                    root.display()
                )
            })?;
    }
    Ok(Some(total))
}

fn file_disk_usage_bytes(metadata: &fs::Metadata) -> u64 {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        metadata.blocks().saturating_mul(512)
    }

    #[cfg(not(unix))]
    {
        metadata.len()
    }
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

fn format_artifact_size(size: Option<u64>) -> String {
    size.map(|size| format_bytes(size).magenta().to_string())
        .unwrap_or_else(|| "(not saved)".dimmed().to_string())
}

fn format_bytes(size: u64) -> String {
    const KIB: u64 = 1024;
    const MIB: u64 = KIB * 1024;
    const GIB: u64 = MIB * 1024;

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

fn command_block_contents(block: &CommandsBlock) -> String {
    if block.commands.is_empty() {
        return "(empty)".to_string();
    }

    block
        .commands
        .iter()
        .map(display_command)
        .collect::<Vec<_>>()
        .join("\n")
}

fn render_command_blocks(run_history: &RunHistory, selected_name: Option<&str>) -> Result<()> {
    if let Some(name) = selected_name {
        let block = run_history.command_block(name).ok_or_else(|| {
            eyre!(
                "Unknown command block '{}'. Available command blocks: {}",
                name,
                run_history
                    .command_blocks
                    .iter()
                    .map(|block| block.name.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        })?;
        info!("\n{}", command_block_contents(block));
        return Ok(());
    }

    if run_history.command_blocks.is_empty() {
        info!("{}", "No active command blocks.".yellow());
        return Ok(());
    }

    info!("{}", "Command blocks".bold().blue());
    info!("\n{}", render_command_blocks_table(run_history));
    Ok(())
}

fn render_command_blocks_table(run_history: &RunHistory) -> String {
    let mut builder = Builder::new();
    builder.push_record([
        "name".bold().blue().to_string(),
        "placeholders".bold().blue().to_string(),
        "commands".bold().blue().to_string(),
    ]);

    for block in &run_history.command_blocks {
        let placeholders = run_history
            .command_block_placeholder_names(&block.name)
            .into_iter()
            .map(|name| format!("{{{name}}}"))
            .collect::<Vec<_>>()
            .join("\n");
        builder.push_record([
            block.name.green().bold().to_string(),
            if placeholders.is_empty() {
                "-".to_string()
            } else {
                placeholders
            },
            command_block_contents(block),
        ]);
    }

    let mut table = builder.build();
    let mut style = Theme::from_style(Style::rounded().remove_horizontals());
    for row in 1..=run_history.command_blocks.len() {
        style.insert_horizontal_line(
            row,
            HorizontalLine::new('─')
                .intersection('┼')
                .left('├')
                .right('┤'),
        );
    }
    table.with(style);
    table.to_string()
}

fn render_named_process_settings(
    state: &State,
    target: &DisplayProcessNamedSettingsArgs,
    setting_kind: NamedProcessSettingKind,
) -> Result<()> {
    let process_id = state.resolve_process_ref(target.process.process.as_ref())?;
    let process = &state.process_list.processes[process_id];
    let integrand_names = selected_integrand_names(process, target.process.integrand_name.as_ref());

    if let Some(name) = target.name.as_deref() {
        return render_named_process_setting_detail(
            state,
            process_id,
            &process.definition.folder_name,
            process.definition.process_id,
            &integrand_names,
            setting_kind,
            name,
        );
    }

    let mut builder = Builder::new();
    let show_integrand = integrand_names.len() > 1;
    let mut header = Vec::new();
    if show_integrand {
        header.push("integrand".bold().blue().to_string());
    }
    header.push("name".bold().blue().to_string());
    header.push("kind".bold().blue().to_string());
    header.push("details".bold().blue().to_string());
    builder.push_record(header);

    let mut num_rows = 0usize;
    for integrand_name in &integrand_names {
        let integrand = state
            .process_list
            .get_integrand(process_id, integrand_name)?
            .require_generated()?;
        match setting_kind {
            NamedProcessSettingKind::Quantity => {
                for (name, settings) in &integrand.get_settings().quantities {
                    num_rows += 1;
                    let mut record = Vec::new();
                    if show_integrand {
                        record.push(integrand_name.clone().cyan().to_string());
                    }
                    record.push(name.clone().green().to_string());
                    record.push(quantity_kind(settings).yellow().to_string());
                    record.push(summarize_quantity(settings));
                    builder.push_record(record);
                }
            }
            NamedProcessSettingKind::Observable => {
                for (name, settings) in &integrand.get_settings().observables {
                    num_rows += 1;
                    let mut record = Vec::new();
                    if show_integrand {
                        record.push(integrand_name.clone().cyan().to_string());
                    }
                    record.push(name.clone().green().to_string());
                    record.push(observable_kind(settings).yellow().to_string());
                    record.push(summarize_observable(settings));
                    builder.push_record(record);
                }
            }
            NamedProcessSettingKind::Selector => {
                for (name, settings) in &integrand.get_settings().selectors {
                    num_rows += 1;
                    let mut record = Vec::new();
                    if show_integrand {
                        record.push(integrand_name.clone().cyan().to_string());
                    }
                    record.push(name.clone().green().to_string());
                    record.push(selector_kind(settings).yellow().to_string());
                    record.push(summarize_selector(settings));
                    builder.push_record(record);
                }
            }
        }
    }

    if num_rows == 0 {
        info!(
            "{}",
            format!(
                "No {} configured for process #{} ({})",
                setting_kind.plural(),
                process.definition.process_id,
                process.definition.folder_name
            )
            .yellow()
        );
        return Ok(());
    }

    let mut table = builder.build();
    table.with(Style::rounded());
    let title = if show_integrand {
        format!(
            "{} for process #{} ({})",
            setting_kind.plural(),
            process.definition.process_id,
            process.definition.folder_name
        )
    } else {
        format!(
            "{} for process #{} ({}) integrand '{}'",
            setting_kind.plural(),
            process.definition.process_id,
            process.definition.folder_name,
            integrand_names[0]
        )
    };
    info!("{}", title.bold().blue());
    info!("\n{table}");
    Ok(())
}

fn render_named_process_setting_detail(
    state: &State,
    process_id: usize,
    process_name: &str,
    process_display_id: usize,
    integrand_names: &[String],
    setting_kind: NamedProcessSettingKind,
    name: &str,
) -> Result<()> {
    let mut selected_roots = Vec::new();
    let mut missing_integrands = Vec::new();

    for integrand_name in integrand_names {
        let integrand = state
            .process_list
            .get_integrand(process_id, integrand_name)?
            .require_generated()?;
        let serialized = serialize_runtime_named_settings(integrand.get_settings())?;
        let map = match setting_kind {
            NamedProcessSettingKind::Quantity => &serialized.quantities,
            NamedProcessSettingKind::Observable => &serialized.observables,
            NamedProcessSettingKind::Selector => &serialized.selectors,
        };
        if let Some(root) = map.get(name) {
            selected_roots.push((integrand_name.clone(), root.clone()));
        } else {
            missing_integrands.push(integrand_name.clone());
        }
    }

    if selected_roots.is_empty() {
        return Err(eyre!(
            "No {} named '{}' found for process #{} ({})",
            setting_kind.singular(),
            name,
            process_display_id,
            process_name
        ));
    }

    if !missing_integrands.is_empty() {
        return Err(eyre!(
            "{} '{}' is missing from integrand(s): {}",
            setting_kind.singular(),
            name,
            missing_integrands.join(", ")
        ));
    }

    for (integrand_name, root) in selected_roots {
        render_settings_table(
            &format!(
                "{} '{}' for process #{} ({}) integrand '{}'",
                setting_kind.singular(),
                name,
                process_display_id,
                process_name,
                integrand_name
            ),
            &root,
            None,
        )?;
    }

    Ok(())
}

fn selected_integrand_names(process: &Process, requested: Option<&String>) -> Vec<String> {
    if let Some(integrand_name) = requested {
        vec![integrand_name.clone()]
    } else {
        process
            .get_integrand_names()
            .into_iter()
            .map(str::to_string)
            .collect()
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
                    let integrand = state
                        .process_list
                        .get_integrand(process_id, name)?
                        .require_generated()?;
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
                            .get_integrand(process_id, integrand_name)?
                            .require_generated()?;
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

    use crate::{
        commands::Commands,
        state::{CommandHistory, CommandsBlock, ProcessRef, RunHistory},
        CLISettings, Repl,
    };
    use gammalooprs::settings::RuntimeSettings;

    use super::{
        command_block_contents, format_bytes, render_command_blocks_table,
        serialize_settings_with_defaults, value_at_path, Display, DisplaySettingsTarget,
        IntegrandDisplayCategory,
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
    fn parse_display_model_target() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "model",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
            "--show-particles",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Model {
                process,
                integrand_name,
                show_particles,
                ..
            }) => {
                assert_eq!(
                    process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(integrand_name, Some("LO".to_string()));
                assert!(show_particles);
            }
            other => panic!("Expected display model command, got {other:?}"),
        }
    }

    #[test]
    fn display_model_requires_process_and_integrand_together() {
        let err = Repl::try_parse_from(["gammaloop", "display", "model", "-p", "epem_a_tth"])
            .unwrap_err();
        assert!(err.to_string().contains("--integrand-name"));
    }

    #[test]
    fn parse_display_quantities_without_name() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "quantities",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Quantities { target }) => {
                assert_eq!(
                    target.process.process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(target.process.integrand_name, Some("LO".to_string()));
                assert_eq!(target.name, None);
            }
            other => panic!("Expected display quantities command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_selector_with_name() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "selectors",
            "-p",
            "epem_a_tth",
            "top_cut",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Selectors { target }) => {
                assert_eq!(
                    target.process.process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(target.process.integrand_name, None);
                assert_eq!(target.name.as_deref(), Some("top_cut"));
            }
            other => panic!("Expected display selectors command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_command_block_without_name() {
        let repl = Repl::try_parse_from(["gammaloop", "display", "command_block"]).unwrap();

        match repl.command {
            Commands::Display(Display::CommandBlock { name }) => {
                assert_eq!(name, None);
            }
            other => panic!("Expected display command_block command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_command_block_with_name() {
        let repl =
            Repl::try_parse_from(["gammaloop", "display", "command_block", "alpha"]).unwrap();

        match repl.command {
            Commands::Display(Display::CommandBlock { name }) => {
                assert_eq!(name.as_deref(), Some("alpha"));
            }
            other => panic!("Expected display command_block command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_command_block_rejects_hyphenated_name() {
        assert!(Repl::try_parse_from(["gammaloop", "display", "command-block"]).is_err());
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
    fn parse_display_integrand_without_process() {
        let repl = Repl::try_parse_from(["gammaloop", "display", "integrand"]).unwrap();

        match repl.command {
            Commands::Display(Display::Integrands {
                process,
                integrand_name,
                graphs,
                categories,
            }) => {
                assert_eq!(process, None);
                assert_eq!(integrand_name, None);
                assert!(graphs.is_empty());
                assert!(categories.is_empty());
            }
            other => panic!("Expected display integrand command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_integrand_detail() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "integrand",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Integrands {
                process,
                integrand_name,
                graphs,
                categories,
            }) => {
                assert_eq!(
                    process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(integrand_name.as_deref(), Some("LO"));
                assert!(graphs.is_empty());
                assert!(categories.is_empty());
            }
            other => panic!("Expected display integrand command, got {other:?}"),
        }
    }

    #[test]
    fn parse_display_integrand_detail_filters() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "display",
            "integrand",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
            "-g",
            "GL0",
            "GL2",
            "--category",
            "generation",
            "orientation",
            "cuts",
        ])
        .unwrap();

        match repl.command {
            Commands::Display(Display::Integrands {
                process,
                integrand_name,
                graphs,
                categories,
            }) => {
                assert_eq!(
                    process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(integrand_name.as_deref(), Some("LO"));
                assert_eq!(graphs, vec!["GL0".to_string(), "GL2".to_string()]);
                assert_eq!(
                    categories,
                    vec![
                        IntegrandDisplayCategory::Generation,
                        IntegrandDisplayCategory::Orientation,
                        IntegrandDisplayCategory::Cuts,
                    ]
                );
            }
            other => panic!("Expected display integrand command, got {other:?}"),
        }
    }

    #[test]
    fn format_bytes_uses_binary_units() {
        assert_eq!(format_bytes(999), "999 B");
        assert_eq!(format_bytes(1024), "1.00 KiB");
        assert_eq!(format_bytes(1024 * 1024), "1.00 MiB");
    }

    #[test]
    fn command_block_contents_preserve_raw_commands() {
        let block = CommandsBlock {
            name: "alpha".to_string(),
            commands: vec![
                CommandHistory::from_raw_string("display processes").unwrap(),
                CommandHistory::from_raw_string("set process -p triangle -i LO defaults").unwrap(),
            ],
        };

        assert_eq!(
            command_block_contents(&block),
            "display processes\nset process -p triangle -i LO defaults"
        );
    }

    #[test]
    fn command_block_contents_marks_empty_blocks() {
        let run_history = RunHistory {
            command_blocks: vec![CommandsBlock {
                name: "empty".to_string(),
                commands: Vec::new(),
            }],
            ..RunHistory::default()
        };

        assert_eq!(
            command_block_contents(&run_history.command_blocks[0]),
            "(empty)"
        );
    }

    #[test]
    fn command_blocks_table_separates_each_block() {
        let run_history = RunHistory {
            command_blocks: vec![
                CommandsBlock {
                    name: "alpha".to_string(),
                    commands: vec![CommandHistory::from_raw_string("display processes").unwrap()],
                },
                CommandsBlock {
                    name: "beta".to_string(),
                    commands: vec![CommandHistory::from_raw_string("display model").unwrap()],
                },
            ],
            ..RunHistory::default()
        };

        let rendered = render_command_blocks_table(&run_history);

        assert!(rendered.contains("alpha"), "{rendered}");
        assert!(rendered.contains("beta"), "{rendered}");
        assert_eq!(
            rendered
                .lines()
                .filter(|line| line.starts_with('├') && line.ends_with('┤'))
                .count(),
            2,
            "{rendered}"
        );
    }

    #[test]
    fn command_blocks_table_lists_command_block_placeholders() {
        let run_history = RunHistory {
            command_blocks: vec![
                CommandsBlock {
                    name: "inner".to_string(),
                    commands: vec![CommandHistory::from_raw_string(
                        "set global kv global.display_directive={level}",
                    )
                    .unwrap()],
                },
                CommandsBlock {
                    name: "outer".to_string(),
                    commands: vec![CommandHistory::from_raw_string(
                        "run inner -D level={outer_level}",
                    )
                    .unwrap()],
                },
            ],
            ..RunHistory::default()
        };

        let rendered = render_command_blocks_table(&run_history);

        assert!(rendered.contains("placeholders"), "{rendered}");
        assert!(rendered.contains("{level}"), "{rendered}");
        assert!(rendered.contains("{outer_level}"), "{rendered}");
        assert_eq!(
            run_history
                .command_block_placeholder_names("outer")
                .into_iter()
                .collect::<Vec<_>>(),
            vec!["outer_level".to_string()]
        );
    }
}
