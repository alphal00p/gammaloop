use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

use clap::{Args, Subcommand, ValueEnum};
use color_eyre::{
    eyre::{eyre, Context},
    Result,
};
use gammalooprs::{
    graph::Graph,
    integrands::process::ProcessIntegrand,
    processes::{Amplitude, CrossSection, Process, ProcessCollection},
    settings::global::UniformNumeratorSamplingScale,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use symbolica::atom::AtomCore;
use three_dimensional_reps::{
    generate_3d_expression, graph_info, render_expression_summary, validate_parsed_graph,
    DisplayOptions, GenerationError, GraphInfo, GraphValidation, NumeratorDisplay,
    NumeratorSamplingScaleMode, OrientationID, RepresentationMode, ThreeDExpression,
    ThreeDGraphSource,
};

use crate::{
    completion::CompletionArgExt,
    state::{ProcessRef, State},
    CLISettings,
};

#[derive(Debug, Subcommand, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum ThreeDRep {
    Validate(Validate),
    Build(Build),
}

impl ThreeDRep {
    pub fn run(&self, state: &State, global_cli_settings: &CLISettings) -> Result<()> {
        match self {
            Self::Validate(command) => command.run(state, global_cli_settings),
            Self::Build(command) => command.run(state, global_cli_settings),
        }
    }
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct GraphSelectorArgs {
    /// Process reference: #<id>, name:<name>, or <id>/<name>.
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// Integrand name inside the selected process.
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// Graph id, graph name, or inspect label such as "#3 : graph_name".
    #[arg(short = 'g', long = "graph", value_name = "GRAPH")]
    pub graph: String,
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Validate {
    #[command(flatten)]
    pub selection: GraphSelectorArgs,

    #[arg(long, value_hint = clap::ValueHint::FilePath)]
    pub json_out: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, ValueEnum, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum CliNumeratorSamplesNormalization {
    #[value(name = "never_M", alias = "never-m", alias = "none")]
    NeverM,
    #[value(name = "M_for_all", alias = "m-for-all", alias = "all")]
    MForAll,
    #[value(
        name = "M_for_beyond_quadratic_only",
        alias = "m-for-beyond-quadratic-only",
        alias = "beyond-quadratic"
    )]
    MForBeyondQuadraticOnly,
}

impl CliNumeratorSamplesNormalization {
    fn resolve(
        override_mode: Option<Self>,
        global_mode: UniformNumeratorSamplingScale,
    ) -> NumeratorSamplingScaleMode {
        match override_mode {
            Some(Self::NeverM) => NumeratorSamplingScaleMode::None,
            Some(Self::MForAll) => NumeratorSamplingScaleMode::All,
            Some(Self::MForBeyondQuadraticOnly) => NumeratorSamplingScaleMode::BeyondQuadratic,
            None => match global_mode {
                UniformNumeratorSamplingScale::None => NumeratorSamplingScaleMode::None,
                UniformNumeratorSamplingScale::BeyondQuadratic => {
                    NumeratorSamplingScaleMode::BeyondQuadratic
                }
                UniformNumeratorSamplingScale::All => NumeratorSamplingScaleMode::All,
            },
        }
    }
}

#[derive(
    Debug, Clone, Copy, ValueEnum, Serialize, Deserialize, JsonSchema, PartialEq, Eq, Default,
)]
#[serde(rename_all = "snake_case")]
pub enum CliRepresentationMode {
    #[default]
    Cff,
    Ltd,
}

impl From<CliRepresentationMode> for RepresentationMode {
    fn from(value: CliRepresentationMode) -> Self {
        match value {
            CliRepresentationMode::Cff => Self::Cff,
            CliRepresentationMode::Ltd => Self::Ltd,
        }
    }
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Build {
    #[command(flatten)]
    pub selection: GraphSelectorArgs,

    /// Three-dimensional representation to build. LTD is reserved but not implemented yet.
    #[serde(default)]
    #[arg(long, alias = "family", value_enum, default_value = "cff")]
    pub representation: CliRepresentationMode,

    /// Override the global numerator sampling normalization for this build.
    #[arg(
        long = "numerator-samples-normalization",
        aliases = ["numerator-sampling-scale-mode", "uniform-numerator-sampling-scale"],
        value_enum
    )]
    pub numerator_samples_normalization: Option<CliNumeratorSamplesNormalization>,

    #[arg(long, default_value_t = false)]
    pub no_save_json: bool,

    #[arg(long, default_value_t = false)]
    pub no_pretty: bool,

    #[arg(long, default_value_t = false)]
    pub no_color: bool,

    #[arg(long)]
    pub show_details_for_orientation: Option<String>,

    #[arg(long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    #[arg(long, value_hint = clap::ValueHint::FilePath)]
    pub json_out: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct ValidateOutput {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    graph: GraphInfo,
    validation: GraphValidation,
}

#[derive(Debug, Serialize)]
struct BuildOutput {
    backend: &'static str,
    family: &'static str,
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    graph: GraphInfo,
    validation: GraphValidation,
    energy_degree_bounds: Option<Vec<(usize, usize)>>,
    numerator_sampling_scale_mode: NumeratorSamplingScaleMode,
    expression: ThreeDExpression<OrientationID>,
}

struct SelectedGraph<'a> {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph: &'a Graph,
}

// Imported graphs and generated integrands store their graph terms in different
// owners. This private borrowed adapter gives Validate and Build one selector path.
enum GraphCatalog<'a> {
    Generated(&'a ProcessIntegrand),
    ImportedAmplitude(&'a Amplitude),
    ImportedCrossSection(&'a CrossSection),
}

impl<'a> GraphCatalog<'a> {
    fn for_integrand(process: &'a Process, integrand_name: &str) -> Result<Self> {
        if let Some(integrand) = process.get_integrand(integrand_name)?.integrand {
            return Ok(Self::Generated(integrand));
        }

        match &process.collection {
            ProcessCollection::Amplitudes(amplitudes) => amplitudes
                .get(integrand_name)
                .map(Self::ImportedAmplitude)
                .ok_or_else(|| eyre!("No amplitude named '{integrand_name}'.")),
            ProcessCollection::CrossSections(cross_sections) => cross_sections
                .get(integrand_name)
                .map(Self::ImportedCrossSection)
                .ok_or_else(|| eyre!("No cross section named '{integrand_name}'.")),
        }
    }

    fn graph_count(&self) -> usize {
        match self {
            Self::Generated(integrand) => integrand.graph_count(),
            Self::ImportedAmplitude(amplitude) => amplitude.graphs.len(),
            Self::ImportedCrossSection(cross_section) => cross_section.supergraphs.len(),
        }
    }

    fn graph_by_id(&self, graph_id: usize) -> Result<&'a Graph> {
        let graph = match self {
            Self::Generated(ProcessIntegrand::Amplitude(integrand)) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| &term.graph),
            Self::Generated(ProcessIntegrand::CrossSection(integrand)) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| &term.graph),
            Self::ImportedAmplitude(amplitude) => {
                amplitude.graphs.get(graph_id).map(|term| &term.graph)
            }
            Self::ImportedCrossSection(cross_section) => cross_section
                .supergraphs
                .get(graph_id)
                .map(|term| &term.graph),
        };

        graph.ok_or_else(|| {
            eyre!(
                "Graph id {graph_id} is out of range; this integrand has {} graphs.",
                self.graph_count()
            )
        })
    }

    fn graph_name_by_id(&self, graph_id: usize) -> Option<&str> {
        match self {
            Self::Generated(ProcessIntegrand::Amplitude(integrand)) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| term.graph.name.as_str()),
            Self::Generated(ProcessIntegrand::CrossSection(integrand)) => integrand
                .data
                .graph_terms
                .get(graph_id)
                .map(|term| term.graph.name.as_str()),
            Self::ImportedAmplitude(amplitude) => amplitude
                .graphs
                .get(graph_id)
                .map(|term| term.graph.name.as_str()),
            Self::ImportedCrossSection(cross_section) => cross_section
                .supergraphs
                .get(graph_id)
                .map(|term| term.graph.name.as_str()),
        }
    }

    fn resolve_graph_id(&self, selector: &str) -> Result<usize> {
        let selector = selector.trim();
        if let Some(id) = parse_graph_id(selector) {
            if id < self.graph_count() {
                return Ok(id);
            }
        }

        let matches = (0..self.graph_count())
            .filter(|&graph_id| {
                self.graph_name_by_id(graph_id).is_some_and(|name| {
                    selector == name || selector == format!("#{graph_id} : {name}")
                })
            })
            .collect::<Vec<_>>();

        match matches.as_slice() {
            [graph_id] => Ok(*graph_id),
            [] => Err(eyre!(
                "Could not resolve graph selector '{selector}'. Use a graph id, graph name, or inspect display name."
            )),
            matches => Err(eyre!(
                "Graph selector '{selector}' is ambiguous and matches graph ids {matches:?}. Use an explicit graph id."
            )),
        }
    }
}

impl Validate {
    fn run(&self, state: &State, global_cli_settings: &CLISettings) -> Result<()> {
        let selected = select_graph(state, &self.selection)?;
        let parsed = selected.graph.to_three_d_parsed_graph()?;
        let output = ValidateOutput {
            process_id: selected.process_id,
            integrand_name: selected.integrand_name,
            graph_id: selected.graph_id,
            graph_name: selected.graph.name.clone(),
            graph: graph_info(&parsed),
            validation: validate_parsed_graph(&parsed),
        };
        if let Some(path) = &self.json_out {
            global_cli_settings
                .ensure_write_target_outside_active_state(path, "write 3Drep validation JSON")?;
        }
        write_or_print(
            self.json_out.as_deref(),
            &serde_json::to_string_pretty(&output)?,
        )
    }
}

impl Build {
    fn run(&self, state: &State, global_cli_settings: &CLISettings) -> Result<()> {
        let representation = RepresentationMode::from(self.representation);
        if representation != RepresentationMode::Cff {
            return Err(GenerationError::NotImplemented {
                mode: representation,
            }
            .into());
        }
        let selected = select_graph(state, &self.selection)?;
        let parsed = selected.graph.to_three_d_parsed_graph()?;
        let numerator_sampling_scale_mode = CliNumeratorSamplesNormalization::resolve(
            self.numerator_samples_normalization,
            global_cli_settings
                .global
                .generation
                .uniform_numerator_sampling_scale,
        );
        let mut options = selected
            .graph
            .cff_3d_expression_options(numerator_sampling_scale_mode)?;
        options.representation = representation;
        let initial_state_cut_edges = selected
            .graph
            .iter_edges_of(&selected.graph.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<BTreeSet<_>>();
        options.preserve_internal_edges_as_four_d_denominators = selected
            .graph
            .iter_edges_of(&selected.graph.tree_edges)
            .filter_map(|(pair, edge_id, edge)| {
                (pair.is_paired()
                    && !edge.data.is_dummy
                    && !initial_state_cut_edges.contains(&edge_id))
                .then_some(usize::from(edge_id))
            })
            .collect();
        let energy_degree_bounds = options.energy_degree_bounds.clone();
        let expression = generate_3d_expression(selected.graph, &options)?.expression;
        let output = BuildOutput {
            backend: "gammaloop-3Drep",
            family: "cff",
            process_id: selected.process_id,
            integrand_name: selected.integrand_name.clone(),
            graph_id: selected.graph_id,
            graph_name: selected.graph.name.clone(),
            graph: graph_info(&parsed),
            validation: validate_parsed_graph(&parsed),
            energy_degree_bounds: energy_degree_bounds.clone(),
            numerator_sampling_scale_mode,
            expression,
        };

        if !self.no_save_json {
            let workspace = self.workspace_path(global_cli_settings);
            let artifact_dir =
                build_artifact_dir(&workspace, &selected, output.numerator_sampling_scale_mode);
            let json_path = self
                .json_out
                .clone()
                .unwrap_or_else(|| artifact_dir.join("oriented_expression.json"));
            global_cli_settings
                .ensure_write_target_outside_active_state(&workspace, "write 3Drep workspace")?;
            global_cli_settings.ensure_write_target_outside_active_state(
                &json_path,
                "write 3Drep oriented expression",
            )?;
            write_path(&json_path, &serde_json::to_string_pretty(&output)?)?;
            write_latest_expression_pointer(&workspace, &json_path)?;
            println!(
                "Saved 3Drep oriented expression to {}",
                relative_display(&json_path)
            );
        }

        if !self.no_pretty || self.show_details_for_orientation.is_some() {
            let numerator = selected.graph.full_numerator_atom().to_canonical_string();
            println!(
                "{}",
                render_expression_summary(
                    &output.expression,
                    &output.graph,
                    energy_degree_bounds.as_deref(),
                    NumeratorDisplay {
                        original: Some(&numerator),
                        simplified: None,
                    },
                    output.numerator_sampling_scale_mode,
                    &DisplayOptions {
                        use_color: !self.no_color,
                        details_for_orientation: self.show_details_for_orientation.clone(),
                    },
                )
            );
        }

        Ok(())
    }

    fn workspace_path(&self, global_cli_settings: &CLISettings) -> PathBuf {
        self.workspace_path
            .clone()
            .unwrap_or_else(|| default_workspace_path(global_cli_settings))
    }
}

fn select_graph<'a>(state: &'a State, selection: &GraphSelectorArgs) -> Result<SelectedGraph<'a>> {
    let (process_id, integrand_name) = state.find_integrand_ref(
        selection.process.as_ref(),
        selection.integrand_name.as_ref(),
    )?;
    let process = &state.process_list.processes[process_id];
    let catalog = GraphCatalog::for_integrand(process, &integrand_name)
        .with_context(|| format!("while resolving integrand '{integrand_name}'"))?;
    let graph_id = catalog.resolve_graph_id(&selection.graph)?;
    let graph = catalog.graph_by_id(graph_id)?;
    Ok(SelectedGraph {
        process_id,
        integrand_name,
        graph_id,
        graph,
    })
}

fn parse_graph_id(selector: &str) -> Option<usize> {
    let stripped = selector.trim().strip_prefix('#').unwrap_or(selector).trim();
    stripped
        .split_once(':')
        .map_or(stripped, |(id, _)| id)
        .trim()
        .parse()
        .ok()
}

fn default_workspace_path(global_cli_settings: &CLISettings) -> PathBuf {
    if global_cli_settings.session.read_only_state {
        global_cli_settings.cwd_output_path_with_state_name("threed_workspace")
    } else {
        global_cli_settings.state.folder.join("threed_workspace")
    }
}

fn graph_workspace_dir(workspace: &Path, selected: &SelectedGraph<'_>) -> PathBuf {
    workspace
        .join(format!(
            "process_{:04}_{}",
            selected.process_id,
            slug(&selected.integrand_name)
        ))
        .join(format!(
            "graph_{:04}_{}",
            selected.graph_id,
            slug(&selected.graph.name)
        ))
}

fn build_artifact_dir(
    workspace: &Path,
    selected: &SelectedGraph<'_>,
    scale_mode: NumeratorSamplingScaleMode,
) -> PathBuf {
    graph_workspace_dir(workspace, selected)
        .join("build")
        .join(format!("cff_{scale_mode:?}").to_lowercase())
}

fn slug(value: &str) -> String {
    let slug = value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                ch
            } else {
                '_'
            }
        })
        .collect::<String>();
    if slug.is_empty() {
        "unnamed".to_string()
    } else {
        slug
    }
}

fn relative_display(path: &Path) -> String {
    std::env::current_dir()
        .ok()
        .and_then(|cwd| path.strip_prefix(cwd).ok().map(Path::to_path_buf))
        .unwrap_or_else(|| path.to_path_buf())
        .display()
        .to_string()
}

fn write_latest_expression_pointer(workspace: &Path, expression_path: &Path) -> Result<()> {
    write_path(
        &workspace.join("latest_oriented_expression_path.txt"),
        &relative_display(expression_path),
    )
}

fn write_or_print(path: Option<&Path>, text: &str) -> Result<()> {
    if let Some(path) = path {
        write_path(path, text)
    } else {
        println!("{text}");
        Ok(())
    }
}

fn write_path(path: &Path, text: &str) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Could not create directory {}", parent.display()))?;
        }
    }
    fs::write(path, text).with_context(|| format!("Could not write {}", path.display()))
}
