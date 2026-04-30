use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use clap::{Args, Subcommand, ValueEnum};
use color_eyre::{
    eyre::{eyre, Context},
    Result,
};
use gammalooprs::{
    cff::expression::GammaLoopThreeDExpression,
    graph::Graph,
    integrands::process::{evaluators::EvaluatorStack, ProcessIntegrand},
    processes::{Amplitude, CrossSection, EvaluatorSettings, Process, ProcessCollection},
    settings::{global::OrientationPattern, RuntimeSettings},
    utils::symbolica_ext::LogPrint,
};
use linnet::half_edge::involution::HedgePair;
use nu_ansi_term::{Color, Style as AnsiStyle};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tabled::{builder::Builder, settings::Style};
use three_dimensional_reps::{
    generate_3d_expression, graph_info, reconstruct_dot_from_expression, render_expression_summary,
    validate_parsed_graph, DisplayOptions, Generate3DExpressionOptions, GraphInfo, GraphValidation,
    OrientationID, ReconstructDotFormat, ReconstructDotOptions, RepresentationMode,
    ThreeDExpression, ThreeDGraphSource,
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
    Evaluate(Evaluate),
    #[command(name = "test-cff-ltd", alias = "test", alias = "compare")]
    TestCffLtd(TestCffLtd),
    #[command(name = "graph-from-signatures")]
    GraphFromSignatures(GraphFromSignatures),
}

impl ThreeDRep {
    pub fn run(
        &self,
        state: &State,
        global_cli_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        match self {
            Self::Validate(command) => command.run(state),
            Self::Build(command) => {
                command.run(state, global_cli_settings, default_runtime_settings)
            }
            Self::Evaluate(command) => command.run(state, global_cli_settings),
            Self::TestCffLtd(command) => {
                command.run(state, global_cli_settings, default_runtime_settings)
            }
            Self::GraphFromSignatures(command) => command.run(),
        }
    }
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct GraphSelectorArgs {
    /// Process reference: #<id>, name:<name>, or <id>/<name>
    #[arg(
        short = 'p',
        long = "process",
        value_name = "PROCESS",
        completion_process_selector(crate::completion::SelectorKind::Any)
    )]
    pub process: Option<ProcessRef>,

    /// The integrand name to use
    #[arg(
        short = 'i',
        long = "integrand-name",
        value_name = "NAME",
        completion_integrand_selector(crate::completion::SelectorKind::Any)
    )]
    pub integrand_name: Option<String>,

    /// Individual graph id, graph name, or inspect display label such as "#3 : graph_name"
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

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Build {
    #[command(flatten)]
    pub selection: GraphSelectorArgs,

    #[arg(long, alias = "family", value_enum, default_value_t = CliRepresentationMode::Cff)]
    pub representation: CliRepresentationMode,

    /// Edge-degree overrides in the form edge:degree,edge:degree.
    #[arg(long)]
    pub energy_degree_bounds: Option<String>,

    #[arg(long, value_enum, default_value_t = CliSamplingScaleMode::Auto)]
    pub numerator_sampling_scale_mode: CliSamplingScaleMode,

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

    #[arg(long, default_value_t = false)]
    pub clean: bool,

    #[arg(long, value_hint = clap::ValueHint::FilePath)]
    pub json_out: Option<PathBuf>,
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Evaluate {
    /// 3Drep artifact workspace containing the oriented expression JSON.
    #[arg(long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    #[arg(long, default_value_t = false)]
    pub clean: bool,
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct TestCffLtd {
    #[command(flatten)]
    pub selection: GraphSelectorArgs,

    /// Edge-degree overrides in the form edge:degree,edge:degree.
    #[arg(long)]
    pub energy_degree_bounds: Option<String>,

    #[arg(long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    #[arg(long, default_value_t = false)]
    pub clean: bool,
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct GraphFromSignatures {
    #[arg(
        long,
        conflicts_with = "signatures_file",
        required_unless_present = "signatures_file"
    )]
    pub signatures: Option<String>,

    #[arg(
        long,
        conflicts_with = "signatures",
        required_unless_present = "signatures"
    )]
    pub signatures_file: Option<PathBuf>,

    #[arg(long, default_value = "-")]
    pub dot_output: PathBuf,

    #[arg(long, default_value = "k")]
    pub loop_prefix: String,

    #[arg(long, default_value = "p,q")]
    pub external_prefixes: String,

    #[arg(long, default_value = "prop(q_,m_)")]
    pub prop_pattern: String,

    #[arg(long)]
    pub num_vertices: Option<usize>,

    #[arg(long)]
    pub allow_disconnected: bool,

    #[arg(long)]
    pub max_degree: Option<usize>,

    #[arg(long)]
    pub minimize_externals: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CliRepresentationMode {
    Ltd,
    Cff,
    PureLtd,
}

impl From<CliRepresentationMode> for RepresentationMode {
    fn from(value: CliRepresentationMode) -> Self {
        match value {
            CliRepresentationMode::Ltd => Self::Ltd,
            CliRepresentationMode::Cff => Self::Cff,
            CliRepresentationMode::PureLtd => Self::PureLtd,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CliSamplingScaleMode {
    Auto,
    None,
    BeyondQuadratic,
    All,
}

impl CliSamplingScaleMode {
    fn resolve(
        self,
        default_runtime_settings: &RuntimeSettings,
    ) -> three_dimensional_reps::NumeratorSamplingScaleMode {
        match self {
            Self::Auto
                if default_runtime_settings
                    .general
                    .numerator_interpolation_scale
                    .is_some() =>
            {
                three_dimensional_reps::NumeratorSamplingScaleMode::BeyondQuadratic
            }
            Self::Auto | Self::None => three_dimensional_reps::NumeratorSamplingScaleMode::None,
            Self::BeyondQuadratic => {
                three_dimensional_reps::NumeratorSamplingScaleMode::BeyondQuadratic
            }
            Self::All => three_dimensional_reps::NumeratorSamplingScaleMode::All,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct BuildOutput {
    backend: String,
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    family: RepresentationMode,
    graph: GraphInfo,
    validation: GraphValidation,
    automatic_energy_degree_bounds: Vec<(usize, usize)>,
    override_energy_degree_bounds: Vec<(usize, usize)>,
    energy_degree_bounds: Vec<(usize, usize)>,
    numerator_interpolation_scale: Option<f64>,
    numerator_sampling_scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    expression: ThreeDExpression<OrientationID>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ValidateOutput {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    graph: GraphInfo,
    validation: GraphValidation,
}

#[derive(Debug, Serialize, Deserialize)]
struct TestCffLtdOutput {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    automatic_energy_degree_bounds: Vec<(usize, usize)>,
    override_energy_degree_bounds: Vec<(usize, usize)>,
    energy_degree_bounds: Vec<(usize, usize)>,
    numerator_interpolation_scale: Option<f64>,
    sampled_m_values: Vec<f64>,
    cases: Vec<TestCffLtdCaseOutput>,
}

#[derive(Debug, Serialize, Deserialize)]
struct TestCffLtdCaseOutput {
    name: String,
    representation: RepresentationMode,
    numerator_sampling_scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    generation_status: String,
    orientation_count: usize,
    unfolded_term_count: usize,
    expression_path: PathBuf,
    symbolica_expression_path: PathBuf,
    evaluator_build_status: String,
    generation_error_path: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct EvaluateOutput {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    expression_path: PathBuf,
    symbolica_expression_path: PathBuf,
    param_builder_path: PathBuf,
}

struct SelectedGraph<'a> {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph: &'a Graph,
}

impl BuildOutput {
    fn matches_request(
        &self,
        selected: &SelectedGraph<'_>,
        family: RepresentationMode,
        energy_degree_bounds: &[(usize, usize)],
        sampling_scale: three_dimensional_reps::NumeratorSamplingScaleMode,
        numerator_interpolation_scale: Option<f64>,
    ) -> bool {
        self.process_id == selected.process_id
            && self.integrand_name == selected.integrand_name
            && self.graph_id == selected.graph_id
            && self.graph_name == selected.graph.name
            && self.family == family
            && self.energy_degree_bounds == energy_degree_bounds
            && self.numerator_sampling_scale_mode == sampling_scale
            && self.numerator_interpolation_scale == numerator_interpolation_scale
    }
}

impl TestCffLtdOutput {
    fn matches_request(
        &self,
        selected: &SelectedGraph<'_>,
        energy_degree_bounds: &[(usize, usize)],
        numerator_interpolation_scale: Option<f64>,
    ) -> bool {
        self.process_id == selected.process_id
            && self.integrand_name == selected.integrand_name
            && self.graph_id == selected.graph_id
            && self.graph_name == selected.graph.name
            && self.energy_degree_bounds == energy_degree_bounds
            && self.numerator_interpolation_scale == numerator_interpolation_scale
    }
}

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

        let mut matches = Vec::new();
        for graph_id in 0..self.graph_count() {
            let Some(name) = self.graph_name_by_id(graph_id) else {
                continue;
            };
            let display = format!("#{graph_id} : {name}");
            if selector == name || selector == display {
                matches.push(graph_id);
            }
        }

        match matches.as_slice() {
            [graph_id] => Ok(*graph_id),
            [] => Err(eyre!(
                "Could not resolve graph selector '{selector}'. Use a graph id, graph name, or inspect display name."
            )),
            matches => Err(eyre!(
                "Graph selector '{selector}' is ambiguous and matches graph ids {:?}. Use an explicit graph id.",
                matches
            )),
        }
    }
}

impl Validate {
    fn run(&self, state: &State) -> Result<()> {
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
        write_or_print(
            self.json_out.as_deref(),
            &serde_json::to_string_pretty(&output)?,
        )
    }
}

impl Build {
    fn run(
        &self,
        state: &State,
        global_cli_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        let selected = select_graph(state, &self.selection)?;
        let parsed = selected.graph.to_three_d_parsed_graph()?;
        let automatic_energy_degree_bounds =
            automatic_energy_degree_bounds_for_graph(selected.graph)?;
        let override_energy_degree_bounds =
            parse_energy_degree_bounds(self.energy_degree_bounds.as_deref())?;
        let energy_degree_bounds = merge_energy_degree_bounds(
            &automatic_energy_degree_bounds,
            &override_energy_degree_bounds,
        );
        let numerator_sampling_scale_mode = self
            .numerator_sampling_scale_mode
            .resolve(default_runtime_settings);
        let representation = RepresentationMode::from(self.representation);
        let workspace = self.workspace_path(global_cli_settings);
        let artifact_dir = build_artifact_dir(
            &workspace,
            &selected,
            representation,
            numerator_sampling_scale_mode,
        );
        let json_path = self
            .json_out
            .clone()
            .unwrap_or_else(|| artifact_dir.join("oriented_expression.json"));
        let numerator_interpolation_scale = default_runtime_settings
            .general
            .numerator_interpolation_scale;
        let cached_output = if !self.clean && json_path.exists() {
            let text = fs::read_to_string(&json_path)
                .with_context(|| format!("Could not read {}", json_path.display()))?;
            let cached = serde_json::from_str::<BuildOutput>(&text)
                .with_context(|| format!("Could not parse {}", json_path.display()))?;
            if cached.matches_request(
                &selected,
                representation,
                &energy_degree_bounds,
                numerator_sampling_scale_mode,
                numerator_interpolation_scale,
            ) {
                println!(
                    "Reusing cached 3Drep oriented expression from {}",
                    relative_display(&json_path)
                );
                Some(cached)
            } else {
                None
            }
        } else {
            None
        };

        let (output, reused_cached_output) = if let Some(output) = cached_output {
            (output, true)
        } else {
            let options = Generate3DExpressionOptions {
                representation,
                energy_degree_bounds: energy_degree_bounds.clone(),
                numerator_sampling_scale: numerator_sampling_scale_mode,
            };
            let expression = generate_3d_expression(selected.graph, &options)?;
            (
                BuildOutput {
                    backend: "gammaloop-3Drep".to_string(),
                    process_id: selected.process_id,
                    integrand_name: selected.integrand_name.clone(),
                    graph_id: selected.graph_id,
                    graph_name: selected.graph.name.clone(),
                    family: representation,
                    graph: graph_info(&parsed),
                    validation: validate_parsed_graph(&parsed),
                    automatic_energy_degree_bounds,
                    override_energy_degree_bounds,
                    energy_degree_bounds,
                    numerator_interpolation_scale,
                    numerator_sampling_scale_mode,
                    expression,
                },
                false,
            )
        };

        if !self.no_save_json {
            write_latest_expression_pointer(&workspace, &json_path)?;
            if !reused_cached_output {
                write_path(&json_path, &serde_json::to_string_pretty(&output)?)?;
                println!(
                    "Saved 3Drep oriented expression to {}",
                    relative_display(&json_path)
                );
            }
        }

        if !self.no_pretty || self.show_details_for_orientation.is_some() {
            let numerator_expr = selected.graph.full_numerator_atom().log_print(Some(120));
            println!(
                "{}",
                render_expression_summary(
                    &output.expression,
                    representation,
                    &output.graph,
                    &output.energy_degree_bounds,
                    Some(&numerator_expr),
                    numerator_sampling_scale_mode,
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

impl Evaluate {
    fn run(&self, state: &State, global_cli_settings: &CLISettings) -> Result<()> {
        let workspace = self
            .workspace_path
            .clone()
            .unwrap_or_else(|| default_workspace_path(global_cli_settings));
        let expression_path = if workspace.join("oriented_expression.json").exists() {
            workspace.join("oriented_expression.json")
        } else {
            read_latest_expression_path(&workspace)?
        };
        println!(
            "Loading 3Drep oriented expression from {}",
            relative_display(&expression_path)
        );
        let expression_text = fs::read_to_string(&expression_path).with_context(|| {
            format!(
                "Could not read 3Drep oriented expression JSON at {}",
                expression_path.display()
            )
        })?;
        let artifact: BuildOutput = serde_json::from_str(&expression_text).with_context(|| {
            format!(
                "Could not parse 3Drep oriented expression JSON at {}",
                expression_path.display()
            )
        })?;
        let selected = select_graph_from_artifact(state, &artifact)?;
        let parametric_atom = artifact
            .expression
            .parametric_atom_gs(selected.graph, &OrientationPattern::default());
        let artifact_dir = expression_path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| workspace.clone());
        let symbolica_expression_path = artifact_dir.join("symbolica_expression.txt");
        let param_builder_path = artifact_dir.join("param_builder.txt");
        let evaluate_manifest_path = artifact_dir.join("evaluate_manifest.json");
        if self.clean || !symbolica_expression_path.exists() {
            write_path(&symbolica_expression_path, &parametric_atom.log_print(None))?;
        } else {
            println!(
                "Reusing cached Symbolica expression from {}",
                relative_display(&symbolica_expression_path)
            );
        }
        if self.clean || !param_builder_path.exists() {
            write_path(
                &param_builder_path,
                &selected.graph.param_builder.to_string(),
            )?;
        } else {
            println!(
                "Reusing cached parameter builder from {}",
                relative_display(&param_builder_path)
            );
        }

        let orientations = artifact
            .expression
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect::<Vec<_>>();
        let evaluator_settings = threedrep_evaluator_settings();
        let _evaluator = EvaluatorStack::new(
            &[&parametric_atom],
            &selected.graph.param_builder,
            &orientations,
            None,
            &evaluator_settings,
        )?;

        let summary = EvaluateOutput {
            process_id: selected.process_id,
            integrand_name: selected.integrand_name,
            graph_id: selected.graph_id,
            graph_name: selected.graph.name.clone(),
            expression_path,
            symbolica_expression_path,
            param_builder_path,
        };
        write_path(
            &evaluate_manifest_path,
            &serde_json::to_string_pretty(&summary)?,
        )?;
        println!("{}", render_evaluate_summary(&summary));
        println!(
            "Saved 3Drep evaluate summary to {}",
            relative_display(&evaluate_manifest_path)
        );
        Ok(())
    }
}

impl TestCffLtd {
    fn run(
        &self,
        state: &State,
        global_cli_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        let selected = select_graph(state, &self.selection)?;
        let workspace = self
            .workspace_path
            .clone()
            .unwrap_or_else(|| default_workspace_path(global_cli_settings));
        let graph_workspace = graph_workspace_dir(&workspace, &selected);
        let manifest_path = graph_workspace.join("test_cff_ltd_manifest.json");
        let automatic_energy_degree_bounds =
            automatic_energy_degree_bounds_for_graph(selected.graph)?;
        let override_energy_degree_bounds =
            parse_energy_degree_bounds(self.energy_degree_bounds.as_deref())?;
        let energy_degree_bounds = merge_energy_degree_bounds(
            &automatic_energy_degree_bounds,
            &override_energy_degree_bounds,
        );
        let sampled_m_values = default_runtime_settings
            .general
            .numerator_interpolation_scale
            .map(|m| vec![m, 2.0 * m, 0.5 * m])
            .unwrap_or_default();

        let scale_modes = if default_runtime_settings
            .general
            .numerator_interpolation_scale
            .is_some()
        {
            vec![
                three_dimensional_reps::NumeratorSamplingScaleMode::None,
                three_dimensional_reps::NumeratorSamplingScaleMode::BeyondQuadratic,
                three_dimensional_reps::NumeratorSamplingScaleMode::All,
            ]
        } else {
            vec![three_dimensional_reps::NumeratorSamplingScaleMode::None]
        };

        let cached_output = if !self.clean && manifest_path.exists() {
            let text = fs::read_to_string(&manifest_path)
                .with_context(|| format!("Could not read {}", manifest_path.display()))?;
            let cached = serde_json::from_str::<TestCffLtdOutput>(&text)
                .with_context(|| format!("Could not parse {}", manifest_path.display()))?;
            if cached.matches_request(
                &selected,
                &energy_degree_bounds,
                default_runtime_settings
                    .general
                    .numerator_interpolation_scale,
            ) {
                println!(
                    "Reusing cached 3Drep comparison from {}",
                    relative_display(&manifest_path)
                );
                Some(cached)
            } else {
                None
            }
        } else {
            None
        };

        let (output, reused_cached_output) = if let Some(output) = cached_output {
            (output, true)
        } else {
            let mut cases = Vec::new();
            for representation in [
                RepresentationMode::Cff,
                RepresentationMode::Ltd,
                RepresentationMode::PureLtd,
            ] {
                for scale_mode in &scale_modes {
                    let case = self.build_case(
                        selected.graph,
                        &graph_workspace,
                        representation,
                        *scale_mode,
                        &energy_degree_bounds,
                    )?;
                    cases.push(case);
                }
            }

            (
                TestCffLtdOutput {
                    process_id: selected.process_id,
                    integrand_name: selected.integrand_name,
                    graph_id: selected.graph_id,
                    graph_name: selected.graph.name.clone(),
                    automatic_energy_degree_bounds,
                    override_energy_degree_bounds,
                    energy_degree_bounds,
                    numerator_interpolation_scale: default_runtime_settings
                        .general
                        .numerator_interpolation_scale,
                    sampled_m_values,
                    cases,
                },
                false,
            )
        };
        if !reused_cached_output {
            write_path(&manifest_path, &serde_json::to_string_pretty(&output)?)?;
        }
        println!("{}", render_test_cff_ltd_summary(&output, &manifest_path));
        if !reused_cached_output {
            println!(
                "Saved 3Drep comparison manifest to {}",
                relative_display(&manifest_path)
            );
        }
        Ok(())
    }

    fn build_case(
        &self,
        graph: &Graph,
        workspace: &Path,
        representation: RepresentationMode,
        scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
        energy_degree_bounds: &[(usize, usize)],
    ) -> Result<TestCffLtdCaseOutput> {
        let name = format!("{representation:?}_{scale_mode:?}").to_lowercase();
        let case_dir = workspace.join("test_cff_ltd").join(&name);
        let expression_path = case_dir.join("oriented_expression.json");
        let symbolica_expression_path = case_dir.join("symbolica_expression.txt");
        let expression = match generate_3d_expression(
            graph,
            &Generate3DExpressionOptions {
                representation,
                energy_degree_bounds: energy_degree_bounds.to_vec(),
                numerator_sampling_scale: scale_mode,
            },
        ) {
            Ok(expression) => expression,
            Err(error) => {
                let generation_error_path = case_dir.join("generation_error.txt");
                write_path(&generation_error_path, &format!("{error:?}"))?;
                return Ok(TestCffLtdCaseOutput {
                    name,
                    representation,
                    numerator_sampling_scale_mode: scale_mode,
                    generation_status: "failed".to_string(),
                    orientation_count: 0,
                    unfolded_term_count: 0,
                    expression_path,
                    symbolica_expression_path,
                    evaluator_build_status: "skipped: generation failed".to_string(),
                    generation_error_path: Some(generation_error_path),
                });
            }
        };
        let parametric_atom = expression.parametric_atom_gs(graph, &OrientationPattern::default());
        write_path(
            &expression_path,
            &serde_json::to_string_pretty(&expression)?,
        )?;
        write_path(&symbolica_expression_path, &parametric_atom.log_print(None))?;

        let orientations = expression
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect::<Vec<_>>();
        let evaluator_settings = threedrep_evaluator_settings();
        let evaluator_build_status = match EvaluatorStack::new(
            &[&parametric_atom],
            &graph.param_builder,
            &orientations,
            None,
            &evaluator_settings,
        ) {
            Ok(_) => "ok".to_string(),
            Err(error) => format!("failed: {error:?}"),
        };

        Ok(TestCffLtdCaseOutput {
            name,
            representation,
            numerator_sampling_scale_mode: scale_mode,
            generation_status: "ok".to_string(),
            orientation_count: expression.orientations.len(),
            unfolded_term_count: expression.num_unfolded_terms(),
            expression_path,
            symbolica_expression_path,
            evaluator_build_status,
            generation_error_path: None,
        })
    }
}

impl GraphFromSignatures {
    fn run(&self) -> Result<()> {
        let expression = if let Some(path) = &self.signatures_file {
            fs::read_to_string(path)
                .with_context(|| format!("Could not read signatures from {}", path.display()))?
        } else {
            self.signatures
                .clone()
                .expect("clap ensures --signatures is present")
        };
        let external_prefixes = self
            .external_prefixes
            .split(',')
            .map(str::trim)
            .filter(|prefix| !prefix.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        let options = ReconstructDotOptions {
            num_vertices: self.num_vertices,
            require_connected: !self.allow_disconnected,
            max_degree: self.max_degree,
            format: ReconstructDotFormat::Gammaloop,
            graph_engine: None,
            minimize_external_legs: true,
        };
        let (dot, _) = reconstruct_dot_from_expression(
            &expression,
            &self.loop_prefix,
            &external_prefixes,
            &self.prop_pattern,
            &options,
        )?;
        if self.dot_output == Path::new("-") {
            print!("{dot}");
            Ok(())
        } else {
            write_path(&self.dot_output, &dot)
        }
    }
}

fn threedrep_evaluator_settings() -> EvaluatorSettings {
    EvaluatorSettings {
        iterative_orientation_optimization: false,
        ..EvaluatorSettings::default()
    }
}

fn automatic_energy_degree_bounds_for_graph(graph: &Graph) -> Result<Vec<(usize, usize)>> {
    let mut bounds = graph
        .underlying
        .iter_edges()
        .filter(|(pair, _, _)| matches!(pair, HedgePair::Paired { .. }))
        .map(|(_, edge, _)| (edge.0, 0usize))
        .collect::<BTreeMap<_, _>>();

    for (edge, degree) in graph.numerator_energy_power_caps()?.iter() {
        bounds.insert(edge.0, degree);
    }

    Ok(bounds.into_iter().collect())
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

fn select_graph_from_artifact<'a>(
    state: &'a State,
    artifact: &BuildOutput,
) -> Result<SelectedGraph<'a>> {
    let process = state
        .process_list
        .processes
        .get(artifact.process_id)
        .ok_or_else(|| {
            eyre!(
                "3Drep artifact references process id {}, but the active state has only {} process(es).",
                artifact.process_id,
                state.process_list.processes.len()
            )
        })?;
    let catalog =
        GraphCatalog::for_integrand(process, &artifact.integrand_name).with_context(|| {
            format!(
                "while resolving integrand '{}' from 3Drep artifact",
                artifact.integrand_name
            )
        })?;
    let graph = catalog.graph_by_id(artifact.graph_id)?;
    Ok(SelectedGraph {
        process_id: artifact.process_id,
        integrand_name: artifact.integrand_name.clone(),
        graph_id: artifact.graph_id,
        graph,
    })
}

fn parse_graph_id(selector: &str) -> Option<usize> {
    let stripped = selector.trim().strip_prefix('#').unwrap_or(selector).trim();
    let id_part = stripped
        .split_once(':')
        .map_or(stripped, |(id, _)| id)
        .trim();
    id_part.parse::<usize>().ok()
}

fn parse_energy_degree_bounds(value: Option<&str>) -> Result<Vec<(usize, usize)>> {
    let Some(value) = value else {
        return Ok(Vec::new());
    };
    if value.trim().is_empty() {
        return Ok(Vec::new());
    }

    value
        .split(',')
        .map(|item| {
            let (edge, degree) = item.split_once(':').ok_or_else(|| {
                eyre!("Invalid energy-degree item '{item}'; expected edge:degree")
            })?;
            let edge = edge
                .trim()
                .parse::<usize>()
                .map_err(|_| eyre!("Invalid edge id in energy-degree item '{item}'"))?;
            let degree = degree
                .trim()
                .parse::<usize>()
                .map_err(|_| eyre!("Invalid degree in energy-degree item '{item}'"))?;
            Ok((edge, degree))
        })
        .collect()
}

fn merge_energy_degree_bounds(
    automatic: &[(usize, usize)],
    overrides: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let mut merged = automatic.iter().copied().collect::<BTreeMap<_, _>>();
    for (edge, degree) in overrides {
        merged.insert(*edge, *degree);
    }
    merged.into_iter().collect()
}

fn default_workspace_path(global_cli_settings: &CLISettings) -> PathBuf {
    if global_cli_settings.session.read_only_state {
        let workspace_name = global_cli_settings
            .state
            .name
            .as_deref()
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .map(|name| format!("threed_workspace_{name}"))
            .unwrap_or_else(|| "threed_workspace".to_string());
        PathBuf::from(".").join(workspace_name)
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
    representation: RepresentationMode,
    scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
) -> PathBuf {
    graph_workspace_dir(workspace, selected)
        .join("build")
        .join(format!("{representation:?}_{scale_mode:?}").to_lowercase())
}

fn latest_expression_pointer(workspace: &Path) -> PathBuf {
    workspace.join("latest_oriented_expression_path.txt")
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
        &latest_expression_pointer(workspace),
        &relative_display(expression_path),
    )
}

fn read_latest_expression_path(workspace: &Path) -> Result<PathBuf> {
    let pointer = latest_expression_pointer(workspace);
    let value = fs::read_to_string(&pointer)
        .with_context(|| format!("Could not read {}", pointer.display()))?;
    let path = PathBuf::from(value.trim());
    if path.is_absolute() {
        Ok(path)
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
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

fn color_text(text: impl AsRef<str>, color: Color) -> String {
    color.paint(text.as_ref()).to_string()
}

fn table_header(text: impl AsRef<str>) -> String {
    AnsiStyle::new()
        .bold()
        .paint(Color::Cyan.paint(text.as_ref()).to_string())
        .to_string()
}

fn status_text(text: &str) -> String {
    let color = if text == "ok" {
        Color::Green
    } else if text.starts_with("failed") {
        Color::Red
    } else {
        Color::Yellow
    };
    color_text(text, color)
}

fn render_evaluate_summary(output: &EvaluateOutput) -> String {
    let mut table = Builder::new();
    table.push_record(vec![table_header("field"), table_header("value")]);
    table.push_record(vec![
        "process".to_string(),
        color_text(output.process_id.to_string(), Color::Blue),
    ]);
    table.push_record(vec![
        "integrand".to_string(),
        color_text(&output.integrand_name, Color::Green),
    ]);
    table.push_record(vec![
        "graph".to_string(),
        format!(
            "{} {}",
            color_text(format!("#{}", output.graph_id), Color::Blue),
            color_text(&output.graph_name, Color::Green)
        ),
    ]);
    table.push_record(vec![
        "oriented expression".to_string(),
        color_text(relative_display(&output.expression_path), Color::Purple),
    ]);
    table.push_record(vec![
        "symbolica expression".to_string(),
        color_text(
            relative_display(&output.symbolica_expression_path),
            Color::Purple,
        ),
    ]);
    table.push_record(vec![
        "param builder".to_string(),
        color_text(relative_display(&output.param_builder_path), Color::Purple),
    ]);
    table.build().with(Style::rounded()).to_string()
}

fn render_test_cff_ltd_summary(output: &TestCffLtdOutput, manifest_path: &Path) -> String {
    let mut table = Builder::new();
    table.push_record(vec![
        table_header("case"),
        table_header("rep"),
        table_header("M mode"),
        table_header("generation"),
        table_header("orientations"),
        table_header("terms"),
        table_header("evaluator"),
    ]);
    for case in &output.cases {
        table.push_record(vec![
            color_text(&case.name, Color::Blue),
            color_text(format!("{:?}", case.representation), Color::Green),
            color_text(
                format!("{:?}", case.numerator_sampling_scale_mode),
                Color::Purple,
            ),
            status_text(&case.generation_status),
            color_text(case.orientation_count.to_string(), Color::Blue),
            color_text(case.unfolded_term_count.to_string(), Color::Blue),
            status_text(&case.evaluator_build_status),
        ]);
    }
    format!(
        "{} for process {}, integrand {}, graph {} {}\n{}: {}\n{}",
        color_text("3Drep comparison", Color::Cyan),
        color_text(format!("#{}", output.process_id), Color::Blue),
        color_text(format!("'{}'", output.integrand_name), Color::Green),
        color_text(format!("#{}", output.graph_id), Color::Blue),
        color_text(&output.graph_name, Color::Green),
        color_text("manifest", Color::Cyan),
        color_text(relative_display(manifest_path), Color::Purple),
        table.build().with(Style::rounded())
    )
}
