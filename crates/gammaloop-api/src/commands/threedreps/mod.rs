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
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
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

    #[arg(long, value_hint = clap::ValueHint::FilePath)]
    pub json_out: Option<PathBuf>,
}

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub struct Evaluate {
    /// 3Drep artifact workspace containing the oriented expression JSON.
    #[arg(long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,
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

#[derive(Debug, Serialize)]
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

#[derive(Debug, Serialize)]
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

struct SelectedGraph<'a> {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph: &'a Graph,
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
        let options = Generate3DExpressionOptions {
            representation,
            energy_degree_bounds: energy_degree_bounds.clone(),
            numerator_sampling_scale: numerator_sampling_scale_mode,
        };
        let expression = generate_3d_expression(selected.graph, &options)?;
        let output = BuildOutput {
            backend: "gammaloop-3Drep".to_string(),
            process_id: selected.process_id,
            integrand_name: selected.integrand_name,
            graph_id: selected.graph_id,
            graph_name: selected.graph.name.clone(),
            family: representation,
            graph: graph_info(&parsed),
            validation: validate_parsed_graph(&parsed),
            automatic_energy_degree_bounds,
            override_energy_degree_bounds,
            energy_degree_bounds,
            numerator_interpolation_scale: default_runtime_settings
                .general
                .numerator_interpolation_scale,
            numerator_sampling_scale_mode,
            expression,
        };

        if !self.no_save_json {
            let json_path = self.json_out.clone().unwrap_or_else(|| {
                self.workspace_path(global_cli_settings)
                    .join("oriented_expression.json")
            });
            write_path(&json_path, &serde_json::to_string_pretty(&output)?)?;
        }

        if !self.no_pretty || self.show_details_for_orientation.is_some() {
            println!(
                "{}",
                render_expression_summary(
                    &output.expression,
                    representation,
                    &output.graph,
                    &output.energy_degree_bounds,
                    None,
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
        let expression_path = workspace.join("oriented_expression.json");
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
        write_path(
            &workspace.join("symbolica_expression.txt"),
            &parametric_atom.log_print(None),
        )?;
        write_path(
            &workspace.join("param_builder.txt"),
            &selected.graph.param_builder.to_string(),
        )?;

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

        println!(
            "Built 3Drep evaluator input for process #{}, integrand '{}', graph #{} in {}",
            selected.process_id,
            selected.integrand_name,
            selected.graph_id,
            workspace.display()
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

        let mut cases = Vec::new();
        for representation in [
            RepresentationMode::Cff,
            RepresentationMode::Ltd,
            RepresentationMode::PureLtd,
        ] {
            for scale_mode in &scale_modes {
                let case = self.build_case(
                    selected.graph,
                    &workspace,
                    representation,
                    *scale_mode,
                    &energy_degree_bounds,
                )?;
                cases.push(case);
            }
        }

        let output = TestCffLtdOutput {
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
        };
        write_path(
            &workspace.join("test_cff_ltd_manifest.json"),
            &serde_json::to_string_pretty(&output)?,
        )?;
        println!(
            "Built {} 3Drep comparison case(s) in {}",
            output.cases.len(),
            workspace.display()
        );
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
            minimize_external_legs: self.minimize_externals,
        };
        let (dot, _) = reconstruct_dot_from_expression(
            &expression,
            &self.loop_prefix,
            &external_prefixes,
            &self.prop_pattern,
            &options,
        )?;
        if self.dot_output == PathBuf::from("-") {
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
