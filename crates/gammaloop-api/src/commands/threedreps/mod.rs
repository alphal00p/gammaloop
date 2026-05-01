use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use clap::{Args, Subcommand, ValueEnum};
use color_eyre::{
    eyre::{eyre, Context},
    Result,
};
use gammalooprs::{
    cff::expression::GammaLoopThreeDExpression,
    graph::{Graph, ThreeDRepMassShift},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            evaluators::{EvaluatorMethod, EvaluatorStack, InputParams, SingleOrAllOrientations},
            param_builder::{ParamBuilder, ParamBuilderInputGroup},
            ProcessIntegrand,
        },
    },
    model::Model,
    processes::{Amplitude, CrossSection, EvaluatorSettings, Process, ProcessCollection},
    settings::{
        global::{FrozenCompilationMode, OrientationPattern, UniformNumeratorSamplingScale},
        runtime::Precision,
        RuntimeSettings,
    },
    utils::{f128, symbolica_ext::LogPrint, ArbPrec, FloatLike, F},
};
use linnet::half_edge::{
    involution::{EdgeVec, HedgePair, Orientation},
    subgraph::subset::SubSet,
};
use nu_ansi_term::{Color, Style as AnsiStyle};
use rand::{rngs::SmallRng, Rng, SeedableRng};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::atom::AtomCore;
use tabled::{builder::Builder, settings::Style};
use three_dimensional_reps::{
    generate_3d_expression, graph_info, reconstruct_dot_from_expression, render_expression_summary,
    validate_parsed_graph, DisplayOptions, Generate3DExpressionOptions, GraphInfo, GraphValidation,
    OrientationID, ReconstructDotFormat, ReconstructDotOptions, RepresentationMode,
    ThreeDExpression, ThreeDGraphSource,
};
use typed_index_collections::TiVec;

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
            Self::Evaluate(command) => {
                command.run(state, global_cli_settings, default_runtime_settings)
            }
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

    #[arg(
        long = "numerator-samples-normalization",
        alias = "numerator-sampling-scale-mode",
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

    #[arg(long, value_enum)]
    pub precision: Option<CliRuntimePrecision>,

    #[arg(long, default_value_t = 1)]
    pub seed: u64,

    #[arg(long, default_value_t = 1.0)]
    pub scale: f64,

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

    #[arg(long, value_enum)]
    pub precision: Option<CliRuntimePrecision>,

    #[arg(long, default_value_t = 1)]
    pub seed: u64,

    #[arg(long, default_value_t = 1.0)]
    pub scale: f64,

    #[arg(long)]
    pub mass_shift: Option<f64>,

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
    fn resolve_from_global(global_cli_settings: &CLISettings) -> Self {
        match global_cli_settings
            .global
            .generation
            .uniform_numerator_sampling_scale
        {
            UniformNumeratorSamplingScale::None => Self::NeverM,
            UniformNumeratorSamplingScale::BeyondQuadratic => Self::MForBeyondQuadraticOnly,
            UniformNumeratorSamplingScale::All => Self::MForAll,
        }
    }

    fn to_generation_mode(self) -> three_dimensional_reps::NumeratorSamplingScaleMode {
        match self {
            Self::NeverM => three_dimensional_reps::NumeratorSamplingScaleMode::None,
            Self::MForAll => three_dimensional_reps::NumeratorSamplingScaleMode::All,
            Self::MForBeyondQuadraticOnly => {
                three_dimensional_reps::NumeratorSamplingScaleMode::BeyondQuadratic
            }
        }
    }

    fn from_generation_mode(mode: three_dimensional_reps::NumeratorSamplingScaleMode) -> Self {
        match mode {
            three_dimensional_reps::NumeratorSamplingScaleMode::None => Self::NeverM,
            three_dimensional_reps::NumeratorSamplingScaleMode::All => Self::MForAll,
            three_dimensional_reps::NumeratorSamplingScaleMode::BeyondQuadratic => {
                Self::MForBeyondQuadraticOnly
            }
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::NeverM => "never_M",
            Self::MForAll => "M_for_all",
            Self::MForBeyondQuadraticOnly => "M_for_beyond_quadratic_only",
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub enum CliRuntimePrecision {
    #[value(name = "Double", alias = "double")]
    Double,
    #[value(name = "Quad", alias = "quad")]
    Quad,
    #[value(name = "ArbPrec", alias = "arb", alias = "arb-prec")]
    ArbPrec,
}

impl CliRuntimePrecision {
    fn resolve(default_runtime_settings: &RuntimeSettings, cli_precision: Option<Self>) -> Self {
        cli_precision.unwrap_or_else(|| {
            default_runtime_settings
                .stability
                .levels
                .first()
                .map(|level| match level.precision {
                    Precision::Double => Self::Double,
                    Precision::Quad => Self::Quad,
                    Precision::Arb => Self::ArbPrec,
                })
                .unwrap_or(Self::Double)
        })
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
    numerator_interpolation_scale: f64,
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
    settings: ThreeDrepRunSettings,
    automatic_energy_degree_bounds: Vec<(usize, usize)>,
    override_energy_degree_bounds: Vec<(usize, usize)>,
    energy_degree_bounds: Vec<(usize, usize)>,
    numerator_interpolation_scale: f64,
    mass_shift_start: f64,
    cases: Vec<TestCffLtdCaseOutput>,
    #[serde(default)]
    verdict: TestCffLtdVerdict,
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
    evaluator_build_timing: Option<String>,
    evaluator_build_timing_seconds: Option<f64>,
    #[serde(default)]
    evaluations: Vec<ThreeDrepEvaluationRecord>,
    generation_error_path: Option<PathBuf>,
}

#[derive(Debug, Serialize, Deserialize)]
struct EvaluateOutput {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph_name: String,
    expression_path: PathBuf,
    symbolica_expression_path: PathBuf,
    param_builder_path: PathBuf,
    settings: ThreeDrepRunSettings,
    evaluation: ThreeDrepEvaluationRecord,
    parameters: Vec<ParameterValueRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ThreeDrepEvaluationRecord {
    #[serde(default)]
    id: usize,
    label: String,
    representation: Option<RepresentationMode>,
    numerator_sampling_scale_mode: Option<three_dimensional_reps::NumeratorSamplingScaleMode>,
    mass_shift: String,
    #[serde(default)]
    mass_shift_values: Vec<MassShiftValueRecord>,
    numerator_interpolation_scale: f64,
    runtime_precision: CliRuntimePrecision,
    evaluator_method: EvaluatorMethod,
    evaluator_backend: String,
    evaluator_build_timing: Option<String>,
    evaluator_build_timing_seconds: Option<f64>,
    value: String,
    value_re: String,
    value_im: String,
    sample_evaluation_timing: String,
    sample_evaluation_timing_seconds: f64,
    status: String,
    error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct MassShiftValueRecord {
    group_index: usize,
    split_index: usize,
    local_edge_id: usize,
    edge_id: usize,
    base_mass: String,
    shifted_mass: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
struct TestCffLtdVerdict {
    status: TestCffLtdStatus,
    direct_relative_tolerance: String,
    direct_absolute_tolerance: String,
    mass_shift_relative_tolerance: String,
    mass_shift_absolute_tolerance: String,
    checks: Vec<TestCffLtdCheckRecord>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
enum TestCffLtdStatus {
    Success,
    Fail,
    #[default]
    Unknown,
}

impl TestCffLtdStatus {
    fn label(self) -> &'static str {
        match self {
            Self::Success => "SUCCESS",
            Self::Fail => "FAIL",
            Self::Unknown => "UNKNOWN",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TestCffLtdCheckRecord {
    name: String,
    evaluation_id: String,
    status: TestCffLtdStatus,
    lhs: String,
    rhs: String,
    abs_diff: String,
    rel_diff: String,
    tolerance: String,
    message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ParameterValueRecord {
    name: String,
    canonical_name: String,
    value: String,
    source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ThreeDrepRunSettings {
    runtime_precision: CliRuntimePrecision,
    evaluator_method: EvaluatorMethod,
    evaluator_backend: String,
    compiled_backend_available: bool,
    numerator_samples_normalization: CliNumeratorSamplesNormalization,
    numerator_interpolation_scale: f64,
    seed: u64,
    scale: f64,
}

struct SelectedGraph<'a> {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph: &'a Graph,
}

struct TestCffLtdCaseBuildRequest<'a> {
    graph: &'a Graph,
    model: &'a Model,
    workspace: &'a Path,
    representation: RepresentationMode,
    scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    energy_degree_bounds: &'a [(usize, usize)],
    run_settings: &'a ThreeDrepRunSettings,
    default_runtime_settings: &'a RuntimeSettings,
    global_cli_settings: &'a CLISettings,
    mass_shift_start: f64,
}

struct StandardEvaluationRequest<'a> {
    name: &'a str,
    graph: &'a Graph,
    expression: &'a ThreeDExpression<OrientationID>,
    parametric_atom: &'a symbolica::atom::Atom,
    default_runtime_settings: &'a RuntimeSettings,
    global_cli_settings: &'a CLISettings,
    workspace: &'a Path,
    representation: RepresentationMode,
    scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    run_settings: &'a ThreeDrepRunSettings,
    orientations: &'a TiVec<OrientationID, EdgeVec<Orientation>>,
}

struct PureLtdMassShiftEvaluationRequest<'a> {
    name: &'a str,
    graph: &'a Graph,
    model: &'a Model,
    energy_degree_bounds: &'a [(usize, usize)],
    default_runtime_settings: &'a RuntimeSettings,
    global_cli_settings: &'a CLISettings,
    workspace: &'a Path,
    scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    run_settings: &'a ThreeDrepRunSettings,
    mass_shift_start: f64,
}

#[derive(Clone, Copy)]
struct EvaluatorBuildContext<'a> {
    workspace: &'a Path,
    global_cli_settings: &'a CLISettings,
    default_runtime_settings: &'a RuntimeSettings,
    run_settings: &'a ThreeDrepRunSettings,
}

struct PreparedEvaluator {
    evaluator: EvaluatorStack,
    info: PreparedEvaluatorInfo,
}

#[derive(Clone)]
struct PreparedEvaluatorInfo {
    build_timing: String,
    build_timing_seconds: f64,
    backend: String,
}

impl BuildOutput {
    fn matches_request(
        &self,
        selected: &SelectedGraph<'_>,
        family: RepresentationMode,
        energy_degree_bounds: &[(usize, usize)],
        sampling_scale: three_dimensional_reps::NumeratorSamplingScaleMode,
        numerator_interpolation_scale: f64,
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
    fn assign_evaluation_ids(&mut self) {
        for (id, evaluation) in self
            .cases
            .iter_mut()
            .flat_map(|case| case.evaluations.iter_mut())
            .enumerate()
        {
            evaluation.id = id;
        }
    }

    fn has_ordered_evaluation_ids(&self) -> bool {
        self.cases
            .iter()
            .flat_map(|case| &case.evaluations)
            .enumerate()
            .all(|(expected, evaluation)| evaluation.id == expected)
    }

    fn matches_request(
        &self,
        selected: &SelectedGraph<'_>,
        energy_degree_bounds: &[(usize, usize)],
        run_settings: &ThreeDrepRunSettings,
        mass_shift_start: f64,
    ) -> bool {
        self.process_id == selected.process_id
            && self.integrand_name == selected.integrand_name
            && self.graph_id == selected.graph_id
            && self.graph_name == selected.graph.name
            && self.energy_degree_bounds == energy_degree_bounds
            && self.numerator_interpolation_scale == run_settings.numerator_interpolation_scale
            && self.mass_shift_start == mass_shift_start
            && self.settings.runtime_precision == run_settings.runtime_precision
            && self.settings.evaluator_method == run_settings.evaluator_method
            && self.settings.evaluator_backend == run_settings.evaluator_backend
            && self.settings.numerator_samples_normalization
                == run_settings.numerator_samples_normalization
            && self.settings.seed == run_settings.seed
            && self.settings.scale == run_settings.scale
            && self.verdict.status != TestCffLtdStatus::Unknown
            && self.has_ordered_evaluation_ids()
            && self.cases.iter().all(|case| {
                case.generation_status != "ok"
                    || !case.evaluator_build_status.starts_with("ok")
                    || !case.evaluations.is_empty()
            })
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
        let numerator_samples_normalization =
            self.numerator_samples_normalization.unwrap_or_else(|| {
                CliNumeratorSamplesNormalization::resolve_from_global(global_cli_settings)
            });
        let numerator_sampling_scale_mode = numerator_samples_normalization.to_generation_mode();
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
        let numerator_interpolation_scale =
            definite_numerator_interpolation_scale(default_runtime_settings);
        let cached_output = if !self.clean && json_path.exists() {
            let text = fs::read_to_string(&json_path)
                .with_context(|| format!("Could not read {}", json_path.display()))?;
            let cached = serde_json::from_str::<BuildOutput>(&text).ok();
            if let Some(cached) = cached.filter(|cached| {
                cached.matches_request(
                    &selected,
                    representation,
                    &energy_degree_bounds,
                    numerator_sampling_scale_mode,
                    numerator_interpolation_scale,
                )
            }) {
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
    fn run(
        &self,
        state: &State,
        global_cli_settings: &CLISettings,
        default_runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
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
            .diagnostic_parametric_atom_gs(selected.graph, &OrientationPattern::default());
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

        let run_settings = threedrep_run_settings(
            global_cli_settings,
            default_runtime_settings,
            self.precision,
            self.seed,
            self.scale,
            CliNumeratorSamplesNormalization::from_generation_mode(
                artifact.numerator_sampling_scale_mode,
            ),
        )?;
        let evaluation = evaluate_threedrep_expression(EvaluationRequest {
            label: "evaluate",
            graph: selected.graph,
            expression: &artifact.expression,
            parametric_atom: &parametric_atom,
            workspace: &artifact_dir,
            global_cli_settings,
            default_runtime_settings,
            representation: Some(artifact.family),
            numerator_sampling_scale_mode: Some(artifact.numerator_sampling_scale_mode),
            run_settings: &run_settings,
            mass_shift: "none",
            mass_shift_values: &[],
        })?;
        let parameters = parameter_records(selected.graph, &run_settings, default_runtime_settings);

        let summary = EvaluateOutput {
            process_id: selected.process_id,
            integrand_name: selected.integrand_name,
            graph_id: selected.graph_id,
            graph_name: selected.graph.name.clone(),
            expression_path,
            symbolica_expression_path,
            param_builder_path,
            settings: run_settings,
            evaluation,
            parameters,
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
        let model =
            state.resolve_model_for_integrand(selected.process_id, &selected.integrand_name)?;
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
        let run_settings = threedrep_run_settings(
            global_cli_settings,
            default_runtime_settings,
            self.precision,
            self.seed,
            self.scale,
            CliNumeratorSamplesNormalization::resolve_from_global(global_cli_settings),
        )?;
        let mass_shift_start = self.mass_shift.unwrap_or(self.scale);

        let cached_output = if !self.clean && manifest_path.exists() {
            let text = fs::read_to_string(&manifest_path)
                .with_context(|| format!("Could not read {}", manifest_path.display()))?;
            let cached = serde_json::from_str::<TestCffLtdOutput>(&text).ok();
            if let Some(cached) = cached.filter(|cached| {
                cached.matches_request(
                    &selected,
                    &energy_degree_bounds,
                    &run_settings,
                    mass_shift_start,
                )
            }) {
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
                let case = self.build_case(TestCffLtdCaseBuildRequest {
                    graph: selected.graph,
                    model: &model,
                    workspace: &graph_workspace,
                    representation,
                    scale_mode: run_settings
                        .numerator_samples_normalization
                        .to_generation_mode(),
                    energy_degree_bounds: &energy_degree_bounds,
                    run_settings: &run_settings,
                    default_runtime_settings,
                    global_cli_settings,
                    mass_shift_start,
                })?;
                cases.push(case);
            }

            let mut output = TestCffLtdOutput {
                process_id: selected.process_id,
                integrand_name: selected.integrand_name,
                graph_id: selected.graph_id,
                graph_name: selected.graph.name.clone(),
                settings: run_settings.clone(),
                automatic_energy_degree_bounds,
                override_energy_degree_bounds,
                energy_degree_bounds,
                numerator_interpolation_scale: run_settings.numerator_interpolation_scale,
                mass_shift_start,
                cases,
                verdict: TestCffLtdVerdict::default(),
            };
            output.assign_evaluation_ids();
            output.verdict = build_test_cff_ltd_verdict(&output);

            (output, false)
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

    fn build_case(&self, request: TestCffLtdCaseBuildRequest<'_>) -> Result<TestCffLtdCaseOutput> {
        let TestCffLtdCaseBuildRequest {
            graph,
            model,
            workspace,
            representation,
            scale_mode,
            energy_degree_bounds,
            run_settings,
            default_runtime_settings,
            global_cli_settings,
            mass_shift_start,
        } = request;
        let name = representation_label(representation).to_lowercase();
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
                    evaluator_build_timing: None,
                    evaluator_build_timing_seconds: None,
                    evaluations: Vec::new(),
                    generation_error_path: Some(generation_error_path),
                });
            }
        };
        let parametric_atom =
            expression.diagnostic_parametric_atom_gs(graph, &OrientationPattern::default());
        write_path(
            &expression_path,
            &serde_json::to_string_pretty(&expression)?,
        )?;
        write_path(&symbolica_expression_path, &parametric_atom.log_print(None))?;

        let orientations = diagnostic_evaluation_orientations(&expression);
        let has_repeated_propagators =
            !three_dimensional_reps::repeated_groups(&graph.to_three_d_parsed_graph()?).is_empty();
        let (evaluator_build_status, evaluations) =
            if representation == RepresentationMode::PureLtd && has_repeated_propagators {
                self.build_pure_ltd_mass_shift_evaluations(PureLtdMassShiftEvaluationRequest {
                    name: &name,
                    graph,
                    model,
                    workspace: &case_dir,
                    energy_degree_bounds,
                    default_runtime_settings,
                    global_cli_settings,
                    scale_mode,
                    run_settings,
                    mass_shift_start,
                })?
            } else {
                self.build_standard_evaluations(StandardEvaluationRequest {
                    name: &name,
                    graph,
                    expression: &expression,
                    parametric_atom: &parametric_atom,
                    default_runtime_settings,
                    global_cli_settings,
                    workspace: &case_dir,
                    representation,
                    scale_mode,
                    run_settings,
                    orientations: &orientations,
                })?
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
            evaluator_build_timing: evaluations
                .first()
                .and_then(|evaluation| evaluation.evaluator_build_timing.clone()),
            evaluator_build_timing_seconds: evaluations
                .first()
                .and_then(|evaluation| evaluation.evaluator_build_timing_seconds),
            evaluations,
            generation_error_path: None,
        })
    }

    fn build_standard_evaluations(
        &self,
        request: StandardEvaluationRequest<'_>,
    ) -> Result<(String, Vec<ThreeDrepEvaluationRecord>)> {
        let prepared = match build_threedrep_evaluator(
            request.name,
            EvaluatorBuildContext {
                workspace: request.workspace,
                global_cli_settings: request.global_cli_settings,
                default_runtime_settings: request.default_runtime_settings,
                run_settings: request.run_settings,
            },
            request.graph,
            &[request.parametric_atom],
            request.orientations,
        ) {
            Ok(prepared) => prepared,
            Err(error) => return Ok((format!("failed: {error:?}"), Vec::new())),
        };
        let evaluator_info = prepared.info.clone();
        let mut evaluator = prepared.evaluator;

        let evaluations = vec![evaluate_with_evaluator(
            &mut evaluator,
            EvaluationRequest {
                label: request.name,
                graph: request.graph,
                expression: request.expression,
                parametric_atom: request.parametric_atom,
                workspace: request.workspace,
                global_cli_settings: request.global_cli_settings,
                default_runtime_settings: request.default_runtime_settings,
                representation: Some(request.representation),
                numerator_sampling_scale_mode: Some(request.scale_mode),
                run_settings: request.run_settings,
                mass_shift: "none",
                mass_shift_values: &[],
            },
            request.orientations,
            &evaluator_info,
        )?];
        Ok(("ok".to_string(), evaluations))
    }

    fn build_pure_ltd_mass_shift_evaluations(
        &self,
        request: PureLtdMassShiftEvaluationRequest<'_>,
    ) -> Result<(String, Vec<ThreeDrepEvaluationRecord>)> {
        let mut evaluations = Vec::new();
        for epsilon in pure_ltd_mass_shift_epsilons(request.mass_shift_start) {
            let (shifted_graph, mass_shifts) = request
                .graph
                .split_repeated_masses_for_three_drep(request.model, epsilon)?;
            let mass_shift_values = mass_shift_value_records(&mass_shifts);
            let mass_shift = mass_shift_label(epsilon, &mass_shift_values);
            let shifted_expression = generate_3d_expression(
                &shifted_graph,
                &Generate3DExpressionOptions {
                    representation: RepresentationMode::Ltd,
                    energy_degree_bounds: request.energy_degree_bounds.to_vec(),
                    numerator_sampling_scale: request.scale_mode,
                },
            )
            .with_context(|| {
                format!(
                    "while generating split-mass LTD expression for pure-LTD mass shift {epsilon}"
                )
            })?;
            let parametric_atom = shifted_expression
                .diagnostic_parametric_atom_gs(&shifted_graph, &OrientationPattern::default());
            let orientations = diagnostic_evaluation_orientations(&shifted_expression);
            let prepared = match build_threedrep_evaluator(
                &format!("{}_{}", request.name, mass_shift_file_label(epsilon)),
                EvaluatorBuildContext {
                    workspace: request.workspace,
                    global_cli_settings: request.global_cli_settings,
                    default_runtime_settings: request.default_runtime_settings,
                    run_settings: request.run_settings,
                },
                &shifted_graph,
                &[&parametric_atom],
                &orientations,
            ) {
                Ok(prepared) => prepared,
                Err(error) => return Ok((format!("failed: {error:?}"), evaluations)),
            };
            let evaluator_info = prepared.info.clone();
            let mut evaluator = prepared.evaluator;
            evaluations.push(evaluate_with_evaluator(
                &mut evaluator,
                EvaluationRequest {
                    label: request.name,
                    graph: &shifted_graph,
                    expression: &shifted_expression,
                    parametric_atom: &parametric_atom,
                    workspace: request.workspace,
                    global_cli_settings: request.global_cli_settings,
                    default_runtime_settings: request.default_runtime_settings,
                    representation: Some(RepresentationMode::PureLtd),
                    numerator_sampling_scale_mode: Some(request.scale_mode),
                    run_settings: request.run_settings,
                    mass_shift: &mass_shift,
                    mass_shift_values: &mass_shift_values,
                },
                &orientations,
                &evaluator_info,
            )?);
        }
        Ok(("ok (mass-shift diagnostics)".to_string(), evaluations))
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

fn threedrep_run_settings(
    global_cli_settings: &CLISettings,
    default_runtime_settings: &RuntimeSettings,
    precision: Option<CliRuntimePrecision>,
    seed: u64,
    scale: f64,
    numerator_samples_normalization: CliNumeratorSamplesNormalization,
) -> Result<ThreeDrepRunSettings> {
    if !scale.is_finite() {
        return Err(eyre!("3Drep evaluation --scale must be finite"));
    }
    let runtime_precision = CliRuntimePrecision::resolve(default_runtime_settings, precision);
    let evaluator_method = default_runtime_settings.general.evaluator_method.clone();
    let evaluator_settings =
        threedrep_evaluator_settings(global_cli_settings, &evaluator_method, runtime_precision);
    let frozen_mode = if runtime_precision == CliRuntimePrecision::Double {
        global_cli_settings
            .global
            .generation
            .compile
            .frozen_mode(&evaluator_settings)
    } else {
        FrozenCompilationMode::Eager
    };

    Ok(ThreeDrepRunSettings {
        runtime_precision,
        evaluator_method,
        evaluator_backend: frozen_mode.active_backend_name().to_string(),
        compiled_backend_available: frozen_mode.compile_enabled(),
        numerator_samples_normalization,
        numerator_interpolation_scale: definite_numerator_interpolation_scale(
            default_runtime_settings,
        ),
        seed,
        scale,
    })
}

fn threedrep_evaluator_settings(
    global_cli_settings: &CLISettings,
    evaluator_method: &EvaluatorMethod,
    runtime_precision: CliRuntimePrecision,
) -> EvaluatorSettings {
    let mut settings = global_cli_settings.global.generation.evaluator;
    settings.iterative_orientation_optimization =
        matches!(evaluator_method, EvaluatorMethod::Iterative);
    settings.summed_function_map = matches!(evaluator_method, EvaluatorMethod::SummedFunctionMap);
    settings.summed = matches!(evaluator_method, EvaluatorMethod::Summed);
    if runtime_precision != CliRuntimePrecision::Double {
        settings.compile = false;
    }
    settings
}

fn threedrep_runtime_settings(
    default_runtime_settings: &RuntimeSettings,
    run_settings: &ThreeDrepRunSettings,
) -> RuntimeSettings {
    let mut settings = default_runtime_settings.clone();
    settings.general.evaluator_method = run_settings.evaluator_method.clone();
    settings
}

fn build_threedrep_evaluator<A: AtomCore>(
    name: &str,
    context: EvaluatorBuildContext<'_>,
    graph: &Graph,
    atoms: &[A],
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
) -> Result<PreparedEvaluator> {
    let evaluator_settings = threedrep_evaluator_settings(
        context.global_cli_settings,
        &context.run_settings.evaluator_method,
        context.run_settings.runtime_precision,
    );
    let frozen_mode = if context.run_settings.runtime_precision == CliRuntimePrecision::Double {
        context
            .global_cli_settings
            .global
            .generation
            .compile
            .frozen_mode(&evaluator_settings)
    } else {
        FrozenCompilationMode::Eager
    };
    let evaluator_dir = context.workspace.join("evaluators");
    fs::create_dir_all(&evaluator_dir)
        .with_context(|| format!("Could not create {}", evaluator_dir.display()))?;

    let start = Instant::now();
    let mut evaluator = EvaluatorStack::new(
        atoms,
        &graph.param_builder,
        &orientations.raw,
        None,
        &evaluator_settings,
    )?;
    evaluator.prepare_f64_backend(name, &evaluator_dir, &frozen_mode)?;
    let build_timing = start.elapsed();

    let runtime_settings =
        threedrep_runtime_settings(context.default_runtime_settings, context.run_settings);
    if runtime_settings.general.evaluator_method != context.run_settings.evaluator_method {
        return Err(eyre!("3Drep evaluator method resolution mismatch"));
    }

    Ok(PreparedEvaluator {
        evaluator,
        info: PreparedEvaluatorInfo {
            build_timing: format_duration_dynamic(build_timing),
            build_timing_seconds: build_timing.as_secs_f64(),
            backend: frozen_mode.active_backend_name().to_string(),
        },
    })
}

fn format_duration_dynamic(duration: Duration) -> String {
    let seconds = duration.as_secs_f64();
    let (value, unit) = if seconds >= 1.0 {
        (seconds, "s")
    } else if seconds >= 1.0e-3 {
        (seconds * 1.0e3, "ms")
    } else {
        (seconds * 1.0e6, "µs")
    };
    let precision = if value >= 100.0 {
        0
    } else if value >= 10.0 {
        1
    } else {
        2
    };
    format!("{value:.precision$} {unit}")
}

fn format_f64_full(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else if value.is_infinite() {
        if value.is_sign_negative() {
            "-inf".to_string()
        } else {
            "+inf".to_string()
        }
    } else {
        format!("{value:+.17e}")
    }
}

fn format_complex_full(value: &Complex<F<f64>>) -> (String, String, String) {
    let re = format_f64_full(value.re.0);
    let im = format_f64_full(value.im.0);
    let separator = if im.starts_with('-') || im.starts_with('+') {
        ""
    } else {
        "+"
    };
    let value = format!("({re}{separator}{im}i)");
    (value, re, im)
}

fn display_decimal_precision<T: FloatLike>(value: &F<T>) -> usize {
    let precision_bits = value.0.get_precision().max(1) as f64;
    (precision_bits * std::f64::consts::LOG10_2).ceil().max(1.0) as usize + 1
}

fn format_complex_full_precise<T: FloatLike>(value: &Complex<F<T>>) -> (String, String, String) {
    let precision = display_decimal_precision(&value.re).max(display_decimal_precision(&value.im));
    let re = format!("{:+.*e}", precision, value.re);
    let im = format!("{:+.*e}", precision, value.im);
    let separator = if im.starts_with('-') || im.starts_with('+') {
        ""
    } else {
        "+"
    };
    let value = format!("({re}{separator}{im}i)");
    (value, re, im)
}

const DIRECT_RELATIVE_TOLERANCE: f64 = 1.0e-7;
const DIRECT_ABSOLUTE_TOLERANCE: f64 = 1.0e-12;
const MASS_SHIFT_RELATIVE_TOLERANCE: f64 = 1.0e-2;
const MASS_SHIFT_ABSOLUTE_TOLERANCE: f64 = 1.0e-8;

#[derive(Clone, Debug)]
struct DecimalValue {
    sign: i8,
    digits: Vec<u8>,
    exponent: i32,
}

impl DecimalValue {
    fn parse(input: &str) -> Option<Self> {
        let mut value = input.trim();
        if value.eq_ignore_ascii_case("nan")
            || value.eq_ignore_ascii_case("inf")
            || value.eq_ignore_ascii_case("+inf")
            || value.eq_ignore_ascii_case("-inf")
        {
            return None;
        }

        let mut sign = 1;
        if let Some(stripped) = value.strip_prefix('-') {
            sign = -1;
            value = stripped;
        } else if let Some(stripped) = value.strip_prefix('+') {
            value = stripped;
        }

        let (mantissa, parsed_exponent) =
            if let Some((mantissa, exponent)) = value.split_once(['e', 'E']) {
                (mantissa, exponent.parse::<i32>().ok()?)
            } else {
                (value, 0)
            };
        let mut exponent = parsed_exponent;
        let mut digits = Vec::new();
        let mut fractional_digits = 0usize;
        let mut past_decimal_point = false;
        for ch in mantissa.chars() {
            match ch {
                '0'..='9' => {
                    digits.push(ch as u8 - b'0');
                    if past_decimal_point {
                        fractional_digits += 1;
                    }
                }
                '.' if !past_decimal_point => past_decimal_point = true,
                _ => return None,
            }
        }
        exponent -= fractional_digits as i32;
        Some(Self::new(sign, digits, exponent))
    }

    fn new(sign: i8, mut digits: Vec<u8>, mut exponent: i32) -> Self {
        let first_nonzero = digits.iter().position(|digit| *digit != 0);
        let Some(first_nonzero) = first_nonzero else {
            return Self {
                sign: 0,
                digits: Vec::new(),
                exponent: 0,
            };
        };
        if first_nonzero > 0 {
            digits.drain(0..first_nonzero);
        }
        while digits.last() == Some(&0) {
            digits.pop();
            exponent += 1;
        }
        Self {
            sign: if sign < 0 { -1 } else { 1 },
            digits,
            exponent,
        }
    }

    fn is_zero(&self) -> bool {
        self.sign == 0
    }

    fn negated(&self) -> Self {
        let mut value = self.clone();
        value.sign = -value.sign;
        value
    }

    fn abs_scientific(&self) -> ScientificValue {
        if self.is_zero() {
            return ScientificValue::zero();
        }
        let take = self.digits.len().min(18);
        let leading = self
            .digits
            .iter()
            .take(take)
            .fold(0_u64, |acc, digit| acc * 10 + u64::from(*digit));
        let mantissa = leading as f64 / 10_f64.powi(take.saturating_sub(1) as i32);
        ScientificValue::new(mantissa, self.exponent + self.digits.len() as i32 - 1)
    }

    fn difference(lhs: &Self, rhs: &Self) -> Self {
        lhs.signed_add(&rhs.negated())
    }

    fn signed_add(&self, rhs: &Self) -> Self {
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }

        let exponent = self.exponent.min(rhs.exponent);
        let lhs_digits = self.aligned_digits(exponent);
        let rhs_digits = rhs.aligned_digits(exponent);
        if self.sign == rhs.sign {
            return Self::new(
                self.sign,
                add_decimal_digits(&lhs_digits, &rhs_digits),
                exponent,
            );
        }

        match compare_decimal_digits(&lhs_digits, &rhs_digits) {
            std::cmp::Ordering::Greater => Self::new(
                self.sign,
                subtract_decimal_digits(&lhs_digits, &rhs_digits),
                exponent,
            ),
            std::cmp::Ordering::Less => Self::new(
                rhs.sign,
                subtract_decimal_digits(&rhs_digits, &lhs_digits),
                exponent,
            ),
            std::cmp::Ordering::Equal => Self::new(0, Vec::new(), 0),
        }
    }

    fn aligned_digits(&self, target_exponent: i32) -> Vec<u8> {
        let mut digits = self.digits.clone();
        digits.extend(std::iter::repeat_n(
            0,
            self.exponent.saturating_sub(target_exponent) as usize,
        ));
        digits
    }
}

fn compare_decimal_digits(lhs: &[u8], rhs: &[u8]) -> std::cmp::Ordering {
    let lhs_first = lhs
        .iter()
        .position(|digit| *digit != 0)
        .unwrap_or(lhs.len());
    let rhs_first = rhs
        .iter()
        .position(|digit| *digit != 0)
        .unwrap_or(rhs.len());
    let lhs = &lhs[lhs_first..];
    let rhs = &rhs[rhs_first..];
    lhs.len().cmp(&rhs.len()).then_with(|| lhs.cmp(rhs))
}

fn add_decimal_digits(lhs: &[u8], rhs: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut carry = 0_u8;
    let mut lhs_iter = lhs.iter().rev();
    let mut rhs_iter = rhs.iter().rev();
    loop {
        let lhs_digit = lhs_iter.next().copied();
        let rhs_digit = rhs_iter.next().copied();
        if lhs_digit.is_none() && rhs_digit.is_none() && carry == 0 {
            break;
        }
        let sum = lhs_digit.unwrap_or(0) + rhs_digit.unwrap_or(0) + carry;
        out.push(sum % 10);
        carry = sum / 10;
    }
    out.reverse();
    out
}

fn subtract_decimal_digits(lhs: &[u8], rhs: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut borrow = 0_i8;
    let mut lhs_iter = lhs.iter().rev();
    let mut rhs_iter = rhs.iter().rev();
    loop {
        let Some(lhs_digit) = lhs_iter.next().copied() else {
            break;
        };
        let mut digit = lhs_digit as i8 - borrow - rhs_iter.next().copied().unwrap_or(0) as i8;
        if digit < 0 {
            digit += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }
        out.push(digit as u8);
    }
    out.reverse();
    out
}

#[derive(Clone, Debug)]
struct EvaluationValueDecimal {
    re: DecimalValue,
    im: DecimalValue,
}

#[derive(Clone, Copy, Debug)]
struct ScientificValue {
    mantissa: f64,
    exponent: i32,
}

impl ScientificValue {
    fn zero() -> Self {
        Self {
            mantissa: 0.0,
            exponent: 0,
        }
    }

    fn new(mantissa: f64, exponent: i32) -> Self {
        if mantissa == 0.0 || !mantissa.is_finite() {
            return Self::zero();
        }
        let mut mantissa = mantissa.abs();
        let mut exponent = exponent;
        while mantissa >= 10.0 {
            mantissa /= 10.0;
            exponent += 1;
        }
        while mantissa < 1.0 {
            mantissa *= 10.0;
            exponent -= 1;
        }
        Self { mantissa, exponent }
    }

    fn from_f64(value: f64) -> Self {
        if value == 0.0 || !value.is_finite() {
            return Self::zero();
        }
        let exponent = value.abs().log10().floor() as i32;
        Self::new(value.abs() / 10_f64.powi(exponent), exponent)
    }

    fn is_zero(self) -> bool {
        self.mantissa == 0.0
    }

    fn hypot(lhs: Self, rhs: Self) -> Self {
        if lhs.is_zero() {
            return rhs;
        }
        if rhs.is_zero() {
            return lhs;
        }
        let exponent = lhs.exponent.max(rhs.exponent);
        let lhs_scaled = lhs.mantissa * 10_f64.powi(lhs.exponent - exponent);
        let rhs_scaled = rhs.mantissa * 10_f64.powi(rhs.exponent - exponent);
        Self::new(lhs_scaled.hypot(rhs_scaled), exponent)
    }

    fn ratio(numerator: Self, denominator: Self) -> Self {
        if numerator.is_zero() {
            return Self::zero();
        }
        if denominator.is_zero() {
            return Self::new(f64::INFINITY, i32::MAX);
        }
        Self::new(
            numerator.mantissa / denominator.mantissa,
            numerator.exponent - denominator.exponent,
        )
    }

    fn max(self, rhs: Self) -> Self {
        if compare_scientific(self, rhs) == std::cmp::Ordering::Less {
            rhs
        } else {
            self
        }
    }

    fn leq_f64(self, rhs: f64) -> bool {
        compare_scientific(self, Self::from_f64(rhs)) != std::cmp::Ordering::Greater
    }
}

fn compare_scientific(lhs: ScientificValue, rhs: ScientificValue) -> std::cmp::Ordering {
    if lhs.is_zero() && rhs.is_zero() {
        return std::cmp::Ordering::Equal;
    }
    if lhs.is_zero() {
        return std::cmp::Ordering::Less;
    }
    if rhs.is_zero() {
        return std::cmp::Ordering::Greater;
    }
    lhs.exponent
        .cmp(&rhs.exponent)
        .then_with(|| lhs.mantissa.total_cmp(&rhs.mantissa))
}

fn format_scientific_significant(
    value: Option<ScientificValue>,
    significant_digits: usize,
) -> String {
    let Some(value) = value else {
        return "-".to_string();
    };
    if value.is_zero() {
        let decimals = significant_digits.saturating_sub(1);
        return format!("{:+.decimals$e}", 0.0);
    }
    let decimals = significant_digits.saturating_sub(1);
    format!("{:+.decimals$}e{}", value.mantissa, value.exponent)
}

#[derive(Clone)]
struct ComparableEvaluation {
    id: usize,
    case_name: String,
    representation: Option<RepresentationMode>,
    scale_mode: Option<three_dimensional_reps::NumeratorSamplingScaleMode>,
    numerator_interpolation_scale: f64,
    mass_shift: String,
    value: EvaluationValueDecimal,
    value_text: String,
}

struct TestCffLtdDistance {
    abs_diff: Option<ScientificValue>,
    rel_diff: Option<ScientificValue>,
    tolerance: Option<ScientificValue>,
}

impl TestCffLtdDistance {
    fn new(abs_diff: ScientificValue, rel_diff: ScientificValue, tolerance: f64) -> Self {
        Self {
            abs_diff: Some(abs_diff),
            rel_diff: Some(rel_diff),
            tolerance: Some(ScientificValue::from_f64(tolerance)),
        }
    }

    fn undefined() -> Self {
        Self {
            abs_diff: None,
            rel_diff: None,
            tolerance: None,
        }
    }
}

impl ComparableEvaluation {
    fn label(&self) -> String {
        let representation = self
            .representation
            .map(|representation| format!("{representation:?}"))
            .unwrap_or_else(|| "-".to_string());
        let scale_mode = self
            .scale_mode
            .map(|scale_mode| format!("{scale_mode:?}"))
            .unwrap_or_else(|| "-".to_string());
        let scale = format_f64_full(self.numerator_interpolation_scale);
        format!(
            "#{} {} {representation} {scale_mode} M={scale} mass_shift={}",
            self.id, self.case_name, self.mass_shift
        )
    }
}

fn build_test_cff_ltd_verdict(output: &TestCffLtdOutput) -> TestCffLtdVerdict {
    let mut checks = Vec::new();

    for case in &output.cases {
        if case.generation_status != "ok" {
            checks.push(test_cff_ltd_check(
                format!("{} generation", case.name),
                None,
                TestCffLtdStatus::Fail,
                "-",
                "-",
                TestCffLtdDistance::undefined(),
                &case.generation_status,
            ));
            continue;
        }

        if !case.evaluator_build_status.starts_with("ok") {
            checks.push(test_cff_ltd_check(
                format!("{} evaluator", case.name),
                None,
                TestCffLtdStatus::Fail,
                "-",
                "-",
                TestCffLtdDistance::undefined(),
                &case.evaluator_build_status,
            ));
        }

        for evaluation in &case.evaluations {
            if evaluation.status != "ok" {
                checks.push(test_cff_ltd_check(
                    format!("{} evaluation {}", case.name, evaluation.mass_shift),
                    Some(evaluation.id),
                    TestCffLtdStatus::Fail,
                    &evaluation.value,
                    "-",
                    TestCffLtdDistance::undefined(),
                    evaluation.error.as_deref().unwrap_or(&evaluation.status),
                ));
            }
        }
    }

    let direct = comparable_evaluations(output, false);
    let mass_shifted = comparable_evaluations(output, true);

    let baseline = direct
        .iter()
        .find(|evaluation| evaluation.representation == Some(RepresentationMode::Cff))
        .cloned();

    let Some(baseline) = baseline else {
        checks.push(test_cff_ltd_check(
            "cff reference".to_string(),
            None,
            TestCffLtdStatus::Fail,
            "-",
            "-",
            TestCffLtdDistance::undefined(),
            "No finite CFF reference evaluation was available.",
        ));
        return finish_test_cff_ltd_verdict(checks);
    };

    for candidate in direct {
        if candidate.label() == baseline.label() {
            continue;
        }
        let (abs_diff, rel_diff) = complex_distance(&candidate.value, &baseline.value);
        let within_tolerance = abs_diff.leq_f64(DIRECT_ABSOLUTE_TOLERANCE)
            || rel_diff.leq_f64(DIRECT_RELATIVE_TOLERANCE);
        let exactly_identical = candidate.representation != baseline.representation
            && candidate.value_text == baseline.value_text;
        let status = if within_tolerance && !exactly_identical {
            TestCffLtdStatus::Success
        } else {
            TestCffLtdStatus::Fail
        };
        let message = if exactly_identical {
            "Different representations produced an exactly identical formatted value."
        } else if within_tolerance {
            "Direct representations agree within the half-precision threshold."
        } else {
            "Direct representations differ beyond the half-precision threshold."
        };
        checks.push(test_cff_ltd_check(
            format!("direct {} vs #{}", candidate.label(), baseline.id),
            Some(candidate.id),
            status,
            &candidate.value_text,
            &baseline.value_text,
            TestCffLtdDistance::new(abs_diff, rel_diff, DIRECT_RELATIVE_TOLERANCE),
            message,
        ));
    }

    for group in mass_shift_groups(&mass_shifted) {
        if !group.is_empty() {
            let best = group
                .iter()
                .map(|evaluation| {
                    let (abs_diff, rel_diff) = complex_distance(&evaluation.value, &baseline.value);
                    (evaluation, abs_diff, rel_diff)
                })
                .min_by(|(_, lhs_abs, lhs_rel), (_, rhs_abs, rhs_rel)| {
                    compare_scientific(*lhs_rel, *rhs_rel)
                        .then_with(|| compare_scientific(*lhs_abs, *rhs_abs))
                });
            let Some((best, abs_diff, rel_diff)) = best else {
                continue;
            };
            let within_tolerance = abs_diff.leq_f64(MASS_SHIFT_ABSOLUTE_TOLERANCE)
                || rel_diff.leq_f64(MASS_SHIFT_RELATIVE_TOLERANCE);
            checks.push(test_cff_ltd_check(
                format!("mass-shift {} vs #{}", best.label(), baseline.id),
                Some(best.id),
                if within_tolerance {
                    TestCffLtdStatus::Success
                } else {
                    TestCffLtdStatus::Fail
                },
                &best.value_text,
                &baseline.value_text,
                TestCffLtdDistance::new(abs_diff, rel_diff, MASS_SHIFT_RELATIVE_TOLERANCE),
                if within_tolerance {
                    "At least one tested mass shift reproduces the CFF reference to the loose diagnostic threshold."
                } else {
                    "No tested mass shift reproduces the CFF reference to the loose diagnostic threshold."
                },
            ));
        }
    }

    finish_test_cff_ltd_verdict(checks)
}

fn comparable_evaluations(
    output: &TestCffLtdOutput,
    require_mass_shift: bool,
) -> Vec<ComparableEvaluation> {
    output
        .cases
        .iter()
        .flat_map(|case| {
            case.evaluations.iter().filter_map(move |evaluation| {
                let has_mass_shift = !evaluation.mass_shift_values.is_empty();
                if has_mass_shift != require_mass_shift || evaluation.status != "ok" {
                    return None;
                }
                let value = parse_evaluation_value(evaluation)?;
                Some(ComparableEvaluation {
                    id: evaluation.id,
                    case_name: case.name.clone(),
                    representation: evaluation.representation,
                    scale_mode: evaluation.numerator_sampling_scale_mode,
                    numerator_interpolation_scale: evaluation.numerator_interpolation_scale,
                    mass_shift: evaluation.mass_shift.clone(),
                    value,
                    value_text: evaluation.value.clone(),
                })
            })
        })
        .collect()
}

fn mass_shift_groups(evaluations: &[ComparableEvaluation]) -> Vec<Vec<ComparableEvaluation>> {
    let mut groups = BTreeMap::<String, Vec<ComparableEvaluation>>::new();
    for evaluation in evaluations {
        let scale = format_f64_full(evaluation.numerator_interpolation_scale);
        let key = format!(
            "{}::{:?}::{:?}::{scale}",
            evaluation.case_name, evaluation.representation, evaluation.scale_mode
        );
        groups.entry(key).or_default().push(evaluation.clone());
    }
    groups.into_values().collect()
}

fn parse_evaluation_value(
    evaluation: &ThreeDrepEvaluationRecord,
) -> Option<EvaluationValueDecimal> {
    Some(EvaluationValueDecimal {
        re: DecimalValue::parse(&evaluation.value_re)?,
        im: DecimalValue::parse(&evaluation.value_im)?,
    })
}

fn complex_distance(
    lhs: &EvaluationValueDecimal,
    rhs: &EvaluationValueDecimal,
) -> (ScientificValue, ScientificValue) {
    let re_diff = DecimalValue::difference(&lhs.re, &rhs.re).abs_scientific();
    let im_diff = DecimalValue::difference(&lhs.im, &rhs.im).abs_scientific();
    let abs_diff = ScientificValue::hypot(re_diff, im_diff);
    let lhs_abs = ScientificValue::hypot(lhs.re.abs_scientific(), lhs.im.abs_scientific());
    let rhs_abs = ScientificValue::hypot(rhs.re.abs_scientific(), rhs.im.abs_scientific());
    let scale = lhs_abs.max(rhs_abs);
    let rel_diff = ScientificValue::ratio(abs_diff, scale);
    (abs_diff, rel_diff)
}

fn test_cff_ltd_check(
    name: String,
    evaluation_id: Option<usize>,
    status: TestCffLtdStatus,
    lhs: impl AsRef<str>,
    rhs: impl AsRef<str>,
    distance: TestCffLtdDistance,
    message: impl AsRef<str>,
) -> TestCffLtdCheckRecord {
    TestCffLtdCheckRecord {
        name,
        evaluation_id: evaluation_id
            .map(|id| format!("#{id}"))
            .unwrap_or_else(|| "-".to_string()),
        status,
        lhs: lhs.as_ref().to_string(),
        rhs: rhs.as_ref().to_string(),
        abs_diff: format_scientific_significant(distance.abs_diff, 3),
        rel_diff: format_scientific_significant(distance.rel_diff, 3),
        tolerance: format_scientific_significant(distance.tolerance, 3),
        message: message.as_ref().to_string(),
    }
}

fn finish_test_cff_ltd_verdict(checks: Vec<TestCffLtdCheckRecord>) -> TestCffLtdVerdict {
    let status = if checks
        .iter()
        .any(|check| check.status != TestCffLtdStatus::Success)
    {
        TestCffLtdStatus::Fail
    } else {
        TestCffLtdStatus::Success
    };
    TestCffLtdVerdict {
        status,
        direct_relative_tolerance: format_f64_full(DIRECT_RELATIVE_TOLERANCE),
        direct_absolute_tolerance: format_f64_full(DIRECT_ABSOLUTE_TOLERANCE),
        mass_shift_relative_tolerance: format_f64_full(MASS_SHIFT_RELATIVE_TOLERANCE),
        mass_shift_absolute_tolerance: format_f64_full(MASS_SHIFT_ABSOLUTE_TOLERANCE),
        checks,
    }
}

fn pure_ltd_mass_shift_epsilons(start: f64) -> [f64; 4] {
    [start, 0.1 * start, 0.01 * start, 0.001 * start]
}

fn mass_shift_file_label(epsilon: f64) -> String {
    format!("eps_{epsilon:.3e}")
        .replace(['+', '-'], "")
        .replace('.', "p")
}

fn definite_numerator_interpolation_scale(default_runtime_settings: &RuntimeSettings) -> f64 {
    default_runtime_settings
        .general
        .numerator_interpolation_scale
        .unwrap_or(default_runtime_settings.general.numerator_sampling_scale)
}

fn representation_label(representation: RepresentationMode) -> &'static str {
    match representation {
        RepresentationMode::Cff => "CFF",
        RepresentationMode::Ltd => "LTD",
        RepresentationMode::PureLtd => "PureLTD",
    }
}

fn mass_shift_value_records(mass_shifts: &[ThreeDRepMassShift]) -> Vec<MassShiftValueRecord> {
    mass_shifts
        .iter()
        .map(|mass_shift| MassShiftValueRecord {
            group_index: mass_shift.group_index,
            split_index: mass_shift.split_index,
            local_edge_id: mass_shift.local_edge_id,
            edge_id: mass_shift.edge_id.0,
            base_mass: format!("{:.3e}", mass_shift.base_mass),
            shifted_mass: format!("{:.3e}", mass_shift.shifted_mass),
        })
        .collect()
}

fn mass_shift_label(epsilon: f64, values: &[MassShiftValueRecord]) -> String {
    let assignments = values
        .iter()
        .map(|value| format!("e{}={}", value.edge_id, value.shifted_mass))
        .collect::<Vec<_>>()
        .join(" ");
    if assignments.is_empty() {
        format!("eps={epsilon:.3e}")
    } else {
        format!("eps={epsilon:.3e} {assignments}")
    }
}

fn diagnostic_evaluation_orientations(
    expression: &ThreeDExpression<OrientationID>,
) -> TiVec<OrientationID, EdgeVec<Orientation>> {
    expression
        .orientations
        .iter()
        .take(1)
        .map(|orientation| orientation.data.orientation.clone())
        .collect()
}

fn set_random_diagnostic_input_values(
    param_builder: &mut ParamBuilder<f64>,
    seed: u64,
    scale: f64,
) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let parameters = param_builder.evaluator_input_parameters();
    for parameter in parameters {
        if matches!(
            parameter.group,
            ParamBuilderInputGroup::ExternalEnergy
                | ParamBuilderInputGroup::ExternalSpatial
                | ParamBuilderInputGroup::LoopMomentumSpatial
                | ParamBuilderInputGroup::Additional
        ) {
            param_builder.values[0][parameter.index] =
                Complex::new_re(F(scale * rng.random::<f64>()));
        }
    }
}

fn parameter_records(
    graph: &Graph,
    run_settings: &ThreeDrepRunSettings,
    default_runtime_settings: &RuntimeSettings,
) -> Vec<ParameterValueRecord> {
    let mut param_builder = graph.param_builder.clone();
    set_random_diagnostic_input_values(&mut param_builder, run_settings.seed, run_settings.scale);
    param_builder.set_runtime_parameter_values(
        Complex::new_re(F(default_runtime_settings.general.m_uv)),
        Complex::new_re(F(
            default_runtime_settings.general.mu_r * default_runtime_settings.general.mu_r
        )),
        Complex::new_re(F(run_settings.numerator_interpolation_scale)),
    );

    param_builder
        .evaluator_input_parameters()
        .into_iter()
        .map(|parameter| {
            let value = &param_builder.values[0][parameter.index];
            let (formatted_value, _, _) = format_complex_full(value);
            ParameterValueRecord {
                name: parameter.atom.log_print(Some(60)),
                canonical_name: parameter.atom.to_canonical_string(),
                value: formatted_value,
                source: parameter_source(parameter.group).to_string(),
            }
        })
        .collect()
}

fn parameter_source(group: ParamBuilderInputGroup) -> &'static str {
    match group {
        ParamBuilderInputGroup::ExternalEnergy
        | ParamBuilderInputGroup::ExternalSpatial
        | ParamBuilderInputGroup::LoopMomentumSpatial
        | ParamBuilderInputGroup::Additional => "automatic diagnostic",
        ParamBuilderInputGroup::Runtime
        | ParamBuilderInputGroup::Model
        | ParamBuilderInputGroup::Polarization
        | ParamBuilderInputGroup::LocalCounterterm => "state/default",
    }
}

struct EvaluationRequest<'a> {
    label: &'a str,
    graph: &'a Graph,
    expression: &'a ThreeDExpression<OrientationID>,
    parametric_atom: &'a symbolica::atom::Atom,
    workspace: &'a Path,
    global_cli_settings: &'a CLISettings,
    default_runtime_settings: &'a RuntimeSettings,
    representation: Option<RepresentationMode>,
    numerator_sampling_scale_mode: Option<three_dimensional_reps::NumeratorSamplingScaleMode>,
    run_settings: &'a ThreeDrepRunSettings,
    mass_shift: &'a str,
    mass_shift_values: &'a [MassShiftValueRecord],
}

fn evaluate_threedrep_expression(
    request: EvaluationRequest<'_>,
) -> Result<ThreeDrepEvaluationRecord> {
    let orientations = diagnostic_evaluation_orientations(request.expression);
    let prepared = build_threedrep_evaluator(
        request.label,
        EvaluatorBuildContext {
            workspace: request.workspace,
            global_cli_settings: request.global_cli_settings,
            default_runtime_settings: request.default_runtime_settings,
            run_settings: request.run_settings,
        },
        request.graph,
        &[request.parametric_atom],
        &orientations,
    )?;
    let evaluator_info = prepared.info.clone();
    let mut evaluator = prepared.evaluator;
    evaluate_with_evaluator(&mut evaluator, request, &orientations, &evaluator_info)
}

fn evaluate_with_evaluator(
    evaluator: &mut EvaluatorStack,
    request: EvaluationRequest<'_>,
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    evaluator_info: &PreparedEvaluatorInfo,
) -> Result<ThreeDrepEvaluationRecord> {
    let mut param_builder = request.graph.param_builder.clone();
    set_random_diagnostic_input_values(
        &mut param_builder,
        request.run_settings.seed,
        request.run_settings.scale,
    );
    param_builder.set_runtime_parameter_values(
        Complex::new_re(F(request.default_runtime_settings.general.m_uv)),
        Complex::new_re(F(request.default_runtime_settings.general.mu_r
            * request.default_runtime_settings.general.mu_r)),
        Complex::new_re(F(request.run_settings.numerator_interpolation_scale)),
    );

    let filter = SubSet::<OrientationID>::full(orientations.len());
    let mut evaluation_metadata = EvaluationMetaData::new_empty();
    let runtime_settings =
        threedrep_runtime_settings(request.default_runtime_settings, request.run_settings);
    let start = Instant::now();
    let result = evaluate_for_precision(
        evaluator,
        &mut param_builder,
        request.run_settings.runtime_precision,
        SingleOrAllOrientations::All {
            all: orientations,
            filter: &filter,
        },
        &runtime_settings,
        &mut evaluation_metadata,
    );
    let timing = start.elapsed();

    match result {
        Ok((formatted_value, value_re, value_im)) => Ok(ThreeDrepEvaluationRecord {
            id: 0,
            label: request.label.to_string(),
            representation: request.representation,
            numerator_sampling_scale_mode: request.numerator_sampling_scale_mode,
            mass_shift: request.mass_shift.to_string(),
            mass_shift_values: request.mass_shift_values.to_vec(),
            numerator_interpolation_scale: request.run_settings.numerator_interpolation_scale,
            runtime_precision: request.run_settings.runtime_precision,
            evaluator_method: request.run_settings.evaluator_method.clone(),
            evaluator_backend: evaluator_info.backend.clone(),
            evaluator_build_timing: Some(evaluator_info.build_timing.clone()),
            evaluator_build_timing_seconds: Some(evaluator_info.build_timing_seconds),
            value: formatted_value,
            value_re,
            value_im,
            sample_evaluation_timing: format_duration_dynamic(timing),
            sample_evaluation_timing_seconds: timing.as_secs_f64(),
            status: "ok".to_string(),
            error: None,
        }),
        Err(error) => Ok(ThreeDrepEvaluationRecord {
            id: 0,
            label: request.label.to_string(),
            representation: request.representation,
            numerator_sampling_scale_mode: request.numerator_sampling_scale_mode,
            mass_shift: request.mass_shift.to_string(),
            mass_shift_values: request.mass_shift_values.to_vec(),
            numerator_interpolation_scale: request.run_settings.numerator_interpolation_scale,
            runtime_precision: request.run_settings.runtime_precision,
            evaluator_method: request.run_settings.evaluator_method.clone(),
            evaluator_backend: evaluator_info.backend.clone(),
            evaluator_build_timing: Some(evaluator_info.build_timing.clone()),
            evaluator_build_timing_seconds: Some(evaluator_info.build_timing_seconds),
            value: "-".to_string(),
            value_re: "-".to_string(),
            value_im: "-".to_string(),
            sample_evaluation_timing: format_duration_dynamic(timing),
            sample_evaluation_timing_seconds: timing.as_secs_f64(),
            status: "failed".to_string(),
            error: Some(format!("{error:?}")),
        }),
    }
}

fn evaluate_for_precision(
    evaluator: &mut EvaluatorStack,
    param_builder: &mut ParamBuilder<f64>,
    precision: CliRuntimePrecision,
    orientations: SingleOrAllOrientations<'_, OrientationID>,
    runtime_settings: &RuntimeSettings,
    evaluation_metadata: &mut EvaluationMetaData,
) -> Result<(String, String, String)> {
    match precision {
        CliRuntimePrecision::Double => {
            let input = param_builder.input_params();
            evaluate_for_precision_impl::<f64>(
                evaluator,
                input,
                orientations,
                runtime_settings,
                evaluation_metadata,
            )
        }
        CliRuntimePrecision::Quad => {
            let input = param_builder.input_params_quad();
            evaluate_for_precision_impl::<f128>(
                evaluator,
                input,
                orientations,
                runtime_settings,
                evaluation_metadata,
            )
        }
        CliRuntimePrecision::ArbPrec => {
            let input = param_builder.input_params_arb();
            evaluate_for_precision_impl::<ArbPrec>(
                evaluator,
                input,
                orientations,
                runtime_settings,
                evaluation_metadata,
            )
        }
    }
}

fn evaluate_for_precision_impl<T: FloatLike>(
    evaluator: &mut EvaluatorStack,
    input: InputParams<'_, T>,
    orientations: SingleOrAllOrientations<'_, OrientationID>,
    runtime_settings: &RuntimeSettings,
    evaluation_metadata: &mut EvaluationMetaData,
) -> Result<(String, String, String)> {
    let mut values = evaluator.evaluate::<T, OrientationID>(
        input,
        orientations,
        runtime_settings,
        evaluation_metadata,
        true,
    )?;
    let value = values
        .pop()
        .ok_or_else(|| eyre!("3Drep evaluator returned no value"))?
        .unwrap_real();
    Ok(format_complex_full_precise(&value))
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
    let color = if text == "ok" || text == "SUCCESS" {
        Color::Green
    } else if text == "FAIL" || text.starts_with("failed") {
        Color::Red
    } else {
        Color::Yellow
    };
    color_text(text, color)
}

fn evaluation_delta_cell(
    output: &TestCffLtdOutput,
    evaluation: &ThreeDrepEvaluationRecord,
) -> String {
    let Some(reference) = output
        .cases
        .iter()
        .flat_map(|case| &case.evaluations)
        .find(|candidate| {
            candidate.representation == Some(RepresentationMode::Cff)
                && candidate.mass_shift_values.is_empty()
                && candidate.status == "ok"
        })
        .and_then(parse_evaluation_value)
    else {
        return color_text("-", Color::Yellow);
    };
    let Some(value) = parse_evaluation_value(evaluation) else {
        return color_text("-", Color::Yellow);
    };
    if evaluation.representation == Some(RepresentationMode::Cff)
        && evaluation.mass_shift_values.is_empty()
    {
        return color_text("ref", Color::Green);
    }
    let (_, delta) = complex_distance(&value, &reference);
    let threshold = if evaluation.mass_shift_values.is_empty() {
        DIRECT_RELATIVE_TOLERANCE
    } else {
        MASS_SHIFT_RELATIVE_TOLERANCE
    };
    color_text(
        format_scientific_significant(Some(delta), 2),
        if delta.leq_f64(threshold) {
            Color::Green
        } else {
            Color::Red
        },
    )
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
    table.push_record(vec![
        "precision".to_string(),
        color_text(
            format!("{:?}", output.settings.runtime_precision),
            Color::Green,
        ),
    ]);
    table.push_record(vec![
        "evaluator method".to_string(),
        color_text(
            format!("{:?}", output.settings.evaluator_method),
            Color::Green,
        ),
    ]);
    table.push_record(vec![
        "evaluator backend".to_string(),
        color_text(&output.settings.evaluator_backend, Color::Green),
    ]);
    table.push_record(vec![
        "numerator samples normalization".to_string(),
        color_text(
            output.settings.numerator_samples_normalization.label(),
            Color::Purple,
        ),
    ]);
    table.push_record(vec![
        "M".to_string(),
        color_text(
            format_f64_full(output.settings.numerator_interpolation_scale),
            Color::Yellow,
        ),
    ]);
    table.push_record(vec![
        "seed".to_string(),
        color_text(output.settings.seed.to_string(), Color::Blue),
    ]);
    table.push_record(vec![
        "input scale".to_string(),
        color_text(format_f64_full(output.settings.scale), Color::Yellow),
    ]);
    table.push_record(vec![
        "evaluator build time".to_string(),
        color_text(
            output
                .evaluation
                .evaluator_build_timing
                .as_deref()
                .unwrap_or("-"),
            Color::Yellow,
        ),
    ]);
    table.push_record(vec![
        "sample evaluation time".to_string(),
        color_text(&output.evaluation.sample_evaluation_timing, Color::Yellow),
    ]);
    table.push_record(vec![
        "value".to_string(),
        color_text(&output.evaluation.value, Color::Green),
    ]);

    let mut parameter_table = Builder::new();
    parameter_table.push_record(vec![
        table_header("parameter"),
        table_header("value"),
        table_header("source"),
    ]);
    for parameter in &output.parameters {
        parameter_table.push_record(vec![
            color_text(&parameter.name, Color::Blue),
            color_text(&parameter.value, Color::Green),
            color_text(&parameter.source, Color::Purple),
        ]);
    }

    format!(
        "{}\n\n{}",
        table.build().with(Style::rounded()),
        parameter_table.build().with(Style::rounded())
    )
}

fn render_test_cff_ltd_summary(output: &TestCffLtdOutput, manifest_path: &Path) -> String {
    let mut settings_table = Builder::new();
    settings_table.push_record(vec![table_header("setting"), table_header("value")]);
    settings_table.push_record(vec![
        "precision".to_string(),
        color_text(
            format!("{:?}", output.settings.runtime_precision),
            Color::Green,
        ),
    ]);
    settings_table.push_record(vec![
        "evaluator method".to_string(),
        color_text(
            format!("{:?}", output.settings.evaluator_method),
            Color::Green,
        ),
    ]);
    settings_table.push_record(vec![
        "evaluator backend".to_string(),
        color_text(&output.settings.evaluator_backend, Color::Green),
    ]);
    settings_table.push_record(vec![
        "numerator samples normalization".to_string(),
        color_text(
            output.settings.numerator_samples_normalization.label(),
            Color::Purple,
        ),
    ]);
    settings_table.push_record(vec![
        "M".to_string(),
        color_text(
            format_f64_full(output.numerator_interpolation_scale),
            Color::Yellow,
        ),
    ]);
    settings_table.push_record(vec![
        "seed".to_string(),
        color_text(output.settings.seed.to_string(), Color::Blue),
    ]);
    settings_table.push_record(vec![
        "input scale".to_string(),
        color_text(format_f64_full(output.settings.scale), Color::Yellow),
    ]);
    settings_table.push_record(vec![
        "mass shift start".to_string(),
        color_text(format!("{:.3e}", output.mass_shift_start), Color::Yellow),
    ]);

    let mut table = Builder::new();
    table.push_record(vec![
        table_header("rep"),
        table_header("generation"),
        table_header("orientations"),
        table_header("#nodes"),
        table_header("evaluator"),
        table_header("evaluator build time"),
    ]);
    for case in &output.cases {
        table.push_record(vec![
            color_text(representation_label(case.representation), Color::Green),
            status_text(&case.generation_status),
            color_text(case.orientation_count.to_string(), Color::Blue),
            color_text(case.unfolded_term_count.to_string(), Color::Blue),
            status_text(&case.evaluator_build_status),
            color_text(
                case.evaluator_build_timing.as_deref().unwrap_or("-"),
                Color::Yellow,
            ),
        ]);
    }
    let mut evaluation_table = Builder::new();
    evaluation_table.push_record(vec![
        table_header("id"),
        table_header("rep"),
        table_header("mass shift"),
        table_header("status"),
        table_header("Δ"),
        table_header("sample evaluation time"),
        table_header("value"),
    ]);
    for case in &output.cases {
        for evaluation in &case.evaluations {
            let delta = evaluation_delta_cell(output, evaluation);
            evaluation_table.push_record(vec![
                color_text(format!("#{}", evaluation.id), Color::Blue),
                color_text(representation_label(case.representation), Color::Green),
                color_text(&evaluation.mass_shift, Color::Purple),
                status_text(&evaluation.status),
                delta,
                color_text(&evaluation.sample_evaluation_timing, Color::Yellow),
                color_text(&evaluation.value, Color::Green),
            ]);
        }
    }

    let mut check_table = Builder::new();
    check_table.push_record(vec![
        table_header("evaluation id"),
        table_header("status"),
        table_header("abs diff"),
        table_header("rel diff"),
        table_header("tolerance"),
        table_header("message"),
    ]);
    for check in &output.verdict.checks {
        check_table.push_record(vec![
            color_text(&check.evaluation_id, Color::Blue),
            status_text(check.status.label()),
            color_text(&check.abs_diff, Color::Yellow),
            color_text(&check.rel_diff, Color::Yellow),
            color_text(&check.tolerance, Color::Purple),
            color_text(&check.message, Color::Green),
        ]);
    }

    format!(
        "{} for process {}, integrand {}, graph {} {}\n{}: {}\n{}: {}\n{}\n\n{}\n\n{}\n\n{}",
        color_text("3Drep comparison", Color::Cyan),
        color_text(format!("#{}", output.process_id), Color::Blue),
        color_text(format!("'{}'", output.integrand_name), Color::Green),
        color_text(format!("#{}", output.graph_id), Color::Blue),
        color_text(&output.graph_name, Color::Green),
        color_text("manifest", Color::Cyan),
        color_text(relative_display(manifest_path), Color::Purple),
        color_text("status", Color::Cyan),
        status_text(output.verdict.status.label()),
        settings_table.build().with(Style::rounded()),
        table.build().with(Style::rounded()),
        evaluation_table.build().with(Style::rounded()),
        check_table.build().with(Style::rounded())
    )
}
