use std::{
    collections::BTreeMap,
    fs,
    io::Cursor,
    ops::Deref,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use clap::{Args, Subcommand, ValueEnum};
use color_eyre::{
    eyre::{eyre, Context},
    Result,
};
use gammalooprs::{
    cff::expression::{
        internal_energy_parameter_atom_gs, numerator_with_internal_energy_parameters_gs,
        GammaLoopGraphOrientation, GammaLoopOrientationExpression,
    },
    cff::surface::GammaLoopSurfaceCache,
    graph::{FeynmanGraph, Graph, ThreeDRepMassShift},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            evaluators::{EvaluatorMethod, EvaluatorStack, InputParams, SingleOrAllOrientations},
            param_builder::{ParamBuilder, ParamBuilderInputGroup, ParamBuilderInputParameter},
            ProcessIntegrand,
        },
    },
    model::Model,
    numerator::symbolica_ext::AtomCoreExt,
    processes::{Amplitude, CrossSection, EvaluatorSettings, Process, ProcessCollection},
    settings::{
        global::{FrozenCompilationMode, OrientationPattern, UniformNumeratorSamplingScale},
        runtime::Precision,
        RuntimeSettings,
    },
    utils::{
        f128,
        symbolica_ext::{CallSymbol, LogPrint},
        ArbPrec, FloatLike, F, FUN_LIB, GS, TENSORLIB, W_,
    },
    DependentMomentaConstructor, GammaLoopContextContainer,
};
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
use indicatif::{ProgressBar, ProgressStyle};
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair, Orientation},
    subgraph::subset::SubSet,
};
use nu_ansi_term::{Color, Style as AnsiStyle};
use rand::{rngs::SmallRng, Rng, SeedableRng};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::{
    algebra::complex::Complex,
    network::{parsing::SPENSO_TAG, ExecutionResult, Sequential, SmallestDegree},
    structure::{
        concrete_index::ExpandedIndex,
        representation::{LibraryRep, Minkowski, RepName},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Indeterminate, Symbol},
    state::State as SymbolicaState,
};
use tabled::{builder::Builder, settings::Style};
use three_dimensional_reps::{
    generate_3d_expression, graph_info, reconstruct_dot_from_expression, render_expression_summary,
    validate_parsed_graph, DisplayOptions, GraphInfo, GraphValidation, NumeratorDisplay,
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

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Default)]
pub struct OptionalGraphSelectorArgs {
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
    pub graph: Option<String>,
}

impl OptionalGraphSelectorArgs {
    fn has_any_selector(&self) -> bool {
        self.process.is_some() || self.integrand_name.is_some() || self.graph.is_some()
    }

    fn require_graph_selection(&self, context: &str) -> Result<GraphSelectorArgs> {
        Ok(GraphSelectorArgs {
            process: self.process.clone(),
            integrand_name: self.integrand_name.clone(),
            graph: self.graph.clone().ok_or_else(|| {
                eyre!("{context} requires --graph/-g when selecting a cached representation")
            })?,
        })
    }
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
    #[command(flatten)]
    pub selection: OptionalGraphSelectorArgs,

    /// Explicit oriented-expression JSON to evaluate. Takes precedence over representation lookup.
    #[arg(long, value_hint = clap::ValueHint::FilePath)]
    pub json_in: Option<PathBuf>,

    /// Cached representation to evaluate when --json-in is not supplied.
    #[arg(long, value_enum)]
    pub representation: Option<CliRepresentationMode>,

    /// Cached numerator-sampling normalization variant to evaluate with --representation.
    #[arg(
        long = "numerator-samples-normalization",
        alias = "numerator-sampling-scale-mode",
        value_enum
    )]
    pub numerator_samples_normalization: Option<CliNumeratorSamplesNormalization>,

    /// 3Drep artifact workspace containing the oriented expression JSON.
    #[arg(long, value_hint = clap::ValueHint::DirPath)]
    pub workspace_path: Option<PathBuf>,

    #[arg(long, value_enum)]
    pub precision: Option<CliRuntimePrecision>,

    #[arg(long, default_value_t = 1)]
    pub seed: u64,

    #[arg(long, default_value_t = 1.0)]
    pub scale: f64,

    /// Profile repeated evaluator calls for approximately this duration.
    #[arg(long)]
    pub profile: Option<String>,

    /// Force the eager evaluator backend for Double precision.
    #[arg(long, default_value_t = false)]
    pub eager: bool,

    /// Evaluate only the numerator, keeping internal edge energies as Q(edge, cind(0)) input parameters.
    #[arg(long, default_value_t = false)]
    pub numerator_only: bool,

    /// Override numerator-only internal energy inputs as edge:value entries.
    #[arg(
        long = "numerator-q0",
        alias = "numerator-energy-components",
        value_name = "EDGE:VALUE",
        action = clap::ArgAction::Append
    )]
    pub numerator_q0: Vec<String>,

    /// File name used for the evaluate summary manifest in the 3Drep workspace.
    #[arg(long, default_value = "evaluate_manifest.json", value_name = "NAME")]
    pub manifest_name: String,

    /// Do not print the input-parameter table in the evaluate command output.
    #[arg(long, default_value_t = false)]
    pub no_show_parameters: bool,

    /// Write the raw Symbolica evaluator replay artifacts without building or running the evaluator.
    #[arg(long, default_value_t = false)]
    pub standalone_rust_only: bool,

    /// Build separate numerator/orientation multi-output evaluators and combine them at runtime.
    #[arg(long, default_value_t = false)]
    pub iterative: bool,

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

    #[arg(long, default_value_t = 5)]
    pub n_epsilon_steps: usize,

    /// Build separate numerator/orientation multi-output evaluators and combine them at runtime.
    #[arg(long, default_value_t = false)]
    pub iterative: bool,

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

fn cli_representation_name(value: CliRepresentationMode) -> &'static str {
    match value {
        CliRepresentationMode::Ltd => "ltd",
        CliRepresentationMode::Cff => "cff",
        CliRepresentationMode::PureLtd => "pure-ltd",
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
    n_epsilon_steps: usize,
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
    #[serde(default)]
    symbolica_expression_pretty_path: PathBuf,
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
    #[serde(default)]
    numerator_only: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    expression_path: Option<PathBuf>,
    symbolica_expression_path: PathBuf,
    symbolica_expression_pretty_path: PathBuf,
    symbolica_expression_raw_path: PathBuf,
    symbolica_expression_raw_script_path: PathBuf,
    param_builder_path: PathBuf,
    settings: ThreeDrepRunSettings,
    evaluation: ThreeDrepEvaluationRecord,
    parameters: Vec<ParameterValueRecord>,
}

const SYMBOLICA_RAW_ARCHIVE_SCHEMA_VERSION: u32 = 2;

#[derive(Debug, Serialize, Deserialize)]
struct SymbolicaEvaluatorInputArchive {
    schema_version: u32,
    description: String,
    evaluator_method: EvaluatorMethod,
    evaluator_backend: String,
    evaluator_settings: SymbolicaEvaluatorSettingsRecord,
    optimization_settings: SymbolicaOptimizationSettingsRecord,
    parameters: Vec<String>,
    function_map_entries: Vec<SymbolicaFunctionMapEntryRecord>,
    representative_input: Vec<SymbolicaComplexValueRecord>,
    calls: Vec<SymbolicaEvaluatorCallRecord>,
}

#[derive(Debug, Serialize, Deserialize)]
struct SymbolicaEvaluatorCallRecord {
    label: String,
    expressions: Vec<String>,
    additional_function_map_entries: Vec<SymbolicaFunctionMapEntryRecord>,
    dual_shape: Option<Vec<Vec<usize>>>,
}

#[derive(Debug, Serialize, Deserialize)]
struct SymbolicaFunctionMapEntryRecord {
    lhs: String,
    rhs: String,
    tags: Vec<String>,
    args: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
struct SymbolicaComplexValueRecord {
    re: String,
    im: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct SymbolicaEvaluatorSettingsRecord {
    do_algebra: bool,
    iterative_orientation_optimization: bool,
    summed: bool,
    summed_function_map: bool,
    compile: bool,
    store_atom: bool,
    do_fn_map_replacements: bool,
    horner_iterations: usize,
    n_cores: usize,
    cpe_iterations: Option<usize>,
    abort_level: usize,
    max_horner_scheme_variables: usize,
    max_common_pair_cache_entries: usize,
    max_common_pair_distance: usize,
    verbose: bool,
}

#[derive(Debug, Serialize, Deserialize)]
struct SymbolicaOptimizationSettingsRecord {
    horner_iterations: usize,
    n_cores: usize,
    cpe_iterations: Option<usize>,
    abort_level: usize,
    max_horner_scheme_variables: usize,
    max_common_pair_cache_entries: usize,
    max_common_pair_distance: usize,
    verbose: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ThreeDrepProfileRecord {
    target_timing: String,
    target_timing_seconds: f64,
    warmup_calls: usize,
    warmup_timing: String,
    warmup_timing_seconds: f64,
    calls: usize,
    total_timing: String,
    total_timing_seconds: f64,
    timing_per_sample: String,
    timing_per_sample_seconds: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum ThreeDrepBuildStrategy {
    Monolithic,
    Iterative,
}

impl ThreeDrepBuildStrategy {
    fn from_iterative_flag(iterative: bool) -> Self {
        if iterative {
            Self::Iterative
        } else {
            Self::Monolithic
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Monolithic => "monolithic",
            Self::Iterative => "iterative",
        }
    }

    fn is_iterative(self) -> bool {
        self == Self::Iterative
    }
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
    #[serde(default = "default_build_strategy")]
    build_strategy: ThreeDrepBuildStrategy,
    evaluator_method: EvaluatorMethod,
    evaluator_backend: String,
    evaluator_build_timing: Option<String>,
    evaluator_build_timing_seconds: Option<f64>,
    value: String,
    value_re: String,
    value_im: String,
    sample_evaluation_timing: String,
    sample_evaluation_timing_seconds: f64,
    #[serde(default)]
    profile: Option<ThreeDrepProfileRecord>,
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
    #[serde(default = "default_build_strategy")]
    build_strategy: ThreeDrepBuildStrategy,
    evaluator_method: EvaluatorMethod,
    evaluator_backend: String,
    compiled_backend_available: bool,
    #[serde(default)]
    force_eager: bool,
    numerator_samples_normalization: CliNumeratorSamplesNormalization,
    numerator_interpolation_scale: f64,
    seed: u64,
    scale: f64,
}

const fn default_build_strategy() -> ThreeDrepBuildStrategy {
    ThreeDrepBuildStrategy::Monolithic
}

struct SelectedGraph<'a> {
    process_id: usize,
    integrand_name: String,
    graph_id: usize,
    graph: &'a Graph,
    integrand_kind: SelectedIntegrandKind,
}

struct EvaluateInput<'a> {
    selected: SelectedGraph<'a>,
    artifact: Option<BuildOutput>,
    expression_path: Option<PathBuf>,
    artifact_dir: PathBuf,
    numerator_sampling_scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
}

#[derive(Clone, Copy)]
enum SelectedIntegrandKind {
    Amplitude,
    CrossSection,
}

struct TestCffLtdCaseBuildRequest<'a> {
    graph: &'a Graph,
    model: &'a Model,
    integrand_kind: SelectedIntegrandKind,
    workspace: &'a Path,
    representation: RepresentationMode,
    scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    energy_degree_bounds: &'a [(usize, usize)],
    run_settings: &'a ThreeDrepRunSettings,
    default_runtime_settings: &'a RuntimeSettings,
    global_cli_settings: &'a CLISettings,
    mass_shift_start: f64,
    n_epsilon_steps: usize,
}

struct StandardEvaluationRequest<'a> {
    name: &'a str,
    graph: &'a Graph,
    model: &'a Model,
    integrand_kind: SelectedIntegrandKind,
    expression: &'a ThreeDExpression<OrientationID>,
    parametric_atom: &'a symbolica::atom::Atom,
    atom_build_timing: Duration,
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
    integrand_kind: SelectedIntegrandKind,
    energy_degree_bounds: &'a [(usize, usize)],
    default_runtime_settings: &'a RuntimeSettings,
    global_cli_settings: &'a CLISettings,
    workspace: &'a Path,
    scale_mode: three_dimensional_reps::NumeratorSamplingScaleMode,
    run_settings: &'a ThreeDrepRunSettings,
    mass_shift_start: f64,
    n_epsilon_steps: usize,
}

#[derive(Clone, Copy)]
struct EvaluatorBuildContext<'a> {
    workspace: &'a Path,
    model: &'a Model,
    global_cli_settings: &'a CLISettings,
    default_runtime_settings: &'a RuntimeSettings,
    run_settings: &'a ThreeDrepRunSettings,
    clean: bool,
}

struct PreparedEvaluator {
    evaluator: EvaluatorStack,
    info: PreparedEvaluatorInfo,
}

struct PreparedComponentEvaluator {
    numerator_evaluator: EvaluatorStack,
    orientation_evaluator: EvaluatorStack,
    info: PreparedEvaluatorInfo,
}

struct DiagnosticComponentAtoms {
    processed_numerator: Atom,
    numerator_atoms: Vec<Atom>,
    orientation_atoms: Vec<Atom>,
}

#[derive(Clone)]
struct PreparedEvaluatorInfo {
    build_timing: Option<String>,
    build_timing_seconds: Option<f64>,
    backend: String,
}

impl PreparedEvaluatorInfo {
    fn include_atom_build_timing(&mut self, atom_build_timing: Duration) {
        let Some(build_timing_seconds) = self.build_timing_seconds else {
            return;
        };
        let combined =
            Duration::from_secs_f64(build_timing_seconds + atom_build_timing.as_secs_f64());
        self.build_timing = Some(format_duration_dynamic(combined));
        self.build_timing_seconds = Some(combined.as_secs_f64());
    }
}

#[derive(Serialize, Deserialize, PartialEq)]
struct EvaluatorCacheManifest {
    cache_version: u32,
    name: String,
    runtime_precision: CliRuntimePrecision,
    evaluator_method: EvaluatorMethod,
    backend: String,
    force_eager: bool,
    atom_hash: String,
    orientation_hash: String,
    parameter_hash: String,
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
        n_epsilon_steps: usize,
    ) -> bool {
        self.process_id == selected.process_id
            && self.integrand_name == selected.integrand_name
            && self.graph_id == selected.graph_id
            && self.graph_name == selected.graph.name
            && self.energy_degree_bounds == energy_degree_bounds
            && self.numerator_interpolation_scale == run_settings.numerator_interpolation_scale
            && self.mass_shift_start == mass_shift_start
            && self.n_epsilon_steps == n_epsilon_steps
            && self.settings.runtime_precision == run_settings.runtime_precision
            && self.settings.build_strategy == run_settings.build_strategy
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

    fn integrand_kind(&self) -> SelectedIntegrandKind {
        match self {
            Self::Generated(ProcessIntegrand::Amplitude(_)) | Self::ImportedAmplitude(_) => {
                SelectedIntegrandKind::Amplitude
            }
            Self::Generated(ProcessIntegrand::CrossSection(_)) | Self::ImportedCrossSection(_) => {
                SelectedIntegrandKind::CrossSection
            }
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
        let numerator_samples_normalization =
            self.numerator_samples_normalization.unwrap_or_else(|| {
                CliNumeratorSamplesNormalization::resolve_from_global(global_cli_settings)
            });
        let numerator_sampling_scale_mode = numerator_samples_normalization.to_generation_mode();
        let representation = RepresentationMode::from(self.representation);
        let mut options = selected
            .graph
            .three_d_expression_options(representation, numerator_sampling_scale_mode)?;
        let automatic_energy_degree_bounds = options.energy_degree_bounds.clone();
        let override_energy_degree_bounds =
            parse_energy_degree_bounds(self.energy_degree_bounds.as_deref())?;
        let energy_degree_bounds = merge_energy_degree_bounds(
            &automatic_energy_degree_bounds,
            &override_energy_degree_bounds,
        );
        options
            .energy_degree_bounds
            .clone_from(&energy_degree_bounds);
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
            let numerator = selected.graph.full_numerator_atom();
            let numerator_expr = numerator.log_print(Some(120));
            let simplified_numerator_expr = numerator
                .simplify_metrics()
                .simplify_color()
                .log_print(Some(120));
            println!(
                "{}",
                render_expression_summary(
                    &output.expression,
                    representation,
                    &output.graph,
                    &output.energy_degree_bounds,
                    NumeratorDisplay {
                        original: Some(&numerator_expr),
                        simplified: Some(&simplified_numerator_expr),
                    },
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
        let input = self.resolve_input(state, global_cli_settings, &workspace)?;
        if let Some(expression_path) = &input.expression_path {
            println!(
                "Loading 3Drep oriented expression from {}",
                relative_display(expression_path)
            );
        } else {
            println!(
                "3Drep evaluate --numerator-only selected graph #{} {} directly; no oriented expression JSON is required.",
                input.selected.graph_id,
                input.selected.graph.name
            );
        }
        let selected = input.selected;
        let model =
            state.resolve_model_for_integrand(selected.process_id, &selected.integrand_name)?;
        let run_settings = threedrep_run_settings(ThreeDrepRunSettingsRequest {
            global_cli_settings,
            default_runtime_settings,
            precision: self.precision,
            seed: self.seed,
            scale: self.scale,
            numerator_samples_normalization: CliNumeratorSamplesNormalization::from_generation_mode(
                input.numerator_sampling_scale_mode,
            ),
            force_eager: self.eager,
            build_strategy: ThreeDrepBuildStrategy::from_iterative_flag(self.iterative),
        })?;
        let input_config = DiagnosticInputConfig {
            numerator_only: self.numerator_only,
            numerator_q0_overrides: parse_numerator_q0_overrides(&self.numerator_q0)?,
        };
        if !self.numerator_only && !input_config.numerator_q0_overrides.is_empty() {
            return Err(eyre!(
                "--numerator-q0 can only be used together with 3Drep evaluate --numerator-only"
            ));
        }
        let evaluator_settings = threedrep_evaluator_settings(
            global_cli_settings,
            &run_settings.evaluator_method,
            run_settings.runtime_precision,
            run_settings.force_eager,
        );
        if self.standalone_rust_only && run_settings.build_strategy.is_iterative() {
            return Err(eyre!(
                "3Drep evaluate --iterative does not support --standalone-rust-only"
            ));
        }
        let numerator_only_orientations = self
            .numerator_only
            .then(|| numerator_only_orientation(selected.graph));
        let expression = input.artifact.as_ref().map(|artifact| &artifact.expression);
        let artifact_dir = input.artifact_dir;
        let symbolica_expression_path = artifact_dir.join("symbolica_expression.txt");
        let symbolica_expression_pretty_path = artifact_dir.join("symbolica_expression_pretty.txt");
        let symbolica_expression_raw_path = artifact_dir.join("symbolica_expression_raw.json");
        let symbolica_expression_raw_script_path = artifact_dir.join("symbolica_expression_raw.rs");
        let param_builder_path = artifact_dir.join("param_builder.txt");
        let evaluate_manifest_path =
            artifact_dir.join(evaluate_manifest_file_name(&self.manifest_name)?);
        let prepared_param_builder = prepare_diagnostic_param_builder(
            selected.graph,
            &model,
            selected.integrand_kind,
            &run_settings,
            default_runtime_settings,
            &input_config,
        )?;
        let mut parametric_atom: Option<Atom> = None;
        let mut component_atoms: Option<DiagnosticComponentAtoms> = None;
        let atom_build_start = Instant::now();
        let atom_build_timing;
        if run_settings.build_strategy.is_iterative() {
            component_atoms = Some(if self.numerator_only {
                numerator_only_component_atoms_for_evaluator(selected.graph)
            } else {
                diagnostic_component_atoms_for_evaluator(
                    expression.ok_or_else(|| {
                        eyre!("3Drep evaluate --iterative requires an oriented expression artifact")
                    })?,
                    selected.graph,
                )?
            });
            atom_build_timing = atom_build_start.elapsed();
            let components = component_atoms
                .as_ref()
                .expect("component atoms were just initialized");
            write_symbolica_component_expression_files(
                components,
                &symbolica_expression_path,
                &symbolica_expression_pretty_path,
            )?;
            write_iterative_raw_stub(
                &symbolica_expression_raw_path,
                &symbolica_expression_raw_script_path,
                components.numerator_atoms.len(),
            )?;
        } else {
            parametric_atom = Some(if self.numerator_only {
                numerator_only_parametric_atom_for_evaluator(selected.graph)
            } else {
                diagnostic_parametric_atom_for_evaluator(
                    expression.ok_or_else(|| {
                        eyre!("3Drep evaluate requires an oriented expression artifact")
                    })?,
                    selected.graph,
                    &evaluator_settings,
                )?
            });
            atom_build_timing = atom_build_start.elapsed();
            let parametric_atom = parametric_atom
                .as_ref()
                .expect("parametric atom was just initialized");
            if self.clean
                || !symbolica_expression_path.exists()
                || !symbolica_expression_pretty_path.exists()
            {
                write_symbolica_expression_files(
                    parametric_atom,
                    &symbolica_expression_path,
                    &symbolica_expression_pretty_path,
                )?;
            } else {
                println!(
                    "Reusing cached Symbolica expression from {}",
                    relative_display(&symbolica_expression_path)
                );
                println!(
                    "Reusing cached pretty Symbolica expression from {}",
                    relative_display(&symbolica_expression_pretty_path)
                );
            }
            if self.clean || !symbolica_raw_archive_cache_is_current(&symbolica_expression_raw_path)
            {
                let processed_atoms =
                    EvaluatorStack::spenso_processed_atoms(&[parametric_atom], &evaluator_settings)
                        .with_context(|| {
                            format!(
                                "Could not build raw Symbolica evaluator input for {}",
                                relative_display(&symbolica_expression_raw_path)
                            )
                        })?;
                let raw_archive = symbolica_evaluator_input_archive(
                    &processed_atoms,
                    &prepared_param_builder.param_builder,
                    diagnostic_evaluation_orientations_for_raw_input(
                        expression,
                        numerator_only_orientations.as_ref(),
                    )?,
                    &evaluator_settings,
                    &run_settings.evaluator_method,
                    &run_settings.evaluator_backend,
                )?;
                write_path(
                    &symbolica_expression_raw_path,
                    &serde_json::to_string_pretty(&raw_archive)?,
                )?;
            } else {
                println!(
                    "Reusing cached raw Symbolica evaluator input from {}",
                    relative_display(&symbolica_expression_raw_path)
                );
            }
            write_path(
                &symbolica_expression_raw_script_path,
                &symbolica_expression_raw_rust_script(),
            )?;
            make_script_executable(&symbolica_expression_raw_script_path)?;
        }
        if self.clean || !param_builder_path.exists() {
            write_path(
                &param_builder_path,
                &prepared_param_builder.param_builder.to_string(),
            )?;
        } else {
            println!(
                "Reusing cached parameter builder from {}",
                relative_display(&param_builder_path)
            );
        }

        if self.standalone_rust_only {
            println!(
                "Saved raw Symbolica evaluator input to {}",
                relative_display(&symbolica_expression_raw_path)
            );
            println!(
                "Saved raw Symbolica evaluator replay script to {}",
                relative_display(&symbolica_expression_raw_script_path)
            );
            println!(
                "3Drep evaluate --standalone-rust-only selected; skipped evaluator build, evaluation, and manifest writing."
            );
            return Ok(());
        }

        let profile = self
            .profile
            .as_deref()
            .map(parse_profile_target_duration)
            .transpose()?;
        let label = if self.numerator_only {
            "numerator_only"
        } else {
            "evaluate"
        };
        let evaluation = if let Some(components) = component_atoms.as_ref() {
            let component_orientations;
            let orientations = if let Some(orientations) = numerator_only_orientations.as_ref() {
                orientations
            } else {
                component_orientations =
                    diagnostic_evaluation_orientations(expression.ok_or_else(|| {
                        eyre!("3Drep evaluate --iterative requires an oriented expression artifact")
                    })?);
                &component_orientations
            };
            evaluate_threedrep_component_expression(ComponentEvaluationRequest {
                label,
                graph: selected.graph,
                model: &model,
                integrand_kind: selected.integrand_kind,
                components,
                orientations,
                workspace: &artifact_dir,
                global_cli_settings,
                default_runtime_settings,
                representation: input.artifact.as_ref().map(|artifact| artifact.family),
                numerator_sampling_scale_mode: (!self.numerator_only)
                    .then_some(input.numerator_sampling_scale_mode),
                run_settings: &run_settings,
                mass_shift: "none",
                mass_shift_values: &[],
                profile_target: profile,
                input_config: &input_config,
                clean: self.clean,
                atom_build_timing,
            })?
        } else {
            evaluate_threedrep_expression(EvaluationRequest {
                label,
                graph: selected.graph,
                model: &model,
                integrand_kind: selected.integrand_kind,
                expression,
                parametric_atom: parametric_atom
                    .as_ref()
                    .expect("parametric atom was initialized for monolithic evaluation"),
                orientations: numerator_only_orientations.as_ref(),
                workspace: &artifact_dir,
                global_cli_settings,
                default_runtime_settings,
                representation: input.artifact.as_ref().map(|artifact| artifact.family),
                numerator_sampling_scale_mode: (!self.numerator_only)
                    .then_some(input.numerator_sampling_scale_mode),
                run_settings: &run_settings,
                mass_shift: "none",
                mass_shift_values: &[],
                profile_target: profile,
                input_config: &input_config,
                clean: self.clean,
                atom_build_timing,
            })?
        };
        let parameters = parameter_records(
            selected.graph,
            &model,
            selected.integrand_kind,
            &run_settings,
            default_runtime_settings,
            &input_config,
        )?;

        let summary = EvaluateOutput {
            process_id: selected.process_id,
            integrand_name: selected.integrand_name,
            graph_id: selected.graph_id,
            graph_name: selected.graph.name.clone(),
            numerator_only: self.numerator_only,
            expression_path: input.expression_path,
            symbolica_expression_path,
            symbolica_expression_pretty_path,
            symbolica_expression_raw_path,
            symbolica_expression_raw_script_path,
            param_builder_path,
            settings: run_settings,
            evaluation,
            parameters,
        };
        write_path(
            &evaluate_manifest_path,
            &serde_json::to_string_pretty(&summary)?,
        )?;
        println!(
            "{}",
            render_evaluate_summary(&summary, !self.no_show_parameters)
        );
        println!(
            "Saved 3Drep evaluate summary to {}",
            relative_display(&evaluate_manifest_path)
        );
        Ok(())
    }

    fn resolve_input<'a>(
        &self,
        state: &'a State,
        global_cli_settings: &CLISettings,
        workspace: &Path,
    ) -> Result<EvaluateInput<'a>> {
        if self.numerator_only && self.selection.has_any_selector() {
            let selector = self
                .selection
                .require_graph_selection("3Drep evaluate --numerator-only")?;
            let selected = select_graph(state, &selector)?;
            if self.representation.is_some() || self.json_in.is_some() {
                println!(
                    "3Drep evaluate --numerator-only selected a graph directly; ignoring oriented-expression lookup options."
                );
            }
            let numerator_sampling_scale_mode = self
                .numerator_samples_normalization
                .unwrap_or_else(|| {
                    CliNumeratorSamplesNormalization::resolve_from_global(global_cli_settings)
                })
                .to_generation_mode();
            let artifact_dir = graph_workspace_dir(workspace, &selected).join("numerator_only");
            return Ok(EvaluateInput {
                selected,
                artifact: None,
                expression_path: None,
                artifact_dir,
                numerator_sampling_scale_mode,
            });
        }

        let (expression_path, artifact) =
            self.load_oriented_expression_artifact(state, global_cli_settings, workspace)?;
        let selected = select_graph_from_artifact(state, &artifact)?;
        let artifact_dir = if self.numerator_only {
            graph_workspace_dir(workspace, &selected).join("numerator_only")
        } else {
            expression_path
                .parent()
                .map(Path::to_path_buf)
                .unwrap_or_else(|| workspace.to_path_buf())
        };
        let numerator_sampling_scale_mode = artifact.numerator_sampling_scale_mode;
        Ok(EvaluateInput {
            selected,
            artifact: Some(artifact),
            expression_path: Some(expression_path),
            artifact_dir,
            numerator_sampling_scale_mode,
        })
    }

    fn load_oriented_expression_artifact(
        &self,
        state: &State,
        global_cli_settings: &CLISettings,
        workspace: &Path,
    ) -> Result<(PathBuf, BuildOutput)> {
        let expression_path =
            self.resolve_expression_path(state, global_cli_settings, workspace)?;
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
        Ok((expression_path, artifact))
    }

    fn resolve_expression_path(
        &self,
        state: &State,
        global_cli_settings: &CLISettings,
        workspace: &Path,
    ) -> Result<PathBuf> {
        if let Some(path) = &self.json_in {
            if self.representation.is_some()
                || self.selection.has_any_selector()
                || self.numerator_samples_normalization.is_some()
            {
                println!(
                    "3Drep evaluate --json-in supplied; using {} and ignoring cached representation lookup options.",
                    relative_display(path)
                );
            }
            return Ok(path.clone());
        }

        if let Some(representation) = self.representation {
            let selector = self
                .selection
                .require_graph_selection("3Drep evaluate --representation")?;
            let selected = select_graph(state, &selector)?;
            let scale_mode = self
                .numerator_samples_normalization
                .unwrap_or_else(|| {
                    CliNumeratorSamplesNormalization::resolve_from_global(global_cli_settings)
                })
                .to_generation_mode();
            let expression_path = build_artifact_dir(
                workspace,
                &selected,
                RepresentationMode::from(representation),
                scale_mode,
            )
            .join("oriented_expression.json");
            if expression_path.exists() {
                Ok(expression_path)
            } else {
                Err(eyre!(
                    "No cached 3Drep {representation:?} oriented expression found at {}. Run `3Drep build -p ... -i ... -g ... --representation {}` first, or pass --json-in.",
                    expression_path.display(),
                    cli_representation_name(representation)
                ))
            }
        } else if self.selection.has_any_selector()
            || self.numerator_samples_normalization.is_some()
        {
            Err(eyre!(
                "3Drep evaluate graph selection and --numerator-samples-normalization require --representation, or pass --json-in"
            ))
        } else if workspace.join("oriented_expression.json").exists() {
            Ok(workspace.join("oriented_expression.json"))
        } else {
            read_latest_expression_path(workspace)
        }
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
        let automatic_energy_degree_bounds = selected
            .graph
            .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;
        let override_energy_degree_bounds =
            parse_energy_degree_bounds(self.energy_degree_bounds.as_deref())?;
        let energy_degree_bounds = merge_energy_degree_bounds(
            &automatic_energy_degree_bounds,
            &override_energy_degree_bounds,
        );
        let run_settings = threedrep_run_settings(ThreeDrepRunSettingsRequest {
            global_cli_settings,
            default_runtime_settings,
            precision: self.precision,
            seed: self.seed,
            scale: self.scale,
            numerator_samples_normalization: CliNumeratorSamplesNormalization::resolve_from_global(
                global_cli_settings,
            ),
            force_eager: false,
            build_strategy: ThreeDrepBuildStrategy::from_iterative_flag(self.iterative),
        })?;
        let mass_shift_start = self.mass_shift.unwrap_or(self.scale);
        if self.n_epsilon_steps == 0 {
            return Err(eyre!(
                "3Drep test-cff-ltd --n-epsilon-steps must be at least 1"
            ));
        }

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
                    self.n_epsilon_steps,
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
                    integrand_kind: selected.integrand_kind,
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
                    n_epsilon_steps: self.n_epsilon_steps,
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
                n_epsilon_steps: self.n_epsilon_steps,
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
            integrand_kind,
            workspace,
            representation,
            scale_mode,
            energy_degree_bounds,
            run_settings,
            default_runtime_settings,
            global_cli_settings,
            mass_shift_start,
            n_epsilon_steps,
        } = request;
        let name = representation_label(representation).to_lowercase();
        let case_dir = workspace.join("test_cff_ltd").join(&name);
        let expression_path = case_dir.join("oriented_expression.json");
        let symbolica_expression_path = case_dir.join("symbolica_expression.txt");
        let symbolica_expression_pretty_path = case_dir.join("symbolica_expression_pretty.txt");
        let mut options = graph.three_d_expression_options(representation, scale_mode)?;
        options.energy_degree_bounds = energy_degree_bounds.to_vec();
        let expression = match generate_3d_expression(graph, &options) {
            Ok(expression) => expression,
            Err(error) => {
                let generation_error_path = case_dir.join("generation_error.txt");
                write_path(&generation_error_path, &format_diagnostic_error(&error))?;
                return Ok(TestCffLtdCaseOutput {
                    name,
                    representation,
                    numerator_sampling_scale_mode: scale_mode,
                    generation_status: "failed".to_string(),
                    orientation_count: 0,
                    unfolded_term_count: 0,
                    expression_path,
                    symbolica_expression_path,
                    symbolica_expression_pretty_path,
                    evaluator_build_status: "skipped: generation failed".to_string(),
                    evaluator_build_timing: None,
                    evaluator_build_timing_seconds: None,
                    evaluations: Vec::new(),
                    generation_error_path: Some(generation_error_path),
                });
            }
        };
        let evaluator_settings = threedrep_evaluator_settings(
            global_cli_settings,
            &run_settings.evaluator_method,
            run_settings.runtime_precision,
            run_settings.force_eager,
        );
        let mut parametric_atom: Option<Atom> = None;
        let mut component_atoms: Option<DiagnosticComponentAtoms> = None;
        let atom_build_start = Instant::now();
        let atom_build_timing;
        if run_settings.build_strategy.is_iterative() {
            component_atoms = Some(diagnostic_component_atoms_for_evaluator(
                &expression,
                graph,
            )?);
        } else {
            parametric_atom = Some(diagnostic_parametric_atom_for_evaluator(
                &expression,
                graph,
                &evaluator_settings,
            )?);
        }
        atom_build_timing = atom_build_start.elapsed();
        write_path(
            &expression_path,
            &serde_json::to_string_pretty(&expression)?,
        )?;
        if let Some(components) = component_atoms.as_ref() {
            write_symbolica_component_expression_files(
                components,
                &symbolica_expression_path,
                &symbolica_expression_pretty_path,
            )?;
        } else {
            write_symbolica_expression_files(
                parametric_atom
                    .as_ref()
                    .expect("parametric atom was initialized for monolithic test"),
                &symbolica_expression_path,
                &symbolica_expression_pretty_path,
            )?;
        }

        let orientations = diagnostic_evaluation_orientations(&expression);
        let has_repeated_propagators =
            !three_dimensional_reps::repeated_groups(&graph.to_three_d_parsed_graph()?).is_empty();
        let (evaluator_build_status, evaluations) = if representation == RepresentationMode::PureLtd
            && has_repeated_propagators
        {
            self.build_pure_ltd_mass_shift_evaluations(PureLtdMassShiftEvaluationRequest {
                name: &name,
                graph,
                model,
                integrand_kind,
                workspace: &case_dir,
                energy_degree_bounds,
                default_runtime_settings,
                global_cli_settings,
                scale_mode,
                run_settings,
                mass_shift_start,
                n_epsilon_steps,
            })?
        } else if let Some(components) = component_atoms.as_ref() {
            let input_config = DiagnosticInputConfig::default();
            let evaluation = evaluate_threedrep_component_expression(ComponentEvaluationRequest {
                label: &name,
                graph,
                model,
                integrand_kind,
                components,
                atom_build_timing,
                orientations: &orientations,
                workspace: &case_dir,
                global_cli_settings,
                default_runtime_settings,
                representation: Some(representation),
                numerator_sampling_scale_mode: Some(scale_mode),
                run_settings,
                mass_shift: "none",
                mass_shift_values: &[],
                profile_target: None,
                input_config: &input_config,
                clean: self.clean,
            });
            match evaluation {
                Ok(evaluation) => ("ok".to_string(), vec![evaluation]),
                Err(error) => (
                    format!("failed: {}", format_diagnostic_error(&error)),
                    Vec::new(),
                ),
            }
        } else {
            self.build_standard_evaluations(StandardEvaluationRequest {
                name: &name,
                graph,
                model,
                integrand_kind,
                expression: &expression,
                parametric_atom: parametric_atom
                    .as_ref()
                    .expect("parametric atom was initialized for monolithic test"),
                atom_build_timing,
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
            symbolica_expression_pretty_path,
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
        let input_config = DiagnosticInputConfig::default();
        let prepared_param_builder = prepare_diagnostic_param_builder(
            request.graph,
            request.model,
            request.integrand_kind,
            request.run_settings,
            request.default_runtime_settings,
            &input_config,
        )?;
        let prepared = match build_threedrep_evaluator(
            request.name,
            EvaluatorBuildContext {
                workspace: request.workspace,
                model: request.model,
                global_cli_settings: request.global_cli_settings,
                default_runtime_settings: request.default_runtime_settings,
                run_settings: request.run_settings,
                clean: self.clean,
            },
            &prepared_param_builder.param_builder,
            &[request.parametric_atom],
            request.orientations,
        ) {
            Ok(prepared) => prepared,
            Err(error) => {
                return Ok((
                    format!("failed: {}", format_diagnostic_error(&error)),
                    Vec::new(),
                ));
            }
        };
        let mut evaluator_info = prepared.info.clone();
        evaluator_info.include_atom_build_timing(request.atom_build_timing);
        let mut evaluator = prepared.evaluator;

        let evaluations = vec![evaluate_with_evaluator(
            &mut evaluator,
            EvaluationRequest {
                label: request.name,
                graph: request.graph,
                model: request.model,
                integrand_kind: request.integrand_kind,
                expression: Some(request.expression),
                parametric_atom: request.parametric_atom,
                atom_build_timing: request.atom_build_timing,
                orientations: None,
                workspace: request.workspace,
                global_cli_settings: request.global_cli_settings,
                default_runtime_settings: request.default_runtime_settings,
                representation: Some(request.representation),
                numerator_sampling_scale_mode: Some(request.scale_mode),
                run_settings: request.run_settings,
                mass_shift: "none",
                mass_shift_values: &[],
                profile_target: None,
                input_config: &input_config,
                clean: self.clean,
            },
            request.orientations,
            &evaluator_info,
            prepared_param_builder,
        )?];
        Ok(("ok".to_string(), evaluations))
    }

    fn build_pure_ltd_mass_shift_evaluations(
        &self,
        request: PureLtdMassShiftEvaluationRequest<'_>,
    ) -> Result<(String, Vec<ThreeDrepEvaluationRecord>)> {
        let mut evaluations = Vec::new();
        for epsilon in
            pure_ltd_mass_shift_epsilons(request.mass_shift_start, request.n_epsilon_steps)
        {
            let (shifted_graph, mass_shifts) = request
                .graph
                .split_repeated_masses_for_three_drep(request.model, epsilon)?;
            let mass_shift_values = mass_shift_value_records(&mass_shifts);
            let mass_shift = mass_shift_label(epsilon, &mass_shift_values);
            let mut shifted_options = shifted_graph
                .three_d_expression_options(RepresentationMode::Ltd, request.scale_mode)?;
            shifted_options.energy_degree_bounds = request.energy_degree_bounds.to_vec();
            let shifted_expression = generate_3d_expression(&shifted_graph, &shifted_options)
                .with_context(|| {
                    format!(
                    "while generating split-mass LTD expression for pure-LTD mass shift {epsilon}"
                )
                })?;
            let evaluator_settings = threedrep_evaluator_settings(
                request.global_cli_settings,
                &request.run_settings.evaluator_method,
                request.run_settings.runtime_precision,
                request.run_settings.force_eager,
            );
            let orientations = diagnostic_evaluation_orientations(&shifted_expression);
            let input_config = DiagnosticInputConfig::default();
            if request.run_settings.build_strategy.is_iterative() {
                let atom_build_start = Instant::now();
                let components =
                    diagnostic_component_atoms_for_evaluator(&shifted_expression, &shifted_graph)?;
                let atom_build_timing = atom_build_start.elapsed();
                evaluations.push(evaluate_threedrep_component_expression(
                    ComponentEvaluationRequest {
                        label: request.name,
                        graph: &shifted_graph,
                        model: request.model,
                        integrand_kind: request.integrand_kind,
                        components: &components,
                        atom_build_timing,
                        orientations: &orientations,
                        workspace: request.workspace,
                        global_cli_settings: request.global_cli_settings,
                        default_runtime_settings: request.default_runtime_settings,
                        representation: Some(RepresentationMode::PureLtd),
                        numerator_sampling_scale_mode: Some(request.scale_mode),
                        run_settings: request.run_settings,
                        mass_shift: &mass_shift,
                        mass_shift_values: &mass_shift_values,
                        profile_target: None,
                        input_config: &input_config,
                        clean: self.clean,
                    },
                )?);
            } else {
                let atom_build_start = Instant::now();
                let parametric_atom = diagnostic_parametric_atom_for_evaluator(
                    &shifted_expression,
                    &shifted_graph,
                    &evaluator_settings,
                )?;
                let atom_build_timing = atom_build_start.elapsed();
                let prepared_param_builder = prepare_diagnostic_param_builder(
                    &shifted_graph,
                    request.model,
                    request.integrand_kind,
                    request.run_settings,
                    request.default_runtime_settings,
                    &input_config,
                )?;
                let prepared = match build_threedrep_evaluator(
                    &format!("{}_{}", request.name, mass_shift_file_label(epsilon)),
                    EvaluatorBuildContext {
                        workspace: request.workspace,
                        model: request.model,
                        global_cli_settings: request.global_cli_settings,
                        default_runtime_settings: request.default_runtime_settings,
                        run_settings: request.run_settings,
                        clean: self.clean,
                    },
                    &prepared_param_builder.param_builder,
                    &[&parametric_atom],
                    &orientations,
                ) {
                    Ok(prepared) => prepared,
                    Err(error) => {
                        return Ok((
                            format!("failed: {}", format_diagnostic_error(&error)),
                            evaluations,
                        ));
                    }
                };
                let mut evaluator_info = prepared.info.clone();
                evaluator_info.include_atom_build_timing(atom_build_timing);
                let mut evaluator = prepared.evaluator;
                evaluations.push(evaluate_with_evaluator(
                    &mut evaluator,
                    EvaluationRequest {
                        label: request.name,
                        graph: &shifted_graph,
                        model: request.model,
                        integrand_kind: request.integrand_kind,
                        expression: Some(&shifted_expression),
                        parametric_atom: &parametric_atom,
                        atom_build_timing,
                        orientations: None,
                        workspace: request.workspace,
                        global_cli_settings: request.global_cli_settings,
                        default_runtime_settings: request.default_runtime_settings,
                        representation: Some(RepresentationMode::PureLtd),
                        numerator_sampling_scale_mode: Some(request.scale_mode),
                        run_settings: request.run_settings,
                        mass_shift: &mass_shift,
                        mass_shift_values: &mass_shift_values,
                        profile_target: None,
                        input_config: &input_config,
                        clean: self.clean,
                    },
                    &orientations,
                    &evaluator_info,
                    prepared_param_builder,
                )?);
            }
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

struct ThreeDrepRunSettingsRequest<'a> {
    global_cli_settings: &'a CLISettings,
    default_runtime_settings: &'a RuntimeSettings,
    precision: Option<CliRuntimePrecision>,
    seed: u64,
    scale: f64,
    numerator_samples_normalization: CliNumeratorSamplesNormalization,
    force_eager: bool,
    build_strategy: ThreeDrepBuildStrategy,
}

fn threedrep_run_settings(
    request: ThreeDrepRunSettingsRequest<'_>,
) -> Result<ThreeDrepRunSettings> {
    let ThreeDrepRunSettingsRequest {
        global_cli_settings,
        default_runtime_settings,
        precision,
        seed,
        scale,
        numerator_samples_normalization,
        force_eager,
        build_strategy,
    } = request;
    if !scale.is_finite() {
        return Err(eyre!("3Drep evaluation --scale must be finite"));
    }
    let runtime_precision = CliRuntimePrecision::resolve(default_runtime_settings, precision);
    let evaluator_method = default_runtime_settings.general.evaluator_method.clone();
    let evaluator_settings = threedrep_evaluator_settings(
        global_cli_settings,
        &evaluator_method,
        runtime_precision,
        force_eager,
    );
    let frozen_mode = if runtime_precision == CliRuntimePrecision::Double && !force_eager {
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
        build_strategy,
        evaluator_method,
        evaluator_backend: frozen_mode.active_backend_name().to_string(),
        compiled_backend_available: frozen_mode.compile_enabled(),
        force_eager,
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
    force_eager: bool,
) -> EvaluatorSettings {
    let mut settings = global_cli_settings.global.generation.evaluator;
    settings.iterative_orientation_optimization =
        matches!(evaluator_method, EvaluatorMethod::Iterative);
    settings.summed_function_map = matches!(evaluator_method, EvaluatorMethod::SummedFunctionMap);
    settings.summed = matches!(evaluator_method, EvaluatorMethod::Summed);
    if runtime_precision != CliRuntimePrecision::Double || force_eager {
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
    param_builder: &ParamBuilder<f64>,
    atoms: &[A],
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
) -> Result<PreparedEvaluator> {
    let evaluator_settings = threedrep_evaluator_settings(
        context.global_cli_settings,
        &context.run_settings.evaluator_method,
        context.run_settings.runtime_precision,
        context.run_settings.force_eager,
    );
    let frozen_mode = if context.run_settings.runtime_precision == CliRuntimePrecision::Double
        && !context.run_settings.force_eager
    {
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
    let cache_manifest = EvaluatorCacheManifest {
        cache_version: 1,
        name: name.to_string(),
        runtime_precision: context.run_settings.runtime_precision,
        evaluator_method: context.run_settings.evaluator_method.clone(),
        backend: frozen_mode.active_backend_name().to_string(),
        force_eager: context.run_settings.force_eager,
        atom_hash: stable_hash_hex(
            atoms
                .iter()
                .map(|atom| atom.as_atom_view().to_canonical_string()),
        ),
        orientation_hash: orientation_hash(orientations),
        parameter_hash: parameter_hash(param_builder),
    };
    if !context.clean {
        if let Some(prepared) =
            try_load_cached_threedrep_evaluator(name, context, &frozen_mode, &cache_manifest)?
        {
            return Ok(prepared);
        }
    }

    let start = Instant::now();
    let mut evaluator = EvaluatorStack::new(
        atoms,
        param_builder,
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
    save_threedrep_evaluator_cache(&evaluator_dir, name, &evaluator, &cache_manifest)?;

    Ok(PreparedEvaluator {
        evaluator,
        info: PreparedEvaluatorInfo {
            build_timing: Some(format_duration_dynamic(build_timing)),
            build_timing_seconds: Some(build_timing.as_secs_f64()),
            backend: frozen_mode.active_backend_name().to_string(),
        },
    })
}

fn build_threedrep_component_evaluator(
    name: &str,
    context: EvaluatorBuildContext<'_>,
    param_builder: &ParamBuilder<f64>,
    components: &DiagnosticComponentAtoms,
) -> Result<PreparedComponentEvaluator> {
    if components.numerator_atoms.len() != components.orientation_atoms.len() {
        return Err(eyre!(
            "Iterative 3Drep component mismatch: {} numerator components but {} orientation components",
            components.numerator_atoms.len(),
            components.orientation_atoms.len()
        ));
    }
    let component_count = components.numerator_atoms.len();
    if component_count == 0 {
        return Err(eyre!(
            "Iterative 3Drep evaluator needs at least one component"
        ));
    }

    let evaluator_settings = threedrep_evaluator_settings(
        context.global_cli_settings,
        &context.run_settings.evaluator_method,
        context.run_settings.runtime_precision,
        context.run_settings.force_eager,
    );
    let frozen_mode = if context.run_settings.runtime_precision == CliRuntimePrecision::Double
        && !context.run_settings.force_eager
    {
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
    let mut numerator_evaluator = build_preprocessed_component_stack(
        &format!("{name} numerator components"),
        &components.numerator_atoms,
        param_builder,
        &evaluator_settings,
    )?;
    let mut orientation_evaluator = build_preprocessed_component_stack(
        &format!("{name} orientation components"),
        &components.orientation_atoms,
        param_builder,
        &evaluator_settings,
    )?;
    numerator_evaluator.prepare_f64_backend(
        format!("{name}_numerator_components"),
        &evaluator_dir,
        &frozen_mode,
    )?;
    orientation_evaluator.prepare_f64_backend(
        format!("{name}_orientation_components"),
        &evaluator_dir,
        &frozen_mode,
    )?;
    let build_timing = start.elapsed();

    Ok(PreparedComponentEvaluator {
        numerator_evaluator,
        orientation_evaluator,
        info: PreparedEvaluatorInfo {
            build_timing: Some(format_duration_dynamic(build_timing)),
            build_timing_seconds: Some(build_timing.as_secs_f64()),
            backend: frozen_mode.active_backend_name().to_string(),
        },
    })
}

fn build_preprocessed_component_stack(
    label: &str,
    atoms: &[Atom],
    param_builder: &ParamBuilder<f64>,
    evaluator_settings: &EvaluatorSettings,
) -> Result<EvaluatorStack> {
    let bar = component_build_progress_bar(label, atoms.len());
    let bar_for_progress = bar.clone();
    let (evaluator, _) = EvaluatorStack::new_preprocessed_components_with_timings(
        atoms,
        param_builder,
        None,
        evaluator_settings,
        move |done, _total| {
            bar_for_progress.set_position(done as u64);
        },
    )?;
    bar.finish_and_clear();
    Ok(evaluator)
}

fn component_build_progress_bar(label: &str, len: usize) -> ProgressBar {
    let bar = ProgressBar::new(len as u64);
    let style = ProgressStyle::with_template(
        "{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len}",
    )
    .unwrap_or_else(|_| ProgressStyle::default_bar());
    bar.set_style(style);
    bar.set_message(label.to_string());
    bar
}

fn threedrep_evaluator_cache_paths(evaluator_dir: &Path, name: &str) -> (PathBuf, PathBuf) {
    (
        evaluator_dir.join(format!("{name}.evaluator.bin")),
        evaluator_dir.join(format!("{name}.evaluator_manifest.json")),
    )
}

fn try_load_cached_threedrep_evaluator(
    name: &str,
    context: EvaluatorBuildContext<'_>,
    frozen_mode: &FrozenCompilationMode,
    expected_manifest: &EvaluatorCacheManifest,
) -> Result<Option<PreparedEvaluator>> {
    let evaluator_dir = context.workspace.join("evaluators");
    let (evaluator_path, manifest_path) = threedrep_evaluator_cache_paths(&evaluator_dir, name);
    if !evaluator_path.exists() || !manifest_path.exists() {
        return Ok(None);
    }

    let manifest_text = fs::read_to_string(&manifest_path).with_context(|| {
        format!(
            "Could not read 3Drep evaluator cache manifest at {}",
            manifest_path.display()
        )
    })?;
    let manifest: EvaluatorCacheManifest =
        serde_json::from_str(&manifest_text).with_context(|| {
            format!(
                "Could not parse 3Drep evaluator cache manifest at {}",
                manifest_path.display()
            )
        })?;
    if &manifest != expected_manifest {
        return Ok(None);
    }

    let binary = fs::read(&evaluator_path).with_context(|| {
        format!(
            "Could not read 3Drep evaluator cache at {}",
            evaluator_path.display()
        )
    })?;
    let mut state_bytes = Vec::new();
    SymbolicaState::export(&mut state_bytes)
        .map_err(|error| eyre!("Could not export current Symbolica state: {error}"))?;
    let mut cursor = Cursor::new(state_bytes);
    let state_map = SymbolicaState::import(&mut cursor, None)
        .map_err(|error| eyre!("Could not import current Symbolica state: {error}"))?;
    let decode_context = GammaLoopContextContainer {
        state_map: &state_map,
        model: context.model,
    };
    let (mut evaluator, _): (EvaluatorStack, usize) = bincode::decode_from_slice_with_context(
        &binary,
        bincode::config::standard(),
        decode_context,
    )
    .with_context(|| {
        format!(
            "Could not decode 3Drep evaluator cache at {}",
            evaluator_path.display()
        )
    })?;
    evaluator
        .activate_cached_f64_backend(frozen_mode)
        .with_context(|| {
            format!(
                "Could not activate cached 3Drep evaluator backend from {}",
                evaluator_path.display()
            )
        })?;
    println!(
        "Reusing cached 3Drep evaluator from {}",
        relative_display(&evaluator_path)
    );
    println!(
        "Reusing cached 3Drep evaluator manifest from {}",
        relative_display(&manifest_path)
    );

    Ok(Some(PreparedEvaluator {
        evaluator,
        info: PreparedEvaluatorInfo {
            build_timing: None,
            build_timing_seconds: None,
            backend: frozen_mode.active_backend_name().to_string(),
        },
    }))
}

fn save_threedrep_evaluator_cache(
    evaluator_dir: &Path,
    name: &str,
    evaluator: &EvaluatorStack,
    manifest: &EvaluatorCacheManifest,
) -> Result<()> {
    let (evaluator_path, manifest_path) = threedrep_evaluator_cache_paths(evaluator_dir, name);
    let binary = bincode::encode_to_vec(evaluator, bincode::config::standard())
        .with_context(|| format!("Could not encode 3Drep evaluator cache for {name}"))?;
    fs::write(&evaluator_path, binary).with_context(|| {
        format!(
            "Could not write 3Drep evaluator cache to {}",
            evaluator_path.display()
        )
    })?;
    write_path(&manifest_path, &serde_json::to_string_pretty(manifest)?)?;
    Ok(())
}

fn stable_hash_hex(parts: impl IntoIterator<Item = String>) -> String {
    let mut hash = 0xcbf29ce484222325_u64;
    for part in parts {
        for byte in part.as_bytes().iter().copied().chain(std::iter::once(0)) {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    format!("{hash:016x}")
}

fn orientation_hash(orientations: &TiVec<OrientationID, EdgeVec<Orientation>>) -> String {
    stable_hash_hex(orientations.iter().map(|orientation| {
        orientation
            .iter()
            .map(|(edge, orientation)| format!("{}:{orientation:?}", edge.0))
            .collect::<Vec<_>>()
            .join(",")
    }))
}

fn parameter_hash(param_builder: &ParamBuilder<f64>) -> String {
    stable_hash_hex(
        param_builder
            .evaluator_input_parameters()
            .into_iter()
            .map(|parameter| {
                format!(
                    "{:?}:{}:{}",
                    parameter.group,
                    parameter.index,
                    parameter.atom.to_canonical_string()
                )
            }),
    )
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

fn evaluate_manifest_file_name(name: &str) -> Result<String> {
    let trimmed = name.trim();
    if trimmed.is_empty() {
        return Err(eyre!("3Drep evaluate --manifest-name cannot be empty"));
    }
    if trimmed.contains('/') || trimmed.contains('\\') {
        return Err(eyre!(
            "3Drep evaluate --manifest-name must be a file name, not a path"
        ));
    }
    if trimmed == "." || trimmed == ".." {
        return Err(eyre!(
            "3Drep evaluate --manifest-name must be a regular file name"
        ));
    }
    if trimmed.ends_with(".json") {
        Ok(trimmed.to_string())
    } else {
        Ok(format!("{trimmed}.json"))
    }
}

const DIAGNOSTIC_ERROR_HEAD_LINES: usize = 200;
const DIAGNOSTIC_ERROR_TAIL_LINES: usize = 200;

fn format_diagnostic_error(error: &impl std::fmt::Debug) -> String {
    truncate_middle_lines(
        &format!("{error:?}"),
        DIAGNOSTIC_ERROR_HEAD_LINES,
        DIAGNOSTIC_ERROR_TAIL_LINES,
    )
}

fn truncate_middle_lines(text: &str, head_lines: usize, tail_lines: usize) -> String {
    let lines = text.lines().collect::<Vec<_>>();
    let retained_lines = head_lines + tail_lines;
    if lines.len() <= retained_lines {
        return text.to_string();
    }

    let omitted = lines.len() - retained_lines;
    let mut truncated = lines[..head_lines].join("\n");
    truncated.push_str(&format!(
        "\n... [truncated {omitted} diagnostic line(s)] ...\n"
    ));
    truncated.push_str(&lines[lines.len() - tail_lines..].join("\n"));
    truncated
}

fn parse_profile_target_duration(input: &str) -> Result<Duration> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(eyre!(
            "3Drep evaluate --profile requires a non-empty duration"
        ));
    }

    let lower = trimmed.to_ascii_lowercase();
    let (number, multiplier) = if let Some(number) = lower.strip_suffix("microseconds") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("microsecond") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("usec") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("us") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("µs") {
        (number, 1.0e-6)
    } else if let Some(number) = lower.strip_suffix("milliseconds") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("millisecond") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("msec") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("ms") {
        (number, 1.0e-3)
    } else if let Some(number) = lower.strip_suffix("seconds") {
        (number, 1.0)
    } else if let Some(number) = lower.strip_suffix("second") {
        (number, 1.0)
    } else if let Some(number) = lower.strip_suffix("sec") {
        (number, 1.0)
    } else if let Some(number) = lower.strip_suffix('s') {
        (number, 1.0)
    } else {
        (trimmed, 1.0)
    };
    let seconds =
        number.trim().parse::<f64>().with_context(|| {
            format!("Could not parse 3Drep evaluate --profile duration '{input}'")
        })? * multiplier;
    if !seconds.is_finite() || seconds <= 0.0 {
        return Err(eyre!(
            "3Drep evaluate --profile duration must be positive and finite"
        ));
    }
    Ok(Duration::from_secs_f64(seconds))
}

fn parse_numerator_q0_overrides(inputs: &[String]) -> Result<BTreeMap<EdgeIndex, f64>> {
    let mut overrides = BTreeMap::new();
    for raw_input in inputs {
        for entry in raw_input
            .split(',')
            .map(str::trim)
            .filter(|entry| !entry.is_empty())
        {
            let (edge, value) = entry
                .split_once(':')
                .or_else(|| entry.split_once('='))
                .ok_or_else(|| {
                    eyre!(
                        "Could not parse --numerator-q0 entry '{entry}'. Expected edge:value or edge=value."
                    )
                })?;
            let edge_id = edge.trim().parse::<usize>().with_context(|| {
                format!("Could not parse numerator-only energy edge id in '{entry}'")
            })?;
            let value = value.trim().parse::<f64>().with_context(|| {
                format!("Could not parse numerator-only energy value in '{entry}'")
            })?;
            if overrides.insert(EdgeIndex(edge_id), value).is_some() {
                return Err(eyre!(
                    "Duplicate --numerator-q0 override for edge {edge_id}"
                ));
            }
        }
    }
    Ok(overrides)
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

const DIRECT_ABSOLUTE_TOLERANCE: f64 = 1.0e-12;
const MASS_SHIFT_ABSOLUTE_TOLERANCE: f64 = 1.0e-8;

fn direct_relative_tolerance(precision: CliRuntimePrecision) -> f64 {
    match precision {
        CliRuntimePrecision::Double => 1.0e-7,
        CliRuntimePrecision::Quad => 1.0e-16,
        CliRuntimePrecision::ArbPrec => 1.0e-150,
    }
}

fn mass_shift_relative_tolerance(precision: CliRuntimePrecision) -> f64 {
    match precision {
        CliRuntimePrecision::Double => 1.0e-2,
        CliRuntimePrecision::Quad => 1.0e-3,
        CliRuntimePrecision::ArbPrec => 1.0e-4,
    }
}

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
    let direct_relative_tolerance = direct_relative_tolerance(output.settings.runtime_precision);
    let mass_shift_relative_tolerance =
        mass_shift_relative_tolerance(output.settings.runtime_precision);

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
        return finish_test_cff_ltd_verdict(
            checks,
            direct_relative_tolerance,
            mass_shift_relative_tolerance,
        );
    };

    for candidate in direct {
        if candidate.label() == baseline.label() {
            continue;
        }
        let (abs_diff, rel_diff) = complex_distance(&candidate.value, &baseline.value);
        let within_tolerance = abs_diff.leq_f64(DIRECT_ABSOLUTE_TOLERANCE)
            || rel_diff.leq_f64(direct_relative_tolerance);
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
            TestCffLtdDistance::new(abs_diff, rel_diff, direct_relative_tolerance),
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
                || rel_diff.leq_f64(mass_shift_relative_tolerance);
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
                TestCffLtdDistance::new(abs_diff, rel_diff, mass_shift_relative_tolerance),
                if within_tolerance {
                    "At least one tested mass shift reproduces the CFF reference to the loose diagnostic threshold."
                } else {
                    "No tested mass shift reproduces the CFF reference to the loose diagnostic threshold."
                },
            ));
        }
    }

    finish_test_cff_ltd_verdict(
        checks,
        direct_relative_tolerance,
        mass_shift_relative_tolerance,
    )
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

fn finish_test_cff_ltd_verdict(
    checks: Vec<TestCffLtdCheckRecord>,
    direct_relative_tolerance: f64,
    mass_shift_relative_tolerance: f64,
) -> TestCffLtdVerdict {
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
        direct_relative_tolerance: format_f64_full(direct_relative_tolerance),
        direct_absolute_tolerance: format_f64_full(DIRECT_ABSOLUTE_TOLERANCE),
        mass_shift_relative_tolerance: format_f64_full(mass_shift_relative_tolerance),
        mass_shift_absolute_tolerance: format_f64_full(MASS_SHIFT_ABSOLUTE_TOLERANCE),
        checks,
    }
}

fn pure_ltd_mass_shift_epsilons(start: f64, n_steps: usize) -> Vec<f64> {
    (0..n_steps)
        .map(|step| start * 10.0_f64.powi(-(step as i32)))
        .collect()
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

fn diagnostic_parametric_atom_for_evaluator(
    expression: &ThreeDExpression<OrientationID>,
    graph: &Graph,
    _evaluator_settings: &EvaluatorSettings,
) -> Result<Atom> {
    let numerator = graph.full_numerator_atom();
    let processed_numerator = spenso_process_numerator(&numerator)?;
    let processed_numerator = expand_simple_minkowski_dots(processed_numerator);

    Ok(gammalooprs::cff::expression::GammaLoopThreeDExpression::diagnostic_parametric_atom_with_numerator_gs(
        expression,
        graph,
        &processed_numerator,
        &OrientationPattern::default(),
    ))
}

fn diagnostic_component_atoms_for_evaluator(
    expression: &ThreeDExpression<OrientationID>,
    graph: &Graph,
) -> Result<DiagnosticComponentAtoms> {
    let numerator = graph.full_numerator_atom();
    let processed_numerator = expand_simple_minkowski_dots(spenso_process_numerator(&numerator)?);
    let mut numerator_atoms = Vec::with_capacity(expression.orientations.len());
    let mut orientation_atoms = Vec::with_capacity(expression.orientations.len());

    for orientation in &expression.orientations {
        let numerator_atom = orientation.numerator_atom_gs(graph, &processed_numerator);
        let orientation_atom = expression
            .surfaces
            .substitute_energies_gs(&orientation.to_atom_gs(), &[]);
        numerator_atoms.push(numerator_atom);
        orientation_atoms.push(orientation_atom);
    }

    Ok(DiagnosticComponentAtoms {
        processed_numerator,
        numerator_atoms,
        orientation_atoms,
    })
}

fn numerator_only_component_atoms_for_evaluator(graph: &Graph) -> DiagnosticComponentAtoms {
    let numerator = numerator_only_parametric_atom_for_evaluator(graph);
    DiagnosticComponentAtoms {
        processed_numerator: numerator.clone(),
        numerator_atoms: vec![numerator],
        orientation_atoms: vec![Atom::num(1)],
    }
}

fn preprocess_numerator(numerator: &Atom) -> Atom {
    numerator.as_atom_view().simplify_color()
}

fn spenso_process_numerator(numerator: &Atom) -> Result<Atom> {
    let preprocessed_numerator = preprocess_numerator(numerator);
    let mut net = preprocessed_numerator
        .as_atom_view()
        .parse_into_net()
        .with_context(|| "Could not parse preprocessed numerator into Spenso network")?;

    net.execute::<Sequential, SmallestDegree, _, _, _>(
        TENSORLIB.read().unwrap().deref(),
        FUN_LIB.deref(),
    )
    .with_context(|| "Could not execute Spenso numerator network")?;

    net.result_scalar()
        .map(|result| match result {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        })
        .map_err(Into::into)
}

fn expand_simple_minkowski_dots(atom: Atom) -> Atom {
    let minkowski_tag = LibraryRep::from(Minkowski {}).to_symbolic([Atom::num(4)]);
    atom.replace_map(|view, _ctx, out| {
        let AtomView::Fun(function) = view else {
            return;
        };
        if function.get_symbol() != SPENSO_TAG.dot || function.get_nargs() != 3 {
            return;
        }

        let args = function.iter().collect::<Vec<_>>();
        let Some(metric_position) = args.iter().position(|arg| *arg == minkowski_tag.as_view())
        else {
            return;
        };
        let vectors = args
            .iter()
            .enumerate()
            .filter_map(|(position, arg)| (position != metric_position).then_some(*arg))
            .collect::<Vec<_>>();
        if vectors.len() != 2 {
            return;
        }

        if let Some(dot) = minkowski_dot_components(vectors[0], vectors[1]) {
            **out = dot;
        }
    })
}

fn minkowski_dot_components(left: AtomView<'_>, right: AtomView<'_>) -> Option<Atom> {
    let mut dot = vector_component(left, 0)? * vector_component(right, 0)?;
    for component in 1..=3 {
        dot -= vector_component(left, component)? * vector_component(right, component)?;
    }
    Some(dot)
}

fn vector_component(vector: AtomView<'_>, component: usize) -> Option<Atom> {
    let component = Atom::from(ExpandedIndex::from_iter([component]));
    match vector {
        AtomView::Fun(function) => {
            let mut args = function
                .iter()
                .map(|arg| arg.to_owned())
                .collect::<Vec<_>>();
            args.push(component);
            Some(
                FunctionBuilder::new(function.get_symbol())
                    .add_args(&args)
                    .finish(),
            )
        }
        AtomView::Var(symbol) => Some(
            FunctionBuilder::new(symbol.get_symbol())
                .add_arg(component.as_view())
                .finish(),
        ),
        _ => None,
    }
}

fn numerator_only_parametric_atom_for_evaluator(graph: &Graph) -> Atom {
    numerator_with_internal_energy_parameters_gs(graph)
        .simplify_metrics()
        .simplify_color()
}

fn numerator_only_orientation(graph: &Graph) -> TiVec<OrientationID, EdgeVec<Orientation>> {
    std::iter::once(graph.underlying.new_edgevec(|_, _, pair| {
        if matches!(pair, HedgePair::Paired { .. }) {
            Orientation::Default
        } else {
            Orientation::Undirected
        }
    }))
    .collect()
}

fn set_random_diagnostic_input_values(
    param_builder: &mut ParamBuilder<f64>,
    seed: u64,
    scale: f64,
) {
    set_random_diagnostic_input_values_with_policy(param_builder, seed, scale, true);
}

fn set_random_diagnostic_input_values_with_policy(
    param_builder: &mut ParamBuilder<f64>,
    seed: u64,
    scale: f64,
    randomize_externals: bool,
) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let parameters = param_builder.evaluator_input_parameters();
    for parameter in parameters {
        let should_randomize_external = randomize_externals
            && matches!(
                parameter.group,
                ParamBuilderInputGroup::ExternalEnergy | ParamBuilderInputGroup::ExternalSpatial
            );
        let should_randomize_non_external = matches!(
            parameter.group,
            ParamBuilderInputGroup::LoopMomentumSpatial | ParamBuilderInputGroup::Additional
        );
        if should_randomize_external || should_randomize_non_external {
            param_builder.values[0][parameter.index] =
                Complex::new_re(F(scale * rng.random::<f64>()));
        }
    }
}

#[derive(Clone, Copy)]
enum DiagnosticExternalSource {
    RuntimeKinematics,
    AutomaticDiagnostic,
}

impl DiagnosticExternalSource {
    fn label(self, group: ParamBuilderInputGroup) -> &'static str {
        match (self, group) {
            (
                Self::RuntimeKinematics,
                ParamBuilderInputGroup::ExternalEnergy
                | ParamBuilderInputGroup::ExternalSpatial
                | ParamBuilderInputGroup::Polarization,
            ) => "runtime kinematics",
            (
                Self::AutomaticDiagnostic,
                ParamBuilderInputGroup::ExternalEnergy | ParamBuilderInputGroup::ExternalSpatial,
            ) => "automatic diagnostic",
            (
                _,
                ParamBuilderInputGroup::LoopMomentumSpatial | ParamBuilderInputGroup::Additional,
            ) => "automatic diagnostic",
            (
                _,
                ParamBuilderInputGroup::Runtime
                | ParamBuilderInputGroup::Model
                | ParamBuilderInputGroup::Polarization
                | ParamBuilderInputGroup::LocalCounterterm,
            ) => "state/default",
        }
    }
}

struct PreparedDiagnosticParamBuilder {
    param_builder: ParamBuilder<f64>,
    external_source: DiagnosticExternalSource,
    user_overridden_parameters: BTreeMap<String, String>,
}

#[derive(Clone, Default)]
struct DiagnosticInputConfig {
    numerator_only: bool,
    numerator_q0_overrides: BTreeMap<EdgeIndex, f64>,
}

fn prepare_diagnostic_param_builder(
    graph: &Graph,
    model: &Model,
    integrand_kind: SelectedIntegrandKind,
    run_settings: &ThreeDrepRunSettings,
    default_runtime_settings: &RuntimeSettings,
    input_config: &DiagnosticInputConfig,
) -> Result<PreparedDiagnosticParamBuilder> {
    let mut param_builder = graph.param_builder.clone();
    let numerator_q0_parameters =
        ensure_numerator_q0_parameters(graph, &mut param_builder, input_config)?;
    let runtime_external_result = set_runtime_external_kinematics(
        &mut param_builder,
        graph,
        integrand_kind,
        default_runtime_settings,
    );
    let external_source = match runtime_external_result {
        Ok(()) => {
            set_random_diagnostic_input_values_with_policy(
                &mut param_builder,
                run_settings.seed,
                run_settings.scale,
                false,
            );
            DiagnosticExternalSource::RuntimeKinematics
        }
        Err(_error) if graph.polarizations.is_empty() => {
            set_random_diagnostic_input_values(
                &mut param_builder,
                run_settings.seed,
                run_settings.scale,
            );
            DiagnosticExternalSource::AutomaticDiagnostic
        }
        Err(error) => {
            return Err(error).with_context(|| {
                format!(
                    "while preparing external kinematics and polarizations for 3Drep graph {}",
                    graph.name
                )
            });
        }
    };

    param_builder.set_runtime_parameter_values(
        Complex::new_re(F(default_runtime_settings.general.m_uv)),
        Complex::new_re(F(
            default_runtime_settings.general.mu_r * default_runtime_settings.general.mu_r
        )),
        Complex::new_re(F(run_settings.numerator_interpolation_scale)),
    );
    param_builder.update_model_values(model);
    let user_overridden_parameters = apply_numerator_q0_overrides(
        &mut param_builder,
        &numerator_q0_parameters,
        &input_config.numerator_q0_overrides,
    )?;

    Ok(PreparedDiagnosticParamBuilder {
        param_builder,
        external_source,
        user_overridden_parameters,
    })
}

fn ensure_numerator_q0_parameters(
    graph: &Graph,
    param_builder: &mut ParamBuilder<f64>,
    input_config: &DiagnosticInputConfig,
) -> Result<BTreeMap<EdgeIndex, ParamBuilderInputParameter>> {
    if !input_config.numerator_only {
        if !input_config.numerator_q0_overrides.is_empty() {
            return Err(eyre!(
                "--numerator-q0 can only be used together with 3Drep evaluate --numerator-only"
            ));
        }
        return Ok(BTreeMap::new());
    }

    let energy_atoms = graph
        .underlying
        .iter_edges()
        .filter_map(|(pair, edge_id, _)| {
            (matches!(pair, HedgePair::Paired { .. })
                && graph.loop_momentum_basis.edge_signatures[edge_id].is_loop_dependent())
            .then(|| (edge_id, internal_energy_parameter_atom_gs(edge_id)))
        })
        .collect::<Vec<_>>();
    let requested_atoms = energy_atoms
        .iter()
        .map(|(_, atom)| atom.clone())
        .collect::<Vec<_>>();
    let ensured = param_builder.ensure_additional_input_parameters(requested_atoms);
    let by_canonical = ensured
        .into_iter()
        .map(|parameter| (parameter.atom.to_canonical_string(), parameter))
        .collect::<BTreeMap<_, _>>();

    let mut by_edge = BTreeMap::new();
    for (edge_id, atom) in energy_atoms {
        let canonical = atom.to_canonical_string();
        let parameter = by_canonical.get(&canonical).cloned().ok_or_else(|| {
            eyre!(
                "Could not register numerator-only energy parameter {}",
                atom.log_print(Some(60))
            )
        })?;
        by_edge.insert(edge_id, parameter);
    }

    for edge_id in input_config.numerator_q0_overrides.keys() {
        if !by_edge.contains_key(edge_id) {
            return Err(eyre!(
                "Cannot override numerator-only energy for edge {} because it is not an internal propagator edge",
                edge_id.0
            ));
        }
    }

    Ok(by_edge)
}

fn apply_numerator_q0_overrides(
    param_builder: &mut ParamBuilder<f64>,
    parameters: &BTreeMap<EdgeIndex, ParamBuilderInputParameter>,
    overrides: &BTreeMap<EdgeIndex, f64>,
) -> Result<BTreeMap<String, String>> {
    let mut user_overridden_parameters = BTreeMap::new();
    for (edge_id, value) in overrides {
        let parameter = parameters.get(edge_id).ok_or_else(|| {
            eyre!(
                "Cannot override numerator-only energy for edge {} because it was not registered",
                edge_id.0
            )
        })?;
        param_builder.values[0][parameter.index] = Complex::new_re(F(*value));
        user_overridden_parameters.insert(
            parameter.atom.to_canonical_string(),
            "user override".to_string(),
        );
    }
    Ok(user_overridden_parameters)
}

fn set_runtime_external_kinematics(
    param_builder: &mut ParamBuilder<f64>,
    graph: &Graph,
    integrand_kind: SelectedIntegrandKind,
    default_runtime_settings: &RuntimeSettings,
) -> Result<()> {
    match integrand_kind {
        SelectedIntegrandKind::Amplitude => {
            let external_signature = graph.get_external_signature();
            param_builder.set_external_kinematics_and_polarizations(
                graph,
                default_runtime_settings,
                DependentMomentaConstructor::Amplitude(&external_signature),
            )
        }
        SelectedIntegrandKind::CrossSection => param_builder
            .set_external_kinematics_and_polarizations(
                graph,
                default_runtime_settings,
                DependentMomentaConstructor::CrossSection,
            ),
    }
}

fn parameter_records(
    graph: &Graph,
    model: &Model,
    integrand_kind: SelectedIntegrandKind,
    run_settings: &ThreeDrepRunSettings,
    default_runtime_settings: &RuntimeSettings,
    input_config: &DiagnosticInputConfig,
) -> Result<Vec<ParameterValueRecord>> {
    let prepared = prepare_diagnostic_param_builder(
        graph,
        model,
        integrand_kind,
        run_settings,
        default_runtime_settings,
        input_config,
    )?;

    Ok(prepared
        .param_builder
        .evaluator_input_parameters()
        .into_iter()
        .map(|parameter| {
            let value = &prepared.param_builder.values[0][parameter.index];
            let (formatted_value, _, _) = format_complex_full(value);
            ParameterValueRecord {
                name: parameter.atom.log_print(Some(60)),
                canonical_name: parameter.atom.to_canonical_string(),
                value: formatted_value,
                source: prepared
                    .user_overridden_parameters
                    .get(&parameter.atom.to_canonical_string())
                    .cloned()
                    .unwrap_or_else(|| prepared.external_source.label(parameter.group).to_string()),
            }
        })
        .collect())
}

struct EvaluationRequest<'a> {
    label: &'a str,
    graph: &'a Graph,
    model: &'a Model,
    integrand_kind: SelectedIntegrandKind,
    expression: Option<&'a ThreeDExpression<OrientationID>>,
    parametric_atom: &'a symbolica::atom::Atom,
    atom_build_timing: Duration,
    orientations: Option<&'a TiVec<OrientationID, EdgeVec<Orientation>>>,
    workspace: &'a Path,
    global_cli_settings: &'a CLISettings,
    default_runtime_settings: &'a RuntimeSettings,
    representation: Option<RepresentationMode>,
    numerator_sampling_scale_mode: Option<three_dimensional_reps::NumeratorSamplingScaleMode>,
    run_settings: &'a ThreeDrepRunSettings,
    mass_shift: &'a str,
    mass_shift_values: &'a [MassShiftValueRecord],
    profile_target: Option<Duration>,
    input_config: &'a DiagnosticInputConfig,
    clean: bool,
}

struct ProfileEvaluationRequest<'a> {
    param_builder: &'a mut ParamBuilder<f64>,
    precision: CliRuntimePrecision,
    orientations: &'a TiVec<OrientationID, EdgeVec<Orientation>>,
    filter: &'a SubSet<OrientationID>,
    runtime_settings: &'a RuntimeSettings,
    evaluation_metadata: &'a mut EvaluationMetaData,
    target: Duration,
}

struct ComponentPrecisionEvaluationRequest<'a> {
    param_builder: &'a mut ParamBuilder<f64>,
    precision: CliRuntimePrecision,
    orientations: &'a TiVec<OrientationID, EdgeVec<Orientation>>,
    filter: &'a SubSet<OrientationID>,
    runtime_settings: &'a RuntimeSettings,
    evaluation_metadata: &'a mut EvaluationMetaData,
}

fn evaluate_threedrep_expression(
    request: EvaluationRequest<'_>,
) -> Result<ThreeDrepEvaluationRecord> {
    let owned_orientations;
    let orientations = if let Some(orientations) = request.orientations {
        orientations
    } else {
        owned_orientations =
            diagnostic_evaluation_orientations(request.expression.ok_or_else(|| {
                eyre!("3Drep evaluation requires an oriented expression unless orientations are supplied")
            })?);
        &owned_orientations
    };
    let prepared_param_builder = prepare_diagnostic_param_builder(
        request.graph,
        request.model,
        request.integrand_kind,
        request.run_settings,
        request.default_runtime_settings,
        request.input_config,
    )?;
    let prepared = build_threedrep_evaluator(
        request.label,
        EvaluatorBuildContext {
            workspace: request.workspace,
            model: request.model,
            global_cli_settings: request.global_cli_settings,
            default_runtime_settings: request.default_runtime_settings,
            run_settings: request.run_settings,
            clean: request.clean,
        },
        &prepared_param_builder.param_builder,
        &[request.parametric_atom],
        orientations,
    )
    .map_err(|error| eyre!("{}", format_diagnostic_error(&error)))?;
    let mut evaluator_info = prepared.info.clone();
    evaluator_info.include_atom_build_timing(request.atom_build_timing);
    let mut evaluator = prepared.evaluator;
    evaluate_with_evaluator(
        &mut evaluator,
        request,
        orientations,
        &evaluator_info,
        prepared_param_builder,
    )
}

struct ComponentEvaluationRequest<'a> {
    label: &'a str,
    graph: &'a Graph,
    model: &'a Model,
    integrand_kind: SelectedIntegrandKind,
    components: &'a DiagnosticComponentAtoms,
    atom_build_timing: Duration,
    orientations: &'a TiVec<OrientationID, EdgeVec<Orientation>>,
    workspace: &'a Path,
    global_cli_settings: &'a CLISettings,
    default_runtime_settings: &'a RuntimeSettings,
    representation: Option<RepresentationMode>,
    numerator_sampling_scale_mode: Option<three_dimensional_reps::NumeratorSamplingScaleMode>,
    run_settings: &'a ThreeDrepRunSettings,
    mass_shift: &'a str,
    mass_shift_values: &'a [MassShiftValueRecord],
    profile_target: Option<Duration>,
    input_config: &'a DiagnosticInputConfig,
    clean: bool,
}

fn evaluate_threedrep_component_expression(
    request: ComponentEvaluationRequest<'_>,
) -> Result<ThreeDrepEvaluationRecord> {
    let prepared_param_builder = prepare_diagnostic_param_builder(
        request.graph,
        request.model,
        request.integrand_kind,
        request.run_settings,
        request.default_runtime_settings,
        request.input_config,
    )?;
    let mut prepared = build_threedrep_component_evaluator(
        request.label,
        EvaluatorBuildContext {
            workspace: request.workspace,
            model: request.model,
            global_cli_settings: request.global_cli_settings,
            default_runtime_settings: request.default_runtime_settings,
            run_settings: request.run_settings,
            clean: request.clean,
        },
        &prepared_param_builder.param_builder,
        request.components,
    )
    .map_err(|error| eyre!("{}", format_diagnostic_error(&error)))?;
    prepared
        .info
        .include_atom_build_timing(request.atom_build_timing);
    evaluate_with_component_evaluator(prepared, request, prepared_param_builder)
}

fn evaluate_with_evaluator(
    evaluator: &mut EvaluatorStack,
    request: EvaluationRequest<'_>,
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    evaluator_info: &PreparedEvaluatorInfo,
    prepared_param_builder: PreparedDiagnosticParamBuilder,
) -> Result<ThreeDrepEvaluationRecord> {
    let mut param_builder = prepared_param_builder.param_builder;

    let filter = SubSet::<OrientationID>::full(orientations.len());
    let mut evaluation_metadata = EvaluationMetaData::new_empty();
    let runtime_settings =
        threedrep_runtime_settings(request.default_runtime_settings, request.run_settings);
    if let Some(profile_target) = request.profile_target {
        let profile = profile_evaluator_for_precision(
            evaluator,
            ProfileEvaluationRequest {
                param_builder: &mut param_builder,
                precision: request.run_settings.runtime_precision,
                orientations,
                filter: &filter,
                runtime_settings: &runtime_settings,
                evaluation_metadata: &mut evaluation_metadata,
                target: profile_target,
            },
        );
        return match profile {
            Ok(profile) => Ok(ThreeDrepEvaluationRecord {
                id: 0,
                label: request.label.to_string(),
                representation: request.representation,
                numerator_sampling_scale_mode: request.numerator_sampling_scale_mode,
                mass_shift: request.mass_shift.to_string(),
                mass_shift_values: request.mass_shift_values.to_vec(),
                numerator_interpolation_scale: request.run_settings.numerator_interpolation_scale,
                runtime_precision: request.run_settings.runtime_precision,
                build_strategy: request.run_settings.build_strategy,
                evaluator_method: request.run_settings.evaluator_method.clone(),
                evaluator_backend: evaluator_info.backend.clone(),
                evaluator_build_timing: evaluator_info.build_timing.clone(),
                evaluator_build_timing_seconds: evaluator_info.build_timing_seconds,
                value: "-".to_string(),
                value_re: "-".to_string(),
                value_im: "-".to_string(),
                sample_evaluation_timing: profile.timing_per_sample.clone(),
                sample_evaluation_timing_seconds: profile.timing_per_sample_seconds,
                profile: Some(profile),
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
                build_strategy: request.run_settings.build_strategy,
                evaluator_method: request.run_settings.evaluator_method.clone(),
                evaluator_backend: evaluator_info.backend.clone(),
                evaluator_build_timing: evaluator_info.build_timing.clone(),
                evaluator_build_timing_seconds: evaluator_info.build_timing_seconds,
                value: "-".to_string(),
                value_re: "-".to_string(),
                value_im: "-".to_string(),
                sample_evaluation_timing: "-".to_string(),
                sample_evaluation_timing_seconds: 0.0,
                profile: None,
                status: "failed".to_string(),
                error: Some(format_diagnostic_error(&error)),
            }),
        };
    }
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
            build_strategy: request.run_settings.build_strategy,
            evaluator_method: request.run_settings.evaluator_method.clone(),
            evaluator_backend: evaluator_info.backend.clone(),
            evaluator_build_timing: evaluator_info.build_timing.clone(),
            evaluator_build_timing_seconds: evaluator_info.build_timing_seconds,
            value: formatted_value,
            value_re,
            value_im,
            sample_evaluation_timing: format_duration_dynamic(timing),
            sample_evaluation_timing_seconds: timing.as_secs_f64(),
            profile: None,
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
            build_strategy: request.run_settings.build_strategy,
            evaluator_method: request.run_settings.evaluator_method.clone(),
            evaluator_backend: evaluator_info.backend.clone(),
            evaluator_build_timing: evaluator_info.build_timing.clone(),
            evaluator_build_timing_seconds: evaluator_info.build_timing_seconds,
            value: "-".to_string(),
            value_re: "-".to_string(),
            value_im: "-".to_string(),
            sample_evaluation_timing: format_duration_dynamic(timing),
            sample_evaluation_timing_seconds: timing.as_secs_f64(),
            profile: None,
            status: "failed".to_string(),
            error: Some(format_diagnostic_error(&error)),
        }),
    }
}

fn evaluate_with_component_evaluator(
    mut prepared: PreparedComponentEvaluator,
    request: ComponentEvaluationRequest<'_>,
    prepared_param_builder: PreparedDiagnosticParamBuilder,
) -> Result<ThreeDrepEvaluationRecord> {
    let mut param_builder = prepared_param_builder.param_builder;
    let filter = SubSet::<OrientationID>::full(request.orientations.len());
    let mut evaluation_metadata = EvaluationMetaData::new_empty();
    let runtime_settings =
        threedrep_runtime_settings(request.default_runtime_settings, request.run_settings);

    if let Some(profile_target) = request.profile_target {
        let profile = profile_component_evaluator_for_precision(
            &mut prepared.numerator_evaluator,
            &mut prepared.orientation_evaluator,
            ProfileEvaluationRequest {
                param_builder: &mut param_builder,
                precision: request.run_settings.runtime_precision,
                orientations: request.orientations,
                filter: &filter,
                runtime_settings: &runtime_settings,
                evaluation_metadata: &mut evaluation_metadata,
                target: profile_target,
            },
        );
        return match profile {
            Ok(profile) => Ok(component_evaluation_record(
                &request,
                &prepared.info,
                "-".to_string(),
                "-".to_string(),
                "-".to_string(),
                profile.timing_per_sample.clone(),
                profile.timing_per_sample_seconds,
                Some(profile),
                "ok",
                None,
            )),
            Err(error) => Ok(component_evaluation_record(
                &request,
                &prepared.info,
                "-".to_string(),
                "-".to_string(),
                "-".to_string(),
                "-".to_string(),
                0.0,
                None,
                "failed",
                Some(format_diagnostic_error(&error)),
            )),
        };
    }

    let start = Instant::now();
    let result = evaluate_component_for_precision(
        &mut prepared.numerator_evaluator,
        &mut prepared.orientation_evaluator,
        ComponentPrecisionEvaluationRequest {
            param_builder: &mut param_builder,
            precision: request.run_settings.runtime_precision,
            orientations: request.orientations,
            filter: &filter,
            runtime_settings: &runtime_settings,
            evaluation_metadata: &mut evaluation_metadata,
        },
    );
    let timing = start.elapsed();

    match result {
        Ok((formatted_value, value_re, value_im)) => Ok(component_evaluation_record(
            &request,
            &prepared.info,
            formatted_value,
            value_re,
            value_im,
            format_duration_dynamic(timing),
            timing.as_secs_f64(),
            None,
            "ok",
            None,
        )),
        Err(error) => Ok(component_evaluation_record(
            &request,
            &prepared.info,
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            format_duration_dynamic(timing),
            timing.as_secs_f64(),
            None,
            "failed",
            Some(format_diagnostic_error(&error)),
        )),
    }
}

#[allow(clippy::too_many_arguments)]
fn component_evaluation_record(
    request: &ComponentEvaluationRequest<'_>,
    evaluator_info: &PreparedEvaluatorInfo,
    value: String,
    value_re: String,
    value_im: String,
    sample_evaluation_timing: String,
    sample_evaluation_timing_seconds: f64,
    profile: Option<ThreeDrepProfileRecord>,
    status: &str,
    error: Option<String>,
) -> ThreeDrepEvaluationRecord {
    ThreeDrepEvaluationRecord {
        id: 0,
        label: request.label.to_string(),
        representation: request.representation,
        numerator_sampling_scale_mode: request.numerator_sampling_scale_mode,
        mass_shift: request.mass_shift.to_string(),
        mass_shift_values: request.mass_shift_values.to_vec(),
        numerator_interpolation_scale: request.run_settings.numerator_interpolation_scale,
        runtime_precision: request.run_settings.runtime_precision,
        build_strategy: request.run_settings.build_strategy,
        evaluator_method: request.run_settings.evaluator_method.clone(),
        evaluator_backend: evaluator_info.backend.clone(),
        evaluator_build_timing: evaluator_info.build_timing.clone(),
        evaluator_build_timing_seconds: evaluator_info.build_timing_seconds,
        value,
        value_re,
        value_im,
        sample_evaluation_timing,
        sample_evaluation_timing_seconds,
        profile,
        status: status.to_string(),
        error,
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
    let value = evaluate_value_for_precision_impl::<T>(
        evaluator,
        input,
        orientations,
        runtime_settings,
        evaluation_metadata,
    )?;
    Ok(format_complex_full_precise(&value))
}

fn evaluate_value_for_precision_impl<T: FloatLike>(
    evaluator: &mut EvaluatorStack,
    input: InputParams<'_, T>,
    orientations: SingleOrAllOrientations<'_, OrientationID>,
    runtime_settings: &RuntimeSettings,
    evaluation_metadata: &mut EvaluationMetaData,
) -> Result<Complex<F<T>>> {
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
    Ok(value)
}

fn evaluate_for_precision_discard(
    evaluator: &mut EvaluatorStack,
    param_builder: &mut ParamBuilder<f64>,
    precision: CliRuntimePrecision,
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    filter: &SubSet<OrientationID>,
    runtime_settings: &RuntimeSettings,
    evaluation_metadata: &mut EvaluationMetaData,
) -> Result<()> {
    match precision {
        CliRuntimePrecision::Double => {
            let input = param_builder.input_params();
            evaluate_value_for_precision_impl::<f64>(
                evaluator,
                input,
                SingleOrAllOrientations::All {
                    all: orientations,
                    filter,
                },
                runtime_settings,
                evaluation_metadata,
            )?;
        }
        CliRuntimePrecision::Quad => {
            let input = param_builder.input_params_quad();
            evaluate_value_for_precision_impl::<f128>(
                evaluator,
                input,
                SingleOrAllOrientations::All {
                    all: orientations,
                    filter,
                },
                runtime_settings,
                evaluation_metadata,
            )?;
        }
        CliRuntimePrecision::ArbPrec => {
            let input = param_builder.input_params_arb();
            evaluate_value_for_precision_impl::<ArbPrec>(
                evaluator,
                input,
                SingleOrAllOrientations::All {
                    all: orientations,
                    filter,
                },
                runtime_settings,
                evaluation_metadata,
            )?;
        }
    }
    Ok(())
}

fn evaluate_component_for_precision(
    numerator_evaluator: &mut EvaluatorStack,
    orientation_evaluator: &mut EvaluatorStack,
    request: ComponentPrecisionEvaluationRequest<'_>,
) -> Result<(String, String, String)> {
    let ComponentPrecisionEvaluationRequest {
        param_builder,
        precision,
        orientations,
        filter,
        runtime_settings,
        evaluation_metadata,
    } = request;
    match precision {
        CliRuntimePrecision::Double => {
            let numerator_values = {
                let input = param_builder.input_params();
                evaluate_component_vector_for_precision_impl::<f64>(
                    numerator_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    true,
                )?
            };
            let orientation_values = {
                let input = param_builder.input_params();
                evaluate_component_vector_for_precision_impl::<f64>(
                    orientation_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    false,
                )?
            };
            let value = dot_component_vectors(numerator_values, orientation_values)?;
            Ok(format_complex_full_precise(&value))
        }
        CliRuntimePrecision::Quad => {
            let numerator_values = {
                let input = param_builder.input_params_quad();
                evaluate_component_vector_for_precision_impl::<f128>(
                    numerator_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    true,
                )?
            };
            let orientation_values = {
                let input = param_builder.input_params_quad();
                evaluate_component_vector_for_precision_impl::<f128>(
                    orientation_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    false,
                )?
            };
            let value = dot_component_vectors(numerator_values, orientation_values)?;
            Ok(format_complex_full_precise(&value))
        }
        CliRuntimePrecision::ArbPrec => {
            let numerator_values = {
                let input = param_builder.input_params_arb();
                evaluate_component_vector_for_precision_impl::<ArbPrec>(
                    numerator_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    true,
                )?
            };
            let orientation_values = {
                let input = param_builder.input_params_arb();
                evaluate_component_vector_for_precision_impl::<ArbPrec>(
                    orientation_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    false,
                )?
            };
            let value = dot_component_vectors(numerator_values, orientation_values)?;
            Ok(format_complex_full_precise(&value))
        }
    }
}

fn dot_component_vectors<T: FloatLike>(
    numerator_values: Vec<Complex<F<T>>>,
    orientation_values: Vec<Complex<F<T>>>,
) -> Result<Complex<F<T>>> {
    if numerator_values.len() != orientation_values.len() {
        return Err(eyre!(
            "Iterative 3Drep evaluator returned {} numerator components and {} orientation components",
            numerator_values.len(),
            orientation_values.len()
        ));
    }
    let mut components = numerator_values.into_iter().zip(orientation_values);
    let (first_numerator, first_orientation) = components
        .next()
        .ok_or_else(|| eyre!("Iterative 3Drep evaluator returned no components"))?;
    let mut value = first_numerator * first_orientation;
    for (numerator, orientation) in components {
        value += numerator * orientation;
    }
    Ok(value)
}

fn evaluate_component_vector_for_precision_impl<T: FloatLike>(
    evaluator: &mut EvaluatorStack,
    input: InputParams<'_, T>,
    orientations: SingleOrAllOrientations<'_, OrientationID>,
    runtime_settings: &RuntimeSettings,
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
) -> Result<Vec<Complex<F<T>>>> {
    evaluator
        .evaluate::<T, OrientationID>(
            input,
            orientations,
            runtime_settings,
            evaluation_metadata,
            record_primary_timing,
        )?
        .into_iter()
        .map(|value| Ok(value.unwrap_real()))
        .collect()
}

fn evaluate_component_dot_for_precision_discard(
    numerator_evaluator: &mut EvaluatorStack,
    orientation_evaluator: &mut EvaluatorStack,
    request: ComponentPrecisionEvaluationRequest<'_>,
) -> Result<()> {
    let ComponentPrecisionEvaluationRequest {
        param_builder,
        precision,
        orientations,
        filter,
        runtime_settings,
        evaluation_metadata,
    } = request;
    match precision {
        CliRuntimePrecision::Double => {
            let numerator_values = {
                let input = param_builder.input_params();
                evaluate_component_vector_for_precision_impl::<f64>(
                    numerator_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    true,
                )?
            };
            let orientation_values = {
                let input = param_builder.input_params();
                evaluate_component_vector_for_precision_impl::<f64>(
                    orientation_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    false,
                )?
            };
            dot_component_vectors(numerator_values, orientation_values)?;
        }
        CliRuntimePrecision::Quad => {
            let numerator_values = {
                let input = param_builder.input_params_quad();
                evaluate_component_vector_for_precision_impl::<f128>(
                    numerator_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    true,
                )?
            };
            let orientation_values = {
                let input = param_builder.input_params_quad();
                evaluate_component_vector_for_precision_impl::<f128>(
                    orientation_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    false,
                )?
            };
            dot_component_vectors(numerator_values, orientation_values)?;
        }
        CliRuntimePrecision::ArbPrec => {
            let numerator_values = {
                let input = param_builder.input_params_arb();
                evaluate_component_vector_for_precision_impl::<ArbPrec>(
                    numerator_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    true,
                )?
            };
            let orientation_values = {
                let input = param_builder.input_params_arb();
                evaluate_component_vector_for_precision_impl::<ArbPrec>(
                    orientation_evaluator,
                    input,
                    SingleOrAllOrientations::All {
                        all: orientations,
                        filter,
                    },
                    runtime_settings,
                    evaluation_metadata,
                    false,
                )?
            };
            dot_component_vectors(numerator_values, orientation_values)?;
        }
    }
    Ok(())
}

fn profile_component_evaluator_for_precision(
    numerator_evaluator: &mut EvaluatorStack,
    orientation_evaluator: &mut EvaluatorStack,
    request: ProfileEvaluationRequest<'_>,
) -> Result<ThreeDrepProfileRecord> {
    const WARMUP_CALLS: usize = 10;
    let warmup_start = Instant::now();
    for _ in 0..WARMUP_CALLS {
        evaluate_component_dot_for_precision_discard(
            numerator_evaluator,
            orientation_evaluator,
            ComponentPrecisionEvaluationRequest {
                param_builder: &mut *request.param_builder,
                precision: request.precision,
                orientations: request.orientations,
                filter: request.filter,
                runtime_settings: request.runtime_settings,
                evaluation_metadata: &mut *request.evaluation_metadata,
            },
        )?;
    }
    let warmup_timing = warmup_start.elapsed();
    let warmup_seconds = warmup_timing.as_secs_f64();
    let target_seconds = request.target.as_secs_f64();
    let calls = if warmup_seconds > 0.0 {
        ((WARMUP_CALLS as f64 * target_seconds / warmup_seconds).ceil() as usize).max(1)
    } else {
        WARMUP_CALLS
    };

    let profile_start = Instant::now();
    for _ in 0..calls {
        evaluate_component_dot_for_precision_discard(
            numerator_evaluator,
            orientation_evaluator,
            ComponentPrecisionEvaluationRequest {
                param_builder: &mut *request.param_builder,
                precision: request.precision,
                orientations: request.orientations,
                filter: request.filter,
                runtime_settings: request.runtime_settings,
                evaluation_metadata: &mut *request.evaluation_metadata,
            },
        )?;
    }
    let total_timing = profile_start.elapsed();
    let timing_per_sample_seconds = total_timing.as_secs_f64() / calls as f64;
    Ok(ThreeDrepProfileRecord {
        target_timing: format_duration_dynamic(request.target),
        target_timing_seconds: target_seconds,
        warmup_calls: WARMUP_CALLS,
        warmup_timing: format_duration_dynamic(warmup_timing),
        warmup_timing_seconds: warmup_seconds,
        calls,
        total_timing: format_duration_dynamic(total_timing),
        total_timing_seconds: total_timing.as_secs_f64(),
        timing_per_sample: format_duration_dynamic(Duration::from_secs_f64(
            timing_per_sample_seconds,
        )),
        timing_per_sample_seconds,
    })
}

fn profile_evaluator_for_precision(
    evaluator: &mut EvaluatorStack,
    request: ProfileEvaluationRequest<'_>,
) -> Result<ThreeDrepProfileRecord> {
    const WARMUP_CALLS: usize = 10;
    let warmup_start = Instant::now();
    for _ in 0..WARMUP_CALLS {
        evaluate_for_precision_discard(
            evaluator,
            &mut *request.param_builder,
            request.precision,
            request.orientations,
            request.filter,
            request.runtime_settings,
            &mut *request.evaluation_metadata,
        )?;
    }
    let warmup_timing = warmup_start.elapsed();
    let warmup_seconds = warmup_timing.as_secs_f64();
    let target_seconds = request.target.as_secs_f64();
    let calls = if warmup_seconds > 0.0 {
        ((WARMUP_CALLS as f64 * target_seconds / warmup_seconds).ceil() as usize).max(1)
    } else {
        WARMUP_CALLS
    };

    let profile_start = Instant::now();
    for _ in 0..calls {
        evaluate_for_precision_discard(
            evaluator,
            &mut *request.param_builder,
            request.precision,
            request.orientations,
            request.filter,
            request.runtime_settings,
            &mut *request.evaluation_metadata,
        )?;
    }
    let total_timing = profile_start.elapsed();
    let timing_per_sample_seconds = total_timing.as_secs_f64() / calls as f64;
    Ok(ThreeDrepProfileRecord {
        target_timing: format_duration_dynamic(request.target),
        target_timing_seconds: target_seconds,
        warmup_calls: WARMUP_CALLS,
        warmup_timing: format_duration_dynamic(warmup_timing),
        warmup_timing_seconds: warmup_seconds,
        calls,
        total_timing: format_duration_dynamic(total_timing),
        total_timing_seconds: total_timing.as_secs_f64(),
        timing_per_sample: format_duration_dynamic(Duration::from_secs_f64(
            timing_per_sample_seconds,
        )),
        timing_per_sample_seconds,
    })
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
        integrand_kind: catalog.integrand_kind(),
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
        integrand_kind: catalog.integrand_kind(),
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

fn atom_string(atom: &Atom) -> String {
    atom.to_canonical_string()
}

fn symbolica_raw_archive_cache_is_current(path: &Path) -> bool {
    let Ok(text) = fs::read_to_string(path) else {
        return false;
    };
    let Ok(value) = serde_json::from_str::<serde_json::Value>(&text) else {
        return false;
    };
    let has_current_schema = value
        .get("schema_version")
        .and_then(serde_json::Value::as_u64)
        .is_some_and(|version| version == u64::from(SYMBOLICA_RAW_ARCHIVE_SCHEMA_VERSION));
    let has_calls = value
        .get("calls")
        .and_then(serde_json::Value::as_array)
        .is_some_and(|calls| !calls.is_empty());
    has_current_schema && has_calls
}

fn write_iterative_raw_stub(
    raw_path: &Path,
    script_path: &Path,
    component_count: usize,
) -> Result<()> {
    #[derive(Serialize)]
    struct IterativeRawStub<'a> {
        schema_version: u32,
        build_strategy: &'a str,
        standalone_rust_replay_supported: bool,
        component_count: usize,
        message: &'a str,
    }

    let stub = IterativeRawStub {
        schema_version: SYMBOLICA_RAW_ARCHIVE_SCHEMA_VERSION,
        build_strategy: "iterative",
        standalone_rust_replay_supported: false,
        component_count,
        message: "3Drep --iterative builds separate numerator/orientation component evaluators; standalone raw replay is not supported for this strategy.",
    };
    write_path(raw_path, &serde_json::to_string_pretty(&stub)?)?;
    write_path(
        script_path,
        "fn main() {\n    panic!(\"3Drep --iterative raw standalone replay is not supported\");\n}\n",
    )?;
    make_script_executable(script_path)?;
    Ok(())
}

fn function_map_entry_record(
    entry: &gammalooprs::integrands::process::param_builder::FnMapEntry,
) -> SymbolicaFunctionMapEntryRecord {
    SymbolicaFunctionMapEntryRecord {
        lhs: atom_string(&entry.lhs),
        rhs: atom_string(&entry.rhs),
        tags: entry.tags.iter().map(atom_string).collect(),
        args: entry
            .args
            .iter()
            .map(|arg| atom_string(&Atom::from(arg.clone())))
            .collect(),
    }
}

fn complex_value_record(value: &Complex<F<f64>>) -> SymbolicaComplexValueRecord {
    let (_, re, im) = format_complex_full(value);
    SymbolicaComplexValueRecord { re, im }
}

fn evaluator_settings_record(settings: &EvaluatorSettings) -> SymbolicaEvaluatorSettingsRecord {
    SymbolicaEvaluatorSettingsRecord {
        do_algebra: settings.do_algebra,
        iterative_orientation_optimization: settings.iterative_orientation_optimization,
        summed: settings.summed,
        summed_function_map: settings.summed_function_map,
        compile: settings.compile,
        store_atom: settings.store_atom,
        do_fn_map_replacements: settings.do_fn_map_replacements,
        horner_iterations: settings.horner_iterations,
        n_cores: settings.n_cores,
        cpe_iterations: settings.cpe_iterations,
        abort_level: settings.abort_level,
        max_horner_scheme_variables: settings.max_horner_scheme_variables,
        max_common_pair_cache_entries: settings.max_common_pair_cache_entries,
        max_common_pair_distance: settings.max_common_pair_distance,
        verbose: settings.verbose,
    }
}

fn optimization_settings_record(
    settings: &EvaluatorSettings,
) -> SymbolicaOptimizationSettingsRecord {
    let optimization_settings = settings.optimization_settings();
    SymbolicaOptimizationSettingsRecord {
        horner_iterations: optimization_settings.horner_iterations,
        n_cores: optimization_settings.n_cores,
        cpe_iterations: optimization_settings.cpe_iterations,
        abort_level: optimization_settings.abort_level,
        max_horner_scheme_variables: optimization_settings.max_horner_scheme_variables,
        max_common_pair_cache_entries: optimization_settings.max_common_pair_cache_entries,
        max_common_pair_distance: optimization_settings.max_common_pair_distance,
        verbose: optimization_settings.verbose,
    }
}

fn param_builder_params(param_builder: &ParamBuilder<f64>) -> Vec<Atom> {
    (&param_builder.pairs)
        .into_iter()
        .flat_map(|pair| pair.params.clone())
        .collect()
}

fn diagnostic_evaluation_orientations_for_raw_input(
    expression: Option<&ThreeDExpression<OrientationID>>,
    explicit_orientations: Option<&TiVec<OrientationID, EdgeVec<Orientation>>>,
) -> Result<TiVec<OrientationID, EdgeVec<Orientation>>> {
    if let Some(orientations) = explicit_orientations {
        Ok(orientations.clone())
    } else {
        Ok(diagnostic_evaluation_orientations(expression.ok_or_else(
            || eyre!("3Drep raw evaluator input requires an oriented expression"),
        )?))
    }
}

fn symbolica_evaluator_input_archive(
    processed_atoms: &[Atom],
    param_builder: &ParamBuilder<f64>,
    orientations: TiVec<OrientationID, EdgeVec<Orientation>>,
    evaluator_settings: &EvaluatorSettings,
    evaluator_method: &EvaluatorMethod,
    evaluator_backend: &str,
) -> Result<SymbolicaEvaluatorInputArchive> {
    let params = param_builder_params(param_builder)
        .into_iter()
        .map(|param| atom_string(&param))
        .collect();
    let base_function_map_entries = param_builder
        .reps
        .iter()
        .map(function_map_entry_record)
        .collect::<Vec<_>>();
    let representative_input = param_builder
        .values
        .first()
        .ok_or_else(|| eyre!("3Drep raw evaluator input requires parameter values"))?
        .iter()
        .map(complex_value_record)
        .collect();
    let calls = symbolica_evaluator_call_records(
        processed_atoms,
        param_builder,
        &orientations,
        evaluator_settings,
        evaluator_method,
    )?;

    Ok(SymbolicaEvaluatorInputArchive {
        schema_version: SYMBOLICA_RAW_ARCHIVE_SCHEMA_VERSION,
        description: "Inputs passed to Symbolica when building the 3Drep diagnostic evaluator after Spenso tensor-network processing.".to_string(),
        evaluator_method: evaluator_method.clone(),
        evaluator_backend: evaluator_backend.to_string(),
        evaluator_settings: evaluator_settings_record(evaluator_settings),
        optimization_settings: optimization_settings_record(evaluator_settings),
        parameters: params,
        function_map_entries: base_function_map_entries,
        representative_input,
        calls,
    })
}

fn symbolica_evaluator_call_records(
    processed_atoms: &[Atom],
    param_builder: &ParamBuilder<f64>,
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    evaluator_settings: &EvaluatorSettings,
    evaluator_method: &EvaluatorMethod,
) -> Result<Vec<SymbolicaEvaluatorCallRecord>> {
    let mut calls = vec![symbolica_evaluator_call_record_for_method(
        processed_atoms,
        param_builder,
        orientations,
        evaluator_settings,
        &EvaluatorMethod::SingleParametric,
    )?];
    if evaluator_method != &EvaluatorMethod::SingleParametric {
        calls.push(symbolica_evaluator_call_record_for_method(
            processed_atoms,
            param_builder,
            orientations,
            evaluator_settings,
            evaluator_method,
        )?);
    }
    Ok(calls)
}

fn symbolica_evaluator_call_record_for_method(
    processed_atoms: &[Atom],
    param_builder: &ParamBuilder<f64>,
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    evaluator_settings: &EvaluatorSettings,
    evaluator_method: &EvaluatorMethod,
) -> Result<SymbolicaEvaluatorCallRecord> {
    match evaluator_method {
        EvaluatorMethod::SingleParametric => Ok(SymbolicaEvaluatorCallRecord {
            label: "3d_rep".to_string(),
            expressions: final_symbolica_call_expressions(
                processed_atoms
                    .iter()
                    .map(|atom| GS.collect_orientation_if(atom.as_atom_view(), false)),
                &param_builder.reps,
                evaluator_settings,
            ),
            additional_function_map_entries: Vec::new(),
            dual_shape: None,
        }),
        EvaluatorMethod::Iterative => Ok(SymbolicaEvaluatorCallRecord {
            label: "iterative".to_string(),
            expressions: final_symbolica_call_expressions(
                processed_atoms.iter().flat_map(|atom| {
                    orientations
                        .iter()
                        .map(|orientation| orientation.select_gs(atom.as_atom_view()))
                }),
                &param_builder.reps,
                evaluator_settings,
            ),
            additional_function_map_entries: Vec::new(),
            dual_shape: None,
        }),
        EvaluatorMethod::Summed => Ok(SymbolicaEvaluatorCallRecord {
            label: "summed".to_string(),
            expressions: final_symbolica_call_expressions(
                summed_symbolica_evaluator_expressions(processed_atoms, orientations),
                &param_builder.reps,
                evaluator_settings,
            ),
            additional_function_map_entries: Vec::new(),
            dual_shape: None,
        }),
        EvaluatorMethod::SummedFunctionMap => {
            let (expressions, entries) =
                summed_function_map_symbolica_evaluator_inputs(processed_atoms, orientations)?;
            let mut all_entries = param_builder.reps.clone();
            all_entries.extend(entries.clone());
            Ok(SymbolicaEvaluatorCallRecord {
                label: "summed_function_map".to_string(),
                expressions: final_symbolica_call_expressions(
                    expressions,
                    &all_entries,
                    evaluator_settings,
                ),
                additional_function_map_entries: entries
                    .iter()
                    .map(function_map_entry_record)
                    .collect(),
                dual_shape: None,
            })
        }
    }
}

fn final_symbolica_call_expressions(
    expressions: impl IntoIterator<Item = Atom>,
    function_map_entries: &[gammalooprs::integrands::process::param_builder::FnMapEntry],
    evaluator_settings: &EvaluatorSettings,
) -> Vec<String> {
    let replacements = if evaluator_settings.do_fn_map_replacements {
        function_map_entries
            .iter()
            .map(|entry| entry.replacement())
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };

    expressions
        .into_iter()
        .map(|expression| {
            let expression = expression
                .replace_multiple(&replacements)
                .replace_multiple(&replacements);
            atom_string(&GS.collect_orientation_if(expression, false))
        })
        .collect()
}

fn summed_symbolica_evaluator_expressions(
    processed_atoms: &[Atom],
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
) -> Vec<Atom> {
    processed_atoms
        .iter()
        .map(|atom| {
            orientations
                .iter()
                .map(|orientation| {
                    GS.collect_orientation_if(
                        orientation.orientation_thetas_gs()
                            * orientation.select_gs(atom.as_atom_view()),
                        true,
                    )
                })
                .fold(Atom::Zero, |acc, term| acc + term)
                .replace(
                    Symbol::IF.f([GS.override_if, W_.b_, W_.c_])
                        + Symbol::IF.f([GS.override_if, W_.d_, W_.e_]),
                )
                .repeat()
                .with(Symbol::IF.f([
                    Atom::var(GS.override_if),
                    Atom::var(W_.b_) + Atom::var(W_.d_),
                    Atom::var(W_.c_) + Atom::var(W_.e_),
                ]))
        })
        .collect()
}

fn summed_function_map_symbolica_evaluator_inputs(
    processed_atoms: &[Atom],
    orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
) -> Result<(
    Vec<Atom>,
    Vec<gammalooprs::integrands::process::param_builder::FnMapEntry>,
)> {
    let first_orientation = orientations
        .first()
        .ok_or_else(|| eyre!("summed function-map raw evaluator input requires orientations"))?;
    let entries = processed_atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            let mut args = Vec::new();
            let mut lhs = FunctionBuilder::new(GS.integrand);
            lhs = lhs.add_arg(i);
            for (edge_id, _) in first_orientation {
                lhs = lhs.add_arg(GS.sign(edge_id));
                args.push(Indeterminate::try_from(GS.sign(edge_id)).unwrap());
            }
            Ok(
                gammalooprs::integrands::process::param_builder::FnMapEntry {
                    lhs: lhs.finish(),
                    rhs: GS.collect_orientation_if(atom.as_atom_view(), false),
                    tags: vec![Atom::num(i)],
                    args,
                },
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let expressions = (0..entries.len())
        .map(|i| {
            orientations
                .iter()
                .map(|orientation| {
                    GS.collect_orientation_if(
                        orientation.orientation_thetas_gs() * GS.integrand(i, orientation),
                        true,
                    )
                })
                .fold(Atom::Zero, |acc, term| acc + term)
                .replace(
                    Symbol::IF.f([GS.override_if, W_.b_, W_.c_])
                        + Symbol::IF.f([GS.override_if, W_.d_, W_.e_]),
                )
                .repeat()
                .with(Symbol::IF.f([
                    Atom::var(GS.override_if),
                    Atom::var(W_.b_) + Atom::var(W_.d_),
                    Atom::var(W_.c_) + Atom::var(W_.e_),
                ]))
        })
        .collect();

    Ok((expressions, entries))
}

#[cfg(unix)]
fn make_script_executable(path: &Path) -> Result<()> {
    use std::os::unix::fs::PermissionsExt;

    let mut permissions = fs::metadata(path)?.permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(path, permissions)?;
    Ok(())
}

#[cfg(not(unix))]
fn make_script_executable(_: &Path) -> Result<()> {
    Ok(())
}

fn symbolica_expression_raw_rust_script() -> String {
    r#"#!/usr/bin/env -S rust-script --debug
//! ```cargo
//! [dependencies]
//! color-eyre = "0.6"
//! serde = { version = "1.0", features = ["derive"] }
//! serde_json = "1"
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "d74554ce882c5511876e4b77cbc63670675db6d3", default-features = false, features = ["bincode", "serde"] }
//!
//! [patch.crates-io]
//! graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "d74554ce882c5511876e4b77cbc63670675db6d3" }
//! numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "d74554ce882c5511876e4b77cbc63670675db6d3" }
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "d74554ce882c5511876e4b77cbc63670675db6d3" }
//! ```

#![allow(dead_code)]

use std::{
    env,
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

use color_eyre::{
    eyre::{eyre, Context},
    Result,
};
use serde::Deserialize;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate},
    domains::{
        float::Complex,
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    evaluate::{ExpressionEvaluator, FunctionMap, OptimizationSettings},
    id::{MatchSettings, Replacement},
    parse_lit, symbol, try_parse,
};

#[derive(Debug, Deserialize)]
struct Archive {
    schema_version: u32,
    evaluator_method: String,
    #[serde(default)]
    evaluator_backend: Option<String>,
    #[serde(default)]
    evaluator_settings: Option<EvaluatorSettingsRecord>,
    optimization_settings: OptimizationSettingsRecord,
    parameters: Vec<String>,
    function_map_entries: Vec<FunctionMapEntryRecord>,
    representative_input: Vec<ComplexValueRecord>,
    calls: Vec<EvaluatorCallRecord>,
}

#[derive(Debug, Deserialize)]
struct EvaluatorCallRecord {
    label: String,
    expressions: Vec<String>,
    additional_function_map_entries: Vec<FunctionMapEntryRecord>,
    dual_shape: Option<Vec<Vec<usize>>>,
}

#[derive(Debug, Deserialize, Clone)]
struct FunctionMapEntryRecord {
    lhs: String,
    rhs: String,
    tags: Vec<String>,
    args: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct ComplexValueRecord {
    re: String,
    im: String,
}

#[derive(Debug, Deserialize)]
struct EvaluatorSettingsRecord {
    do_algebra: bool,
    iterative_orientation_optimization: bool,
    summed: bool,
    summed_function_map: bool,
    compile: bool,
    store_atom: bool,
    do_fn_map_replacements: bool,
    horner_iterations: usize,
    n_cores: usize,
    cpe_iterations: Option<usize>,
    abort_level: usize,
    max_horner_scheme_variables: usize,
    max_common_pair_cache_entries: usize,
    max_common_pair_distance: usize,
    verbose: bool,
}

#[derive(Debug, Deserialize)]
struct OptimizationSettingsRecord {
    horner_iterations: usize,
    n_cores: usize,
    cpe_iterations: Option<usize>,
    abort_level: usize,
    max_horner_scheme_variables: usize,
    max_common_pair_cache_entries: usize,
    max_common_pair_distance: usize,
    verbose: bool,
}

type ParsedFnMapEntry = (Atom, Atom, Vec<Atom>, Vec<Indeterminate>);

fn default_input_path() -> Result<PathBuf> {
    if let Some(path) = env::args_os().nth(1) {
        return Ok(PathBuf::from(path));
    }
    Ok(PathBuf::from(file!())
        .canonicalize()
        .with_context(|| format!("Could not canonicalize script path {}", file!()))?
        .parent()
        .ok_or_else(|| eyre!("Could not determine script directory"))?
        .join("symbolica_expression_raw.json"))
}

fn parse_atom(value: &str) -> Result<Atom> {
    try_parse!(value).map_err(|error| eyre!(error))
}

fn parse_complex(value: &ComplexValueRecord) -> Result<Complex<f64>> {
    Ok(Complex::new(value.re.parse()?, value.im.parse()?))
}

fn parse_fn_map_entries(entries: &[FunctionMapEntryRecord]) -> Result<Vec<ParsedFnMapEntry>> {
    entries
        .iter()
        .map(|entry| {
            let lhs = parse_atom(&entry.lhs)?;
            let rhs = parse_atom(&entry.rhs)?;
            let tags = entry
                .tags
                .iter()
                .map(|tag| parse_atom(tag))
                .collect::<Result<Vec<_>>>()?;
            let args = entry
                .args
                .iter()
                .map(|arg| {
                    let atom = parse_atom(arg)?;
                    let atom_display = atom.to_string();
                    atom.try_into().map_err(|_| {
                        eyre!("Expected indeterminate function-map argument, got {atom_display}")
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            Ok((lhs, rhs, tags, args))
        })
        .collect()
}

fn build_function_map(entries: Vec<ParsedFnMapEntry>) -> Result<FunctionMap> {
    let mut fn_map = FunctionMap::new();
    fn_map.add_constant(
        parse_lit!(gammalooprs::x),
        Complex::<Rational>::try_from(Atom::Zero.as_view()).unwrap(),
    );

    for (lhs, rhs, tags, args) in entries {
        if let AtomView::Var(_) = lhs.as_view() {
            if let Ok(value) = Complex::<Rational>::try_from(rhs.as_view()) {
                fn_map.add_constant(lhs, value);
            }
        } else if let AtomView::Fun(function) = lhs.as_view() {
            if tags.is_empty() {
                fn_map
                    .add_function(
                        function.get_symbol(),
                        function.get_symbol().get_name().into(),
                        args.clone(),
                        rhs.clone(),
                    )
                    .map_err(|error| eyre!(error))?;

                let wildcards = args
                    .iter()
                    .enumerate()
                    .map(|(i, arg)| {
                        Replacement::new(
                            Atom::from(arg.clone()).to_pattern(),
                            Atom::var(symbol!(format!("x{i}_"))),
                        )
                        .with_settings(MatchSettings {
                            allow_new_wildcards_on_rhs: true,
                            ..Default::default()
                        })
                    })
                    .collect::<Vec<_>>();
                let _ = lhs.replace_multiple(&wildcards);
            } else {
                fn_map
                    .add_tagged_function(
                        function.get_symbol(),
                        tags,
                        function.get_symbol().get_name().into(),
                        args,
                        rhs,
                    )
                    .map_err(|error| eyre!(error))?;
            }
        }
    }

    Ok(fn_map)
}

fn optimization_settings(record: &OptimizationSettingsRecord) -> OptimizationSettings {
    OptimizationSettings {
        horner_iterations: record.horner_iterations,
        n_cores: record.n_cores,
        cpe_iterations: record.cpe_iterations,
        abort_check: Some(Box::new((|| false) as fn() -> bool)),
        abort_level: record.abort_level,
        max_horner_scheme_variables: record.max_horner_scheme_variables,
        max_common_pair_cache_entries: record.max_common_pair_cache_entries,
        max_common_pair_distance: record.max_common_pair_distance,
        verbose: record.verbose,
        ..OptimizationSettings::default()
    }
}

fn print_replay_options(archive: &Archive) {
    let archived_backend = archive.evaluator_backend.as_deref().unwrap_or_else(|| {
        if archive
            .evaluator_settings
            .as_ref()
            .is_some_and(|settings| settings.compile)
        {
            "compiled backend requested; exact backend not archived"
        } else {
            "eager"
        }
    });
    println!("archived GammaLoop backend: {archived_backend}");
    println!("standalone replay backend: eager Symbolica expression evaluator");
    if let Some(settings) = &archive.evaluator_settings {
        println!(
            "archived evaluator options: compile={}, do_algebra={}, fn_map_replacements={}, store_atom={}",
            settings.compile,
            settings.do_algebra,
            settings.do_fn_map_replacements,
            settings.store_atom
        );
    }
    println!(
        "optimization options: horner_iterations={}, cpe_iterations={:?}, n_cores={}, max_horner_scheme_variables={}",
        archive.optimization_settings.horner_iterations,
        archive.optimization_settings.cpe_iterations,
        archive.optimization_settings.n_cores,
        archive.optimization_settings.max_horner_scheme_variables
    );
}

fn build_evaluator(
    expressions: &[Atom],
    fn_map: &FunctionMap,
    params: &[Atom],
    settings: OptimizationSettings,
) -> Result<ExpressionEvaluator<Complex<f64>>> {
    let mut tree: Option<ExpressionEvaluator<Complex<Fraction<IntegerRing>>>> = None;
    for expression in expressions {
        let evaluator = expression
            .evaluator(fn_map, params, settings.clone())
            .map_err(|error| eyre!("Failed to build evaluator for {expression}: {error}"))?;
        tree = Some(if let Some(mut tree) = tree {
            tree.merge(evaluator, settings.cpe_iterations)
                .map_err(|error| eyre!("Failed to merge evaluator trees: {error}"))?;
            tree
        } else {
            evaluator
        });
    }
    let tree = tree.ok_or_else(|| eyre!("No expressions found in raw evaluator input"))?;
    Ok(tree.map_coeff(&|value| Complex {
        re: value.re.to_f64(),
        im: value.im.to_f64(),
    }))
}

fn run(path: &Path) -> Result<()> {
    let archive: Archive = serde_json::from_str(
        &fs::read_to_string(path)
            .with_context(|| format!("Could not read {}", path.display()))?,
    )
    .with_context(|| format!("Could not parse {}", path.display()))?;
    if archive.schema_version != 1 && archive.schema_version != 2 {
        return Err(eyre!(
            "Unsupported raw evaluator input schema {}; expected 1 or 2",
            archive.schema_version
        ));
    }
    let params = archive
        .parameters
        .iter()
        .map(|param| parse_atom(param))
        .collect::<Result<Vec<_>>>()?;
    let input = archive
        .representative_input
        .iter()
        .map(parse_complex)
        .collect::<Result<Vec<_>>>()?;

    println!("Loaded {}", path.display());
    print_replay_options(&archive);
    println!("parameter count: {}", params.len());
    println!("representative input length: {}", input.len());

    for call in &archive.calls {
        if call.dual_shape.is_some() {
            return Err(eyre!(
                "Raw evaluator replay script does not support dual_shape for call {}",
                call.label
            ));
        }
        let expressions = call
            .expressions
            .iter()
            .map(|expression| parse_atom(expression))
            .collect::<Result<Vec<_>>>()?;
        let mut entries = archive.function_map_entries.clone();
        entries.extend(call.additional_function_map_entries.clone());
        let fn_map = build_function_map(parse_fn_map_entries(&entries)?)?;
        let started = Instant::now();
        let mut evaluator = build_evaluator(
            &expressions,
            &fn_map,
            &params,
            optimization_settings(&archive.optimization_settings),
        )?;
        println!(
            "built call '{}' with {} expression(s) in {:?}",
            call.label,
            expressions.len(),
            started.elapsed()
        );
        let mut output = vec![Complex::new(0.0, 0.0); expressions.len()];
        evaluator.evaluate(&input, &mut output);
        for (index, value) in output.iter().enumerate() {
            println!("{}[{index}] = {} + {} i", call.label, value.re, value.im);
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    color_eyre::install()?;
    run(&default_input_path()?)
}
"#
    .to_string()
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

fn write_symbolica_expression_files(
    atom: &Atom,
    canonical_path: &Path,
    pretty_path: &Path,
) -> Result<()> {
    write_path(canonical_path, &atom.to_canonical_string())?;
    println!(
        "Saved Symbolica expression to {}",
        relative_display(canonical_path)
    );
    write_path(pretty_path, &atom.log_print(Some(80)))?;
    println!(
        "Saved pretty Symbolica expression to {}",
        relative_display(pretty_path)
    );
    Ok(())
}

fn write_symbolica_component_expression_files(
    components: &DiagnosticComponentAtoms,
    canonical_path: &Path,
    pretty_path: &Path,
) -> Result<()> {
    let canonical = render_symbolica_component_expression(components, false);
    write_path(canonical_path, &canonical)?;
    println!(
        "Saved Symbolica expression to {}",
        relative_display(canonical_path)
    );
    let pretty = render_symbolica_component_expression(components, true);
    write_path(pretty_path, &pretty)?;
    println!(
        "Saved pretty Symbolica expression to {}",
        relative_display(pretty_path)
    );
    Ok(())
}

fn render_symbolica_component_expression(
    components: &DiagnosticComponentAtoms,
    pretty: bool,
) -> String {
    let atom_text = |atom: &Atom| {
        if pretty {
            atom.log_print(Some(80))
        } else {
            atom.to_canonical_string()
        }
    };
    let mut text = String::new();
    text.push_str("# 3Drep iterative component expression\n");
    text.push_str("# value = sum_i numerator_component[i] * orientation_component[i]\n");
    text.push_str("# numerator_component[i] is obtained by applying orientation i energy replacements to the processed numerator\n");
    text.push_str(&format!(
        "component_count = {}\n\n",
        components.numerator_atoms.len()
    ));
    text.push_str("[processed_numerator]\n");
    text.push_str(&atom_text(&components.processed_numerator));
    text.push_str("\n\n");
    for (index, orientation) in components.orientation_atoms.iter().enumerate() {
        text.push_str(&format!("[component {index} orientation]\n"));
        text.push_str(&atom_text(orientation));
        text.push_str("\n\n");
    }
    text
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
        direct_relative_tolerance(output.settings.runtime_precision)
    } else {
        mass_shift_relative_tolerance(output.settings.runtime_precision)
    };
    color_text(
        format_scientific_significant(Some(delta), 3),
        if delta.leq_f64(threshold) {
            Color::Green
        } else {
            Color::Red
        },
    )
}

fn render_evaluate_summary(output: &EvaluateOutput, show_parameters: bool) -> String {
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
        "evaluation mode".to_string(),
        color_text(
            if output.numerator_only {
                "numerator-only"
            } else {
                "3D expression"
            },
            Color::Green,
        ),
    ]);
    table.push_record(vec![
        "oriented expression".to_string(),
        output.expression_path.as_ref().map_or_else(
            || color_text("not used for numerator-only", Color::Purple),
            |path| color_text(relative_display(path), Color::Purple),
        ),
    ]);
    table.push_record(vec![
        "symbolica expression".to_string(),
        color_text(
            relative_display(&output.symbolica_expression_path),
            Color::Purple,
        ),
    ]);
    table.push_record(vec![
        "pretty symbolica expression".to_string(),
        color_text(
            relative_display(&output.symbolica_expression_pretty_path),
            Color::Purple,
        ),
    ]);
    table.push_record(vec![
        "raw symbolica expression".to_string(),
        color_text(
            relative_display(&output.symbolica_expression_raw_path),
            Color::Purple,
        ),
    ]);
    table.push_record(vec![
        "raw symbolica replay script".to_string(),
        color_text(
            relative_display(&output.symbolica_expression_raw_script_path),
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
        "build strategy".to_string(),
        color_text(output.settings.build_strategy.label(), Color::Green),
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
        "force eager".to_string(),
        color_text(output.settings.force_eager.to_string(), Color::Blue),
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
    if let Some(profile) = &output.evaluation.profile {
        table.push_record(vec![
            "profile target".to_string(),
            color_text(&profile.target_timing, Color::Yellow),
        ]);
        table.push_record(vec![
            "warmup calls".to_string(),
            color_text(profile.warmup_calls.to_string(), Color::Blue),
        ]);
        table.push_record(vec![
            "warmup time".to_string(),
            color_text(&profile.warmup_timing, Color::Yellow),
        ]);
        table.push_record(vec![
            "profile calls".to_string(),
            color_text(profile.calls.to_string(), Color::Blue),
        ]);
        table.push_record(vec![
            "profile total time".to_string(),
            color_text(&profile.total_timing, Color::Yellow),
        ]);
        table.push_record(vec![
            "profile time per sample".to_string(),
            color_text(&profile.timing_per_sample, Color::Yellow),
        ]);
    } else {
        table.push_record(vec![
            "value".to_string(),
            color_text(&output.evaluation.value, Color::Green),
        ]);
    }

    if !show_parameters {
        return table.build().with(Style::rounded()).to_string();
    }

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
    settings_table.push_record(vec![
        "epsilon steps".to_string(),
        color_text(output.n_epsilon_steps.to_string(), Color::Blue),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimal_complex_distance_keeps_sub_double_precision_differences() {
        let reference = EvaluationValueDecimal {
            re: DecimalValue::parse("+0").unwrap(),
            im: DecimalValue::parse(
                "+2.60633369798997730231172326228159376258650967087478475749183960285052030575615217922069468716473401485881434704408195115849206644935348832565092389878772266122111709668676254415561085167195287053414949997300736045316267041815023568348660819283931785067626512952923333570583338281250466259009884923439956e-5",
            )
            .unwrap(),
        };
        let candidate = EvaluationValueDecimal {
            re: DecimalValue::parse("+0").unwrap(),
            im: DecimalValue::parse(
                "+2.60633369798997730231172326228159376258650967087478475749183960285052030575615217922069468716473401485881434704408195115849206644935348832565092389878772266122111709668676254415561085167195287053414949997300736045316267041815023568348660819283931785067626512952923333570583338281250466259009885027507109e-5",
            )
            .unwrap(),
        };

        let (abs_diff, rel_diff) = complex_distance(&candidate, &reference);

        assert_eq!(
            format_scientific_significant(Some(abs_diff), 3),
            "+1.04e-299"
        );
        assert_eq!(
            format_scientific_significant(Some(rel_diff), 3),
            "+3.99e-295"
        );
    }

    #[test]
    fn profile_duration_parser_accepts_common_units() {
        assert_eq!(
            parse_profile_target_duration("5").unwrap(),
            Duration::from_secs(5)
        );
        assert_eq!(
            parse_profile_target_duration("250ms").unwrap(),
            Duration::from_millis(250)
        );
        assert_eq!(
            parse_profile_target_duration("10 us").unwrap(),
            Duration::from_micros(10)
        );
    }

    #[test]
    fn diagnostic_error_truncation_keeps_first_and_last_200_lines() {
        let text = (0..450)
            .map(|line| format!("line {line}"))
            .collect::<Vec<_>>()
            .join("\n");
        let truncated = truncate_middle_lines(&text, 200, 200);
        let lines = truncated.lines().collect::<Vec<_>>();

        assert_eq!(lines.len(), 401);
        assert_eq!(lines[0], "line 0");
        assert_eq!(lines[199], "line 199");
        assert_eq!(lines[200], "... [truncated 50 diagnostic line(s)] ...");
        assert_eq!(lines[201], "line 250");
        assert_eq!(lines[400], "line 449");
    }

    #[test]
    fn threedrep_precision_dependent_comparison_tolerances_are_configured() {
        assert_eq!(
            direct_relative_tolerance(CliRuntimePrecision::Double),
            1.0e-7
        );
        assert_eq!(
            direct_relative_tolerance(CliRuntimePrecision::Quad),
            1.0e-16
        );
        assert_eq!(
            direct_relative_tolerance(CliRuntimePrecision::ArbPrec),
            1.0e-150
        );
        assert_eq!(
            mass_shift_relative_tolerance(CliRuntimePrecision::Double),
            1.0e-2
        );
        assert_eq!(
            mass_shift_relative_tolerance(CliRuntimePrecision::Quad),
            1.0e-3
        );
        assert_eq!(
            mass_shift_relative_tolerance(CliRuntimePrecision::ArbPrec),
            1.0e-4
        );
    }
}
