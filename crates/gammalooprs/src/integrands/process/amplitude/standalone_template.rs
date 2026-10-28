#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { version = "=3.0.1", default-features = false, features = ["bincode", "float-mpfr", "integer-gmp", "native_code_generation", "serde"] }
//! ```

#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::io::Cursor;
use std::ops::Neg;
use std::path::{Path, PathBuf};
use std::{env, process, slice};

use crate::integrands::process::standalone;
use bincode::config;
use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::evaluate::{
    FunctionRegistrationOptions, InliningPolicy, JITCompiledEvaluator, ParameterRowEvaluator,
    ParameterRowPlan, ParameterRowWorkspace, ParameterRows,
};
use symbolica::prelude::*;
use symbolica::state::StateMap;

const STANDALONE_EVALUATORS_VERSION: u32 = 18;
const ARB_PRECISION_BITS: u32 = 1000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Encode, Decode)]
#[serde(rename_all = "snake_case")]
enum StandaloneNumericTarget {
    Double,
    Quad,
    Arb,
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
struct StandaloneComplexInput {
    re: String,
    im: String,
}

impl StandaloneComplexInput {
    fn parse<T: StandaloneNumber>(&self) -> Result<Complex<T>> {
        Ok(Complex::new(
            T::parse_standalone_input(&self.re)?,
            T::parse_standalone_input(&self.im)?,
        ))
    }
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneEvaluatorArchive<S = Vec<u8>, T = Vec<u8>> {
    version: u32,
    three_d_representations: Vec<String>,
    numeric_target: StandaloneNumericTarget,
    symbolica_state: S,
    graph_terms: Vec<StandaloneGraphTermArchive<T>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneGraphTermArchive<A = Vec<u8>> {
    graph_name: String,
    orientations: Vec<Vec<i8>>,
    param_builder_params: Vec<A>,
    fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    original_integrand: BTreeMap<String, StandaloneEvaluatorStackArchive<A>>,
    threshold_counterterms: Vec<BTreeMap<String, Vec<StandaloneIndexedEvaluatorStackArchive<A>>>>,
    threshold_counterterms_are_variants: bool,
    threshold_variants: Vec<StandaloneAmplitudeThresholdVariant>,
    threshold_multipliers: Option<StandaloneThresholdMultiplierCollectionArchive<A>>,
    metadata_registry: Option<ThresholdCountertermMetadataRegistry>,
}

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode, Serialize, Deserialize,
)]
struct StandaloneCutCFFIndex {
    left_threshold_order: Option<usize>,
    right_threshold_order: Option<usize>,
    lu_cut_order: Option<usize>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneIndexedEvaluatorStackArchive<A = Vec<u8>> {
    cut_cff_index: StandaloneCutCFFIndex,
    evaluator_stack: StandaloneEvaluatorStackArchive<A>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneEvaluatorStackArchive<A = Vec<u8>> {
    explicit_orientation_sum_only: bool,
    production_orientation_ids: Vec<usize>,
    parametric_rows: Option<ParameterRows<A>>,
    single_parametric: StandaloneGenericEvaluatorArchive<A>,
    iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    representative_input: Vec<StandaloneComplexInput>,
    start: usize,
    residue_map_id_start: usize,
    mult_offset: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneGenericEvaluatorArchive<A = Vec<u8>> {
    parameter_rows: Option<ParameterRows<A>>,
    exprs: Vec<A>,
    additional_fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    dual_shape: Option<Vec<Vec<usize>>>,
    #[serde(default)]
    zero_components: Vec<(usize, usize)>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneAmplitudeThresholdVariant {
    variant_id: usize,
    name: String,
    raised_esurface_id: usize,
    generated: bool,
    active: bool,
    requested_subspace: Option<Vec<usize>>,
    requested_parent_lmb: Option<Vec<usize>>,
    resolved_parent_lmb: Vec<usize>,
    resolved_subspace: Vec<usize>,
    subspace_loop_count: usize,
    multiplier_expression: Option<String>,
    multiplier_symmetrize: bool,
    multiplier_opaque_derivatives: bool,
    threshold_edge_sets: Vec<Vec<usize>>,
    explicit_associations: Vec<bool>,
}

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode, Serialize, Deserialize,
)]
enum StandaloneThresholdMultiplierPoint {
    Effective,
    Star,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode, Serialize, Deserialize)]
enum StandaloneThresholdMultiplierInput {
    ModelParameter {
        index: usize,
    },
    AdditionalParameter {
        index: usize,
    },
    ExternalMomentum {
        position: usize,
        component: usize,
    },
    EdgeMomentum {
        point: StandaloneThresholdMultiplierPoint,
        edge: usize,
        component: usize,
    },
    EdgeEnergy {
        point: StandaloneThresholdMultiplierPoint,
        edge: usize,
    },
    Esurface {
        point: StandaloneThresholdMultiplierPoint,
        esurface: usize,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode, Serialize, Deserialize)]
struct StandaloneThresholdMultiplierEsurface {
    edges: Vec<usize>,
    external_shift: Vec<(usize, i64)>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneThresholdMultiplierLayoutArchive<A = Vec<u8>> {
    model_parameter_count: usize,
    additional_parameters: Vec<A>,
    external_count: usize,
    edges: Vec<usize>,
    esurfaces: Vec<StandaloneThresholdMultiplierEsurface>,
    inputs: Vec<StandaloneThresholdMultiplierInput>,
    parameters: Vec<A>,
}

#[derive(Clone, Copy, Debug, Encode, Decode, Serialize, Deserialize)]
struct StandaloneThresholdMultiplierVariantReference {
    variant_id: usize,
    evaluator_id: Option<usize>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneThresholdMultiplierCollectionArchive<A = Vec<u8>> {
    layout: StandaloneThresholdMultiplierLayoutArchive<A>,
    evaluators: Vec<StandaloneGenericEvaluatorArchive<A>>,
    left_variants: Vec<StandaloneThresholdMultiplierVariantReference>,
    right_variants: Vec<StandaloneThresholdMultiplierVariantReference>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum ThresholdCountertermSide {
    Amplitude,
    Left,
    Right,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum ThresholdCountertermOrigin {
    Explicit,
    Autogenerated,
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermAssociationMetadata {
    cut_id: Option<usize>,
    cut_edges: Vec<usize>,
    threshold_edges: Vec<usize>,
    esurface_id: usize,
    eligible: bool,
    origin: ThresholdCountertermOrigin,
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermMultiplierMetadata {
    expression: String,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    function_map: BTreeMap<String, String>,
    symmetrize: bool,
    opaque_derivatives: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermVariantMetadata {
    variant_id: usize,
    name: String,
    #[serde(default)]
    group_id: Option<usize>,
    cut_group_id: Option<usize>,
    associations: Vec<ThresholdCountertermAssociationMetadata>,
    side: ThresholdCountertermSide,
    threshold_esurface_ids: Vec<usize>,
    requested_subspace: Option<Vec<usize>>,
    resolved_subspace: Vec<usize>,
    requested_parent_lmb: Option<Vec<usize>>,
    resolved_parent_lmb: Vec<usize>,
    subspace_loop_count: usize,
    multiplier: Option<ThresholdCountertermMultiplierMetadata>,
    generated: bool,
    active: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermEvaluatorMetadata {
    evaluator_id: usize,
    cut_group_id: Option<usize>,
    collection_evaluator_id: usize,
    expression: String,
    variant_ids: Vec<usize>,
}

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode, Serialize, Deserialize,
)]
#[serde(rename_all = "snake_case")]
enum ThresholdCountertermComponentKind {
    Local,
    Integrated,
    LocalLocal,
    LocalIntegrated,
    IntegratedLocal,
    IntegratedIntegrated,
}

impl ThresholdCountertermComponentKind {
    fn variant_count(self) -> usize {
        match self {
            Self::Local | Self::Integrated => 1,
            Self::LocalLocal
            | Self::LocalIntegrated
            | Self::IntegratedLocal
            | Self::IntegratedIntegrated => 2,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermComponentMetadata {
    component_id: usize,
    cut_group_id: Option<usize>,
    kind: ThresholdCountertermComponentKind,
    variant_ids: Vec<usize>,
    evaluator_ids: Vec<Option<usize>>,
    sign: i8,
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermMetadataRegistry {
    graph_name: String,
    variants: Vec<ThresholdCountertermVariantMetadata>,
    evaluators: Vec<ThresholdCountertermEvaluatorMetadata>,
    components: Vec<ThresholdCountertermComponentMetadata>,
}

type SerializedFnMapEntry<A> = (A, A, Vec<A>, Vec<A>, InliningPolicy, bool);
type ParsedFnMapEntry = (
    Atom,
    Atom,
    Vec<Atom>,
    Vec<Indeterminate>,
    InliningPolicy,
    bool,
);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StandaloneBackend {
    Eager,
    Cpp,
    Assembly,
    Symjit,
}

impl StandaloneBackend {
    fn parse(value: &str) -> Result<Self> {
        match value {
            "eager" => Ok(Self::Eager),
            "c++" | "cpp" => Ok(Self::Cpp),
            "assembly" => Ok(Self::Assembly),
            "symjit" => Ok(Self::Symjit),
            _ => Err(eyre!(
                "Unsupported backend '{value}', expected eager, c++, assembly, or symjit",
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Eager => "eager",
            Self::Cpp => "c++",
            Self::Assembly => "assembly",
            Self::Symjit => "symjit",
        }
    }

    fn inline_asm(self) -> InlineASM {
        match self {
            Self::Assembly => InlineASM::default(),
            Self::Eager | Self::Cpp | Self::Symjit => InlineASM::None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StandaloneMethod {
    SingleParametric,
    Iterative,
    SummedFunctionMap,
    Summed,
}

impl StandaloneMethod {
    fn parse(value: &str) -> Result<Self> {
        match value {
            "single_parametric" | "parametric" => Ok(Self::SingleParametric),
            "iterative" => Ok(Self::Iterative),
            "summed_function_map" | "summed_fnmap" => Ok(Self::SummedFunctionMap),
            "summed" => Ok(Self::Summed),
            _ => Err(eyre!(
                "Unsupported method '{value}', expected single_parametric, iterative, summed_function_map, or summed",
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::SingleParametric => "single_parametric",
            Self::Iterative => "iterative",
            Self::SummedFunctionMap => "summed_function_map",
            Self::Summed => "summed",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum StandaloneStackSelection {
    Original,
    ThresholdCounterterm((usize, usize)),
}

impl StandaloneStackSelection {
    fn parse(value: &str) -> Result<Self> {
        if value == "original" {
            return Ok(Self::Original);
        }

        if let Some(rest) = value.strip_prefix("ct:") {
            let mut parts = rest.split(',');
            let first = parts
                .next()
                .ok_or_else(|| eyre!("Invalid threshold counterterm format"))?
                .parse::<usize>()?;
            let second = parts
                .next()
                .ok_or_else(|| eyre!("Invalid threshold counterterm format"))?
                .parse::<usize>()?;
            return Ok(Self::ThresholdCounterterm((first, second)));
        }

        Err(eyre!(
            "Unsupported stack '{value}', expected original or ct:<first>,<second>",
        ))
    }

    fn label(&self) -> String {
        match self {
            Self::Original => "original".to_string(),
            Self::ThresholdCounterterm((first, second)) => format!("ct_{first}_{second}"),
        }
    }
}

#[derive(Debug)]
struct StandaloneCliOptions {
    input: PathBuf,
    input_json: Option<PathBuf>,
    graph_index: usize,
    graph_name: Option<String>,
    representation: Option<String>,
    stack: StandaloneStackSelection,
    method: StandaloneMethod,
    orientation_index: Option<usize>,
    backend: StandaloneBackend,
    compare_backends: Vec<StandaloneBackend>,
    artifact_dir: Option<PathBuf>,
    print_input: bool,
}

impl Default for StandaloneCliOptions {
    fn default() -> Self {
        Self {
            input: PathBuf::from("standalone_evaluators.bin"),
            input_json: None,
            graph_index: 0,
            graph_name: None,
            representation: None,
            stack: StandaloneStackSelection::Original,
            method: StandaloneMethod::SingleParametric,
            orientation_index: None,
            backend: StandaloneBackend::Eager,
            compare_backends: Vec::new(),
            artifact_dir: None,
            print_input: false,
        }
    }
}

enum StandaloneRuntimeEvaluator<'a> {
    Eager(&'a mut ExpressionEvaluator<Complex<f64>>),
    Compiled(CompiledComplexEvaluator),
    Symjit(Box<JITCompiledEvaluator<Complex<f64>>>),
}

impl<'a> StandaloneRuntimeEvaluator<'a> {
    fn build(
        evaluator: &'a mut ExpressionEvaluator<Complex<f64>>,
        backend: StandaloneBackend,
        artifact_root: &Path,
        label: &str,
    ) -> Result<Self> {
        match backend {
            StandaloneBackend::Eager => Ok(Self::Eager(evaluator)),
            StandaloneBackend::Cpp | StandaloneBackend::Assembly => {
                fs::create_dir_all(artifact_root)?;
                let function_name = format!(
                    "standalone_{}_{}",
                    sanitize_label(label),
                    sanitize_label(backend.as_str())
                );
                let source_path = artifact_root.join(format!("{function_name}.cpp"));
                let library_path = artifact_root.join(format!("{function_name}.so"));
                let exported = evaluator
                    .export_cpp::<Complex<f64>>(
                        &source_path,
                        &function_name,
                        ExportSettings::new()
                            .include_header(true)
                            .inline_asm(backend.inline_asm())
                            .custom_header(None),
                    )
                    .map_err(|error| eyre!(error))?;

                let compiled = exported
                    .compile(&library_path, CompileOptions::default())
                    .map_err(|error| eyre!(error))?
                    .load()
                    .map_err(|error| eyre!(error))?;
                Ok(Self::Compiled(compiled))
            }
            StandaloneBackend::Symjit => Ok(Self::Symjit(Box::new(
                evaluator
                    .jit_compile(
                        JITCompilationSettings::new()
                            .direct_translation(true)
                            .optimization_level(2)
                            .with_option("compact", "false")
                            .with_option("cse", "false"),
                    )
                    .map_err(|error| eyre!(error))?,
            ))),
        }
    }

    fn evaluate(&mut self, args: &[Complex<f64>], out: &mut [Complex<f64>]) {
        match self {
            Self::Eager(evaluator) => evaluator.evaluate(args, out),
            Self::Compiled(evaluator) => evaluator.evaluate(args, out),
            Self::Symjit(evaluator) => evaluator.evaluate(args, out),
        }
    }
}

trait ImportWithMap {
    fn import_with_map(&self, state_map: &StateMap) -> Result<Atom>;
}

impl ImportWithMap for Vec<u8> {
    fn import_with_map(&self, state_map: &StateMap) -> Result<Atom> {
        let mut cursor = Cursor::new(self);
        Atom::import_with_map(&mut cursor, state_map).map_err(|error| eyre!(error))
    }
}

impl ImportWithMap for String {
    fn import_with_map(&self, _: &StateMap) -> Result<Atom> {
        try_parse!(self).map_err(|error| eyre!(error))
    }
}

trait StandaloneNumber:
    Real + RealLike + SingleFloat + Clone + std::fmt::Display + std::fmt::LowerExp
{
    fn parse_standalone_input(value: &str) -> Result<Self>;
    fn exact_from_rational(value: &Rational) -> Self;
    fn zero_value() -> Self;
    fn one_value() -> Self;
}

impl StandaloneNumber for f64 {
    fn parse_standalone_input(value: &str) -> Result<Self> {
        Ok(value.parse()?)
    }

    fn exact_from_rational(value: &Rational) -> Self {
        0.0f64.from_rational(value)
    }

    fn zero_value() -> Self {
        0.0
    }

    fn one_value() -> Self {
        1.0
    }
}

impl StandaloneNumber for DoubleFloat {
    fn parse_standalone_input(value: &str) -> Result<Self> {
        Ok(
            Float::parse(value, Some(DoubleFloat::default().get_precision()))
                .map_err(|error| eyre!(error))?
                .to_double_float(),
        )
    }

    fn exact_from_rational(value: &Rational) -> Self {
        value.into()
    }

    fn zero_value() -> Self {
        DoubleFloat::from(0.0)
    }

    fn one_value() -> Self {
        DoubleFloat::from(1.0)
    }
}

impl StandaloneNumber for Float {
    fn parse_standalone_input(value: &str) -> Result<Self> {
        Float::parse(value, Some(ARB_PRECISION_BITS)).map_err(|error| eyre!(error))
    }

    fn exact_from_rational(value: &Rational) -> Self {
        value.to_multi_prec_float(ARB_PRECISION_BITS)
    }

    fn zero_value() -> Self {
        Float::new(ARB_PRECISION_BITS)
    }

    fn one_value() -> Self {
        Float::new(ARB_PRECISION_BITS).one()
    }
}

fn parse_fn_map_entries<A: ImportWithMap>(
    entries: &[SerializedFnMapEntry<A>],
    state_map: &StateMap,
) -> Result<Vec<ParsedFnMapEntry>> {
    entries
        .iter()
        .map(|(lhs, rhs, tags, args, inlining, is_alias)| {
            let lhs_atom = lhs.import_with_map(state_map)?;
            let rhs_atom = rhs.import_with_map(state_map)?;
            let tags = tags
                .iter()
                .map(|tag| tag.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            let args = args
                .iter()
                .map(|arg| {
                    let arg = arg.import_with_map(state_map)?;
                    if let Ok(indeterminate) = arg.clone().try_into() {
                        Ok(indeterminate)
                    } else {
                        Err(eyre!(
                            "Expected indeterminate in function argument, got {}",
                            arg
                        ))
                    }
                })
                .collect::<Result<Vec<_>>>()?;

            Ok((lhs_atom, rhs_atom, tags, args, *inlining, *is_alias))
        })
        .collect()
}

fn apply_fn_map_entries(
    parsed_entries: Vec<ParsedFnMapEntry>,
) -> Result<(Vec<Replacement>, FunctionMap)> {
    let mut fn_map = FunctionMap::new();
    let mut replacements = Vec::new();

    fn_map
        .add_aliases([(parse_lit!(gammalooprs::x), Atom::Zero)])
        .map_err(|error| eyre!(error))?;
    // Graph and evaluator archives can share definitions. Policy and alias scope
    // are part of their identity; conflicting definitions remain errors.
    let mut seen = HashSet::new();
    for (lhs, rhs, tags, args, inlining, is_alias) in parsed_entries
        .into_iter()
        .filter(|entry| seen.insert(entry.clone()))
    {
        if is_alias {
            if !args.is_empty() || inlining != InliningPolicy::Always {
                return Err(eyre!(
                    "Caller-scope aliases require no formal arguments and Always inlining"
                ));
            }
            fn_map
                .add_aliases([(lhs.clone(), rhs.clone())])
                .map_err(|error| eyre!(error))?;

            continue;
        }
        let options = FunctionRegistrationOptions::new()
            .inlining(inlining)
            .prune_unused_arguments(true);
        let name = match lhs.as_view() {
            AtomView::Fun(function) => Some(function.get_symbol()),
            AtomView::Var(variable) if tags.is_empty() => Some(variable.get_symbol()),
            _ => None,
        };
        if let Some(name) = name {
            fn_map
                .add_tagged_function_with_options(
                    name,
                    tags.clone(),
                    args.clone(),
                    rhs.clone(),
                    options,
                )
                .map_err(|error| eyre!(error))?;
            // Keep ordinary calls in their registered scope: substituting their
            // bodies here would bind retained aliases to the caller's parameters.
        } else {
            if inlining == InliningPolicy::Never {
                return Err(eyre!(
                    "Non-inlined definitions require a symbol or function call"
                ));
            }
            let replacement = Replacement::new(lhs.to_pattern(), rhs);
            replacements.push(replacement);
        }
    }
    Ok((replacements, fn_map))
}

// Backend dispatch consumes the body/output size and optional native plan/prelude.
#[allow(clippy::type_complexity)]
fn build_evaluator<T: StandaloneNumber, A: ImportWithMap>(
    payload: &StandaloneGenericEvaluatorArchive<A>,
    params: &[Atom],
    mut fn_map_entries: Vec<ParsedFnMapEntry>,
    state_map: &StateMap,
    iterate: bool,
    staged_rows: Option<&ParameterRows<Atom>>,
) -> Result<(
    (ExpressionEvaluator<Complex<T>>, usize),
    Option<(ParameterRowPlan, ExpressionEvaluator<Complex<T>>)>,
)>
where
    Complex<T>: EvaluationDomain,
{
    if params.iter().collect::<HashSet<_>>().len() != params.len() {
        return Err(eyre!("Standalone evaluator parameters must be distinct"));
    }
    let optimization_settings = OptimizationSettings::new()
        .horner_iterations(10)
        .cores(1)
        .abort_check(None);
    let cpe_iterations = None;
    let exprs = payload
        .exprs
        .iter()
        .map(|expr| expr.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;

    let additional_fn_map_entries =
        parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?;
    fn_map_entries.extend(additional_fn_map_entries);

    let (replacements, fn_map) = apply_fn_map_entries(fn_map_entries)?;
    let parameter_rows = payload
        .parameter_rows
        .as_ref()
        .map(|rows| {
            let parameters = rows
                .parameters()
                .iter()
                .map(|parameter| parameter.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            ParameterRows::new(parameters, rows.rows().to_vec()).map_err(|error| eyre!(error))
        })
        .transpose()?;
    if parameter_rows.is_some() && staged_rows.is_some() {
        return Err(eyre!(
            "A standalone evaluator cannot have both summed and staged parameter rows"
        ));
    }
    let source_parameters = parameter_rows.as_ref().or(staged_rows).map(|rows| {
        let mut source = params.to_vec();
        source.extend_from_slice(rows.parameters());
        source
    });
    let source_params = source_parameters.as_deref().unwrap_or(params);
    if staged_rows.is_some_and(|rows| rows.rows().is_empty()) {
        return Err(eyre!(
            "A physical residue catalog must contain at least one row"
        ));
    }
    let bound_positions = staged_rows
        .map(|rows| rows.input_positions(source_params))
        .transpose()
        .map_err(|error| eyre!(error))?;
    let build = |roots: &[Atom]| -> Result<ExpressionEvaluator<Complex<Rational>>> {
        let expressions = roots
            .iter()
            .map(|expr| expr.replace_multiple(&replacements))
            .collect::<Vec<_>>();
        let mut builder = Atom::evaluator_multiple(&expressions, source_params)
            .function_map(fn_map.clone())
            .optimization_settings(optimization_settings.clone());
        if let Some(rows) = &parameter_rows {
            builder = builder
                .sum_over_parameters(rows.parameters(), rows.rows())
                .map_err(|error| eyre!(error))?;
        }
        builder.build().map_err(|error| eyre!(error))
    };
    let mut evaluator = if iterate && payload.dual_shape.is_none() {
        let mut tree: Option<ExpressionEvaluator<Complex<Rational>>> = None;
        for expr in &exprs {
            let next = build(slice::from_ref(expr))?;
            tree = Some(if let Some(mut merged) = tree {
                merged
                    .merge(next, cpe_iterations)
                    .map_err(|error| eyre!(error))?;
                merged
            } else {
                next
            });
        }
        tree.ok_or_else(|| eyre!("No expressions in evaluator payload"))?
    } else {
        build(&exprs)?
    };
    let shared = if let Some(rows) = staged_rows {
        let mut staged = ParameterRowEvaluator::new(
            evaluator,
            bound_positions.expect("validated staged row positions"),
            rows.rows().to_vec(),
        )
        .map_err(|error| eyre!(error))?;
        if let Some(shape) = &payload.dual_shape {
            staged = staged
                .vectorize(&Dualizer::new(
                    HyperDual::<Complex<Rational>>::new(shape.clone()),
                    payload.zero_components.clone(),
                ))
                .map_err(|error| eyre!(error))?;
        }
        let (plan, prelude, body) = staged.into_parts();
        evaluator = body;
        Some((
            plan,
            prelude.map_coeff_with_prec(
                &|value| {
                    Complex::new(
                        T::exact_from_rational(&value.re),
                        T::exact_from_rational(&value.im),
                    )
                },
                T::one_value().get_precision(),
            ),
        ))
    } else {
        if let Some(shape) = &payload.dual_shape {
            evaluator = evaluator
                .vectorize(&Dualizer::new(
                    HyperDual::<Complex<Rational>>::new(shape.clone()),
                    payload.zero_components.clone(),
                ))
                .map_err(|error| eyre!(error))?;
        }
        None
    };
    let multiplicity = payload.dual_shape.as_ref().map_or(1, Vec::len);
    let output_len = evaluator.get_output_len();
    if shared.as_ref().map_or_else(
        || evaluator.get_input_len(),
        |(plan, _)| plan.get_input_len(),
    ) != params.len() * multiplicity
        || output_len != exprs.len() * multiplicity
    {
        return Err(eyre!(
            "Standalone evaluator has inconsistent input or output dimensions"
        ));
    }
    Ok((
        (
            evaluator.map_coeff_with_prec(
                &|value| {
                    Complex::new(
                        T::exact_from_rational(&value.re),
                        T::exact_from_rational(&value.im),
                    )
                },
                T::one_value().get_precision(),
            ),
            output_len,
        ),
        shared,
    ))
}

fn set_orientation_values_impl<A: Clone + Neg<Output = A>>(
    values: &mut [A],
    one: A,
    zero: A,
    mult_offset: usize,
    start: usize,
    orientation: &[i8],
) {
    let minus_one = -(one.clone());
    let mut orientation_start = start * mult_offset;

    for value in orientation {
        match value {
            1 => values[orientation_start] = one.clone(),
            -1 => values[orientation_start] = minus_one.clone(),
            0 => values[orientation_start] = zero.clone(),
            _ => panic!("Orientation values must be -1, 0, or 1"),
        }
        orientation_start += mult_offset;
    }
}

fn expected_threshold_multiplier_inputs<A>(
    layout: &StandaloneThresholdMultiplierLayoutArchive<A>,
) -> Vec<StandaloneThresholdMultiplierInput> {
    let mut inputs = (0..layout.model_parameter_count)
        .map(|index| StandaloneThresholdMultiplierInput::ModelParameter { index })
        .collect::<Vec<_>>();
    inputs.extend(
        (0..layout.additional_parameters.len())
            .map(|index| StandaloneThresholdMultiplierInput::AdditionalParameter { index }),
    );
    for position in 0..layout.external_count {
        for component in 0..4 {
            inputs.push(StandaloneThresholdMultiplierInput::ExternalMomentum {
                position,
                component,
            });
        }
    }
    for point in [
        StandaloneThresholdMultiplierPoint::Effective,
        StandaloneThresholdMultiplierPoint::Star,
    ] {
        for &edge in &layout.edges {
            for component in 1..4 {
                inputs.push(StandaloneThresholdMultiplierInput::EdgeMomentum {
                    point,
                    edge,
                    component,
                });
            }
            inputs.push(StandaloneThresholdMultiplierInput::EdgeEnergy { point, edge });
        }
    }
    for point in [
        StandaloneThresholdMultiplierPoint::Effective,
        StandaloneThresholdMultiplierPoint::Star,
    ] {
        for esurface in 0..layout.esurfaces.len() {
            inputs.push(StandaloneThresholdMultiplierInput::Esurface { point, esurface });
        }
    }
    inputs
}

impl<A: PartialEq> StandaloneThresholdMultiplierCollectionArchive<A> {
    fn validate(&self, expected_variants: usize) -> Result<()> {
        let layout = &self.layout;
        if layout.edges.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(eyre!(
                "threshold-multiplier layout graph edges are not strictly ordered"
            ));
        }
        let edge_set = layout.edges.iter().copied().collect::<BTreeSet<_>>();
        for (index, esurface) in layout.esurfaces.iter().enumerate() {
            if esurface.edges.is_empty()
                || esurface.edges.windows(2).any(|pair| pair[0] >= pair[1])
                || esurface
                    .external_shift
                    .windows(2)
                    .any(|pair| pair[0].0 >= pair[1].0)
                || esurface
                    .external_shift
                    .iter()
                    .any(|(_, coefficient)| *coefficient == 0)
            {
                return Err(eyre!(
                    "threshold-multiplier E-surface {index} is not canonical"
                ));
            }
            if let Some(edge) = esurface
                .edges
                .iter()
                .chain(esurface.external_shift.iter().map(|(edge, _)| edge))
                .find(|edge| !edge_set.contains(edge))
            {
                return Err(eyre!(
                    "threshold-multiplier E-surface {index} refers to missing graph edge {edge}"
                ));
            }
        }
        if layout.esurfaces.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(eyre!(
                "threshold-multiplier E-surface catalog is not strictly ordered"
            ));
        }
        let expected_inputs = expected_threshold_multiplier_inputs(layout);
        if layout.inputs != expected_inputs || layout.parameters.len() != expected_inputs.len() {
            return Err(eyre!(
                "threshold-multiplier binding layout is inconsistent: inputs={}, parameters={}, expected={}",
                layout.inputs.len(),
                layout.parameters.len(),
                expected_inputs.len(),
            ));
        }
        let additional_start = layout.model_parameter_count;
        let additional_end = additional_start + layout.additional_parameters.len();
        if layout.parameters.get(additional_start..additional_end)
            != Some(layout.additional_parameters.as_slice())
        {
            return Err(eyre!(
                "threshold-multiplier parameter list does not preserve its declared additional-parameter atoms and order"
            ));
        }
        if self.left_variants.len() != expected_variants || !self.right_variants.is_empty() {
            return Err(eyre!(
                "amplitude threshold-multiplier variant dimensions are {}x{}, expected {expected_variants}x0",
                self.left_variants.len(),
                self.right_variants.len(),
            ));
        }
        if self.evaluators.is_empty() {
            return Err(eyre!(
                "empty threshold-multiplier archive must be represented by null"
            ));
        }
        let mut referenced = vec![false; self.evaluators.len()];
        for (variant_id, reference) in self.left_variants.iter().enumerate() {
            if reference.variant_id != variant_id {
                return Err(eyre!(
                    "amplitude threshold-multiplier references are not ordered by contiguous variant ID"
                ));
            }
            if let Some(evaluator_id) = reference.evaluator_id {
                let slot = referenced.get_mut(evaluator_id).ok_or_else(|| {
                    eyre!(
                        "threshold variant {variant_id} references missing multiplier evaluator {evaluator_id}"
                    )
                })?;
                *slot = true;
            }
        }
        if let Some(evaluator_id) = referenced.iter().position(|referenced| !referenced) {
            return Err(eyre!(
                "threshold-multiplier evaluator {evaluator_id} is not referenced by any variant"
            ));
        }
        for (index, evaluator) in self.evaluators.iter().enumerate() {
            if evaluator.exprs.len() != 1 || evaluator.dual_shape.is_some() {
                return Err(eyre!(
                    "threshold-multiplier evaluator {index} must contain one scalar, non-dual expression"
                ));
            }
        }
        Ok(())
    }
}

impl ThresholdCountertermMetadataRegistry {
    fn validate(&self) -> Result<()> {
        for (variant_id, variant) in self.variants.iter().enumerate() {
            if variant.variant_id != variant_id {
                return Err(eyre!(
                    "threshold variant metadata ID {} is stored at index {variant_id}",
                    variant.variant_id,
                ));
            }
        }
        for (evaluator_id, evaluator) in self.evaluators.iter().enumerate() {
            if evaluator.evaluator_id != evaluator_id
                || evaluator
                    .variant_ids
                    .iter()
                    .any(|variant_id| *variant_id >= self.variants.len())
            {
                return Err(eyre!(
                    "threshold evaluator metadata {evaluator_id} has inconsistent IDs"
                ));
            }
        }
        for (component_id, component) in self.components.iter().enumerate() {
            if component.component_id != component_id
                || component.variant_ids.len() != component.kind.variant_count()
                || component.evaluator_ids.len() != component.kind.variant_count()
            {
                return Err(eyre!(
                    "threshold component metadata {component_id} has inconsistent IDs or dimensions"
                ));
            }
            for (&variant_id, &evaluator_id) in
                component.variant_ids.iter().zip(&component.evaluator_ids)
            {
                let variant = self.variants.get(variant_id).ok_or_else(|| {
                    eyre!(
                        "threshold component {component_id} refers to missing variant {variant_id}"
                    )
                })?;
                if variant.cut_group_id != component.cut_group_id {
                    return Err(eyre!(
                        "threshold component {component_id} and variant {variant_id} have different cut groups"
                    ));
                }
                match evaluator_id {
                    Some(evaluator_id) => {
                        let evaluator = self.evaluators.get(evaluator_id).ok_or_else(|| {
                            eyre!(
                                "threshold component {component_id} refers to missing evaluator {evaluator_id}"
                            )
                        })?;
                        if evaluator.cut_group_id != component.cut_group_id
                            || !evaluator.variant_ids.contains(&variant_id)
                        {
                            return Err(eyre!(
                                "threshold evaluator {evaluator_id} is not registered for component {component_id} variant {variant_id}"
                            ));
                        }
                    }
                    None if variant.multiplier.is_none() => {}
                    None => {
                        return Err(eyre!(
                            "threshold component {component_id} omits the evaluator for multiplier variant {variant_id}"
                        ));
                    }
                }
            }
            let expected_sign = if component.kind.variant_count() == 1 {
                -1
            } else {
                1
            };
            if component.sign != expected_sign {
                return Err(eyre!(
                    "threshold component {component_id} has sign {}, expected {expected_sign}",
                    component.sign,
                ));
            }
        }
        Ok(())
    }

    fn validate_amplitude_payload<A: PartialEq>(
        &self,
        graph_name: &str,
        generalized: bool,
        evaluator_slots: usize,
        variants: &[StandaloneAmplitudeThresholdVariant],
        multipliers: Option<&StandaloneThresholdMultiplierCollectionArchive<A>>,
    ) -> Result<()> {
        self.validate()?;
        if self.graph_name != graph_name {
            return Err(eyre!(
                "amplitude graph '{graph_name}' metadata registry belongs to '{}'",
                self.graph_name,
            ));
        }
        let expected_variant_count = if generalized {
            variants.len()
        } else {
            evaluator_slots
        };
        if self.variants.len() != expected_variant_count {
            return Err(eyre!(
                "amplitude graph '{graph_name}' metadata registry has {} variants but its runtime payload has {expected_variant_count}",
                self.variants.len(),
            ));
        }

        let mut raised_to_esurfaces = BTreeMap::<usize, Vec<usize>>::new();
        let mut esurfaces_to_raised = BTreeMap::<Vec<usize>, usize>::new();
        for (variant_id, metadata) in self.variants.iter().enumerate() {
            if metadata.cut_group_id.is_some()
                || metadata.side != ThresholdCountertermSide::Amplitude
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} has a cut group or non-amplitude side",
                ));
            }
            if metadata.name.trim().is_empty()
                || metadata.subspace_loop_count == 0
                || metadata.resolved_subspace.len() != metadata.subspace_loop_count
                || metadata.resolved_parent_lmb.is_empty()
                || metadata
                    .resolved_parent_lmb
                    .iter()
                    .collect::<BTreeSet<_>>()
                    .len()
                    != metadata.resolved_parent_lmb.len()
                || metadata
                    .resolved_subspace
                    .iter()
                    .collect::<BTreeSet<_>>()
                    .len()
                    != metadata.resolved_subspace.len()
                || metadata
                    .resolved_subspace
                    .iter()
                    .any(|edge| !metadata.resolved_parent_lmb.contains(edge))
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} has invalid name or resolved subspace metadata",
                ));
            }
            for (label, edges) in [
                ("requested subspace", metadata.requested_subspace.as_deref()),
                (
                    "requested parent LMB",
                    metadata.requested_parent_lmb.as_deref(),
                ),
            ] {
                if edges.is_some_and(|edges| {
                    edges.is_empty() || edges.iter().collect::<BTreeSet<_>>().len() != edges.len()
                }) {
                    return Err(eyre!(
                        "amplitude graph '{graph_name}' metadata variant {variant_id} has an invalid {label}",
                    ));
                }
            }
            if metadata.associations.is_empty()
                || metadata.threshold_esurface_ids
                    != metadata
                        .associations
                        .iter()
                        .map(|association| association.esurface_id)
                        .collect::<Vec<_>>()
                || metadata
                    .threshold_esurface_ids
                    .iter()
                    .collect::<BTreeSet<_>>()
                    .len()
                    != metadata.threshold_esurface_ids.len()
                || metadata.associations.iter().any(|association| {
                    association.cut_id.is_some()
                        || !association.cut_edges.is_empty()
                        || association.threshold_edges.is_empty()
                        || association
                            .threshold_edges
                            .windows(2)
                            .any(|pair| pair[0] >= pair[1])
                        || !association.eligible
                })
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} has invalid amplitude associations",
                ));
            }
            if metadata.multiplier.as_ref().is_some_and(|multiplier| {
                multiplier.expression.trim().is_empty()
                    || multiplier.symmetrize
                    || !multiplier.opaque_derivatives
            }) {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} has unsupported multiplier metadata",
                ));
            }
            if !generalized {
                if metadata.multiplier.is_some() {
                    return Err(eyre!(
                        "legacy amplitude graph '{graph_name}' metadata variant {variant_id} cannot have a multiplier",
                    ));
                }
                continue;
            }

            let archived = &variants[variant_id];
            if metadata.name != archived.name
                || metadata.requested_subspace != archived.requested_subspace
                || metadata.requested_parent_lmb != archived.requested_parent_lmb
                || metadata.resolved_parent_lmb != archived.resolved_parent_lmb
                || metadata.resolved_subspace != archived.resolved_subspace
                || metadata.subspace_loop_count != archived.subspace_loop_count
                || metadata.generated != archived.generated
                || metadata.active != archived.active
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} disagrees with its runtime record",
                ));
            }
            if metadata.associations.len() != archived.threshold_edge_sets.len()
                || metadata
                    .associations
                    .iter()
                    .zip(&archived.threshold_edge_sets)
                    .zip(&archived.explicit_associations)
                    .any(|((association, threshold_edges), explicit)| {
                        association.threshold_edges != *threshold_edges
                            || (association.origin == ThresholdCountertermOrigin::Explicit)
                                != *explicit
                    })
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} association provenance disagrees with its runtime record",
                ));
            }
            match (&metadata.multiplier, &archived.multiplier_expression) {
                (None, None) => {}
                (Some(multiplier), Some(expression))
                    if multiplier.expression == *expression
                        && multiplier.symmetrize == archived.multiplier_symmetrize
                        && multiplier.opaque_derivatives
                            == archived.multiplier_opaque_derivatives => {}
                _ => {
                    return Err(eyre!(
                        "amplitude graph '{graph_name}' metadata variant {variant_id} multiplier disagrees with its runtime record",
                    ));
                }
            }
            if let Some(previous) = raised_to_esurfaces.insert(
                archived.raised_esurface_id,
                metadata.threshold_esurface_ids.clone(),
            ) && previous != metadata.threshold_esurface_ids
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' raised E-surface {} has inconsistent threshold E-surface groups",
                    archived.raised_esurface_id,
                ));
            }
            if let Some(previous) = esurfaces_to_raised.insert(
                metadata.threshold_esurface_ids.clone(),
                archived.raised_esurface_id,
            ) && previous != archived.raised_esurface_id
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' threshold E-surface group {:?} refers to inconsistent raised E-surfaces",
                    metadata.threshold_esurface_ids,
                ));
            }
        }

        let collection_references = multipliers.map(|collection| &collection.left_variants);
        let mut variant_evaluators = vec![None; self.variants.len()];
        let mut collection_evaluators = BTreeSet::new();
        for (evaluator_id, evaluator) in self.evaluators.iter().enumerate() {
            if evaluator.cut_group_id.is_some() || evaluator.expression.trim().is_empty() {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata evaluator {evaluator_id} has a cut group or blank expression",
                ));
            }
            if !collection_evaluators.insert(evaluator.collection_evaluator_id) {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata refers to cut-local evaluator {} more than once",
                    evaluator.collection_evaluator_id,
                ));
            }
            let collection = multipliers.ok_or_else(|| {
                eyre!(
                    "amplitude graph '{graph_name}' metadata evaluator {evaluator_id} has no multiplier archive",
                )
            })?;
            if evaluator.collection_evaluator_id >= collection.evaluators.len() {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata evaluator {evaluator_id} refers to missing cut-local evaluator {}",
                    evaluator.collection_evaluator_id,
                ));
            }
            let archived_variant_ids = collection
                .left_variants
                .iter()
                .filter_map(|reference| {
                    (reference.evaluator_id == Some(evaluator.collection_evaluator_id))
                        .then_some(reference.variant_id)
                })
                .collect::<Vec<_>>();
            if evaluator.variant_ids != archived_variant_ids {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata evaluator {evaluator_id} variant links disagree with the multiplier archive",
                ));
            }
            for &variant_id in &evaluator.variant_ids {
                let variant = self.variants.get(variant_id).ok_or_else(|| {
                    eyre!(
                        "amplitude graph '{graph_name}' metadata evaluator {evaluator_id} refers to missing variant {variant_id}",
                    )
                })?;
                if variant.multiplier.is_none()
                    || variant_evaluators[variant_id]
                        .replace(evaluator_id)
                        .is_some()
                {
                    return Err(eyre!(
                        "amplitude graph '{graph_name}' metadata variant {variant_id} has an invalid evaluator assignment",
                    ));
                }
            }
        }
        if collection_evaluators.len()
            != multipliers.map_or(0, |collection| collection.evaluators.len())
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' metadata does not cover every multiplier evaluator",
            ));
        }
        for (variant_id, variant) in self.variants.iter().enumerate() {
            let archived_evaluator = collection_references
                .and_then(|references| references.get(variant_id))
                .and_then(|reference| reference.evaluator_id);
            let metadata_evaluator = variant_evaluators[variant_id]
                .map(|evaluator_id| self.evaluators[evaluator_id].collection_evaluator_id);
            if archived_evaluator != metadata_evaluator
                || variant.multiplier.is_some() != metadata_evaluator.is_some()
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata variant {variant_id} multiplier/evaluator links disagree",
                ));
            }
        }

        let mut actual_components = BTreeSet::new();
        for (component_id, component) in self.components.iter().enumerate() {
            if component.cut_group_id.is_some()
                || component.variant_ids.len() != 1
                || component.evaluator_ids.len() != 1
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata component {component_id} is not a single-variant amplitude component",
                ));
            }
            let variant_id = component.variant_ids[0];
            if component.evaluator_ids[0] != variant_evaluators.get(variant_id).copied().flatten() {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' metadata component {component_id} has an inconsistent evaluator link",
                ));
            }
            if !actual_components.insert((component.kind, variant_id)) {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' has duplicate threshold component metadata for variant {variant_id}",
                ));
            }
        }
        let expected_components = (0..self.variants.len())
            .flat_map(|variant_id| {
                [
                    (ThresholdCountertermComponentKind::Local, variant_id),
                    (ThresholdCountertermComponentKind::Integrated, variant_id),
                ]
            })
            .collect::<BTreeSet<_>>();
        if actual_components != expected_components {
            return Err(eyre!(
                "amplitude graph '{graph_name}' component registry does not contain exactly one local and integrated component per variant",
            ));
        }
        Ok(())
    }
}

impl<A: PartialEq> StandaloneGraphTermArchive<A> {
    fn validate_threshold_payload(&self) -> Result<()> {
        if self.threshold_counterterms_are_variants {
            if self.threshold_counterterms.len() != self.threshold_variants.len() {
                return Err(eyre!(
                    "amplitude graph '{}' has {} variant evaluator slots but {} variant metadata entries",
                    self.graph_name,
                    self.threshold_counterterms.len(),
                    self.threshold_variants.len(),
                ));
            }
        } else if !self.threshold_variants.is_empty() || self.threshold_multipliers.is_some() {
            return Err(eyre!(
                "legacy amplitude graph '{}' cannot contain generalized variant payloads",
                self.graph_name,
            ));
        }
        for (variant_id, variant) in self.threshold_variants.iter().enumerate() {
            if variant.variant_id != variant_id
                || variant.name.trim().is_empty()
                || variant.subspace_loop_count == 0
                || variant.resolved_subspace.len() != variant.subspace_loop_count
                || variant.threshold_edge_sets.len() != variant.explicit_associations.len()
                || variant.threshold_edge_sets.is_empty()
                || variant.threshold_edge_sets.iter().any(|edges| {
                    edges.is_empty() || edges.windows(2).any(|pair| pair[0] >= pair[1])
                })
                || variant.multiplier_symmetrize
                || !variant.multiplier_opaque_derivatives
            {
                return Err(eyre!(
                    "amplitude graph '{}' threshold variant {variant_id} has invalid metadata",
                    self.graph_name,
                ));
            }
        }
        let has_multiplier = self
            .threshold_variants
            .iter()
            .any(|variant| variant.multiplier_expression.is_some());
        if has_multiplier != self.threshold_multipliers.is_some() {
            return Err(eyre!(
                "amplitude graph '{}' multiplier metadata and evaluator payload disagree",
                self.graph_name,
            ));
        }
        if let Some(multipliers) = &self.threshold_multipliers {
            multipliers.validate(self.threshold_variants.len())?;
            if multipliers
                .left_variants
                .iter()
                .enumerate()
                .any(|(variant_id, reference)| {
                    reference.evaluator_id.is_some()
                        != self.threshold_variants[variant_id]
                            .multiplier_expression
                            .is_some()
                })
            {
                return Err(eyre!(
                    "amplitude graph '{}' multiplier references disagree with its variants",
                    self.graph_name,
                ));
            }
        }
        if let Some(registry) = &self.metadata_registry {
            registry.validate_amplitude_payload(
                &self.graph_name,
                self.threshold_counterterms_are_variants,
                self.threshold_counterterms.len(),
                &self.threshold_variants,
                self.threshold_multipliers.as_ref(),
            )?;
        } else if self.threshold_counterterms_are_variants {
            return Err(eyre!(
                "generalized amplitude graph '{}' is missing threshold metadata",
                self.graph_name,
            ));
        }
        Ok(())
    }
}

impl<S, A: PartialEq> StandaloneEvaluatorArchive<S, A> {
    fn validate(&self) -> Result<()> {
        if self.version != STANDALONE_EVALUATORS_VERSION {
            return Err(eyre!(
                "Unsupported version {} (expected {})",
                self.version,
                STANDALONE_EVALUATORS_VERSION,
            ));
        }
        let representations = self.three_d_representations.iter().collect::<BTreeSet<_>>();
        if representations.is_empty()
            || representations.len() != self.three_d_representations.len()
            || representations
                .iter()
                .any(|representation| !matches!(representation.as_str(), "cff" | "ltd"))
        {
            return Err(eyre!(
                "Standalone archive has an invalid generated representation list"
            ));
        }
        for graph in &self.graph_terms {
            if graph.original_integrand.keys().collect::<BTreeSet<_>>() != representations
                || graph
                    .threshold_counterterms
                    .iter()
                    .any(|slot| slot.keys().collect::<BTreeSet<_>>() != representations)
            {
                return Err(eyre!(
                    "Standalone graph '{}' has evaluator representations inconsistent with its generated list",
                    graph.graph_name
                ));
            }
            graph.validate_threshold_payload()?;
            for stack in graph.original_integrand.values().chain(
                graph
                    .threshold_counterterms
                    .iter()
                    .flat_map(|slot| slot.values())
                    .flatten()
                    .map(|entry| &entry.evaluator_stack),
            ) {
                standalone::validate_standalone_input_layout(
                    graph.param_builder_params.len(),
                    Some(stack.representative_input.len()),
                    stack.mult_offset,
                    [stack.start, stack.residue_map_id_start],
                    std::iter::once(&stack.single_parametric)
                        .chain(stack.iterative.as_ref())
                        .chain(stack.summed_function_map.as_ref())
                        .chain(stack.summed.as_ref())
                        .map(|payload| payload.dual_shape.as_deref()),
                )?;
            }
        }
        Ok(())
    }
}

impl<A> StandaloneEvaluatorStackArchive<A> {
    fn selected_payload(
        &self,
        method: StandaloneMethod,
    ) -> Result<(&StandaloneGenericEvaluatorArchive<A>, bool)> {
        match method {
            StandaloneMethod::SingleParametric => Ok((&self.single_parametric, false)),
            StandaloneMethod::Iterative => self
                .iterative
                .as_ref()
                .map(|payload| (payload, true))
                .ok_or_else(|| eyre!("Missing iterative evaluator in standalone archive")),
            StandaloneMethod::SummedFunctionMap => self
                .summed_function_map
                .as_ref()
                .map(|payload| (payload, false))
                .ok_or_else(|| {
                    eyre!("Missing summed_function_map evaluator in standalone archive")
                }),
            StandaloneMethod::Summed => self
                .summed
                .as_ref()
                .map(|payload| (payload, false))
                .ok_or_else(|| eyre!("Missing summed evaluator in standalone archive")),
        }
    }

    fn representative_input<T: StandaloneNumber>(&self) -> Result<Vec<Complex<T>>> {
        self.representative_input
            .iter()
            .map(StandaloneComplexInput::parse::<T>)
            .collect()
    }

    fn set_orientation<T: StandaloneNumber>(
        &self,
        orientation_index: usize,
        orientation: &[i8],
    ) -> Result<Vec<Complex<T>>> {
        let mut input = self.representative_input::<T>()?;
        let production_orientation_id = self
            .production_orientation_ids
            .get(orientation_index)
            .ok_or_else(|| {
                eyre!(
                    "Missing production residue-map ID for runtime orientation {}",
                    orientation_index
                )
            })?;
        input[self.residue_map_id_start * self.mult_offset] = Complex::new(
            T::exact_from_rational(&Rational::from(*production_orientation_id)),
            T::zero_value(),
        );
        let zero = T::zero_value();
        let one = T::one_value();
        set_orientation_values_impl(
            &mut input,
            Complex::new(one, T::zero_value()),
            Complex::new(zero, T::zero_value()),
            self.mult_offset,
            self.start,
            orientation,
        );
        Ok(input)
    }
}

impl<A> StandaloneGraphTermArchive<A> {
    fn stack(
        &self,
        representation: &str,
        selection: &StandaloneStackSelection,
    ) -> Result<&StandaloneEvaluatorStackArchive<A>> {
        match selection {
            StandaloneStackSelection::Original => self
                .original_integrand
                .get(representation)
                .ok_or_else(|| eyre!("Representation {representation} was not generated")),
            StandaloneStackSelection::ThresholdCounterterm((first, second)) => self
                .threshold_counterterms
                .get(*first)
                .and_then(|representations| representations.get(representation))
                .and_then(|orders| orders.get(*second))
                .map(|entry| &entry.evaluator_stack)
                .ok_or_else(|| {
                    eyre!(
                        "Threshold counterterm index {},{} is out of range for graph {}",
                        first,
                        second,
                        self.graph_name
                    )
                }),
        }
    }
}

impl<S, A> StandaloneEvaluatorArchive<S, A> {
    fn graph_term(
        &self,
        graph_index: usize,
        graph_name: Option<&str>,
    ) -> Result<&StandaloneGraphTermArchive<A>> {
        if let Some(graph_name) = graph_name {
            return self
                .graph_terms
                .iter()
                .find(|term| term.graph_name == graph_name)
                .ok_or_else(|| eyre!("Unknown graph '{graph_name}'"));
        }

        self.graph_terms.get(graph_index).ok_or_else(|| {
            eyre!(
                "Graph index {} is out of range for {} graph terms",
                graph_index,
                self.graph_terms.len()
            )
        })
    }
}

fn evaluate_eager<T: StandaloneNumber>(
    evaluator: &mut ExpressionEvaluator<Complex<T>>,
    mut shared: Option<&mut (ParameterRowPlan, ExpressionEvaluator<Complex<T>>)>,
    output_len: usize,
    inputs: &[Vec<Complex<T>>],
) -> Result<Vec<Complex<T>>>
where
    Complex<T>: EvaluationDomain,
{
    let zero = Complex::new(T::zero_value(), T::zero_value());
    let mut accumulated = vec![zero.clone(); output_len];
    let mut current = accumulated.clone();
    let mut workspace = ParameterRowWorkspace::default();
    for input in inputs {
        if let Some((plan, prelude)) = shared.as_mut() {
            plan.evaluate_with(
                input,
                &mut current,
                &mut workspace,
                |value| Complex::new(T::exact_from_rational(value), T::zero_value()),
                |input, output| prelude.try_evaluate(input, output),
                |input, output| evaluator.try_evaluate(input, output),
            )
            .map_err(|error| eyre!(error))?;
        } else {
            evaluator
                .try_evaluate(input, &mut current)
                .map_err(|error| eyre!(error))?;
        }
        for (sum, value) in accumulated.iter_mut().zip(&current) {
            *sum += value;
        }
    }
    Ok(accumulated)
}

#[allow(clippy::too_many_arguments)]
fn evaluate_with_backend_f64(
    evaluator: &mut ExpressionEvaluator<Complex<f64>>,
    shared: Option<&mut (ParameterRowPlan, ExpressionEvaluator<Complex<f64>>)>,
    backend: StandaloneBackend,
    inputs: &[Vec<Complex<f64>>],
    artifact_root: &Path,
    label: &str,
) -> Result<Vec<Complex<f64>>> {
    let output_len = evaluator.get_output_len();
    let input_len = evaluator.get_input_len();
    let mut runtime = StandaloneRuntimeEvaluator::build(evaluator, backend, artifact_root, label)?;
    let mut shared = shared
        .map(|(plan, prelude)| -> Result<_> {
            Ok((
                plan,
                (prelude.get_output_len() != 0)
                    .then(|| {
                        StandaloneRuntimeEvaluator::build(
                            prelude,
                            backend,
                            artifact_root,
                            &format!("{label}_auxiliary"),
                        )
                    })
                    .transpose()?,
            ))
        })
        .transpose()?;
    let mut accumulated = vec![Complex::new(0.0, 0.0); output_len];
    let mut current = accumulated.clone();
    let mut workspace = ParameterRowWorkspace::default();
    for input in inputs {
        if let Some((plan, prelude)) = &mut shared {
            plan.evaluate_with(
                input,
                &mut current,
                &mut workspace,
                |value| Complex::new(value.to_f64(), 0.0),
                |input, output| {
                    prelude
                        .as_mut()
                        .expect("the library plan skips an empty prelude")
                        .evaluate(input, output);
                    Ok(())
                },
                |input, output| {
                    runtime.evaluate(input, output);
                    Ok(())
                },
            )
            .map_err(|error| eyre!(error))?;
        } else {
            if input.len() != input_len {
                return Err(eyre!(
                    "Physical input length does not match standalone evaluator"
                ));
            }
            runtime.evaluate(input, &mut current);
        }
        for (sum, value) in accumulated.iter_mut().zip(&current) {
            *sum += value;
        }
    }
    Ok(accumulated)
}

fn parse_backend_list(value: &str) -> Result<Vec<StandaloneBackend>> {
    value
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(StandaloneBackend::parse)
        .collect()
}

fn sanitize_label(value: &str) -> String {
    value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() {
                character.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect()
}

fn print_usage(program: &str) {
    eprintln!(
        "Usage: {program} [standalone_evaluators.bin|json] [options]\n\
         \n\
         Options:\n\
           --backend <eager|c++|assembly|symjit>\n\
           --compare-backends <backend[,backend...]>\n\
           --input-json <path>\n\
           --graph-index <usize>\n\
           --graph-name <name>\n\
           --three-dimensional-representation <cff|ltd>\n\
           --stack <original|ct:N,M>\n\
           --method <single_parametric|iterative|summed_function_map|summed>\n\
           --orientation-index <usize> (single_parametric only)\n\
           --artifact-dir <path>\n\
           --print-input\n\
           --help"
    );
}

fn parse_cli_options() -> Result<StandaloneCliOptions> {
    let mut options = StandaloneCliOptions::default();
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--help" | "-h" => {
                let program = env::args()
                    .next()
                    .unwrap_or_else(|| "standalone_evaluators_rust.rs".to_string());
                print_usage(&program);
                process::exit(0);
            }
            "--backend" => {
                options.backend = StandaloneBackend::parse(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --backend"))?,
                )?;
            }
            "--compare-backends" => {
                options.compare_backends = parse_backend_list(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --compare-backends"))?,
                )?;
            }
            "--input-json" => {
                options.input_json = Some(PathBuf::from(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --input-json"))?,
                ));
            }
            "--graph-index" => {
                options.graph_index = args
                    .next()
                    .ok_or_else(|| eyre!("Missing value for --graph-index"))?
                    .parse()?;
            }
            "--three-dimensional-representation" => {
                options.representation = Some(args.next().ok_or_else(|| {
                    eyre!("Missing value for --three-dimensional-representation")
                })?);
            }
            "--graph-name" => {
                options.graph_name = Some(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --graph-name"))?,
                );
            }
            "--stack" => {
                options.stack = StandaloneStackSelection::parse(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --stack"))?,
                )?;
            }
            "--method" => {
                options.method = StandaloneMethod::parse(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --method"))?,
                )?;
            }
            "--orientation-index" => {
                options.orientation_index = Some(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --orientation-index"))?
                        .parse()?,
                );
            }
            "--artifact-dir" => {
                options.artifact_dir = Some(PathBuf::from(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --artifact-dir"))?,
                ));
            }
            "--print-input" => {
                options.print_input = true;
            }
            _ if arg.starts_with("--") => {
                return Err(eyre!("Unsupported option '{arg}'"));
            }
            _ => {
                options.input = PathBuf::from(arg);
            }
        }
    }

    if options.orientation_index.is_some() && options.method != StandaloneMethod::SingleParametric {
        return Err(eyre!(
            "--orientation-index can only be used with --method single_parametric"
        ));
    }

    Ok(options)
}

#[derive(Deserialize)]
#[serde(untagged)]
enum RawStandaloneInputValue {
    StringComponents { re: String, im: String },
    NumericPair([f64; 2]),
}

impl From<RawStandaloneInputValue> for StandaloneComplexInput {
    fn from(value: RawStandaloneInputValue) -> Self {
        match value {
            RawStandaloneInputValue::StringComponents { re, im } => Self { re, im },
            RawStandaloneInputValue::NumericPair([re, im]) => Self {
                re: re.to_string(),
                im: im.to_string(),
            },
        }
    }
}

fn load_custom_input(path: impl AsRef<Path>) -> Result<Vec<StandaloneComplexInput>> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let raw: Vec<RawStandaloneInputValue> = serde_json::from_slice(&binary)
        .with_context(|| format!("Failed to parse {}", path.as_ref().display()))?;
    Ok(raw.into_iter().map(StandaloneComplexInput::from).collect())
}

fn load_bin(path: impl AsRef<Path>) -> Result<(StandaloneEvaluatorArchive, StateMap)> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (version, _): (u32, _) = bincode::decode_from_slice(&binary, config::standard())?;
    if version != STANDALONE_EVALUATORS_VERSION {
        return Err(eyre!(
            "Unsupported version {version} (expected {}); regenerate the standalone archive",
            STANDALONE_EVALUATORS_VERSION
        ));
    }
    let (archive, _): (StandaloneEvaluatorArchive, _) =
        bincode::decode_from_slice(&binary, config::standard())?;
    archive.validate()?;
    let state_map = State::import(&mut Cursor::new(&archive.symbolica_state), None)
        .map_err(|error| eyre!(error))?;
    Ok((archive, state_map))
}

fn load_json(path: impl AsRef<Path>) -> Result<StandaloneEvaluatorArchive<(), String>> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let archive: StandaloneEvaluatorArchive<(), String> = serde_json::from_slice(&binary)?;
    archive.validate()?;
    Ok(archive)
}

fn current_state_map() -> Result<StateMap> {
    let mut symbolica_state = Vec::new();
    State::export(&mut symbolica_state)
        .with_context(|| "Failed to export Symbolica state for standalone evaluators")?;

    let mut state_cursor = Cursor::new(&symbolica_state);
    State::import(&mut state_cursor, None).map_err(|error| eyre!(error))
}

fn validate_backends(
    numeric_target: StandaloneNumericTarget,
    compare_backends: &[StandaloneBackend],
) -> Result<()> {
    if numeric_target == StandaloneNumericTarget::Double {
        return Ok(());
    }

    if compare_backends.len() > 1 {
        return Err(eyre!(
            "Only a single eager backend is available for {} standalone exports",
            match numeric_target {
                StandaloneNumericTarget::Double => "double",
                StandaloneNumericTarget::Quad => "quad",
                StandaloneNumericTarget::Arb => "arb",
            }
        ));
    }

    if compare_backends
        .iter()
        .any(|backend| *backend != StandaloneBackend::Eager)
    {
        return Err(eyre!(
            "Only the eager backend is available for {} standalone exports",
            match numeric_target {
                StandaloneNumericTarget::Double => "double",
                StandaloneNumericTarget::Quad => "quad",
                StandaloneNumericTarget::Arb => "arb",
            }
        ));
    }

    Ok(())
}

fn diff_ratio(lhs: Complex<f64>, rhs: Complex<f64>) -> Option<Complex<f64>> {
    if rhs.re == 0.0 && rhs.im == 0.0 {
        None
    } else {
        Some(lhs / rhs)
    }
}

fn evaluate_double_archive<A: ImportWithMap, S>(
    archive: StandaloneEvaluatorArchive<S, A>,
    state_map: &StateMap,
    options: &StandaloneCliOptions,
    custom_input: Option<&[StandaloneComplexInput]>,
) -> Result<()> {
    let compare_backends = if options.compare_backends.is_empty() {
        vec![options.backend]
    } else {
        options.compare_backends.clone()
    };
    validate_backends(StandaloneNumericTarget::Double, &compare_backends)?;

    let representation = options
        .representation
        .as_ref()
        .or_else(|| archive.three_d_representations.first())
        .ok_or_else(|| eyre!("No generated 3D representation"))?;
    let graph = archive.graph_term(options.graph_index, options.graph_name.as_deref())?;
    let stack = graph.stack(representation, &options.stack)?;
    let (payload, iterate) = stack.selected_payload(options.method)?;

    if stack.explicit_orientation_sum_only && options.orientation_index.is_some() {
        return Err(eyre!(
            "`--orientation-index` is invalid for an explicit orientation-sum evaluator because its single-parametric expression already contains the complete orientation sum"
        ));
    }

    let params = graph
        .param_builder_params
        .iter()
        .map(|param| param.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;
    let fn_map_entries = parse_fn_map_entries(&graph.fn_map_entries, state_map)?;
    let rows = stack
        .parametric_rows
        .as_ref()
        .map(|rows| {
            let parameters = rows
                .parameters()
                .iter()
                .map(|parameter| parameter.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            ParameterRows::new(parameters, rows.rows().to_vec()).map_err(|error| eyre!(error))
        })
        .transpose()?;
    let staged_rows = if options.method == StandaloneMethod::SingleParametric {
        rows.as_ref()
    } else {
        None
    };
    let ((mut evaluator, _), mut shared) = build_evaluator::<f64, _>(
        payload,
        &params,
        fn_map_entries,
        state_map,
        iterate,
        staged_rows,
    )?;

    let inputs = if let Some(custom_input) = custom_input {
        if custom_input.len() != stack.representative_input.len() {
            return Err(eyre!(
                "Custom input length differs from the physical parameter count"
            ));
        }
        vec![
            custom_input
                .iter()
                .map(StandaloneComplexInput::parse::<f64>)
                .collect::<Result<Vec<_>>>()?,
        ]
    } else {
        match options.method {
            StandaloneMethod::SingleParametric if stack.explicit_orientation_sum_only => {
                vec![stack.representative_input::<f64>()?]
            }
            StandaloneMethod::SingleParametric => {
                if let Some(index) = options.orientation_index {
                    let orientation = graph.orientations.get(index).ok_or_else(|| {
                        eyre!(
                            "Orientation index {} is out of range for {} orientations",
                            index,
                            graph.orientations.len()
                        )
                    })?;
                    vec![stack.set_orientation::<f64>(index, orientation)?]
                } else {
                    graph
                        .orientations
                        .iter()
                        .enumerate()
                        .map(|(index, orientation)| {
                            stack.set_orientation::<f64>(index, orientation)
                        })
                        .collect::<Result<Vec<_>>>()?
                }
            }
            StandaloneMethod::Iterative
            | StandaloneMethod::SummedFunctionMap
            | StandaloneMethod::Summed => vec![stack.representative_input::<f64>()?],
        }
    };

    let artifact_root = options.artifact_dir.clone().unwrap_or_else(|| {
        options
            .input
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join("standalone_backend_artifacts")
    });

    println!(
        "precision=double graph={} stack={} method={} orientation={} artifact_dir={}",
        graph.graph_name,
        options.stack.label(),
        options.method.as_str(),
        options
            .orientation_index
            .map(|index| index.to_string())
            .unwrap_or_else(|| "all".to_string()),
        artifact_root.display()
    );

    if options.print_input {
        for (index, input) in inputs.iter().enumerate() {
            println!("input[{index}]={input:?}");
        }
    }

    let mut results = Vec::new();
    for backend in compare_backends {
        let values = evaluate_with_backend_f64(
            &mut evaluator,
            shared.as_mut(),
            backend,
            &inputs,
            &artifact_root.join(sanitize_label(&graph.graph_name)),
            &format!(
                "{}_{representation}_{}_{}",
                graph.graph_name,
                options.stack.label(),
                options.method.as_str(),
            ),
        )?;
        println!("backend={}", backend.as_str());
        for (index, value) in values.iter().enumerate() {
            println!("  result[{index}] = {value}");
        }
        results.push((backend, values));
    }

    if results.len() == 2 {
        let (lhs_backend, lhs_values) = &results[0];
        let (rhs_backend, rhs_values) = &results[1];
        println!(
            "comparison {} -> {}",
            lhs_backend.as_str(),
            rhs_backend.as_str()
        );
        for (index, (lhs_value, rhs_value)) in lhs_values.iter().zip(rhs_values).enumerate() {
            println!("  diff[{index}] = {}", *rhs_value - *lhs_value);
            match diff_ratio(*rhs_value, *lhs_value) {
                Some(ratio) => println!("  ratio[{index}] = {ratio}"),
                None => println!("  ratio[{index}] = undefined"),
            }
        }
    }

    Ok(())
}

fn evaluate_higher_precision_archive<T: StandaloneNumber, A: ImportWithMap, S>(
    archive: StandaloneEvaluatorArchive<S, A>,
    state_map: &StateMap,
    options: &StandaloneCliOptions,
    custom_input: Option<&[StandaloneComplexInput]>,
    precision_label: &str,
    numeric_target: StandaloneNumericTarget,
) -> Result<()>
where
    Complex<T>: EvaluationDomain,
{
    let compare_backends = if options.compare_backends.is_empty() {
        vec![options.backend]
    } else {
        options.compare_backends.clone()
    };
    validate_backends(numeric_target, &compare_backends)?;

    let representation = options
        .representation
        .as_ref()
        .or_else(|| archive.three_d_representations.first())
        .ok_or_else(|| eyre!("No generated 3D representation"))?;
    let graph = archive.graph_term(options.graph_index, options.graph_name.as_deref())?;
    let stack = graph.stack(representation, &options.stack)?;
    let (payload, iterate) = stack.selected_payload(options.method)?;

    if stack.explicit_orientation_sum_only && options.orientation_index.is_some() {
        return Err(eyre!(
            "`--orientation-index` is invalid for an explicit orientation-sum evaluator because its single-parametric expression already contains the complete orientation sum"
        ));
    }
    let params = graph
        .param_builder_params
        .iter()
        .map(|param| param.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;
    let fn_map_entries = parse_fn_map_entries(&graph.fn_map_entries, state_map)?;
    let rows = stack
        .parametric_rows
        .as_ref()
        .map(|rows| {
            let parameters = rows
                .parameters()
                .iter()
                .map(|parameter| parameter.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            ParameterRows::new(parameters, rows.rows().to_vec()).map_err(|error| eyre!(error))
        })
        .transpose()?;
    let staged_rows = if options.method == StandaloneMethod::SingleParametric {
        rows.as_ref()
    } else {
        None
    };
    let ((mut evaluator, output_len), mut shared) = build_evaluator::<T, _>(
        payload,
        &params,
        fn_map_entries,
        state_map,
        iterate,
        staged_rows,
    )?;

    let inputs = if let Some(custom_input) = custom_input {
        if custom_input.len() != stack.representative_input.len() {
            return Err(eyre!(
                "Custom input length differs from the physical parameter count"
            ));
        }
        vec![
            custom_input
                .iter()
                .map(StandaloneComplexInput::parse::<T>)
                .collect::<Result<Vec<_>>>()?,
        ]
    } else {
        match options.method {
            StandaloneMethod::SingleParametric if stack.explicit_orientation_sum_only => {
                vec![stack.representative_input::<T>()?]
            }
            StandaloneMethod::SingleParametric => {
                if let Some(index) = options.orientation_index {
                    let orientation = graph.orientations.get(index).ok_or_else(|| {
                        eyre!(
                            "Orientation index {} is out of range for {} orientations",
                            index,
                            graph.orientations.len()
                        )
                    })?;
                    vec![stack.set_orientation::<T>(index, orientation)?]
                } else {
                    graph
                        .orientations
                        .iter()
                        .enumerate()
                        .map(|(index, orientation)| stack.set_orientation::<T>(index, orientation))
                        .collect::<Result<Vec<_>>>()?
                }
            }
            StandaloneMethod::Iterative
            | StandaloneMethod::SummedFunctionMap
            | StandaloneMethod::Summed => vec![stack.representative_input::<T>()?],
        }
    };

    println!(
        "precision={} graph={} stack={} method={} orientation={} backend=eager",
        precision_label,
        graph.graph_name,
        options.stack.label(),
        options.method.as_str(),
        options
            .orientation_index
            .map(|index| index.to_string())
            .unwrap_or_else(|| "all".to_string()),
    );

    if options.print_input {
        for (index, input) in inputs.iter().enumerate() {
            println!("input[{index}]={input:?}");
        }
    }

    let values = evaluate_eager(&mut evaluator, shared.as_mut(), output_len, &inputs)?;
    for (index, value) in values.iter().enumerate() {
        println!("  result[{index}] = {value}");
    }

    Ok(())
}

fn main() -> Result<()> {
    let options = parse_cli_options()?;
    let input = options.input.clone();

    let Some(extension) = input.extension() else {
        return Err(eyre!("No extension, expected .bin or .json"));
    };
    let custom_input = options
        .input_json
        .as_ref()
        .map(load_custom_input)
        .transpose()?;

    match extension.to_string_lossy().as_ref() {
        "bin" => {
            let (archive, state_map) = load_bin(&input)?;
            match archive.numeric_target {
                StandaloneNumericTarget::Double => {
                    evaluate_double_archive(archive, &state_map, &options, custom_input.as_deref())
                }
                StandaloneNumericTarget::Quad => {
                    evaluate_higher_precision_archive::<DoubleFloat, _, _>(
                        archive,
                        &state_map,
                        &options,
                        custom_input.as_deref(),
                        "quad",
                        StandaloneNumericTarget::Quad,
                    )
                }
                StandaloneNumericTarget::Arb => evaluate_higher_precision_archive::<Float, _, _>(
                    archive,
                    &state_map,
                    &options,
                    custom_input.as_deref(),
                    "arb",
                    StandaloneNumericTarget::Arb,
                ),
            }
        }
        "json" => {
            let archive = load_json(&input)?;
            let state_map = current_state_map()?;
            match archive.numeric_target {
                StandaloneNumericTarget::Double => {
                    evaluate_double_archive(archive, &state_map, &options, custom_input.as_deref())
                }
                StandaloneNumericTarget::Quad => {
                    evaluate_higher_precision_archive::<DoubleFloat, _, _>(
                        archive,
                        &state_map,
                        &options,
                        custom_input.as_deref(),
                        "quad",
                        StandaloneNumericTarget::Quad,
                    )
                }
                StandaloneNumericTarget::Arb => evaluate_higher_precision_archive::<Float, _, _>(
                    archive,
                    &state_map,
                    &options,
                    custom_input.as_deref(),
                    "arb",
                    StandaloneNumericTarget::Arb,
                ),
            }
        }
        _ => Err(eyre!(
            "Unsupported file extension {}, expected .bin or .json",
            extension.to_string_lossy()
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::printer::PrintOptions;

    #[test]
    fn standalone_exact_row_payload_restores_sum_and_retained_derivatives() -> Result<()> {
        let params = [
            parse_lit!(standalone_exact_rows::x),
            parse_lit!(standalone_exact_rows::m),
        ];
        let coefficient = parse_lit!(standalone_exact_rows::c);
        let u = parse_lit!(standalone_exact_rows::u);
        let v = parse_lit!(standalone_exact_rows::v);
        let numerator = symbol!("standalone_exact_rows::numerator");
        let print = |atom: &Atom| atom.printer(PrintOptions::file()).to_string();
        let root = Symbol::IF.call_args([
            coefficient.clone(),
            &coefficient * numerator.call_args(params.clone()),
            Atom::Zero,
        ]);
        let rows = ParameterRows::new(
            vec![print(&coefficient)],
            vec![
                vec![Rational::zero()],
                vec![Rational::new(1, 4)],
                vec![Rational::new(1, 2)],
                vec![Rational::new(1, 2)],
            ],
        )
        .map_err(|error| eyre!(error))?;
        for inlining in [InliningPolicy::Always, InliningPolicy::Never] {
            let payload = StandaloneGenericEvaluatorArchive {
                parameter_rows: Some(rows.clone()),
                exprs: vec![print(&root), print(&params[0])],

                additional_fn_map_entries: vec![(
                    print(&numerator.call_args([u.clone(), v.clone()])),
                    print(&(u.pow(2) + v.pow(Atom::num(Rational::new(1, 2))))),
                    Vec::new(),
                    vec![print(&u), print(&v)],
                    inlining,
                    false,
                )],
                dual_shape: Some(vec![vec![0], vec![1]]),
                zero_components: vec![(1, 1)],
            };
            let json: StandaloneGenericEvaluatorArchive<String> =
                serde_json::from_slice(&serde_json::to_vec(&payload)?)?;
            let bytes = bincode::encode_to_vec(&payload, config::standard())?;
            let (binary, consumed): (StandaloneGenericEvaluatorArchive<String>, _) =
                bincode::decode_from_slice(&bytes, config::standard())?;
            assert_eq!(consumed, bytes.len());
            for restored in [json, binary] {
                assert_eq!(restored.parameter_rows.as_ref(), Some(&rows));
                for dual in [false, true] {
                    let mut payload = restored.clone();
                    if !dual {
                        payload.dual_shape = None;
                        payload.zero_components.clear();
                    }
                    for iterate in [false, true] {
                        let (mut evaluator, length) = build_evaluator::<f64, _>(
                            &payload,
                            &params,
                            Vec::new(),
                            &StateMap::default(),
                            iterate,
                            None,
                        )
                        .map(|(evaluator, _)| evaluator)?;
                        let mut output = vec![Complex::new(0.0, 0.0); length];
                        let input = if dual {
                            vec![2.0, 1.0, 0.0, 0.0]
                        } else {
                            vec![2.0, 0.0]
                        };
                        let input = input
                            .into_iter()
                            .map(|value| Complex::new(value, 0.0))
                            .collect::<Vec<_>>();
                        assert_eq!(
                            evaluator.get_input_len(),
                            input.len(),
                            "bound coefficient columns must not survive loading"
                        );
                        evaluator.evaluate(&input, &mut output);
                        let expected = if dual {
                            vec![5.0, 5.0, 8.0, 4.0]
                        } else {
                            vec![5.0, 8.0]
                        };
                        assert_eq!(
                            output,
                            expected
                                .into_iter()
                                .map(|value| Complex::new(value, 0.0))
                                .collect::<Vec<_>>()
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn standalone_template_cpp_complex_constants_match_eager_and_assembly() -> Result<()> {
        let root = env::temp_dir().join(format!(
            "gammaloop_standalone_template_complex_export_{}",
            process::id()
        ));
        let x = parse!("standalone_template_complex_export::x");
        let q = parse!("standalone_template_complex_export::q");
        for imaginary in [true, false] {
            for retained in [false, true] {
                let coefficient = Atom::num(Rational::new(1, 8))
                    * if imaginary {
                        Atom::one() + Atom::i()
                    } else {
                        Atom::one()
                    };
                let offset = Atom::num(Rational::new(1, 16))
                    * if imaginary { Atom::i() } else { Atom::one() };
                let mut functions = FunctionMap::new();
                let product = if retained {
                    let inner = symbol!("standalone_template_complex_export::inner");
                    let outer = symbol!("standalone_template_complex_export::outer");
                    for (function, body) in [
                        (inner, &coefficient * &q),
                        (outer, inner.call_args([q.clone()])),
                    ] {
                        functions.add_function_with_options(
                            function,
                            vec![Indeterminate::try_from(q.clone()).map_err(|error| eyre!(error))?],
                            body,
                            FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
                        )?;
                    }
                    outer.call_args([x.clone()])
                } else {
                    &coefficient * &x
                };
                let expressions = [
                    product + &x * Atom::num(Rational::new(1, 2)) + offset.clone(),
                    offset,
                    Atom::num(Rational::new(3, 8)),
                ];
                let mut trees = expressions.iter().map(|expression| {
                    expression
                        .evaluator(slice::from_ref(&x))
                        .function_map(functions.clone())
                        .optimization_settings(OptimizationSettings::new().cores(1))
                        .build()
                        .map_err(|error| eyre!(error))
                });
                let mut tree = trees.next().unwrap()?;
                for next in trees {
                    tree.merge(next?, Some(1)).map_err(|error| eyre!(error))?;
                }
                let mut evaluator =
                    tree.map_coeff(&|value| Complex::new(value.re.to_f64(), value.im.to_f64()));
                let label = format!("complex_{imaginary}_retained_{retained}");
                for backend in [
                    StandaloneBackend::Eager,
                    StandaloneBackend::Cpp,
                    StandaloneBackend::Assembly,
                ] {
                    let mut runtime =
                        StandaloneRuntimeEvaluator::build(&mut evaluator, backend, &root, &label)?;
                    for input in [2.0, -2.0] {
                        let mut output = [Complex::new(0.0, 0.0); 3];
                        runtime.evaluate(&[Complex::new(input, 0.0)], &mut output);
                        assert_eq!(
                            output,
                            [
                                Complex::new(
                                    input * 0.625 + if imaginary { 0.0 } else { 0.0625 },
                                    if imaginary {
                                        input * 0.125 + 0.0625
                                    } else {
                                        0.0
                                    },
                                ),
                                if imaginary {
                                    Complex::new(0.0, 0.0625)
                                } else {
                                    Complex::new(0.0625, 0.0)
                                },
                                Complex::new(0.375, 0.0),
                            ]
                        );
                    }
                    if matches!(backend, StandaloneBackend::Cpp) {
                        let path = root.join(format!("standalone_{label}_c__.cpp"));
                        let source = fs::read_to_string(&path)?;
                        assert!(source.contains("T(3.75e-1)"));
                        assert!(!path.with_extension("symbolica.cpp").exists());
                    }
                }
            }
        }
        fs::remove_dir_all(root)?;
        Ok(())
    }

    #[test]
    fn standalone_symjit_reserves_retained_call_capacity_after_roundtrip() -> Result<()> {
        for arity in [65, 2049] {
            let formals = (0..arity)
                .map(|i| format!("standalone_wide::q{i}"))
                .collect::<Vec<_>>();
            let function = "standalone_wide::function";
            let entries: Vec<SerializedFnMapEntry<String>> = vec![(
                format!("{function}({})", formals.join(",")),
                formals
                    .iter()
                    .map(|q| format!("{q}^2"))
                    .collect::<Vec<_>>()
                    .join("+"),
                vec![],
                formals,
                InliningPolicy::Never,
                false,
            )];
            let json = serde_json::to_vec(&entries)?;
            let parsed: Vec<SerializedFnMapEntry<String>> = serde_json::from_slice(&json)?;
            let bytes = bincode::encode_to_vec(&parsed, config::standard())?;
            let (restored, consumed): (Vec<SerializedFnMapEntry<String>>, _) =
                bincode::decode_from_slice(&bytes, config::standard())?;
            assert_eq!(consumed, bytes.len());
            assert_eq!(restored, entries);
            let parsed = parse_fn_map_entries(&restored, &StateMap::default())?;
            let (_, fn_map) = apply_fn_map_entries(parsed)?;
            let parameter = parse_lit!(standalone_wide::x);
            let root = symbol!(function).call_args((0..arity).map(|i| &parameter + i));
            let mut evaluator = root
                .evaluator(slice::from_ref(&parameter))
                .function_map(fn_map)
                .build()
                .map_err(|error| eyre!(error))?
                .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
            let result = StandaloneRuntimeEvaluator::build(
                &mut evaluator,
                StandaloneBackend::Symjit,
                Path::new("."),
                "wide_call_capacity",
            );
            if arity > 2048 {
                let error = match result {
                    Err(error) => error,
                    Ok(_) => panic!("oversized call must be rejected before invoking SymJIT"),
                };
                assert!(
                    error.to_string().contains("call arguments")
                        && error.to_string().contains("2048"),
                    "{error}"
                );
                continue;
            }
            let StandaloneRuntimeEvaluator::Symjit(compiled) = result? else {
                unreachable!()
            };
            let json = serde_json::to_vec(&compiled)?;
            let mut restored_json: JITCompiledEvaluator<Complex<f64>> =
                serde_json::from_slice(&json)?;
            let saved = bincode::encode_to_vec(&compiled, config::standard())?;
            let (mut restored, consumed): (JITCompiledEvaluator<Complex<f64>>, _) =
                bincode::decode_from_slice(&saved, config::standard())?;
            assert_eq!(consumed, saved.len());
            let input = [Complex::new(2.0, 1.0)];
            let mut output = [Complex::new(0.0, 0.0)];
            restored.evaluate(&input, &mut output);
            let expected = (0..arity).fold(Complex::new(0.0, 0.0), |sum, i| {
                let a = (i + 2) as f64;
                sum + Complex::new(a * a - 1.0, 2.0 * a)
            });
            assert_eq!(output[0], expected);
            restored_json.evaluate(&input, &mut output);
            assert_eq!(output[0], expected);
            assert_eq!(evaluator.evaluate_single(&input), expected);
        }
        Ok(())
    }

    #[test]
    fn standalone_function_policy_roundtrip_preserves_shared_bodies_and_alias_scope() -> Result<()>
    {
        let x = parse!("standalone_policy::x");
        let global = parse!("standalone_policy::global");
        let alias = parse!("standalone_policy::alias(0)");
        let function = symbol!("standalone_policy::function");
        let print = |atom: &Atom| atom.printer(PrintOptions::file()).to_string();
        for tags in [Vec::new(), vec![Atom::num(7)]] {
            let call = |argument: Atom| {
                FunctionBuilder::new(function)
                    .add_args(&tags)
                    .add_arg(argument)
                    .finish()
            };
            for inlining in [InliningPolicy::Always, InliningPolicy::Never] {
                let entries: Vec<SerializedFnMapEntry<String>> = vec![
                    (
                        print(&alias),
                        print(&x.pow(2)),
                        vec!["0".into()],
                        vec![],
                        InliningPolicy::Always,
                        true,
                    ),
                    (
                        print(&call(x.clone())),
                        print(&(&alias + &global)),
                        tags.iter().map(print).collect(),
                        vec![print(&x)],
                        inlining,
                        false,
                    ),
                ];
                let json = serde_json::to_vec(&entries)?;
                let parsed: Vec<SerializedFnMapEntry<String>> = serde_json::from_slice(&json)?;
                let bytes = bincode::encode_to_vec(&parsed, config::standard())?;
                let (restored, consumed): (Vec<SerializedFnMapEntry<String>>, _) =
                    bincode::decode_from_slice(&bytes, config::standard())?;
                assert_eq!(consumed, bytes.len());
                assert_eq!(restored, entries);
                let mut parsed = parse_fn_map_entries(&restored, &StateMap::default())?;
                parsed.extend(parsed.clone()); // Graph and evaluator metadata overlap.
                let (replacements, fn_map) = apply_fn_map_entries(parsed)?;
                let root = call(Atom::num(2)) + call(Atom::num(3));
                assert_eq!(root.replace_multiple(&replacements), root);

                let mut evaluator = root
                    .replace_multiple(&replacements)
                    .evaluator(&[global.clone(), x.clone()])
                    .function_map(fn_map)
                    .horner_iterations(1)
                    .build()
                    .map_err(|error| eyre!(error))?
                    .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
                assert_eq!(
                    evaluator.export_instructions().sub_evaluators.len(),
                    usize::from(inlining == InliningPolicy::Never)
                );
                assert_eq!(
                    evaluator.evaluate_single(&[Complex::new(10.0, 0.0), Complex::new(11.0, 0.0)]),
                    Complex::new(33.0, 0.0)
                );
            }
        }
        Ok(())
    }

    #[test]
    fn standalone_template_compiled_artifacts_keep_representation_and_method() -> Result<()> {
        let stack = |value: i64| {
            let payload = StandaloneGenericEvaluatorArchive {
                parameter_rows: None,
                exprs: vec![format!(
                    "{value}*y^2*standalone_template_joint_numerator(x+c)"
                )],
                additional_fn_map_entries: vec![(
                    "standalone_template_joint_numerator(u)".to_owned(),
                    "u^2".to_owned(),
                    Vec::new(),
                    vec!["u".to_owned()],
                    InliningPolicy::Never,
                    false,
                )],
                dual_shape: Some(vec![vec![0], vec![1]]),
                zero_components: vec![(1, 1)],
            };
            let mut summed = payload.clone();
            summed.exprs = vec![format!(
                "{value}*y^2*(standalone_template_joint_numerator(x+1)+standalone_template_joint_numerator(x+3))"
            )];
            StandaloneEvaluatorStackArchive {
                explicit_orientation_sum_only: true,
                production_orientation_ids: vec![0],
                parametric_rows: Some(
                    ParameterRows::new(
                        vec!["c".to_owned()],
                        vec![vec![Rational::from(1)], vec![Rational::from(3)]],
                    )
                    .unwrap(),
                ),
                single_parametric: payload,
                iterative: None,
                summed_function_map: Some(summed.clone()),
                summed: Some(summed),
                representative_input: vec![
                    StandaloneComplexInput {
                        re: "2".to_owned(),
                        im: "0".to_owned(),
                    },
                    StandaloneComplexInput {
                        re: "1".to_owned(),
                        im: "0".to_owned(),
                    },
                    StandaloneComplexInput {
                        re: "1".to_owned(),
                        im: "0".to_owned(),
                    },
                    StandaloneComplexInput {
                        re: "0".to_owned(),
                        im: "0".to_owned(),
                    },
                ],
                start: 0,
                residue_map_id_start: 0,
                mult_offset: 2,
            }
        };
        let archive = StandaloneEvaluatorArchive {
            version: STANDALONE_EVALUATORS_VERSION,
            three_d_representations: vec!["cff".to_owned(), "ltd".to_owned()],
            numeric_target: StandaloneNumericTarget::Double,
            symbolica_state: (),
            graph_terms: vec![StandaloneGraphTermArchive {
                graph_name: "artifact_dispatch".to_owned(),
                orientations: Vec::new(),
                param_builder_params: vec!["x".to_owned(), "y".to_owned()],
                fn_map_entries: Vec::new(),
                original_integrand: BTreeMap::from([
                    ("cff".to_owned(), stack(3)),
                    ("ltd".to_owned(), stack(7)),
                ]),
                threshold_counterterms: Vec::new(),
                threshold_counterterms_are_variants: false,
                threshold_variants: Vec::new(),
                threshold_multipliers: None,
                metadata_registry: None,
            }],
        };
        archive.validate()?;
        for corruption in 0..6 {
            let mut malformed = archive.clone();
            let stack = malformed.graph_terms[0]
                .original_integrand
                .get_mut("cff")
                .unwrap();
            match corruption {
                0 => stack.mult_offset = 1,
                1 => {
                    stack.representative_input.pop();
                }
                2 => stack.summed_function_map.as_mut().unwrap().dual_shape = None,
                3 => stack.residue_map_id_start = 3,
                4 => stack.summed.as_mut().unwrap().dual_shape = None,
                5 => stack.single_parametric.dual_shape = Some(Vec::new()),
                _ => unreachable!(),
            }
            let bytes = bincode::encode_to_vec(&malformed, config::standard())?;
            let (binary, _): (StandaloneEvaluatorArchive<(), String>, _) =
                bincode::decode_from_slice(&bytes, config::standard())?;
            let json: StandaloneEvaluatorArchive<(), String> =
                serde_json::from_slice(&serde_json::to_vec(&malformed)?)?;
            for restored in [binary, json] {
                assert!(
                    restored
                        .validate()
                        .unwrap_err()
                        .to_string()
                        .contains("physical input or dual layout")
                );
            }
        }
        let artifact_dir = env::temp_dir().join(format!(
            "gammaloop-standalone-artifact-dispatch-{}",
            process::id()
        ));
        let state_map = current_state_map()?;
        for (representation, coefficient) in [("cff", 3.0), ("ltd", 7.0)] {
            let stack = &archive.graph_terms[0].original_integrand[representation];
            for method in [
                StandaloneMethod::SingleParametric,
                StandaloneMethod::SummedFunctionMap,
            ] {
                let (payload, iterate) = stack.selected_payload(method)?;
                let duplicate_error = build_evaluator::<f64, _>(
                    payload,
                    &[parse_lit!(x), parse_lit!(x)],
                    Vec::new(),
                    &state_map,
                    iterate,
                    None,
                )
                .map(|(evaluator, _)| evaluator)
                .unwrap_err();
                assert!(
                    duplicate_error
                        .to_string()
                        .contains("parameters must be distinct")
                );
                let params = vec![parse_lit!(x), parse_lit!(y)];
                let rows = if method == StandaloneMethod::SingleParametric {
                    stack
                        .parametric_rows
                        .clone()
                        .map(|rows| {
                            rows.try_map_parameters(|parameter| {
                                parameter.import_with_map(&state_map)
                            })
                        })
                        .transpose()?
                } else {
                    None
                };
                let ((mut evaluator, _), mut shared) = build_evaluator::<f64, _>(
                    payload,
                    &params,
                    Vec::new(),
                    &state_map,
                    iterate,
                    rows.as_ref(),
                )?;
                if method == StandaloneMethod::SingleParametric {
                    // y^2 is shared across every residue row and has a zero
                    // derivative seed, while the numerator still depends on c.
                    assert!(shared.as_ref().unwrap().1.get_output_len() > 0);
                }
                for (x, expected) in [(2.0, [34.0, 16.0]), (3.0, [52.0, 20.0])] {
                    let inputs = vec![vec![
                        Complex::new(x, 0.0),
                        Complex::new(1.0, 0.0),
                        Complex::new(1.0, 0.0),
                        Complex::new(0.0, 0.0),
                    ]];
                    for backend in [
                        StandaloneBackend::Eager,
                        StandaloneBackend::Symjit,
                        StandaloneBackend::Cpp,
                        StandaloneBackend::Assembly,
                    ] {
                        assert_eq!(
                            evaluate_with_backend_f64(
                                &mut evaluator,
                                shared.as_mut(),
                                backend,
                                &inputs,
                                &artifact_dir,
                                &format!("row_values_{representation}_{}_{x}", method.as_str())
                            )?,
                            expected.map(|value| Complex::new(coefficient * value, 0.0))
                        );
                    }
                }
            }
        }

        for representation in [None, Some("ltd".to_owned())] {
            for method in [StandaloneMethod::SingleParametric, StandaloneMethod::Summed] {
                evaluate_double_archive(
                    archive.clone(),
                    &state_map,
                    &StandaloneCliOptions {
                        representation: representation.clone(),
                        method,
                        compare_backends: vec![StandaloneBackend::Cpp, StandaloneBackend::Assembly],
                        artifact_dir: Some(artifact_dir.clone()),
                        ..Default::default()
                    },
                    None,
                )?;
            }
        }
        // All files must coexist after both representations and methods have
        // run: a shared basename would silently overwrite an earlier program.
        let artifacts = artifact_dir.join("artifact_dispatch");
        assert_eq!(fs::read_dir(&artifacts)?.count(), 24);
        for representation in ["cff", "ltd"] {
            for method in ["single_parametric", "summed"] {
                for backend in ["c__", "assembly"] {
                    for extension in ["cpp", "so"] {
                        assert!(artifacts.join(format!(
                            "standalone_artifact_dispatch_{representation}_original_{method}_{backend}.{extension}"
                        )).is_file());
                    }
                }
            }
        }
        for representation in ["cff", "ltd"] {
            let label = format!(
                "standalone_artifact_dispatch_{representation}_original_single_parametric_auxiliary"
            );
            for backend in ["c__", "assembly"] {
                for extension in ["cpp", "so"] {
                    assert!(
                        artifacts
                            .join(format!("{label}_{backend}.{extension}"))
                            .is_file()
                    );
                }
            }
        }
        fs::remove_dir_all(artifact_dir)?;
        Ok(())
    }

    #[test]
    fn standalone_template_joint_rows_keep_exact_constants_and_zero_dual_seeds() -> Result<()> {
        let precise = Rational::from((9_007_199_254_740_993_i64, 9_007_199_254_740_992_i64));
        let catalog = ParameterRows::new(
            vec!["c".to_owned(), "q".to_owned(), "active".to_owned()],
            vec![
                vec![precise.clone(), precise.clone(), 1.into()],
                vec![(-1).into(), precise.clone(), 1.into()],
                vec![7.into(), precise, 0.into()],
            ],
        )?;
        let catalog: ParameterRows<String> =
            serde_json::from_slice(&serde_json::to_vec(&catalog)?)?;
        let state = current_state_map()?;
        let rows = catalog.try_map_parameters(|parameter| parameter.import_with_map(&state))?;
        let payload = StandaloneGenericEvaluatorArchive {
            parameter_rows: None,
            exprs: vec!["if(active,c*x*(q-1),0)".to_owned()],
            additional_fn_map_entries: Vec::new(),
            dual_shape: Some(vec![vec![0], vec![1]]),
            zero_components: Vec::new(),
        };
        let ((mut evaluator, output_len), mut shared) = build_evaluator::<Float, _>(
            &payload,
            &[parse_lit!(x)],
            Vec::new(),
            &state,
            false,
            Some(&rows),
        )?;
        let inputs = vec![vec![
            Complex::new(Float::parse_standalone_input("3")?, Float::zero_value()),
            Complex::new(Float::one_value(), Float::zero_value()),
        ]];
        let result = evaluate_eager(&mut evaluator, shared.as_mut(), output_len, &inputs)?;
        // Both exact coefficient roles matter: rounding either c or q through
        // f64 erases this value and derivative. The third row remains inactive.
        assert_eq!(
            result,
            [3_i64, 1].map(|numerator| Complex::new(
                Float::exact_from_rational(
                    &(Rational::from((numerator, 9_007_199_254_740_992_i64))
                        / Rational::from(9_007_199_254_740_992_i64))
                ),
                Float::zero_value(),
            ))
        );
        Ok(())
    }

    #[test]
    fn standalone_template_arb_external_constants_use_declared_precision() -> Result<()> {
        let x = parse_lit!(standalone_arb_constants::x);
        let coefficient = parse_lit!(standalone_arb_constants::c);
        let print = |atom: &Atom| atom.printer(PrintOptions::file()).to_string();
        let rows = ParameterRows::new(
            vec![coefficient.clone()],
            vec![vec![1.into()], vec![2.into()]],
        )?;
        let state = current_state_map()?;
        let inputs = vec![vec![Complex::new(
            Float::parse_standalone_input("1.25")?,
            Float::zero_value(),
        )]];
        let tolerance = Float::parse_standalone_input("1e-250")?;
        for staged in [false, true] {
            // The reciprocal stays one row-independent operation, with π in
            // its input, instead of flattening into the row-dependent product.
            let mut root = (Atom::var(Symbol::PI) + &x).pow(-1);
            if staged {
                root *= &coefficient;
            }
            let payload = StandaloneGenericEvaluatorArchive {
                parameter_rows: None,
                exprs: vec![print(&root)],
                additional_fn_map_entries: Vec::new(),
                dual_shape: None,
                zero_components: Vec::new(),
            };
            let ((mut evaluator, output_len), mut shared) = build_evaluator::<Float, _>(
                &payload,
                slice::from_ref(&x),
                Vec::new(),
                &state,
                false,
                staged.then_some(&rows),
            )?;
            if staged {
                assert!(shared.as_ref().unwrap().1.get_output_len() > 0);
            }
            let result = evaluate_eager(&mut evaluator, shared.as_mut(), output_len, &inputs)?;
            let scale = if staged { "3" } else { "1" };
            let expected = Float::parse_standalone_input(scale)?
                / (Float::new(ARB_PRECISION_BITS).pi() + Float::parse_standalone_input("1.25")?);
            assert_eq!(result.len(), 1);
            assert_eq!(result[0].re.get_precision(), ARB_PRECISION_BITS);
            assert!((result[0].re.clone() - expected).norm() < tolerance);
            assert_eq!(result[0].im, Float::zero_value());
        }
        Ok(())
    }

    #[test]
    fn standalone_template_symjit_preserves_independent_guards() -> Result<()> {
        let q = parse_lit!(standalone_template_guards::q);
        let s = parse_lit!(standalone_template_guards::s);
        let t = parse_lit!(standalone_template_guards::t);
        let power = q.pow(parse_lit!(-1 / 2));
        let expressions = [s.clone(), t.clone()]
            .map(|guard| Symbol::IF.call_args([guard, power.clone(), Atom::Zero]));
        let parameters = [q, s, t];
        let mut tree = expressions[0]
            .evaluator(&parameters)
            .optimization_settings(OptimizationSettings::new().cores(1))
            .build()
            .map_err(|error| eyre!(error))?;
        let second = expressions[1]
            .evaluator(&parameters)
            .optimization_settings(OptimizationSettings::new().cores(1))
            .build()
            .map_err(|error| eyre!(error))?;
        tree.merge(second, Some(2)).map_err(|error| eyre!(error))?;
        let mut evaluator = tree.map_coeff(&|coefficient| {
            Complex::new(coefficient.re.to_f64(), coefficient.im.to_f64())
        });
        let mut runtime = StandaloneRuntimeEvaluator::build(
            &mut evaluator,
            StandaloneBackend::Symjit,
            &env::temp_dir(),
            "independent_guards",
        )?;
        for guards in [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]] {
            let mut result = [Complex::new(0.0, 0.0); 2];
            runtime.evaluate(
                &[
                    Complex::new(4.0, 0.0),
                    Complex::new(guards[0], 0.0),
                    Complex::new(guards[1], 0.0),
                ],
                &mut result,
            );
            assert_eq!(result, guards.map(|guard| Complex::new(guard * 0.5, 0.0)));
        }
        Ok(())
    }

    #[test]
    fn standalone_template_deduplicates_shared_function_definitions() -> Result<()> {
        let x = parse_lit!(standalone_template_duplicate::x);
        let call = function!(symbol!("standalone_template_duplicate::f"), &x);
        let entry = (
            call.clone(),
            x.clone().pow(2),
            Vec::new(),
            vec![Indeterminate::try_from(x.clone()).unwrap()],
            InliningPolicy::Never,
            false,
        );
        let (_, function_map) = apply_fn_map_entries(vec![entry.clone(), entry])?;
        let mut evaluator = call
            .evaluator(&[x])
            .function_map(function_map)
            .optimization_settings(OptimizationSettings::new().cores(1))
            .build()
            .map_err(|error| eyre!(error))?
            .map_coeff(&|coefficient| coefficient.re.to_f64());
        let mut result = [0.0];
        evaluator.evaluate(&[3.0], &mut result);
        assert_eq!(result, [9.0]);
        Ok(())
    }

    #[test]
    fn standalone_template_binary_load_restores_archived_symbol_map() -> Result<()> {
        let atom = parse_lit!(standalone_template_binary::archived_parameter);
        let mut encoded_atom = Vec::new();
        atom.as_view().write(&mut encoded_atom)?;
        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)?;
        let archive = StandaloneEvaluatorArchive::<Vec<u8>, Vec<u8>> {
            version: STANDALONE_EVALUATORS_VERSION,
            three_d_representations: vec!["cff".to_owned(), "ltd".to_owned()],
            numeric_target: StandaloneNumericTarget::Quad,
            symbolica_state,
            graph_terms: Vec::new(),
        };
        let path = env::temp_dir().join(format!(
            "gammaloop-amplitude-template-symbol-map-{}.bin",
            process::id()
        ));
        fs::write(&path, bincode::encode_to_vec(&archive, config::standard())?)?;
        let (loaded, state_map) = load_bin(&path)?;
        assert_eq!(loaded.three_d_representations, ["cff", "ltd"]);
        assert_eq!(encoded_atom.import_with_map(&state_map)?, atom);
        fs::remove_file(path)?;
        Ok(())
    }
}
