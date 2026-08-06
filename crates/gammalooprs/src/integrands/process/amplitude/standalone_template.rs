#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "0441bd7a511209dce2ca99925fe87f8b18e4bf03", default-features = false, features = ["bincode", "gmp", "native_code_generation", "serde"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "0441bd7a511209dce2ca99925fe87f8b18e4bf03" }
//! graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "0441bd7a511209dce2ca99925fe87f8b18e4bf03" }
//! ```

#![allow(dead_code)]

use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Cursor,
    ops::Neg,
    path::{Path, PathBuf},
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    domains::rational::Fraction, evaluate::JITCompiledEvaluator, prelude::*, state::StateMap,
};

const STANDALONE_EVALUATORS_VERSION: u32 = 8;
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
    original_integrand: StandaloneEvaluatorStackArchive<A>,
    threshold_counterterms: Vec<Vec<StandaloneIndexedEvaluatorStackArchive<A>>>,
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
    single_parametric: StandaloneGenericEvaluatorArchive<A>,
    iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    representative_input: Vec<StandaloneComplexInput>,
    start: usize,
    mult_offset: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneGenericEvaluatorArchive<A = Vec<u8>> {
    exprs: Vec<A>,
    additional_fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    dual_shape: Option<Vec<Vec<usize>>>,
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
    symmetrize: bool,
    opaque_derivatives: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode, Serialize, Deserialize)]
struct ThresholdCountertermVariantMetadata {
    variant_id: usize,
    name: String,
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

type SerializedFnMapEntry<A> = (A, A, Vec<A>, Vec<A>);
type ParsedFnMapEntry = (Atom, Atom, Vec<Atom>, Vec<Indeterminate>);

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
    Symjit(JITCompiledEvaluator<Complex<f64>>),
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
                let compiled = evaluator
                    .export_cpp::<Complex<f64>>(
                        &source_path,
                        &function_name,
                        ExportSettings::new()
                            .include_header(true)
                            .inline_asm(backend.inline_asm())
                            .custom_header(None),
                    )
                    .map_err(|error| eyre!(error))?
                    .compile(&library_path, CompileOptions::default())
                    .map_err(|error| eyre!(error))?
                    .load()
                    .map_err(|error| eyre!(error))?;
                Ok(Self::Compiled(compiled))
            }
            StandaloneBackend::Symjit => Ok(Self::Symjit(
                evaluator
                    // SymJIT 2.21 supports optimization levels up to O2 and cannot compact some
                    // complex temporary layouts.
                    .jit_compile(
                        JITCompilationSettings::new()
                            .optimization_level(2)
                            .with_option("compact", "false"),
                    )
                    .map_err(|error| eyre!(error))?,
            )),
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
        .map(|(lhs, rhs, tags, args)| {
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

            Ok((lhs_atom, rhs_atom, tags, args))
        })
        .collect()
}

fn apply_fn_map_entries(
    parsed_entries: Vec<ParsedFnMapEntry>,
) -> Result<(Vec<Replacement>, FunctionMap)> {
    let mut fn_map = FunctionMap::new();
    let mut replacements: Vec<Replacement> = vec![];
    fn_map
        .add_aliases([(parse_lit!(gammalooprs::x), Atom::Zero)])
        .map_err(|error| eyre!(error))?;

    for (lhs, rhs, tags, args) in parsed_entries {
        if let AtomView::Var(_) = lhs.as_view() {
            if let Ok(constant) = Complex::<Rational>::try_from(rhs.as_view()) {
                fn_map
                    .add_aliases([(lhs.clone(), Atom::num(constant))])
                    .map_err(|error| eyre!(error))?;
            } else {
                replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            }
        } else if let AtomView::Fun(function) = lhs.as_view() {
            if tags.is_empty() {
                let mut wildcards = Vec::new();
                for (index, arg) in args.iter().enumerate() {
                    let atom: Atom = arg.clone().into();
                    wildcards.push(
                        Replacement::new(
                            atom.to_pattern(),
                            Atom::var(symbol!(format!("x{index}_"))),
                        )
                        .allow_new_wildcards_on_rhs(true),
                    );
                }

                fn_map
                    .add_function(function.get_symbol(), args, rhs.clone())
                    .map_err(|error| eyre!(error))?;

                replacements.push(Replacement::new(
                    lhs.replace_multiple(&wildcards).to_pattern(),
                    rhs.replace_multiple(&wildcards),
                ));
            } else {
                fn_map
                    .add_tagged_function(function.get_symbol(), tags, args, rhs.clone())
                    .map_err(|error| eyre!(error))?;
            }
        } else {
            replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
        }
    }

    Ok((replacements, fn_map))
}

fn build_evaluator<T: StandaloneNumber, A: ImportWithMap>(
    payload: &StandaloneGenericEvaluatorArchive<A>,
    params: &[Atom],
    mut fn_map_entries: Vec<ParsedFnMapEntry>,
    state_map: &StateMap,
    iterate: bool,
) -> Result<(ExpressionEvaluator<Complex<T>>, usize)>
where
    Complex<T>: EvaluationDomain,
{
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

    if iterate {
        let mut tree: Option<ExpressionEvaluator<Complex<Fraction<IntegerRing>>>> = None;

        for expr in &exprs {
            let eval = expr
                .replace_multiple(&replacements)
                .evaluator(params)
                .function_map(fn_map.clone())
                .optimization_settings(optimization_settings.clone())
                .build()
                .map_err(|error| {
                    eyre!(
                        "{error} while building iterative evaluator for {}",
                        expr.replace_multiple(&replacements)
                    )
                })?;

            tree = Some(if let Some(mut merged) = tree {
                merged
                    .merge(eval, cpe_iterations)
                    .map_err(|error| eyre!(error))?;
                merged
            } else {
                eval
            });
        }

        tree.map(|eval| {
            (
                eval.map_coeff(&|value| {
                    Complex::new(
                        T::exact_from_rational(&value.re),
                        T::exact_from_rational(&value.im),
                    )
                }),
                exprs.len(),
            )
        })
        .ok_or_else(|| eyre!("No expressions in evaluator payload"))
    } else {
        let replaced_exprs = exprs
            .iter()
            .map(|expr| expr.replace_multiple(&replacements))
            .collect::<Vec<_>>();

        Atom::evaluator_multiple(&replaced_exprs, params)
            .function_map(fn_map)
            .optimization_settings(optimization_settings)
            .build()
            .map(|eval| {
                (
                    eval.map_coeff(&|value| {
                        Complex::new(
                            T::exact_from_rational(&value.re),
                            T::exact_from_rational(&value.im),
                        )
                    }),
                    exprs.len(),
                )
            })
            .map_err(|error| eyre!("{error}"))
    }
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
            if evaluator.exprs.len() != 1
                || evaluator.dual_shape.is_some()
                || !evaluator.additional_fn_map_entries.is_empty()
            {
                return Err(eyre!(
                    "threshold-multiplier evaluator {index} must contain one scalar, non-dual expression without function-map entries"
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
        for graph in &self.graph_terms {
            graph.validate_threshold_payload()?;
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

    fn set_orientation<T: StandaloneNumber>(&self, orientation: &[i8]) -> Result<Vec<Complex<T>>> {
        let mut input = self.representative_input::<T>()?;
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
        selection: &StandaloneStackSelection,
    ) -> Result<&StandaloneEvaluatorStackArchive<A>> {
        match selection {
            StandaloneStackSelection::Original => Ok(&self.original_integrand),
            StandaloneStackSelection::ThresholdCounterterm((first, second)) => self
                .threshold_counterterms
                .get(*first)
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
    output_len: usize,
    inputs: &[Vec<Complex<T>>],
) -> Vec<Complex<T>> {
    let mut accumulated = vec![Complex::new(T::zero_value(), T::zero_value()); output_len];

    for input in inputs {
        let mut current = vec![Complex::new(T::zero_value(), T::zero_value()); output_len];
        evaluator.evaluate(input, &mut current);
        for (accumulated_value, current_value) in accumulated.iter_mut().zip(current) {
            *accumulated_value += current_value;
        }
    }

    accumulated
}

fn evaluate_with_backend_f64(
    evaluator: &mut ExpressionEvaluator<Complex<f64>>,
    backend: StandaloneBackend,
    output_len: usize,
    inputs: &[Vec<Complex<f64>>],
    artifact_root: &Path,
    label: &str,
) -> Result<Vec<Complex<f64>>> {
    let mut runtime = StandaloneRuntimeEvaluator::build(evaluator, backend, artifact_root, label)?;
    let mut accumulated = vec![Complex::new(0.0, 0.0); output_len];

    for input in inputs {
        let mut current = vec![Complex::new(0.0, 0.0); output_len];
        runtime.evaluate(input, &mut current);
        for (accumulated_value, current_value) in accumulated.iter_mut().zip(current) {
            *accumulated_value += current_value;
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
    let mut args = std::env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--help" | "-h" => {
                let program = std::env::args()
                    .next()
                    .unwrap_or_else(|| "standalone_evaluators_rust.rs".to_string());
                print_usage(&program);
                std::process::exit(0);
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

fn load_bin(path: impl AsRef<Path>) -> Result<StandaloneEvaluatorArchive> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (archive, _): (StandaloneEvaluatorArchive, _) =
        bincode::decode_from_slice(&binary, bincode::config::standard())?;
    archive.validate()?;
    Ok(archive)
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

    let graph = archive.graph_term(options.graph_index, options.graph_name.as_deref())?;
    let stack = graph.stack(&options.stack)?;
    let (payload, iterate) = stack.selected_payload(options.method)?;

    let params = graph
        .param_builder_params
        .iter()
        .map(|param| param.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;
    let fn_map_entries = parse_fn_map_entries(&graph.fn_map_entries, state_map)?;
    let (mut evaluator, output_len) =
        build_evaluator::<f64, _>(payload, &params, fn_map_entries, state_map, iterate)?;

    let inputs = if let Some(custom_input) = custom_input {
        vec![
            custom_input
                .iter()
                .map(StandaloneComplexInput::parse::<f64>)
                .collect::<Result<Vec<_>>>()?,
        ]
    } else {
        match options.method {
            StandaloneMethod::SingleParametric => {
                if let Some(index) = options.orientation_index {
                    let orientation = graph.orientations.get(index).ok_or_else(|| {
                        eyre!(
                            "Orientation index {} is out of range for {} orientations",
                            index,
                            graph.orientations.len()
                        )
                    })?;
                    vec![stack.set_orientation::<f64>(orientation)?]
                } else {
                    graph
                        .orientations
                        .iter()
                        .map(|orientation| stack.set_orientation::<f64>(orientation))
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
            backend,
            output_len,
            &inputs,
            &artifact_root.join(sanitize_label(&graph.graph_name)),
            &format!("{}_{}", graph.graph_name, options.stack.label()),
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

    let graph = archive.graph_term(options.graph_index, options.graph_name.as_deref())?;
    let stack = graph.stack(&options.stack)?;
    let (payload, iterate) = stack.selected_payload(options.method)?;
    let params = graph
        .param_builder_params
        .iter()
        .map(|param| param.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;
    let fn_map_entries = parse_fn_map_entries(&graph.fn_map_entries, state_map)?;
    let (mut evaluator, output_len) =
        build_evaluator::<T, _>(payload, &params, fn_map_entries, state_map, iterate)?;

    let inputs = if let Some(custom_input) = custom_input {
        vec![
            custom_input
                .iter()
                .map(StandaloneComplexInput::parse::<T>)
                .collect::<Result<Vec<_>>>()?,
        ]
    } else {
        match options.method {
            StandaloneMethod::SingleParametric => {
                if let Some(index) = options.orientation_index {
                    let orientation = graph.orientations.get(index).ok_or_else(|| {
                        eyre!(
                            "Orientation index {} is out of range for {} orientations",
                            index,
                            graph.orientations.len()
                        )
                    })?;
                    vec![stack.set_orientation::<T>(orientation)?]
                } else {
                    graph
                        .orientations
                        .iter()
                        .map(|orientation| stack.set_orientation::<T>(orientation))
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

    let values = evaluate_eager(&mut evaluator, output_len, &inputs);
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
    let state_map = current_state_map()?;
    let custom_input = options
        .input_json
        .as_ref()
        .map(load_custom_input)
        .transpose()?;

    match extension.to_string_lossy().as_ref() {
        "bin" => {
            let archive = load_bin(&input)?;
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
