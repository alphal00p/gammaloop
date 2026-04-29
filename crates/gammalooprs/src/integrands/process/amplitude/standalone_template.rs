#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { git = "https://github.com/benruijl/symbolica", rev = "aa51532febc5f88e5b3443123a075523faa2883e", default-features = false, features = ["bincode", "serde"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/benruijl/symbolica", rev = "aa51532febc5f88e5b3443123a075523faa2883e" }
//! graphica = { git = "https://github.com/benruijl/symbolica", rev = "aa51532febc5f88e5b3443123a075523faa2883e" }
//! ```

#![allow(dead_code)]

use std::{
    fs,
    io::Cursor,
    ops::Neg,
    path::{Path, PathBuf},
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate},
    domains::{
        float::{
            Complex, DoubleFloat, Float as ArbFloat, FloatLike as SymbolicaFloatLike, Real,
            RealLike, SingleFloat,
        },
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    evaluate::{
        CompileOptions, CompiledComplexEvaluator, ExportSettings, ExpressionEvaluator, FunctionMap,
        InlineASM, JITCompiledEvaluator, OptimizationSettings,
    },
    id::{MatchSettings, Replacement},
    parse_lit,
    state::{State, StateMap},
    symbol, try_parse,
};

const STANDALONE_EVALUATORS_VERSION: u32 = 3;
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
    threshold_counterterms: Vec<Vec<StandaloneEvaluatorStackArchive<A>>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneEvaluatorStackArchive<A = Vec<u8>> {
    single_parametric: StandaloneGenericEvaluatorArchive<A>,
    iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    representative_input: Vec<StandaloneComplexInput>,
    start: usize,
    override_pos: usize,
    mult_offset: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
struct StandaloneGenericEvaluatorArchive<A = Vec<u8>> {
    exprs: Vec<A>,
    additional_fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    dual_shape: Option<Vec<Vec<usize>>>,
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
                        ExportSettings {
                            include_header: true,
                            inline_asm: backend.inline_asm(),
                            custom_header: None,
                        },
                    )
                    .map_err(|error| eyre!(error))?
                    .compile(&library_path, CompileOptions::default())
                    .map_err(|error| eyre!(error))?
                    .load()
                    .map_err(|error| eyre!(error))?;
                Ok(Self::Compiled(compiled))
            }
            StandaloneBackend::Symjit => Ok(Self::Symjit(
                evaluator.jit_compile().map_err(|error| eyre!(error))?,
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
            ArbFloat::parse(value, Some(DoubleFloat::default().get_precision()))
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

impl StandaloneNumber for ArbFloat {
    fn parse_standalone_input(value: &str) -> Result<Self> {
        ArbFloat::parse(value, Some(ARB_PRECISION_BITS)).map_err(|error| eyre!(error))
    }

    fn exact_from_rational(value: &Rational) -> Self {
        value.to_multi_prec_float(ARB_PRECISION_BITS)
    }

    fn zero_value() -> Self {
        ArbFloat::new(ARB_PRECISION_BITS)
    }

    fn one_value() -> Self {
        ArbFloat::new(ARB_PRECISION_BITS).one()
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
    fn_map.add_constant(
        parse_lit!(gammalooprs::x),
        Complex::<Rational>::try_from(Atom::Zero.as_view()).unwrap(),
    );

    for (lhs, rhs, tags, args) in parsed_entries {
        if let AtomView::Var(_) = lhs.as_view() {
            if let Ok(constant) = Complex::<Rational>::try_from(rhs.as_view()) {
                fn_map.add_constant(lhs.clone(), constant);
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
                        .with_settings(MatchSettings {
                            allow_new_wildcards_on_rhs: true,
                            ..Default::default()
                        }),
                    );
                }

                fn_map
                    .add_function(
                        function.get_symbol(),
                        function.get_symbol().get_name().into(),
                        args,
                        rhs.clone(),
                    )
                    .map_err(|error| eyre!(error))?;

                replacements.push(Replacement::new(
                    lhs.replace_multiple(&wildcards).to_pattern(),
                    rhs.replace_multiple(&wildcards),
                ));
            } else {
                fn_map
                    .add_tagged_function(
                        function.get_symbol(),
                        tags,
                        function.get_symbol().get_name().into(),
                        args,
                        rhs.clone(),
                    )
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
) -> Result<(ExpressionEvaluator<Complex<T>>, usize)> {
    let optimization_settings = OptimizationSettings {
        horner_iterations: 10,
        n_cores: 1,
        abort_check: None,
        ..Default::default()
    };
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
                .evaluator(&fn_map, params, optimization_settings.clone())
                .map_err(|error| {
                    eyre!(
                        "{error} while building iterative evaluator for {}",
                        expr.replace_multiple(&replacements)
                    )
                })?;

            tree = Some(if let Some(mut merged) = tree {
                merged
                    .merge(eval, optimization_settings.cpe_iterations)
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

        Atom::evaluator_multiple(
            &replaced_exprs,
            &fn_map,
            params,
            optimization_settings.clone(),
        )
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

fn set_override_if<A: Clone>(
    values: &mut [A],
    one: A,
    zero: A,
    over_ride: bool,
    start: usize,
    multiplicative_offset: usize,
) {
    let override_start = start * multiplicative_offset;
    values[override_start] = if over_ride { one } else { zero };
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

    fn override_input<T: StandaloneNumber>(&self) -> Result<Vec<Complex<T>>> {
        let mut input = self.representative_input::<T>()?;
        let zero = T::zero_value();
        let one = T::one_value();
        set_override_if(
            &mut input,
            Complex::new(one, T::zero_value()),
            Complex::new(zero, T::zero_value()),
            true,
            self.override_pos,
            self.mult_offset,
        );
        Ok(input)
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
           --orientation-index <usize>\n\
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

    Ok(archive)
}

fn load_json(path: impl AsRef<Path>) -> Result<StandaloneEvaluatorArchive<(), String>> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let archive: StandaloneEvaluatorArchive<(), String> = serde_json::from_slice(&binary)?;

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
            | StandaloneMethod::Summed => vec![stack.override_input::<f64>()?],
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
) -> Result<()> {
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
            | StandaloneMethod::Summed => vec![stack.override_input::<T>()?],
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
                StandaloneNumericTarget::Arb => {
                    evaluate_higher_precision_archive::<ArbFloat, _, _>(
                        archive,
                        &state_map,
                        &options,
                        custom_input.as_deref(),
                        "arb",
                        StandaloneNumericTarget::Arb,
                    )
                }
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
                StandaloneNumericTarget::Arb => {
                    evaluate_higher_precision_archive::<ArbFloat, _, _>(
                        archive,
                        &state_map,
                        &options,
                        custom_input.as_deref(),
                        "arb",
                        StandaloneNumericTarget::Arb,
                    )
                }
            }
        }
        _ => Err(eyre!(
            "Unsupported file extension {}, expected .bin or .json",
            extension.to_string_lossy()
        )),
    }
}
