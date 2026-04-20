//#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! rand = "0.9"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { git = "https://github.com/benruijl/symbolica", branch = "dev", default-features = false, features = ["bincode", "serde"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! graphica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! ```

#![allow(dead_code)]
use std::{
    fs,
    io::Cursor,
    ops::Neg,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use rand::Rng;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate},
    domains::{
        float::Complex,
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

use crate::processes::StandaloneNumericTarget;

pub const STANDALONE_EVALUATORS_VERSION: u32 = 3;
pub const STANDALONE_MODE_RUST: u8 = 0;

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneEvaluatorArchive<S = Vec<u8>, T = Vec<u8>> {
    pub(crate) version: u32,
    pub(crate) numeric_target: StandaloneNumericTarget,
    pub(crate) symbolica_state: S,
    pub(crate) graph_terms: Vec<StandaloneGraphTermArchive<T>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneComplexInput {
    pub(crate) re: String,
    pub(crate) im: String,
}

impl StandaloneComplexInput {
    pub(crate) fn to_f64(&self) -> Result<Complex<f64>> {
        Ok(Complex::new(self.re.parse()?, self.im.parse()?))
    }
}
impl StandaloneEvaluatorArchive<(), String> {
    pub fn load(self) -> Result<LoadedStandaloneEvaluators> {
        if self.version != STANDALONE_EVALUATORS_VERSION {
            return Err(eyre!(
                "Unsupported version {} (expected {})",
                self.version,
                STANDALONE_EVALUATORS_VERSION
            ));
        }

        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)
            .with_context(|| "Failed to export Symbolica state for standalone evaluators")?;

        let mut state_cursor = Cursor::new(&symbolica_state);
        let state_map = State::import(&mut state_cursor, None)?;

        self.load_impl(&state_map)
    }
}

impl<S, A: ImportWithMap> StandaloneEvaluatorArchive<S, A> {
    #[allow(clippy::type_complexity)]
    pub fn load_impl(self, state_map: &StateMap) -> Result<LoadedStandaloneEvaluators> {
        let mut graph_terms = Vec::new();

        for graph in self.graph_terms {
            let graph_name = graph.graph_name.as_str();
            let params = graph
                .param_builder_params
                .iter()
                .map(|b| b.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;

            for p in params.iter() {
                println!("Loaded param builder param: {}", p);
            }
            let replacements = parse_fn_map_entries(&graph.fn_map_entries, state_map)?;
            let timed_build = |label: &str,
                               payload: StandaloneGenericEvaluatorArchive<A>,
                               iterate: bool|
             -> Result<LoadedGenericEvaluator> {
                let started = Instant::now();
                let evaluator =
                    build_evaluator(payload, &params, replacements.clone(), state_map, iterate)?;
                println!(
                    "[timing] build_evaluator {}::{} took {:?}",
                    graph_name,
                    label,
                    started.elapsed()
                );
                Ok(evaluator)
            };

            let (exprs, all_reps, single, result) = timed_build(
                "original.parametric",
                graph.original_integrand.single_parametric,
                false,
            )?;

            let iterative = graph
                .original_integrand
                .iterative
                .map(|payload| timed_build("original.iterative", payload, true))
                .transpose()?;

            let summed = graph
                .original_integrand
                .summed
                .map(|payload| timed_build("original.summed", payload, false))
                .transpose()?;

            let mut fnmap_integrand = None;

            let summed_fnmap = graph
                .original_integrand
                .summed_function_map
                .map(|payload| {
                    let (_, rhs, _, _) =
                        &parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?[0];
                    fnmap_integrand = Some(rhs.clone());
                    // parse_lit!(gammaloop::integrand(1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1))
                    timed_build("original.summed_fnmap", payload, false)
                })
                .transpose()?;

            if let Some(a) = fnmap_integrand {
                println!("Comparing fnmap summed fn and parametric epression");

                if a != exprs[0] {
                    println!("They are the different:\n {}!", (&a - &exprs[0]).expand());
                } else {
                    println!("They are the same!")
                }
            }

            let original_integrand = LoadedStandaloneEvaluatorStack {
                parametric: (exprs, all_reps, single, result),
                orientation_start: graph.original_integrand.start,
                override_pos: graph.original_integrand.override_pos,
                mult_offset: graph.original_integrand.mult_offset,
                representative_input: graph
                    .original_integrand
                    .representative_input
                    .iter()
                    .map(StandaloneComplexInput::to_f64)
                    .collect::<Result<Vec<_>>>()?,
                iterative,
                summed,
                summed_fnmap,
            };
            let mut threshold_counterterms = Vec::new();
            for (ct_idx, ct) in graph.threshold_counterterms.into_iter().enumerate() {
                let ct_parametric_label = format!("ct[{ct_idx}].parametric");
                let parametric = timed_build(&ct_parametric_label, ct.single_parametric, false)?;
                let iterative = ct
                    .iterative
                    .map(|payload| {
                        let ct_iterative_label = format!("ct[{ct_idx}].iterative");
                        timed_build(&ct_iterative_label, payload, true)
                    })
                    .transpose()?;
                let ct_evaluator = LoadedStandaloneEvaluatorStack {
                    orientation_start: ct.start,
                    mult_offset: ct.mult_offset,
                    representative_input: ct
                        .representative_input
                        .iter()
                        .map(StandaloneComplexInput::to_f64)
                        .collect::<Result<Vec<_>>>()?,
                    override_pos: ct.override_pos,
                    parametric,
                    iterative,
                    summed: None,
                    summed_fnmap: None,
                };
                threshold_counterterms.push(ct_evaluator);
            }

            println!("Loaded evaluators for graph {}", graph.graph_name);
            graph_terms.push(LoadedStandaloneGraphTerm {
                orientations: graph.orientations,
                graph_name: graph.graph_name,
                param_builder_params: params,
                original_integrand,
                threshold_counterterms,
            });
        }

        Ok(LoadedStandaloneEvaluators { graph_terms })
    }
}

impl StandaloneEvaluatorArchive {
    pub fn load(self) -> Result<LoadedStandaloneEvaluators> {
        if self.version != STANDALONE_EVALUATORS_VERSION {
            return Err(eyre!(
                "Unsupported version {} (expected {})",
                self.version,
                STANDALONE_EVALUATORS_VERSION
            ));
        }

        let mut state_cursor = Cursor::new(&self.symbolica_state);
        let state_map = State::import(&mut state_cursor, None)?;

        self.load_impl(&state_map)
    }
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneGraphTermArchive<A = Vec<u8>> {
    pub(crate) graph_name: String,
    pub(crate) orientations: Vec<Vec<i8>>,
    pub(crate) param_builder_params: Vec<A>,
    pub(crate) fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    pub(crate) original_integrand: StandaloneEvaluatorStackArchive<A>,
    pub(crate) threshold_counterterms: Vec<StandaloneEvaluatorStackArchive<A>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneEvaluatorStackArchive<A = Vec<u8>> {
    pub(crate) single_parametric: StandaloneGenericEvaluatorArchive<A>,
    pub(crate) iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) representative_input: Vec<StandaloneComplexInput>,
    pub(crate) start: usize,
    pub(crate) override_pos: usize,
    pub(crate) mult_offset: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneGenericEvaluatorArchive<A = Vec<u8>> {
    pub(crate) exprs: Vec<A>,
    pub(crate) additional_fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    pub(crate) dual_shape: Option<Vec<Vec<usize>>>,
}

type SerializedFnMapEntry<A> = (A, A, Vec<A>, Vec<A>);
type ParsedFnMapEntry = (Atom, Atom, Vec<Atom>, Vec<Indeterminate>);
type LoadedGenericEvaluator = (
    Vec<Atom>,
    Vec<Replacement>,
    ExpressionEvaluator<Complex<f64>>,
    Vec<Complex<f64>>,
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
                "Unsupported backend '{value}', expected eager, c++, assembly, or symjit"
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
                "Unsupported method '{value}', expected single_parametric, iterative, summed_function_map, or summed"
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
    ThresholdCounterterm(usize),
}

impl StandaloneStackSelection {
    fn parse(value: &str) -> Result<Self> {
        if value == "original" {
            return Ok(Self::Original);
        }

        if let Some(rest) = value.strip_prefix("ct:") {
            return Ok(Self::ThresholdCounterterm(rest.parse::<usize>()?));
        }

        Err(eyre!(
            "Unsupported stack '{value}', expected original or ct:<index>"
        ))
    }

    fn label(&self) -> String {
        match self {
            Self::Original => "original".to_string(),
            Self::ThresholdCounterterm(index) => format!("ct_{index}"),
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
                    .map_err(|err| eyre!(err))?
                    .compile(&library_path, CompileOptions::default())
                    .map_err(|err| eyre!(err))?
                    .load()
                    .map_err(|err| eyre!(err))?;
                Ok(Self::Compiled(compiled))
            }
            StandaloneBackend::Symjit => Ok(Self::Symjit(
                evaluator.jit_compile().map_err(|err| eyre!(err))?,
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

pub trait ImportWithMap {
    fn import_with_map(&self, state_map: &StateMap) -> Result<Atom>;
}

impl ImportWithMap for Vec<u8> {
    fn import_with_map(&self, state_map: &StateMap) -> Result<Atom> {
        let mut cursor = Cursor::new(self);
        Atom::import_with_map(&mut cursor, state_map).map_err(|e| eyre!(e))
    }
}

impl ImportWithMap for String {
    fn import_with_map(&self, _: &StateMap) -> Result<Atom> {
        try_parse!(self).map_err(|e| eyre!(e))
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
                .map(|t| t.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            let args = args
                .iter()
                .map(|a| {
                    let a = a.import_with_map(state_map)?;

                    if let Ok(s) = a.clone().try_into() {
                        Ok(s)
                    } else {
                        Err(eyre!(
                            "Expected indeterminate in function argument, got {}",
                            a
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
) -> Result<(Vec<Replacement>, Vec<Replacement>, FunctionMap)> {
    let mut all_replacements: Vec<Replacement> = vec![];
    let mut fn_map: FunctionMap = FunctionMap::new();
    let mut replacements: Vec<Replacement> = vec![];
    fn_map.add_constant(
        parse_lit!(gammalooprs::x),
        Complex::<Rational>::try_from(Atom::Zero.as_view()).unwrap(),
    );
    for (lhs, rhs, tags, args) in parsed_entries {
        if let AtomView::Var(_) = lhs.as_view() {
            if let Ok(t) = Complex::<Rational>::try_from(rhs.as_view()) {
                fn_map.add_constant(lhs.clone(), t);

                all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            } else {
                replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            }
        } else if let AtomView::Fun(f) = lhs.as_view() {
            if tags.is_empty() {
                let mut wildcards = Vec::new();
                for (i, a) in args.iter().enumerate() {
                    let atom: Atom = a.clone().into();
                    wildcards.push(
                        Replacement::new(atom.to_pattern(), Atom::var(symbol!(format!("x{i}_"))))
                            .with_settings(MatchSettings {
                                allow_new_wildcards_on_rhs: true,
                                ..Default::default()
                            }),
                    )
                }

                fn_map
                    .add_function(
                        f.get_symbol(),
                        f.get_symbol().get_name().into(),
                        args,
                        rhs.clone(),
                    )
                    .map_err(|e| eyre!(e))?;

                all_replacements.push(Replacement::new(
                    lhs.replace_multiple(&wildcards).to_pattern(),
                    rhs.replace_multiple(&wildcards),
                ));
            } else {
                fn_map
                    .add_tagged_function(
                        f.get_symbol(),
                        tags,
                        f.get_symbol().get_name().into(),
                        args,
                        rhs.clone(),
                    )
                    .map_err(|e| eyre!(e))?;
                all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            }
        } else {
            all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
        }
    }

    Ok((replacements, all_replacements, fn_map))
}

#[allow(clippy::type_complexity)]
fn build_evaluator<A: ImportWithMap>(
    payload: StandaloneGenericEvaluatorArchive<A>,
    params: &[Atom],
    mut fn_map_entries: Vec<ParsedFnMapEntry>,
    state_map: &StateMap,
    iterate: bool,
) -> Result<LoadedGenericEvaluator> {
    let optimization_settings = OptimizationSettings {
        horner_iterations: 10,
        n_cores: 10,
        abort_check: Some(Box::new(crate::is_interrupt_requested as fn() -> bool)),
        ..Default::default()
    };
    let exprs = payload
        .exprs
        .iter()
        .map(|b| b.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;

    let additional_reps = parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?;
    fn_map_entries.extend(additional_reps);

    // for (a, _, _, _) in &fn_map_entries {
    //     exprs.push(a.clone());
    // }

    let result = vec![Complex::new(0.0, 0.); exprs.len()];

    let (replacements, all_replacements, fn_map) = apply_fn_map_entries(fn_map_entries)?;
    // for e in exprs.iter() {
    //     println!("Loaded expression: {}", e);
    // }

    if iterate {
        let mut tree: Option<(
            Vec<Atom>,
            ExpressionEvaluator<Complex<Fraction<IntegerRing>>>,
        )> = None;

        for expr in &exprs {
            let eval = expr
                .replace_multiple(&replacements)
                .evaluator(&fn_map, params, optimization_settings.clone())
                .map_err(|e| eyre!("{e} for {expr}:{}", expr.replace_multiple(&replacements)))?;

            tree = Some(if let Some((e, mut t)) = tree {
                t.merge(eval, optimization_settings.cpe_iterations)
                    .map_err(|e| eyre!(e))?;
                (e, t)
            } else {
                (exprs.clone(), eval)
            });
        }
        tree.map(|(a, eval)| {
            (
                a,
                all_replacements,
                eval.map_coeff(&|r| Complex {
                    re: r.re.to_f64(),
                    im: r.im.to_f64(),
                }),
                result,
            )
        })
        .ok_or_else(|| eyre!("No expressions in evaluator payload"))
    } else {
        let mut replaced_exprs = vec![];
        for expr in &exprs {
            replaced_exprs.push(expr.replace_multiple(&replacements))
        }

        Atom::evaluator_multiple(
            &replaced_exprs,
            &fn_map,
            params,
            optimization_settings.clone(),
        )
        .map(|eval| {
            (
                exprs,
                all_replacements,
                eval.map_coeff(&|r| Complex {
                    re: r.re.to_f64(),
                    im: r.im.to_f64(),
                }),
                result,
            )
        })
        .map_err(|e| eyre!("{e}"))
    }
}

pub struct LoadedStandaloneEvaluators {
    pub graph_terms: Vec<LoadedStandaloneGraphTerm>,
}

pub struct LoadedStandaloneGraphTerm {
    pub graph_name: String,
    pub orientations: Vec<Vec<i8>>,
    pub param_builder_params: Vec<Atom>,
    pub original_integrand: LoadedStandaloneEvaluatorStack,
    pub threshold_counterterms: Vec<LoadedStandaloneEvaluatorStack>,
}

#[allow(clippy::type_complexity)]
pub struct LoadedStandaloneEvaluatorStack {
    pub(crate) representative_input: Vec<Complex<f64>>,
    pub(crate) orientation_start: usize,
    pub(crate) override_pos: usize,
    pub(crate) mult_offset: usize,
    pub parametric: LoadedGenericEvaluator,
    pub iterative: Option<LoadedGenericEvaluator>,
    pub summed: Option<LoadedGenericEvaluator>,
    pub summed_fnmap: Option<LoadedGenericEvaluator>,
}

pub(crate) fn set_override_if<A: Clone>(
    values: &mut [A],
    one: A,
    zero: A,
    over_ride: bool,
    start: usize,
    multiplicative_offset: usize,
) {
    let o_start = start * multiplicative_offset;

    if over_ride {
        values[o_start] = one;
    } else {
        values[o_start] = zero;
    }
}

pub(crate) fn set_orientation_values_impl<A: Clone + Neg<Output = A>>(
    values: &mut [A],
    one: A,
    zero: A,
    mult_offset: usize,
    start: usize,
    orientation: &[i8],
) {
    let minusone = -(one.clone());
    let mut o_start = start * mult_offset;

    for i in orientation {
        match i {
            1 => {
                values[o_start] = one.clone();
                o_start += mult_offset;
            }
            -1 => {
                values[o_start] = minusone.clone();
                o_start += mult_offset;
            }
            0 => {
                values[o_start] = zero.clone();
                o_start += mult_offset;
            }
            _ => panic!("Should be -1,0,1"),
        }
    }
}

impl LoadedStandaloneEvaluators {
    fn select_graph_term_mut(
        &mut self,
        graph_index: usize,
        graph_name: Option<&str>,
    ) -> Result<&mut LoadedStandaloneGraphTerm> {
        if let Some(graph_name) = graph_name {
            return self
                .graph_terms
                .iter_mut()
                .find(|term| term.graph_name == graph_name)
                .ok_or_else(|| eyre!("Unknown graph '{graph_name}'"));
        }

        let graph_count = self.graph_terms.len();
        self.graph_terms.get_mut(graph_index).ok_or_else(|| {
            eyre!(
                "Graph index {} is out of range for {} graph terms",
                graph_index,
                graph_count
            )
        })
    }
}

impl LoadedStandaloneGraphTerm {
    fn stack_mut(
        &mut self,
        selection: &StandaloneStackSelection,
    ) -> Result<&mut LoadedStandaloneEvaluatorStack> {
        match selection {
            StandaloneStackSelection::Original => Ok(&mut self.original_integrand),
            StandaloneStackSelection::ThresholdCounterterm(index) => {
                self.threshold_counterterms.get_mut(*index).ok_or_else(|| {
                    eyre!(
                        "Threshold counterterm index {} is out of range for graph {}",
                        index,
                        self.graph_name
                    )
                })
            }
        }
    }
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
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
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
           --stack <original|ct:N>\n\
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

fn diff_ratio(lhs: Complex<f64>, rhs: Complex<f64>) -> Option<Complex<f64>> {
    if rhs.re == 0.0 && rhs.im == 0.0 {
        None
    } else {
        Some(lhs / rhs)
    }
}

struct StandaloneEvaluationRequest<'a> {
    backend: StandaloneBackend,
    method: StandaloneMethod,
    orientations: &'a [Vec<i8>],
    orientation_index: Option<usize>,
    custom_input: Option<&'a [Complex<f64>]>,
    artifact_root: &'a Path,
    label: &'a str,
}

impl LoadedStandaloneEvaluatorStack {
    fn selected_evaluator_mut(
        &mut self,
        method: StandaloneMethod,
    ) -> Result<&mut LoadedGenericEvaluator> {
        match method {
            StandaloneMethod::SingleParametric => Ok(&mut self.parametric),
            StandaloneMethod::Iterative => self
                .iterative
                .as_mut()
                .ok_or_else(|| eyre!("Missing iterative evaluator in standalone archive")),
            StandaloneMethod::SummedFunctionMap => self.summed_fnmap.as_mut().ok_or_else(|| {
                eyre!("Missing summed_function_map evaluator in standalone archive")
            }),
            StandaloneMethod::Summed => self
                .summed
                .as_mut()
                .ok_or_else(|| eyre!("Missing summed evaluator in standalone archive")),
        }
    }

    fn override_input(&self) -> Vec<Complex<f64>> {
        let mut new_input = self.representative_input.clone();
        set_override_if(
            &mut new_input,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            true,
            self.override_pos,
            self.mult_offset,
        );
        new_input
    }

    fn evaluate_with_backend(
        &mut self,
        request: StandaloneEvaluationRequest<'_>,
    ) -> Result<Vec<Complex<f64>>> {
        let inputs = if let Some(custom_input) = request.custom_input {
            if custom_input.len() != self.representative_input.len() {
                return Err(eyre!(
                    "Custom input length {} does not match evaluator input length {}",
                    custom_input.len(),
                    self.representative_input.len()
                ));
            }
            vec![custom_input.to_vec()]
        } else {
            match request.method {
                StandaloneMethod::SingleParametric => {
                    if let Some(index) = request.orientation_index {
                        let orientation = request.orientations.get(index).ok_or_else(|| {
                            eyre!(
                                "Orientation index {} is out of range for {} orientations",
                                index,
                                request.orientations.len()
                            )
                        })?;
                        vec![self.set_orientation(orientation)]
                    } else {
                        request
                            .orientations
                            .iter()
                            .map(|orientation| self.set_orientation(orientation))
                            .collect()
                    }
                }
                StandaloneMethod::Iterative
                | StandaloneMethod::SummedFunctionMap
                | StandaloneMethod::Summed => vec![self.override_input()],
            }
        };

        let evaluator = self.selected_evaluator_mut(request.method)?;
        let (_, _, eval, result_template) = evaluator;
        let mut runtime = StandaloneRuntimeEvaluator::build(
            eval,
            request.backend,
            request.artifact_root,
            &format!("{}_{}", request.label, request.method.as_str()),
        )?;
        let mut accumulated = vec![Complex::new(0.0, 0.0); result_template.len()];

        for input in inputs {
            let mut current = vec![Complex::new(0.0, 0.0); result_template.len()];
            runtime.evaluate(&input, &mut current);
            for (accumulated_value, current_value) in accumulated.iter_mut().zip(current) {
                *accumulated_value += current_value;
            }
        }

        Ok(accumulated)
    }

    // fn benchmark_parametric(&self){
    //     self.parametric.evaluate_single(params)
    // }
    //
    fn benchmark_summed<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
        n_samples: usize,
        compile: bool,
    ) -> Option<(Duration, Duration)> {
        let samples: Vec<_> = (0..n_samples).map(|_| self.scramble_input(rng)).collect();
        let Some((_e, _r, eval, result)) = &mut self.summed else {
            return None;
        };

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;
        if compile {
            let compile_started = Instant::now();
            let e = eval
                .export_cpp::<Complex<f64>>(
                    "bench_summed.cpp",
                    "bench_summed",
                    ExportSettings {
                        include_header: true,
                        inline_asm: symbolica::evaluate::InlineASM::AArch64,
                        custom_header: None,
                    },
                )
                .unwrap()
                .compile(
                    "bench_summed.so",
                    CompileOptions {
                        ..Default::default()
                    },
                )
                .unwrap();
            println!(
                "[timing] benchmark_summed compile took {:?}",
                compile_started.elapsed()
            );

            let mut eval = e.load().unwrap();

            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        } else {
            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        }
        Some((sum / (n_samples as u32), max))
    }

    fn benchmark_iterative<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
        n_samples: usize,
        compile: bool,
    ) -> Option<(Duration, Duration, Complex<f64>)> {
        let samples: Vec<_> = (0..n_samples).map(|_| self.scramble_input(rng)).collect();
        let Some((_e, _r, eval, result)) = &mut self.iterative else {
            return None;
        };

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;
        let mut orientation_sum = Complex::new(0.0, 0.0);
        if compile {
            let compile_started = Instant::now();
            let e = eval
                .export_cpp::<Complex<f64>>(
                    "bench_summed.cpp",
                    "bench_summed",
                    ExportSettings {
                        include_header: true,
                        inline_asm: symbolica::evaluate::InlineASM::AArch64,
                        custom_header: None,
                    },
                )
                .unwrap()
                .compile(
                    "bench_summed.so",
                    CompileOptions {
                        ..Default::default()
                    },
                )
                .unwrap();
            println!(
                "[timing] benchmark_iterative compile took {:?}",
                compile_started.elapsed()
            );

            let mut eval = e.load().unwrap();

            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                orientation_sum = Complex::new(0.0, 0.0);
                for _r in result.iter() {
                    orientation_sum += Complex::new(0.0, 0.0);
                }
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        } else {
            let mut orientation_sum;
            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                orientation_sum = Complex::new(0.0, 0.0);
                for _r in result.iter() {
                    orientation_sum += Complex::new(0.0, 0.0);
                }
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        }
        Some((sum / (n_samples as u32), max, orientation_sum))
    }

    #[allow(clippy::type_complexity)]
    fn benchmark_parametric<R: Rng + ?Sized>(
        &mut self,
        orientations: &[Vec<i8>],
        rng: &mut R,
        n_samples: usize,
    ) -> (
        Vec<Vec<Complex<f64>>>,
        Vec<Vec<Complex<f64>>>,
        Duration,
        Duration,
    ) {
        let samples: Vec<_> = (0..n_samples)
            .map(|_| {
                let mut samples = vec![];
                for o in orientations {
                    samples.push(self.scramble_input_with_orientation(o, rng))
                }
                samples
            })
            .collect();
        let (_, _, eval, result) = &mut self.parametric;
        let mut result_per_orientation =
            vec![vec![Complex::new_zero(); result.len()]; orientations.len()];

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;
        for s in &samples {
            for r in result.iter_mut() {
                *r = Complex::new(0.0, 0.0);
            }
            let instant = Instant::now();
            for (i, o) in s.iter().enumerate() {
                eval.evaluate(o, &mut result_per_orientation[i]);

                for (r, a) in result.iter_mut().zip(&result_per_orientation[i]) {
                    *r += a
                }
            }
            let duration = instant.elapsed();
            if max < duration {
                max = duration;
            }
            sum += duration;
        }

        (
            result_per_orientation,
            samples.last().unwrap().clone(),
            sum / (n_samples as u32),
            max,
        )
    }

    fn benchmark_summed_fnmap<R: Rng + ?Sized>(
        &mut self,

        rng: &mut R,
        n_samples: usize,
    ) -> Option<(Vec<Complex<f64>>, Duration, Duration)> {
        let samples: Vec<_> = (0..n_samples).map(|_| self.scramble_input(rng)).collect();
        let Some((_e, _r, eval, result)) = &mut self.summed_fnmap else {
            return None;
        };

        // for e in e {
        //     println!("{:120}", e);
        // }

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;

        for s in &samples {
            let instant = Instant::now();
            eval.evaluate(s, result);
            let duration = instant.elapsed();
            if max < duration {
                max = duration;
            }
            sum += duration;
        }

        Some((
            samples.last().unwrap().clone(),
            sum / (n_samples as u32),
            max,
        ))
    }

    fn scramble_input_with_orientation<R: Rng + ?Sized>(
        &self,
        orientation: &[i8],
        _rng: &mut R,
    ) -> Vec<Complex<f64>> {
        let mut new_input = self.representative_input.clone();
        // for n in &mut new_input {
        //     *n = *n * rng.random_range(0.8..1.2);
        // }
        set_orientation_values_impl(
            &mut new_input,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            self.mult_offset,
            self.orientation_start,
            orientation,
        );
        set_override_if(
            &mut new_input,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            false,
            self.override_pos,
            self.mult_offset,
        );
        new_input
    }

    fn scramble_input<R: Rng + ?Sized>(&self, _rng: &mut R) -> Vec<Complex<f64>> {
        self.override_input()
    }

    fn set_orientation(&self, orientation: &[i8]) -> Vec<Complex<f64>> {
        let mut sample = self.representative_input.clone();
        set_orientation_values_impl(
            &mut sample,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            self.mult_offset,
            self.orientation_start,
            orientation,
        );
        sample
    }
}

fn load_bin(path: impl AsRef<Path>) -> Result<LoadedStandaloneEvaluators> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (archive, _): (StandaloneEvaluatorArchive, _) =
        bincode::decode_from_slice(&binary, bincode::config::standard())?;

    archive.load()
}

fn load_json(path: impl AsRef<Path>) -> Result<LoadedStandaloneEvaluators> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let archive: StandaloneEvaluatorArchive<(), String> = serde_json::from_slice(&binary)?;

    archive.load()
}

fn load_custom_input(path: impl AsRef<Path>) -> Result<Vec<Complex<f64>>> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let raw: Vec<[f64; 2]> = serde_json::from_slice(&binary)
        .with_context(|| format!("Failed to parse {}", path.as_ref().display()))?;
    Ok(raw
        .into_iter()
        .map(|[re, im]| Complex::new(re, im))
        .collect())
}

fn print_backend_result(backend: StandaloneBackend, values: &[Complex<f64>]) {
    println!("backend={}", backend.as_str());
    for (index, value) in values.iter().enumerate() {
        println!("  result[{index}] = {value}");
    }
}

fn main() -> Result<()> {
    let options = parse_cli_options()?;
    let input = options.input.clone();

    let Some(ext) = input.extension() else {
        return Err(eyre!("No extension, expected .bin or .json"));
    };
    let mut loaded = match ext.to_string_lossy().as_ref() {
        "bin" => load_bin(&input)?,
        "json" => load_json(&input)?,
        _ => {
            return Err(eyre!(
                "Unsupported file extension {}, expected .bin or .json",
                ext.to_string_lossy()
            ));
        }
    };

    let compare_backends = if options.compare_backends.is_empty() {
        vec![options.backend]
    } else {
        options.compare_backends.clone()
    };
    let custom_input = options
        .input_json
        .as_ref()
        .map(load_custom_input)
        .transpose()?;
    let artifact_root = options.artifact_dir.clone().unwrap_or_else(|| {
        input
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join("standalone_backend_artifacts")
    });

    let graph = loaded.select_graph_term_mut(options.graph_index, options.graph_name.as_deref())?;
    let graph_name = graph.graph_name.clone();
    let orientations = graph.orientations.clone();
    let stack_label = options.stack.label();
    let stack = graph.stack_mut(&options.stack)?;

    println!(
        "graph={} stack={} method={} orientation={} artifact_dir={}",
        graph_name,
        stack_label,
        options.method.as_str(),
        options
            .orientation_index
            .map(|index| index.to_string())
            .unwrap_or_else(|| "all".to_string()),
        artifact_root.display()
    );

    if options.print_input {
        match options.method {
            StandaloneMethod::SingleParametric => {
                if let Some(custom_input) = custom_input.as_ref() {
                    println!("input={custom_input:?}");
                } else if let Some(index) = options.orientation_index {
                    let orientation = orientations.get(index).ok_or_else(|| {
                        eyre!(
                            "Orientation index {} is out of range for {} orientations",
                            index,
                            orientations.len()
                        )
                    })?;
                    println!("input={:?}", stack.set_orientation(orientation));
                } else {
                    for (index, orientation) in orientations.iter().enumerate() {
                        println!("orientation[{index}]={orientation:?}");
                        println!("input[{index}]={:?}", stack.set_orientation(orientation));
                    }
                }
            }
            StandaloneMethod::Iterative
            | StandaloneMethod::SummedFunctionMap
            | StandaloneMethod::Summed => {
                println!(
                    "input={:?}",
                    custom_input
                        .as_ref()
                        .cloned()
                        .unwrap_or_else(|| stack.override_input())
                );
            }
        }
    }

    let mut results = Vec::new();
    for backend in compare_backends {
        let values = stack.evaluate_with_backend(StandaloneEvaluationRequest {
            backend,
            method: options.method,
            orientations: &orientations,
            orientation_index: options.orientation_index,
            custom_input: custom_input.as_deref(),
            artifact_root: &artifact_root.join(sanitize_label(&graph_name)),
            label: &format!("{}_{}", graph_name, stack_label),
        })?;
        print_backend_result(backend, &values);
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
