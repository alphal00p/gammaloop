//#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
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
    path::{Path, PathBuf},
    time::Instant,
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate},
    domains::{
        float::Complex,
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    evaluate::{ExpressionEvaluator, FunctionMap, OptimizationSettings},
    id::{MatchSettings, Replacement},
    parse_lit,
    state::{State, StateMap},
    try_parse,
};

pub const STANDALONE_EVALUATORS_VERSION: u32 = 3;

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneCrossSectionArchive<S = Vec<u8>, T = Vec<u8>> {
    pub(crate) version: u32,
    pub(crate) symbolica_state: S,
    pub(crate) graph_terms: Vec<StandaloneCrossSectionGraphTermArchive<T>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneCrossSectionGraphTermArchive<A = Vec<u8>> {
    pub(crate) graph_name: String,
    pub(crate) orientations: Vec<Vec<i8>>,
    pub(crate) param_builder_params: Vec<A>,
    pub(crate) fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    pub(crate) raised_cut_integrands: Vec<Vec<StandaloneEvaluatorStackArchive<A>>>,
    pub(crate) counterterms: Vec<StandaloneCountertermArchive<A>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneCountertermArchive<A = Vec<u8>> {
    pub(crate) left_thresholds_evaluator: Vec<Vec<StandaloneEvaluatorStackArchive<A>>>,
    pub(crate) right_thresholds_evaluator: Vec<Vec<StandaloneEvaluatorStackArchive<A>>>,
    pub(crate) iterated_evaluator:
        StandaloneIteratedCollectionArchive<Vec<StandaloneEvaluatorStackArchive<A>>>,
    pub(crate) pass_two_evaluator: StandaloneGenericEvaluatorArchive<A>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneIteratedCollectionArchive<T> {
    pub(crate) data: Vec<T>,
    pub(crate) num_left_thresholds: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneEvaluatorStackArchive<A = Vec<u8>> {
    pub(crate) single_parametric: StandaloneGenericEvaluatorArchive<A>,
    pub(crate) iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) representative_input: Vec<Complex<f64>>,
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
                .map(|tag| tag.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            let args = args
                .iter()
                .map(|arg| {
                    let arg = arg.import_with_map(state_map)?;
                    arg.clone().try_into().map_err(|_| {
                        eyre!("Expected indeterminate in function argument, got {arg}")
                    })
                })
                .collect::<Result<Vec<_>>>()?;

            Ok((lhs_atom, rhs_atom, tags, args))
        })
        .collect()
}

fn apply_fn_map_entries(
    parsed_entries: Vec<ParsedFnMapEntry>,
) -> Result<(Vec<Replacement>, Vec<Replacement>, FunctionMap)> {
    let mut all_replacements = Vec::new();
    let mut fn_map = FunctionMap::new();
    let mut replacements = Vec::new();
    fn_map.add_constant(
        parse_lit!(gammalooprs::x),
        Complex::<Rational>::try_from(Atom::Zero.as_view()).unwrap(),
    );

    for (lhs, rhs, tags, args) in parsed_entries {
        if let AtomView::Var(_) = lhs.as_view() {
            if let Ok(value) = Complex::<Rational>::try_from(rhs.as_view()) {
                fn_map.add_constant(lhs.clone(), value);
                all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            } else {
                replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            }
        } else if let AtomView::Fun(f) = lhs.as_view() {
            if tags.is_empty() {
                let wildcards = args
                    .iter()
                    .enumerate()
                    .map(|(i, arg)| {
                        let atom: Atom = arg.clone().into();
                        Replacement::new(
                            atom.to_pattern(),
                            Atom::var(symbolica::symbol!(format!("x{i}_"))),
                        )
                        .with_settings(MatchSettings {
                            allow_new_wildcards_on_rhs: true,
                            ..Default::default()
                        })
                    })
                    .collect::<Vec<_>>();

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
        ..Default::default()
    };
    let exprs = payload
        .exprs
        .iter()
        .map(|expr| expr.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;

    let additional_reps = parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?;
    fn_map_entries.extend(additional_reps);

    let result = vec![Complex::new(0.0, 0.0); exprs.len()];
    let (replacements, all_replacements, fn_map) = apply_fn_map_entries(fn_map_entries)?;

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

            tree = Some(if let Some((atoms, mut existing)) = tree {
                existing
                    .merge(eval, optimization_settings.cpe_iterations)
                    .map_err(|e| eyre!(e))?;
                (atoms, existing)
            } else {
                (exprs.clone(), eval)
            });
        }

        tree.map(|(atoms, eval)| {
            (
                atoms,
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

fn build_stack<A: ImportWithMap>(
    stack: StandaloneEvaluatorStackArchive<A>,
    params: &[Atom],
    fn_map_entries: &[SerializedFnMapEntry<A>],
    state_map: &StateMap,
    graph_name: &str,
    label: &str,
) -> Result<LoadedStandaloneEvaluatorStack> {
    let parsed_fn_map_entries = parse_fn_map_entries(fn_map_entries, state_map)?;
    let timed_build = |payload: StandaloneGenericEvaluatorArchive<A>,
                       iterate: bool,
                       component: &str|
     -> Result<LoadedGenericEvaluator> {
        let started = Instant::now();
        let evaluator = build_evaluator(
            payload,
            params,
            parsed_fn_map_entries.clone(),
            state_map,
            iterate,
        )?;
        println!(
            "[timing] build_evaluator {graph_name}::{label}::{component} took {:?}",
            started.elapsed()
        );
        Ok(evaluator)
    };

    Ok(LoadedStandaloneEvaluatorStack {
        single_parametric: timed_build(stack.single_parametric, false, "single_parametric")?,
        iterative: stack
            .iterative
            .map(|payload| timed_build(payload, true, "iterative"))
            .transpose()?,
        summed_function_map: stack
            .summed_function_map
            .map(|payload| timed_build(payload, false, "summed_function_map"))
            .transpose()?,
        summed: stack
            .summed
            .map(|payload| timed_build(payload, false, "summed"))
            .transpose()?,
    })
}

#[derive(Default)]
pub struct LoadedStandaloneCrossSection {
    pub graph_terms: Vec<LoadedStandaloneCrossSectionGraphTerm>,
}

pub struct LoadedStandaloneCrossSectionGraphTerm {
    pub graph_name: String,
    pub orientations: Vec<Vec<i8>>,
    pub param_builder_params: Vec<Atom>,
    pub raised_cut_integrands: Vec<Vec<LoadedStandaloneEvaluatorStack>>,
    pub counterterms: Vec<LoadedStandaloneCounterterm>,
}

pub struct LoadedStandaloneCounterterm {
    pub left_thresholds_evaluator: Vec<Vec<LoadedStandaloneEvaluatorStack>>,
    pub right_thresholds_evaluator: Vec<Vec<LoadedStandaloneEvaluatorStack>>,
    pub iterated_evaluator: LoadedStandaloneIteratedCollection<Vec<LoadedStandaloneEvaluatorStack>>,
    pub pass_two_evaluator: LoadedGenericEvaluator,
}

pub struct LoadedStandaloneIteratedCollection<T> {
    pub data: Vec<T>,
    pub num_left_thresholds: usize,
}

pub struct LoadedStandaloneEvaluatorStack {
    pub single_parametric: LoadedGenericEvaluator,
    pub iterative: Option<LoadedGenericEvaluator>,
    pub summed_function_map: Option<LoadedGenericEvaluator>,
    pub summed: Option<LoadedGenericEvaluator>,
}

impl StandaloneCrossSectionArchive<(), String> {
    pub fn load(self) -> Result<LoadedStandaloneCrossSection> {
        if self.version != STANDALONE_EVALUATORS_VERSION {
            return Err(eyre!(
                "Unsupported version {} (expected {})",
                self.version,
                STANDALONE_EVALUATORS_VERSION
            ));
        }

        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)
            .with_context(|| "Failed to export Symbolica state for standalone cross section")?;
        let mut state_cursor = Cursor::new(&symbolica_state);
        let state_map = State::import(&mut state_cursor, None)?;

        self.load_impl(&state_map)
    }
}

impl StandaloneCrossSectionArchive {
    pub fn load(self) -> Result<LoadedStandaloneCrossSection> {
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

impl<S, A: ImportWithMap + Clone> StandaloneCrossSectionArchive<S, A> {
    pub fn load_impl(self, state_map: &StateMap) -> Result<LoadedStandaloneCrossSection> {
        let mut graph_terms = Vec::new();

        for graph in self.graph_terms {
            let params = graph
                .param_builder_params
                .iter()
                .map(|param| param.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;

            let raised_cut_integrands = graph
                .raised_cut_integrands
                .into_iter()
                .enumerate()
                .map(|(raised_cut_id, derivative_stacks)| {
                    derivative_stacks
                        .into_iter()
                        .enumerate()
                        .map(|(derivative_order, stack)| {
                            build_stack(
                                stack,
                                &params,
                                &graph.fn_map_entries,
                                state_map,
                                &graph.graph_name,
                                &format!(
                                    "raised_cut[{raised_cut_id}].derivative[{derivative_order}]"
                                ),
                            )
                        })
                        .collect::<Result<Vec<_>>>()
                })
                .collect::<Result<Vec<_>>>()?;

            let counterterms = graph
                .counterterms
                .into_iter()
                .enumerate()
                .map(|(cut_id, counterterm)| {
                    let build_stack_collection =
                        |payloads: Vec<StandaloneEvaluatorStackArchive<A>>,
                         label: &str|
                         -> Result<Vec<LoadedStandaloneEvaluatorStack>> {
                            payloads
                                .into_iter()
                                .enumerate()
                                .map(|(i, payload)| {
                                    build_stack(
                                        payload,
                                        &params,
                                        &graph.fn_map_entries,
                                        state_map,
                                        &graph.graph_name,
                                        &format!("counterterm[{cut_id}]::{label}[{i}]"),
                                    )
                                })
                                .collect::<Result<Vec<_>>>()
                        };

                    let build_stack_iterated =
                        |payloads: StandaloneIteratedCollectionArchive<_>,
                         label: &str|
                         -> Result<
                            LoadedStandaloneIteratedCollection<Vec<LoadedStandaloneEvaluatorStack>>,
                        > {
                            let data = payloads
                                .data
                                .into_iter()
                                .enumerate()
                                .map(|(i, payload)| {
                                    build_stack_collection(payload, &format!("{label}[{i}]"))
                                })
                                .collect::<Result<Vec<_>>>()?;

                            Ok(LoadedStandaloneIteratedCollection {
                                data,
                                num_left_thresholds: payloads.num_left_thresholds,
                            })
                        };

                    let pass_two_started = Instant::now();
                    let pass_two_evaluator = build_evaluator(
                        counterterm.pass_two_evaluator,
                        &params,
                        parse_fn_map_entries(&graph.fn_map_entries, state_map)?,
                        state_map,
                        false,
                    )?;
                    println!(
                        "[timing] build_evaluator {}::counterterm[{cut_id}]::pass_two_evaluator took {:?}",
                        graph.graph_name,
                        pass_two_started.elapsed()
                    );

                    Ok(LoadedStandaloneCounterterm {
                        left_thresholds_evaluator: counterterm
                            .left_thresholds_evaluator
                            .into_iter()
                            .enumerate()
                            .map(|(i, payload)| {
                                build_stack_collection(
                                    payload,
                                    &format!("left_thresholds_evaluator[{i}]"),
                                )
                            })
                            .collect::<Result<Vec<_>>>()?,
                        right_thresholds_evaluator: counterterm
                            .right_thresholds_evaluator
                            .into_iter()
                            .enumerate()
                            .map(|(i, payload)| {
                                build_stack_collection(
                                    payload,
                                    &format!("right_thresholds_evaluator[{i}]"),
                                )
                            })
                            .collect::<Result<Vec<_>>>()?,
                        iterated_evaluator: build_stack_iterated(
                            counterterm.iterated_evaluator,
                            "iterated_evaluator",
                        )?,
                        pass_two_evaluator,
                    })
                })
                .collect::<Result<Vec<_>>>()?;

            graph_terms.push(LoadedStandaloneCrossSectionGraphTerm {
                graph_name: graph.graph_name,
                orientations: graph.orientations,
                param_builder_params: params,
                raised_cut_integrands,
                counterterms,
            });
        }

        Ok(LoadedStandaloneCrossSection { graph_terms })
    }
}

fn load_bin(path: impl AsRef<Path>) -> Result<LoadedStandaloneCrossSection> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (archive, _): (StandaloneCrossSectionArchive, _) =
        bincode::decode_from_slice(&binary, bincode::config::standard())?;
    archive.load()
}

fn load_json(path: impl AsRef<Path>) -> Result<LoadedStandaloneCrossSection> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let archive: StandaloneCrossSectionArchive<(), String> = serde_json::from_slice(&binary)?;
    archive.load()
}

fn main() -> Result<()> {
    let input = PathBuf::from(
        std::env::args()
            .nth(1)
            .unwrap_or_else(|| "standalone_cross_section.bin".to_string()),
    );

    let Some(ext) = input.extension() else {
        return Err(eyre!("No extension, expected .bin or .json"));
    };

    let loaded = match ext.to_string_lossy().as_ref() {
        "bin" => load_bin(&input)?,
        "json" => load_json(&input)?,
        _ => {
            return Err(eyre!(
                "Unsupported file extension {}, expected .bin or .json",
                ext.to_string_lossy()
            ));
        }
    };

    println!("Loaded {} graph terms", loaded.graph_terms.len());
    for graph in &loaded.graph_terms {
        println!(
            "graph={} orientations={} raised_cuts={} counterterms={}",
            graph.graph_name,
            graph.orientations.len(),
            graph.raised_cut_integrands.len(),
            graph.counterterms.len()
        );
        for (raised_cut_id, derivative_stacks) in graph.raised_cut_integrands.iter().enumerate() {
            println!(
                "  raised_cut[{raised_cut_id}] derivative_evaluators={}",
                derivative_stacks.len()
            );
        }
        for (cut_id, counterterm) in graph.counterterms.iter().enumerate() {
            println!(
                "  counterterm[{cut_id}] left={} right={} iterated={} pass_two_exprs={}",
                counterterm.left_thresholds_evaluator.len(),
                counterterm.right_thresholds_evaluator.len(),
                counterterm.iterated_evaluator.data.len(),
                counterterm.pass_two_evaluator.0.len()
            );
        }
    }

    Ok(())
}
