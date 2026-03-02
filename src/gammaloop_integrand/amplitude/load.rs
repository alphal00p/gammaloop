//#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! rand = "0.9"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { git = "https://github.com/benruijl/symbolica", branch = "dev", features = ["bincode", "serde"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! graphica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! ```

use std::{
    fs,
    io::Cursor,
    ops::Neg,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{eyre, Context, Result};
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
        CompileOptions, ExportSettings, ExpressionEvaluator, FunctionMap, OptimizationSettings,
    },
    id::{MatchSettings, Replacement},
    parse_lit,
    state::{State, StateMap},
    symbol, try_parse,
};

pub const STANDALONE_EVALUATORS_VERSION: u32 = 2;
pub const STANDALONE_MODE_RUST: u8 = 0;

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneEvaluatorArchive<S = Vec<u8>, T = Vec<u8>> {
    pub(crate) version: u32,
    pub(crate) symbolica_state: S,
    pub(crate) graph_terms: Vec<StandaloneGraphTermArchive<T>>,
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
    pub fn load_impl(self, state_map: &StateMap) -> Result<LoadedStandaloneEvaluators> {
        let mut graph_terms = Vec::new();

        for graph in self.graph_terms {
            let graph_name = graph.graph_name.as_str();
            let params = graph
                .param_builder_params
                .iter()
                .map(|b| b.import_with_map(&state_map))
                .collect::<Result<Vec<_>>>()?;

            for p in params.iter() {
                println!("Loaded param builder param: {}", p);
            }
            let replacements = parse_fn_map_entries(&graph.fn_map_entries, &state_map)?;
            let timed_build = |label: &str,
                               payload: StandaloneGenericEvaluatorArchive<A>,
                               iterate: bool|
             -> Result<(
                Vec<Atom>,
                Vec<Replacement>,
                ExpressionEvaluator<Complex<f64>>,
                Vec<Complex<f64>>,
            )> {
                let started = Instant::now();
                let evaluator =
                    build_evaluator(payload, &params, replacements.clone(), &state_map, iterate)?;
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
                        &parse_fn_map_entries(&payload.additional_fn_map_entries, &state_map)?[0];
                    fnmap_integrand = Some(rhs.clone());
                    // parse_lit!(gammaloop::integrand(1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1))
                    timed_build("original.summed_fnmap", payload, false)
                })
                .transpose()?;

            if let Some(a) = fnmap_integrand {
                println!("Comparing fnmap summed fn and parametric epression");

                if &a != &exprs[0] {
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
                representative_input: graph.original_integrand.representative_input,
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
                    representative_input: ct.representative_input,
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

trait ImportWithMap {
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

fn build_evaluator<A: ImportWithMap>(
    payload: StandaloneGenericEvaluatorArchive<A>,
    params: &[Atom],
    mut fn_map_entries: Vec<ParsedFnMapEntry>,
    state_map: &StateMap,
    iterate: bool,
) -> Result<(
    Vec<Atom>,
    Vec<Replacement>,
    ExpressionEvaluator<Complex<f64>>,
    Vec<Complex<f64>>,
)> {
    let optimization_settings = OptimizationSettings {
        horner_iterations: 10,
        n_cores: 10,
        ..Default::default()
    };
    let mut exprs = payload
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

pub(crate) struct LoadedStandaloneGraphTerm {
    pub graph_name: String,
    pub orientations: Vec<Vec<i8>>,
    pub param_builder_params: Vec<Atom>,
    pub original_integrand: LoadedStandaloneEvaluatorStack,
    pub threshold_counterterms: Vec<LoadedStandaloneEvaluatorStack>,
}

pub(crate) struct LoadedStandaloneEvaluatorStack {
    pub(crate) representative_input: Vec<Complex<f64>>,
    pub(crate) orientation_start: usize,
    pub(crate) override_pos: usize,
    pub(crate) mult_offset: usize,
    pub parametric: (
        Vec<Atom>,
        Vec<Replacement>,
        ExpressionEvaluator<Complex<f64>>,
        Vec<Complex<f64>>,
    ),
    pub iterative: Option<(
        Vec<Atom>,
        Vec<Replacement>,
        ExpressionEvaluator<Complex<f64>>,
        Vec<Complex<f64>>,
    )>,
    pub summed: Option<(
        Vec<Atom>,
        Vec<Replacement>,
        ExpressionEvaluator<Complex<f64>>,
        Vec<Complex<f64>>,
    )>,
    pub summed_fnmap: Option<(
        Vec<Atom>,
        Vec<Replacement>,
        ExpressionEvaluator<Complex<f64>>,
        Vec<Complex<f64>>,
    )>,
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

impl LoadedStandaloneEvaluatorStack {
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
        let samples: Vec<_> = (0..n_samples)
            .into_iter()
            .map(|_| self.scramble_input(rng))
            .collect();
        let Some((e, r, eval, result)) = &mut self.summed else {
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
        let samples: Vec<_> = (0..n_samples)
            .into_iter()
            .map(|_| self.scramble_input(rng))
            .collect();
        let Some((e, r, eval, result)) = &mut self.iterative else {
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
                for r in result.iter() {
                    orientation_sum += Complex::new(0.0, 0.0);
                }
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        } else {
            let mut orientation_sum = Complex::new(0.0, 0.0);
            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                orientation_sum = Complex::new(0.0, 0.0);
                for r in result.iter() {
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
            .into_iter()
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
        let samples: Vec<_> = (0..n_samples)
            .into_iter()
            .map(|_| self.scramble_input(rng))
            .collect();
        let Some((e, r, eval, result)) = &mut self.summed_fnmap else {
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
        rng: &mut R,
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

    fn scramble_input<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec<Complex<f64>> {
        let mut new_input = self.representative_input.clone();
        // for n in &mut new_input {
        //     *n = *n * rng.random_range(0.8..1.2);
        // }
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

fn main() -> Result<()> {
    let input = PathBuf::from(
        std::env::args()
            .nth(1)
            .unwrap_or_else(|| "standalone_evaluators_rust.bin".to_string()),
    );

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

    let mut rng = rand::rng();

    let orientations = loaded.graph_terms[0].orientations.clone();

    let (result_per_orientation, samples, average, max) = loaded.graph_terms[0]
        .original_integrand
        .benchmark_parametric(&orientations, &mut rng, 1000);

    let (exprs, all_reps, single, result) = &loaded.graph_terms[0].original_integrand.parametric;

    println!(" average {average:?} max {max:?}");
    for (e, value) in exprs.iter().zip(result) {
        println!("for Last value {value},");
    }
    // for (i, o) in orientations.iter().enumerate() {
    //     println!("{:?}", o);
    // }
    //
    if let Some((exprs, all_reps, single, result)) =
        &loaded.graph_terms[0].original_integrand.summed_fnmap
    {
        for t in single.export_instructions().0 {
            println!("{t}");
        }
    }

    if let Some((samples, average, max)) = loaded.graph_terms[0]
        .original_integrand
        .benchmark_summed_fnmap(&mut rng, 1000)
    {
        // for (p, v) in loaded.graph_terms[0]
        //     .param_builder_params
        //     .iter()
        //     .zip(samples)
        // {
        //     println!("{:>20} = {:<}", p.to_string(), v);
        // }

        if let Some((e, r, _, res)) = &loaded.graph_terms[0].original_integrand.summed_fnmap {
            // println!(
            //     "{:120}",
            //     parse_lit!(gammalooprs::integrand(
            //         1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, -1
            //     ))
            //     .replace_multiple(r)
            //     .rationalize(&(1, 100).into())
            //     .collect_num()
            // );
            for (e, value) in e.iter().zip(res) {
                println!("for {e:120}Last value {value}, average {average:?} max {max:?}",);
            }
        }
    }
    println!("Loaded {} graph evaluator sets", loaded.graph_terms.len());
    Ok(())
}
