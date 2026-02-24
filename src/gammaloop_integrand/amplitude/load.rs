//#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! rand = "0.9"
//! symbolica = { git = "https://github.com/benruijl/symbolica", branch = "dev", features = ["bincode"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! graphica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! ```

use std::{
    fs,
    io::Cursor,
    ops::Neg,
    path::Path,
    time::{Duration, Instant},
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use rand::Rng;
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
    symbol,
};

pub const STANDALONE_EVALUATORS_VERSION: u32 = 2;
pub const STANDALONE_MODE_RUST: u8 = 0;

#[derive(Clone, Encode, Decode)]
pub struct StandaloneEvaluatorArchive {
    pub(crate) version: u32,
    pub(crate) mode: u8,
    pub(crate) symbolica_state: Vec<u8>,
    pub(crate) graph_terms: Vec<StandaloneGraphTermArchive>,
}

#[derive(Clone, Encode, Decode)]
pub struct StandaloneGraphTermArchive {
    pub(crate) graph_name: String,
    pub(crate) orientations: Vec<Vec<i8>>,
    pub(crate) param_builder_params: Vec<Vec<u8>>,
    pub(crate) fn_map_entries: Vec<SerializedFnMapEntry>,
    pub(crate) original_integrand: StandaloneEvaluatorStackArchive,
    pub(crate) threshold_counterterms: Vec<StandaloneEvaluatorStackArchive>,
}

#[derive(Clone, Encode, Decode)]
pub struct StandaloneEvaluatorStackArchive {
    pub(crate) single_parametric: StandaloneGenericEvaluatorArchive,
    pub(crate) iterative: Option<StandaloneGenericEvaluatorArchive>,
    pub(crate) summed_function_map: Option<StandaloneGenericEvaluatorArchive>,
    pub(crate) summed: Option<StandaloneGenericEvaluatorArchive>,
    pub(crate) representative_input: Vec<Complex<f64>>,
    pub(crate) start: usize,
    pub(crate) override_pos: usize,
    pub(crate) mult_offset: usize,
}

#[derive(Clone, Encode, Decode)]
pub struct StandaloneGenericEvaluatorArchive {
    pub(crate) exprs: Vec<Vec<u8>>,
    pub(crate) additional_fn_map_entries: Vec<SerializedFnMapEntry>,
    pub(crate) dual_shape: Option<Vec<Vec<usize>>>,
}

type SerializedFnMapEntry = (Vec<u8>, Vec<u8>, Vec<Vec<u8>>, Vec<Vec<u8>>);
type ParsedFnMapEntry = (Atom, Atom, Vec<Atom>, Vec<Indeterminate>);

fn import_atom(bytes: &[u8], state_map: &StateMap) -> Result<Atom> {
    let mut cursor = Cursor::new(bytes);
    Atom::import_with_map(&mut cursor, state_map).map_err(|e| eyre!(e))
}

fn parse_fn_map_entries(
    entries: &[SerializedFnMapEntry],
    state_map: &StateMap,
) -> Result<Vec<ParsedFnMapEntry>> {
    entries
        .iter()
        .map(|(lhs, rhs, tags, args)| {
            let lhs_atom = import_atom(lhs, state_map)?;
            let rhs_atom = import_atom(rhs, state_map)?;
            let tags = tags
                .iter()
                .map(|t| import_atom(t, state_map))
                .collect::<Result<Vec<_>>>()?;
            let args = args
                .iter()
                .map(|a| {
                    let a = import_atom(a, state_map)?;

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

fn build_evaluator(
    payload: StandaloneGenericEvaluatorArchive,
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
        .map(|b| import_atom(b, state_map))
        .collect::<Result<Vec<_>>>()?;

    let additional_reps = parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?;
    fn_map_entries.extend(additional_reps);

    for (a, _, _, _) in &fn_map_entries {
        exprs.push(a.clone());
    }

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
        for s in &samples {
            let instant = Instant::now();
            eval.evaluate(s, result);
            let duration = instant.elapsed();
            if max < duration {
                max = duration;
            }
            sum += duration;
        }

        Some((sum / (n_samples as u32), max))
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

fn load(path: impl AsRef<Path>) -> Result<LoadedStandaloneEvaluators> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (archive, _): (StandaloneEvaluatorArchive, _) =
        bincode::decode_from_slice(&binary, bincode::config::standard())?;

    if archive.version != STANDALONE_EVALUATORS_VERSION {
        return Err(eyre!(
            "Unsupported version {} (expected {})",
            archive.version,
            STANDALONE_EVALUATORS_VERSION
        ));
    }
    if archive.mode != STANDALONE_MODE_RUST {
        return Err(eyre!("This loader expects rust-mode payloads"));
    }

    let mut state_cursor = Cursor::new(&archive.symbolica_state);
    let state_map = State::import(&mut state_cursor, None)?;

    let mut graph_terms = Vec::new();

    for graph in archive.graph_terms {
        let params = graph
            .param_builder_params
            .iter()
            .map(|b| import_atom(b, &state_map))
            .collect::<Result<Vec<_>>>()?;

        for p in params.iter() {
            println!("Loaded param builder param: {}", p);
        }
        let replacements = parse_fn_map_entries(&graph.fn_map_entries, &state_map)?;

        let (exprs, all_reps, single, result) = build_evaluator(
            graph.original_integrand.single_parametric,
            &params,
            replacements.clone(),
            &state_map,
            false,
        )?;

        let iterative = graph
            .original_integrand
            .iterative
            .map(|payload| {
                build_evaluator(payload, &params, replacements.clone(), &state_map, true)
            })
            .transpose()?;

        let summed = graph
            .original_integrand
            .summed
            .map(|payload| {
                build_evaluator(payload, &params, replacements.clone(), &state_map, false)
            })
            .transpose()?;

        let summed_fnmap = graph
            .original_integrand
            .summed_function_map
            .map(|payload| {
                // parse_lit!(gammaloop::integrand(1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1))
                build_evaluator(payload, &params, replacements.clone(), &state_map, false)
            })
            .transpose()?;

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
        for ct in graph.threshold_counterterms {
            let parametric = build_evaluator(
                ct.single_parametric,
                &params,
                replacements.clone(),
                &state_map,
                false,
            )?;
            let iterative = ct
                .iterative
                .map(|payload| {
                    build_evaluator(payload, &params, replacements.clone(), &state_map, true)
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

fn main() -> Result<()> {
    let input = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "standalone_evaluators_rust.bin".to_string());
    let mut loaded = load(&input)?;
    let mut rng = rand::rng();

    let orientations = loaded.graph_terms[0].orientations.clone();

    let (result_per_orientation, samples, average, max) = loaded.graph_terms[0]
        .original_integrand
        .benchmark_parametric(&orientations, &mut rng, 10);

    for (i, o) in orientations.iter().enumerate() {
        println!("{:?}", o);
        for (e, value) in loaded.graph_terms[0]
            .original_integrand
            .parametric
            .0
            .iter()
            .zip(&result_per_orientation[i])
        {
            println!("for Last value {value}, average {average:?} max {max:?}",);
        }
    }

    if let Some((samples, average, max)) = loaded.graph_terms[0]
        .original_integrand
        .benchmark_summed_fnmap(&mut rng, 1000)
    {
        for (p, v) in loaded.graph_terms[0]
            .param_builder_params
            .iter()
            .zip(samples)
        {
            println!("{:>20} = {:<}", p.to_string(), v);
        }

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
