use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use eyre::{Context, eyre};
use idenso::{
    color::{ColorSimplifier, ColorSimplifySettings},
    dirac::GammaSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubSetIter, SubSetLike, subset::SubSet},
    typed_vec::IndexLike,
};

use color_eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::{
    algebra::{
        algebraic_traits::RefOne,
        complex::{Complex, symbolica_traits::CompiledComplexEvaluatorSpenso},
    },
    network::{
        DEFAULT_EXACT_JOIN_LIMIT, ExecutionResult, MinResultRank, MinResultRankWith,
        PAIR_SCORE_ATOM_AWARE, PAIR_SCORE_ENTRY_AWARE, PAIR_SCORE_RESULT_RANK_ONLY, Sequential,
        SequentialExtract, SequentialRef, SmallestDegree,
    },
    shadowing::symbolica_utils::{LogPrint, SpensoPrintSettings},
};
use std::ops::Deref;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{mem::transmute, ops::Neg, path::Path};
use symbolica::{
    domains::{dual::HyperDual, float::Complex as SymComplex, rational::Fraction},
    evaluate::JITCompiledEvaluator,
    prelude::*,
};
use tracing::{debug, instrument};
use typed_index_collections::TiVec;

use crate::{
    GammaLoopContext,
    cff::orientations::GraphOrientation,
    graph::Graph,
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            amplitude::load::set_override_if,
            param_builder::{FnMapEntry, LUParams},
        },
    },
    momentum::{Helicity, sample::MomentumSample},
    numerator::symbolica_ext::NumeratorAtomExt,
    processes::{
        ContractionMode, EvaluatorBuildTimings, EvaluatorSettings, ExecutionMode,
        TensorNetworkContractionOrder,
    },
    settings::{RuntimeSettings, global::FrozenCompilationMode},
    utils::{
        ArbPrec, F, FUN_LIB, FloatLike, GS, Length, TENSORLIB, W_, f128,
        hyperdual_utils::{DualOrNot, new_from_values},
    },
};

use super::{
    ParamBuilder, RuntimeCache,
    param_builder::{ThresholdParams, UpdateAndGetParams},
};

const NETWORK_SCALAR_ALIAS_MIN_BYTES: usize = 4096;
const DUMP_EVALUATOR_PRE_NETWORK_PARSE_ENV: &str = "GAMMALOOP_DUMP_EVALUATOR_PRE_NETWORK_PARSE";
const STOP_AFTER_EVALUATOR_PRE_NETWORK_PARSE_ENV: &str =
    "GAMMALOOP_STOP_AFTER_EVALUATOR_PRE_NETWORK_PARSE";
const TRACE_PARAMETRIC_NONFINITE_ENV: &str = "GAMMALOOP_TRACE_PARAMETRIC_NONFINITE";
const DUMP_PARAMETRIC_NONFINITE_DIR_ENV: &str = "GAMMALOOP_DUMP_PARAMETRIC_NONFINITE_DIR";
static PARAMETRIC_NONFINITE_DUMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone, Copy)]
pub enum SingleOrAllOrientations<'a, OID> {
    Single {
        orientation: &'a EdgeVec<Orientation>,
        id: OID,
    },
    All {
        all: &'a TiVec<OID, EdgeVec<Orientation>>,
        filter: &'a SubSet<OID>,
    },
}
impl<'a, OID: IndexLike> SingleOrAllOrientations<'a, OID> {
    pub fn is_all(&self) -> bool {
        let SingleOrAllOrientations::All { filter, .. } = self else {
            return false;
        };
        (**filter).is_full()
    }
    pub fn iter(&self) -> SingleOrAllOrientationsIterator<'_, OID> {
        match self {
            SingleOrAllOrientations::All { all, filter } => SingleOrAllOrientationsIterator::All {
                all: *all,
                filter: (*filter).included_iter(),
            },
            SingleOrAllOrientations::Single { orientation, id } => {
                SingleOrAllOrientationsIterator::Single {
                    orientation: Some(orientation),
                    id: *id,
                }
            }
        }
    }
}

impl<OID: IndexLike> Length for SingleOrAllOrientations<'_, OID> {
    fn len(&self) -> usize {
        match self {
            SingleOrAllOrientations::Single { .. } => 1,
            SingleOrAllOrientations::All { all, filter } => {
                if filter.is_full() {
                    all.len()
                } else {
                    filter.included_iter().count()
                }
            }
        }
    }
}

pub(crate) fn evaluate_evaluator_single<T: FloatLike + GenericEvaluatorFloat>(
    generic_evaluator: &mut GenericEvaluator,
    params: &[Complex<F<T>>],
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
) -> Complex<F<T>> {
    if !record_primary_timing {
        return <T as GenericEvaluatorFloat>::get_evaluator_single(generic_evaluator)(params);
    }

    let start = std::time::Instant::now();
    let result = <T as GenericEvaluatorFloat>::get_evaluator_single(generic_evaluator)(params);
    evaluation_metadata.evaluator_evaluation_time = evaluation_metadata
        .evaluator_evaluation_time
        .saturating_add(start.elapsed());
    result
}

pub(crate) fn evaluate_evaluator<T: FloatLike + GenericEvaluatorFloat>(
    generic_evaluator: &mut GenericEvaluator,
    params: &[Complex<F<T>>],
    evaluation_metadata: &mut EvaluationMetaData,
    record_primary_timing: bool,
) -> Vec<DualOrNot<Complex<F<T>>>> {
    if !record_primary_timing {
        return <T as GenericEvaluatorFloat>::get_evaluator(generic_evaluator)(params);
    }

    let start = std::time::Instant::now();
    let result = <T as GenericEvaluatorFloat>::get_evaluator(generic_evaluator)(params);
    evaluation_metadata.evaluator_evaluation_time = evaluation_metadata
        .evaluator_evaluation_time
        .saturating_add(start.elapsed());
    result
}

#[derive(Clone)]
pub enum SingleOrAllOrientationsIterator<'a, OID> {
    Single {
        orientation: Option<&'a EdgeVec<Orientation>>,
        id: OID,
    },
    All {
        all: &'a TiVec<OID, EdgeVec<Orientation>>,
        filter: SubSetIter<'a, OID>,
    },
}

impl<'a, OID: IndexLike> Iterator for SingleOrAllOrientationsIterator<'a, OID>
where
    usize: From<OID>,
{
    type Item = (OID, &'a EdgeVec<Orientation>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SingleOrAllOrientationsIterator::Single { orientation, id } => {
                orientation.take().map(|a| (*id, a))
            }
            SingleOrAllOrientationsIterator::All { all, filter } => {
                let a = filter.next()?;
                Some((a, &all[a]))
            }
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ActiveF64Backend {
    Eager,
    Cpp,
    Assembly,
    Symjit,
}

impl ActiveF64Backend {
    pub fn as_str(self) -> &'static str {
        match self {
            ActiveF64Backend::Eager => "eager",
            ActiveF64Backend::Cpp => "c++",
            ActiveF64Backend::Assembly => "assembly",
            ActiveF64Backend::Symjit => "symjit",
        }
    }

    pub fn from_frozen_mode(mode: &FrozenCompilationMode) -> Self {
        match mode {
            FrozenCompilationMode::Eager => ActiveF64Backend::Eager,
            FrozenCompilationMode::Symjit => ActiveF64Backend::Symjit,
            FrozenCompilationMode::Cpp(_) => ActiveF64Backend::Cpp,
            FrozenCompilationMode::Assembly(_) => ActiveF64Backend::Assembly,
        }
    }
}

impl std::fmt::Display for ActiveF64Backend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, Copy, Default, Encode, Decode, PartialEq, Eq)]
pub enum EvaluatorBackendPolicy {
    #[default]
    FollowIntegrand,
    EagerOnly,
}

#[derive(Clone)]
pub struct SymjitComplexEvaluatorGL(JITCompiledEvaluator<SymComplex<f64>>);

impl std::fmt::Debug for SymjitComplexEvaluatorGL {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("SymjitComplexEvaluatorGL(..)")
    }
}

impl SymjitComplexEvaluatorGL {
    pub fn evaluate(&mut self, args: &[Complex<F<f64>>], out: &mut [Complex<F<f64>>]) {
        unsafe {
            self.0.evaluate(
                transmute::<&[Complex<F<f64>>], &[SymComplex<f64>]>(args),
                transmute::<&mut [Complex<F<f64>>], &mut [SymComplex<f64>]>(out),
            );
        }
    }
}
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema)]
pub enum EvaluatorMethod {
    SingleParametric,
    Iterative,
    SummedFunctionMap,
    Summed,
}

#[allow(clippy::derivable_impls)]
impl Default for EvaluatorMethod {
    fn default() -> Self {
        EvaluatorMethod::SingleParametric
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct EvaluatorStack {
    pub single_parametric: GenericEvaluator,
    pub iterative: Option<(GenericEvaluator, usize)>,
    // pub iterative_function_map: Option<GenericEvaluator>,
    pub summed_function_map: Option<GenericEvaluator>,
    pub summed: Option<GenericEvaluator>,
}

fn sum_iterative_outputs_for_selected_orientations<V>(
    output: Vec<V>,
    num_generated_orientations: usize,
    selected_orientation_ids: &[usize],
) -> Vec<V>
where
    V: Clone + std::ops::AddAssign,
{
    assert!(
        num_generated_orientations > 0,
        "iterative evaluators must be generated with at least one orientation"
    );
    assert!(
        !selected_orientation_ids.is_empty(),
        "iterative subset evaluation requires at least one selected orientation"
    );
    assert!(
        selected_orientation_ids
            .iter()
            .all(|selected_id| *selected_id < num_generated_orientations),
        "selected orientation ids {:?} exceed generated orientation count {}",
        selected_orientation_ids,
        num_generated_orientations
    );
    assert!(
        output.len().is_multiple_of(num_generated_orientations),
        "iterative evaluator output length {} is not divisible by the generated orientation count {}",
        output.len(),
        num_generated_orientations
    );

    output
        .chunks(num_generated_orientations)
        .map(|chunk| {
            let mut selected = selected_orientation_ids.iter();
            let first = *selected
                .next()
                .expect("selected orientation ids must be non-empty");
            let mut sum = chunk[first].clone();
            for selected_id in selected {
                sum += chunk[*selected_id].clone();
            }
            sum
        })
        .collect()
}

impl EvaluatorStack {
    pub(crate) fn generic_evaluator_count(&self) -> usize {
        let mut count = 1;
        if self.iterative.is_some() {
            count += 1;
        }
        if self.summed_function_map.is_some() {
            count += 1;
        }
        if self.summed.is_some() {
            count += 1;
        }
        count
    }

    #[instrument(skip_all)]
    fn new_single_parametric<A: AtomCore>(
        parametric_atom: &[A],
        param_builder: &ParamBuilder,
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard = crate::processes::enter_detailed_progress_span(
            "Generating Single Parametric Evaluator",
        );
        let opt_settings = settings.optimization_settings();

        GenericEvaluator::new_from_builder(
            parametric_atom
                .iter()
                .map(|atom| GS.collect_orientation_if(atom.as_atom_view(), false)),
            param_builder,
            dual_shape.clone(),
            opt_settings.clone(),
            settings,
        )
    }
    #[instrument(skip_all)]
    fn new_iterative<A: AtomCore>(
        parametric_atom: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(GenericEvaluator, usize)> {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Generating Iterative Evaluator");
        // let  n_orientations=;

        Ok((
            GenericEvaluator::new_from_builder(
                parametric_atom.iter().flat_map(|atom| {
                    orientations.iter().map(|a| {
                        let selected =
                            GS.collect_orientation_if(a.select(atom.as_atom_view()), false);
                        debug!(selected_expr = %selected.log_print(None), "Iterative");
                        selected
                    })
                }),
                param_builder,
                dual_shape.clone(),
                settings.optimization_settings(),
                settings,
            )?,
            orientations.len(),
        ))
    }

    #[instrument(skip_all)]
    fn new_summed_function_map<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard = crate::processes::enter_detailed_progress_span(
            "Generating Summed Function Map Evaluator",
        );
        let params: Vec<Atom> = (&param_builder.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();
        let mut fn_map = param_builder.fn_map.clone();

        //I(sign(1), sign(2), sign(3),...) -> I(σ1, σ2, σ3,...)

        let entries: Vec<FnMapEntry> = atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let mut args = vec![];
                let mut lhs = FunctionBuilder::new(GS.integrand);
                lhs = lhs.add_arg(i);
                for (e, _) in &orientations[0] {
                    lhs = lhs.add_arg(GS.sign(e));
                    args.push(Indeterminate::try_from(GS.sign(e)).unwrap());
                }
                let param_integrand = GS.collect_orientation_if(a.as_atom_view(), false);
                fn_map
                    .add_tagged_function(
                        GS.integrand,
                        vec![Atom::num(i)],
                        args.clone(),
                        param_integrand.clone(),
                    )
                    .map_err(|a| eyre!(a))?;
                Ok(FnMapEntry {
                    lhs: lhs.finish(),
                    rhs: param_integrand.clone(),
                    tags: vec![Atom::num(i)],
                    args,
                })
            })
            .collect::<Result<_>>()?;

        // fn_map.add_conditional(name)

        let sum = (0..entries.len()).map(|i| {
            orientations
                .iter()
                .map(|a| {
                    GS.collect_orientation_if(a.orientation_thetas() * GS.integrand(i, a), true)
                    // GS.integrand(a)
                })
                .fold(Atom::Zero, |acc, n| acc + n)
                .replace(
                    Symbol::IF.call_args([GS.override_if, W_.b_, W_.c_])
                        + Symbol::IF.call_args([GS.override_if, W_.d_, W_.e_]),
                )
                .repeat()
                .with(Symbol::IF.call_args([
                    Atom::var(GS.override_if),
                    Atom::var(W_.b_) + Atom::var(W_.d_),
                    Atom::var(W_.c_) + Atom::var(W_.e_),
                ]))
        });

        GenericEvaluator::new_from_raw_params(
            sum,
            &params,
            &fn_map,
            entries,
            settings.optimization_settings(),
            dual_shape.clone(),
            settings,
        )
    }

    #[instrument(skip_all)]
    fn new_summed<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Generating Summed Evaluator");
        let sum = atoms.iter().map(|atom| {
            orientations
                .iter()
                .map(|a| {
                    let selected = GS.collect_orientation_if(
                        a.orientation_thetas() * a.select(atom.as_atom_view()),
                        true,
                    );
                    debug!(selected_expr = %selected.log_print(None), "Iterative");
                    selected
                })
                .fold(Atom::Zero, |acc, n| acc + n)
                .replace(
                    Symbol::IF.call_args([GS.override_if, W_.b_, W_.c_])
                        + Symbol::IF.call_args([GS.override_if, W_.d_, W_.e_]),
                )
                .repeat()
                .with(Symbol::IF.call_args([
                    Atom::var(GS.override_if),
                    Atom::var(W_.b_) + Atom::var(W_.d_),
                    Atom::var(W_.c_) + Atom::var(W_.e_),
                ]))
        });

        GenericEvaluator::new_from_builder(
            sum,
            param_builder,
            dual_shape.clone(),
            settings.optimization_settings(),
            settings,
        )
    }
    pub fn new<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<Self> {
        Ok(Self::new_with_timings(atoms, param_builder, orientations, dual_shape, settings)?.0)
    }

    #[instrument(skip_all, err)]
    pub fn new_with_timings<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Building Evaluator Stack");
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_stack_new_start",
            atom_count = atoms.len(),
            orientation_count = orientations.len(),
            iterative_orientation_optimization = settings.iterative_orientation_optimization,
            summed_function_map = settings.summed_function_map,
            summed = settings.summed,
            do_algebra = settings.do_algebra,
            "Evaluator timing milestone"
        );
        let mut timings = EvaluatorBuildTimings::default();
        let spenso_started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_stack_parse_atoms_start",
            atom_count = atoms.len(),
            do_algebra = settings.do_algebra,
            "Evaluator timing milestone"
        );
        let parsed_atoms = atoms
            .iter()
            .enumerate()
            .map(|(atom_index, a)| {
                let atom_started = std::time::Instant::now();
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_start",
                    atom_index,
                    do_algebra = settings.do_algebra,
                    "Evaluator timing milestone"
                );
                // println!("Parsing {}", a.as_atom_view().log_print(Some(120)));
                let instant = std::time::Instant::now();
                let network_input = if settings.do_algebra {
                    let color_simplified = a.as_atom_view().simplify_color_with(
                        ColorSimplifySettings::default().with_cof_dimension_invariants(),
                    );
                    let gamma_simplified = color_simplified.simplify_gamma();
                    crate::debug_tags!(#generation, #profile, #compile, #term, #dump;
                        stage = "evaluator_stack_parse_atom_after_simplify_gamma",
                        atom_index,
                        log.after_gamma = gamma_simplified,
                        "Evaluator atom after gamma simplification"
                    );
                    let simplified = gamma_simplified.simplify_metrics().to_dots();
                    crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                        stage = "evaluator_stack_parse_atom_simplify_done",
                        atom_index,
                        elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                        "Evaluator timing milestone"
                    );
                    simplified
                } else {
                    crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                        stage = "evaluator_stack_parse_atom_simplify_skipped",
                        atom_index,
                        elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                        "Evaluator timing milestone"
                    );
                    a.as_atom_view().to_cof_dimension_invariants()
                };
                crate::debug_tags!(#generation, #profile, #compile, #term, #dump;
                    stage = "evaluator_stack_parse_atom_before_network_parse",
                    atom_index,
                    log.atom = network_input,
                    "Evaluator atom before network parsing"
                );
                if let Some(path) =
                    std::env::var_os(DUMP_EVALUATOR_PRE_NETWORK_PARSE_ENV).map(std::path::PathBuf::from)
                {
                    std::fs::write(&path, network_input.to_plain_string()).with_context(|| {
                        format!("failed to write evaluator pre-network atom to {}", path.display())
                    })?;
                }
                if std::env::var_os(STOP_AFTER_EVALUATOR_PRE_NETWORK_PARSE_ENV).is_some() {
                    return Err(eyre!(
                        "stopped after evaluator pre-network parse dump because {STOP_AFTER_EVALUATOR_PRE_NETWORK_PARSE_ENV} is set"
                    ));
                }
                let mut net = network_input.parse_into_net()?;
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_net_done",
                    atom_index,
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );

                // println!("Net: {}", net.dot_pretty());
                let scalar_aliases = net.alias_scalar_refs(|_, scalar| {
                    scalar.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES
                });
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_scalar_aliases_done",
                    atom_index,
                    threshold_bytes = NETWORK_SCALAR_ALIAS_MIN_BYTES,
                    aliases_created = scalar_aliases.aliases_created(),
                    aliased_terms = scalar_aliases.aliased_terms(),
                    aliased_bytes = scalar_aliases.aliased_bytes(),
                    max_aliased_bytes = scalar_aliases.max_aliased_bytes(),
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                crate::debug_tags!(#generation, #compile, #term, #dump;
                    stage = "evaluator_stack_parse_atom_network_dump",
                    atom_index,
                    file.atom = %a.as_atom_view().to_canonical_string(),
                    file.network = %net.dot_pretty(),
                    "Parsed evaluator network dump"
                );

                let parse_elapsed = instant.elapsed();
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_parse_elapsed",
                    atom_index,
                    elapsed_ms = parse_elapsed.as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                let instant = std::time::Instant::now();

                macro_rules! execute_min_result_rank {
                    ($execution_strategy:ty) => {
                        match settings.tensor_network_contraction_order {
                            TensorNetworkContractionOrder::SparseAtomAware => {
                                net.execute::<$execution_strategy, MinResultRank, _, _, _>(
                                    TENSORLIB.read().unwrap().deref(),
                                    FUN_LIB.deref(),
                                )
                            }
                            TensorNetworkContractionOrder::AtomAware => {
                                net.execute::<
                                    $execution_strategy,
                                    MinResultRankWith<
                                        { PAIR_SCORE_ATOM_AWARE },
                                        { DEFAULT_EXACT_JOIN_LIMIT },
                                    >,
                                    _,
                                    _,
                                    _,
                                >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref())
                            }
                            TensorNetworkContractionOrder::ResultRankOnly => {
                                net.execute::<
                                    $execution_strategy,
                                    MinResultRankWith<
                                        { PAIR_SCORE_RESULT_RANK_ONLY },
                                        { DEFAULT_EXACT_JOIN_LIMIT },
                                    >,
                                    _,
                                    _,
                                    _,
                                >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref())
                            }
                            TensorNetworkContractionOrder::EntryAware => {
                                net.execute::<
                                    $execution_strategy,
                                    MinResultRankWith<
                                        { PAIR_SCORE_ENTRY_AWARE },
                                        { DEFAULT_EXACT_JOIN_LIMIT },
                                    >,
                                    _,
                                    _,
                                    _,
                                >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref())
                            }
                        }
                    };
                }

                match settings.spenso_execution_mode {
                    (ExecutionMode::Sequential, ContractionMode::SmallestDegree) => {
                        net.execute::<Sequential, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )?;
                    }
                    (ExecutionMode::Sequential, ContractionMode::MinResultRank) => {
                        execute_min_result_rank!(Sequential)?;
                    }
                    (ExecutionMode::SequentialRef, ContractionMode::SmallestDegree) => {
                        net.execute::<SequentialRef, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )?;
                    }
                    (ExecutionMode::SequentialRef, ContractionMode::MinResultRank) => {
                        execute_min_result_rank!(SequentialRef)?;
                    }
                    (ExecutionMode::SequentialExtract, ContractionMode::SmallestDegree) => {
                        net.execute::<SequentialExtract, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )?;
                    }
                    (ExecutionMode::SequentialExtract, ContractionMode::MinResultRank) => {
                        execute_min_result_rank!(SequentialExtract)?;
                    }
                    _ => {
                        net.execute::<Sequential, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )?;
                    }
                }

                // println!("Executing ", net.dot_pretty());
                net.execute::<SequentialRef, SmallestDegree, _, _, _>(
                    TENSORLIB.read().unwrap().deref(),
                    FUN_LIB.deref(),
                )?;

                let execute_elapsed = instant.elapsed();
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_execute_elapsed",
                    atom_index,
                    elapsed_ms = execute_elapsed.as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_execute_done",
                    atom_index,
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );

                let result = net
                    .result_scalar()
                    .map(|a| match a {
                        ExecutionResult::One => Atom::num(1),
                        ExecutionResult::Zero => Atom::Zero,
                        ExecutionResult::Val(v) => v.into_owned(),
                    })
                    .map(|root| net.resolve_scalar_aliases(&scalar_aliases, root))
                    .map_err(|a| {
                        Report::from(a)
                            .with_note(|| format!("Network looks like: {}", net.dot_pretty()))
                    });
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_done",
                    atom_index,
                    success = result.is_ok(),
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                result
            })
            .collect::<Result<Vec<_>>>()?;
        timings.spenso_time += spenso_started.elapsed();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_stack_parse_atoms_done",
            atom_count = parsed_atoms.len(),
            orientation_count = orientations.len(),
            elapsed_ms = timings.spenso_time.as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );

        let symbolica_started = std::time::Instant::now();
        let iterative_started = std::time::Instant::now();
        let iterative = if settings.iterative_orientation_optimization {
            Some(
                Self::new_iterative(
                    &parsed_atoms,
                    param_builder,
                    orientations,
                    &dual_shape,
                    settings,
                )
                .with_context(|| "Failed to create iterative evaluator")?,
            )
        } else {
            None
        };
        if settings.iterative_orientation_optimization {
            crate::debug_tags!(#generation, #profile, #compile, #summary;
                stage = "evaluator_stack_new_iterative_done",
                orientation_count = orientations.len(),
                elapsed_ms = iterative_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Evaluator timing milestone"
            );
        }

        let summed_function_map_started = std::time::Instant::now();
        let summed_function_map = if settings.summed_function_map {
            Some(
                Self::new_summed_function_map(
                    &parsed_atoms,
                    param_builder,
                    orientations,
                    &dual_shape,
                    settings,
                )
                .with_context(|| "Failed to create summed function map")?,
            )
        } else {
            None
        };
        if settings.summed_function_map {
            crate::debug_tags!(#generation, #profile, #compile, #summary;
                stage = "evaluator_stack_new_summed_function_map_done",
                orientation_count = orientations.len(),
                elapsed_ms = summed_function_map_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Evaluator timing milestone"
            );
        }

        let summed_started = std::time::Instant::now();
        let summed = if settings.summed {
            Some(
                Self::new_summed(
                    &parsed_atoms,
                    param_builder,
                    orientations,
                    &dual_shape,
                    settings,
                )
                .with_context(|| "Failed to create summed ")?,
            )
        } else {
            None
        };
        if settings.summed {
            crate::debug_tags!(#generation, #profile, #compile, #summary;
                stage = "evaluator_stack_new_summed_done",
                orientation_count = orientations.len(),
                elapsed_ms = summed_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Evaluator timing milestone"
            );
        }

        let single_started = std::time::Instant::now();
        let single_parametric =
            Self::new_single_parametric(&parsed_atoms, param_builder, &dual_shape, settings)
                .with_context(|| "Failed to create parametric")?;
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_stack_new_single_parametric_done",
            orientation_count = orientations.len(),
            elapsed_ms = single_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );
        timings.symbolica_time += symbolica_started.elapsed();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_stack_new_done",
            atom_count = parsed_atoms.len(),
            orientation_count = orientations.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            spenso_ms = timings.spenso_time.as_secs_f64() * 1000.0,
            symbolica_ms = timings.symbolica_time.as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );

        Ok((
            EvaluatorStack {
                single_parametric,
                iterative,
                summed_function_map,
                summed,
            },
            timings,
        ))
    }

    fn evaluate_parametric<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Vec<DualOrNot<Complex<F<T>>>>
    where
        usize: From<OID>,
    {
        let mut result: Option<Vec<DualOrNot<Complex<F<T>>>>> = None;
        for (orientation_id, e) in orientations.iter() {
            input.set_orientation_values(e);
            let output = evaluate_evaluator(
                &mut self.single_parametric,
                input.as_slice(),
                evaluation_metadata,
                record_primary_timing,
            );
            if std::env::var_os(TRACE_PARAMETRIC_NONFINITE_ENV).is_some() {
                let output_nonfinite = output.iter().any(|entry| match entry {
                    DualOrNot::Dual(dual_result) => dual_result.values.iter().any(|value| {
                        value.re.is_nan()
                            || value.re.is_infinite()
                            || value.im.is_nan()
                            || value.im.is_infinite()
                    }),
                    DualOrNot::NonDual(value) => {
                        value.re.is_nan()
                            || value.re.is_infinite()
                            || value.im.is_nan()
                            || value.im.is_infinite()
                    }
                });

                if output_nonfinite {
                    let f128_params = input
                        .as_slice()
                        .iter()
                        .map(|value| {
                            Complex::new(
                                F::<f128>::from_ff64(value.re.into_ff64()),
                                F::<f128>::from_ff64(value.im.into_ff64()),
                            )
                        })
                        .collect::<Vec<_>>();
                    let mut f128_out =
                        vec![Complex::default(); self.single_parametric.compute_out_size()];
                    self.single_parametric
                        .f128
                        .evaluate(&f128_params, &mut f128_out);

                    let arb_params = input
                        .as_slice()
                        .iter()
                        .map(|value| {
                            Complex::new(
                                F::<ArbPrec>::from_ff64(value.re.into_ff64()),
                                F::<ArbPrec>::from_ff64(value.im.into_ff64()),
                            )
                        })
                        .collect::<Vec<_>>();
                    let mut arb_out =
                        vec![Complex::default(); self.single_parametric.compute_out_size()];
                    self.single_parametric
                        .arb
                        .evaluate(&arb_params, &mut arb_out);

                    let f128_nonfinite_count = f128_out
                        .iter()
                        .filter(|value| {
                            value.re.is_nan()
                                || value.re.is_infinite()
                                || value.im.is_nan()
                                || value.im.is_infinite()
                        })
                        .count();
                    let arb_nonfinite_count = arb_out
                        .iter()
                        .filter(|value| {
                            value.re.is_nan()
                                || value.re.is_infinite()
                                || value.im.is_nan()
                                || value.im.is_infinite()
                        })
                        .count();
                    let f128_dump = f128_out
                        .iter()
                        .enumerate()
                        .map(|(index, value)| format!("{index}: {value:+16e}"))
                        .collect::<Vec<_>>()
                        .join("\n");
                    let arb_dump = arb_out
                        .iter()
                        .enumerate()
                        .map(|(index, value)| format!("{index}: {value:+16e}"))
                        .collect::<Vec<_>>()
                        .join("\n");
                    let params_dump = input
                        .as_slice()
                        .iter()
                        .enumerate()
                        .map(|(index, value)| format!("{index:04}\t{value:+16e}"))
                        .collect::<Vec<_>>()
                        .join("\n");

                    if let Some(dump_dir) = std::env::var_os(DUMP_PARAMETRIC_NONFINITE_DIR_ENV) {
                        let dump_dir = std::path::PathBuf::from(dump_dir);
                        let _ = std::fs::create_dir_all(&dump_dir);
                        let dump_index =
                            PARAMETRIC_NONFINITE_DUMP_COUNTER.fetch_add(1, Ordering::Relaxed);
                        let stem = format!(
                            "parametric_nonfinite_{dump_index:04}_orientation{}",
                            usize::from(orientation_id)
                        );
                        let _ = std::fs::write(
                            dump_dir.join(format!("{stem}.txt")),
                            format!(
                                "orientation_id={}\norientations_start={}\noverride_pos={}\nmultiplicative_offset={}\nf128_nonfinite_count={}\narb_nonfinite_count={}\n\n# f128_out\n{}\n\n# arb_out\n{}\n\n# params_after_orientation\n{}\n",
                                usize::from(orientation_id),
                                input.orientations_start,
                                input.override_pos,
                                input.multiplicative_offset,
                                f128_nonfinite_count,
                                arb_nonfinite_count,
                                f128_dump,
                                arb_dump,
                                params_dump
                            ),
                        );
                    }
                    crate::debug_tags!(#integration, #inspect, #dump;
                        stage = "parametric_orientation_nonfinite",
                        orientation_id = usize::from(orientation_id),
                        output_nonfinite = output_nonfinite,
                        f128_nonfinite_count = f128_nonfinite_count,
                        arb_nonfinite_count = arb_nonfinite_count,
                        f128 = %f128_dump,
                        arb = %arb_dump,
                        "single-parametric orientation produced nonfinite output"
                    );
                }
            }
            if let Some(result) = &mut result {
                for (r, v) in result.iter_mut().zip(output) {
                    *r += v;
                }
            } else {
                result = Some(output)
            }
        }
        result.unwrap()
    }

    fn evaluate_iterative<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>>
    where
        usize: From<OID>,
    {
        let Some((iterative, len)) = &mut self.iterative else {
            return Err(eyre!(
                "Iterative evaluator not available. Regenerate with iterative set to true."
            ));
        };

        let mut result = vec![];
        if orientations.is_all() {
            let mut push = 0;
            let mut val = None;
            for r in evaluate_evaluator(
                iterative,
                input.as_slice(),
                evaluation_metadata,
                record_primary_timing,
            ) {
                push += 1;
                if let Some(mut value) = val.take() {
                    value += r;
                    if push == *len {
                        push = 0;
                        result.push(value);
                    } else {
                        val = Some(value)
                    }
                } else {
                    val = Some(r);
                }
            }
            Ok(result)
        } else {
            let selected_orientation_ids = orientations
                .iter()
                .map(|(id, _)| usize::from(id))
                .collect::<Vec<_>>();
            let output = evaluate_evaluator(
                iterative,
                input.as_slice(),
                evaluation_metadata,
                record_primary_timing,
            );
            Ok(sum_iterative_outputs_for_selected_orientations(
                output,
                *len,
                &selected_orientation_ids,
            ))
        }
    }

    fn evaluate_summed_fnmap<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>>
    where
        usize: From<OID>,
    {
        let Some(summed_function_map) = &mut self.summed_function_map else {
            return Err(eyre!(
                "Runtime requested evaluator_method=SummedFunctionMap, but this integrand was generated without a summed function-map evaluator. Regenerate with global.generation.evaluator.summed_function_map=true, or set process runtime general.evaluator_method=SingleParametric."
            ));
        };

        if orientations.is_all() {
            input.override_if(true);
            // if let Some(exprs) = &summed_function_map.exprs {
            //     for e in exprs {
            //         debug!(expr=%e.log_print(None),"Summed evaluator");
            //     }
            // }

            Ok(evaluate_evaluator(
                summed_function_map,
                input.as_slice(),
                evaluation_metadata,
                record_primary_timing,
            ))
        } else {
            let mut result: Option<Vec<DualOrNot<Complex<F<T>>>>> = None;
            for (_, e) in orientations.iter() {
                input.set_orientation_values(e);
                let output = evaluate_evaluator(
                    summed_function_map,
                    input.as_slice(),
                    evaluation_metadata,
                    record_primary_timing,
                );
                if let Some(result) = &mut result {
                    for (r, v) in result.iter_mut().zip(output) {
                        *r += v;
                    }
                } else {
                    result = Some(output)
                }
            }
            Ok(result.unwrap())
        }
    }

    fn evaluate_summed<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        mut input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>>
    where
        usize: From<OID>,
    {
        let Some(summed) = &mut self.summed else {
            return Err(eyre!(
                "Summed evaluator not available. Regenerate with summed set to true."
            ));
        };

        if orientations.is_all() {
            input.override_if(true);

            Ok(evaluate_evaluator(
                summed,
                input.as_slice(),
                evaluation_metadata,
                record_primary_timing,
            ))
        } else {
            let mut result: Option<Vec<DualOrNot<Complex<F<T>>>>> = None;
            for (_, e) in orientations.iter() {
                input.set_orientation_values(e);
                let output = evaluate_evaluator(
                    summed,
                    input.as_slice(),
                    evaluation_metadata,
                    record_primary_timing,
                );
                if let Some(result) = &mut result {
                    for (r, v) in result.iter_mut().zip(output) {
                        *r += v;
                    }
                } else {
                    result = Some(output)
                }
            }
            Ok(result.unwrap())
        }
    }
    #[instrument(
        name = "evaluate",
        level = "debug",
        skip(
            self,
            input,
            orientations,
            settings,
            evaluation_metadata,
            record_primary_timing
        ),
        fields(
            num_orientations = orientations.len(),
            method = ?settings.general.evaluator_method,
        )
    )]
    pub fn evaluate<'a, T: FloatLike, OID: IndexLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        orientations: SingleOrAllOrientations<'a, OID>,
        settings: &RuntimeSettings,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>>
    where
        usize: From<OID>,
    {
        match settings.general.evaluator_method {
            EvaluatorMethod::SingleParametric => Ok(self.evaluate_parametric(
                input,
                orientations,
                evaluation_metadata,
                record_primary_timing,
            )),
            EvaluatorMethod::Iterative => self.evaluate_iterative(
                input,
                orientations,
                evaluation_metadata,
                record_primary_timing,
            ),
            EvaluatorMethod::SummedFunctionMap => self.evaluate_summed_fnmap(
                input,
                orientations,
                evaluation_metadata,
                record_primary_timing,
            ),
            EvaluatorMethod::Summed => self.evaluate_summed(
                input,
                orientations,
                evaluation_metadata,
                record_primary_timing,
            ),
        }
    }

    #[instrument(
          name = "compile",
          level = "info",
          skip(self, path, name, frozen_mode),
          fields(
              name = %name.as_ref(),
              path = %path.as_ref().display(),
          )
      )]
    pub fn compile(
        &mut self,
        name: impl AsRef<str>,
        path: impl AsRef<Path>,
        frozen_mode: &FrozenCompilationMode,
    ) -> Result<()> {
        let name = name.as_ref();
        self.single_parametric.compile_external(
            path.as_ref().join(name).with_extension("cpp"),
            name,
            path.as_ref().join(name).with_extension("so"),
            frozen_mode,
        )?;

        if let Some((iterative, _)) = &mut self.iterative {
            iterative.compile_external(
                path.as_ref()
                    .join(format!("{}_iterative", name))
                    .with_extension("cpp"),
                format!("{}_iterative", name),
                path.as_ref()
                    .join(format!("{}_iterative", name))
                    .with_extension("so"),
                frozen_mode,
            )?;
        }

        if let Some(summed_function_map) = &mut self.summed_function_map {
            summed_function_map.compile_external(
                path.as_ref()
                    .join(format!("{}_summed_function_map", name))
                    .with_extension("cpp"),
                format!("{}_summed_function_map", name),
                path.as_ref()
                    .join(format!("{}_summed_function_map", name))
                    .with_extension("so"),
                frozen_mode,
            )?;
        }

        if let Some(summed) = &mut self.summed {
            summed.compile_external(
                path.as_ref()
                    .join(format!("{}_summed", name))
                    .with_extension("cpp"),
                format!("{}_summed", name),
                path.as_ref()
                    .join(format!("{}_summed", name))
                    .with_extension("so"),
                frozen_mode,
            )?;
        }
        Ok(())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        f(&mut self.single_parametric)?;

        if let Some((iterative, _)) = &mut self.iterative {
            f(iterative)?;
        }

        if let Some(summed_function_map) = &mut self.summed_function_map {
            f(summed_function_map)?;
        }

        if let Some(summed) = &mut self.summed {
            f(summed)?;
        }

        Ok(())
    }
}

#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenericEvaluator {
    pub exprs: Option<Vec<Atom>>,
    pub fn_map_entries: Vec<FnMapEntry>,
    pub exprs_len: usize,
    pub backend_policy: EvaluatorBackendPolicy,
    pub rational: Option<ExpressionEvaluator<symbolica::domains::float::Complex<Rational>>>,
    pub f64_compiled: Option<CompiledCode<Complex<f64>>>,
    pub f64_eager: ExpressionEvaluator<Complex<F<f64>>>,
    pub f128: ExpressionEvaluator<Complex<F<f128>>>,
    pub dual_shape: Option<Vec<Vec<usize>>>,
    pub arb: ExpressionEvaluator<Complex<F<ArbPrec>>>,
    pub(crate) loaded_f64_compiled: RuntimeCache<CompiledComplexEvaluatorSpenso>,
    pub(crate) symjit_f64: RuntimeCache<SymjitComplexEvaluatorGL>,
    pub(crate) active_f64_backend: RuntimeCache<ActiveF64Backend>,
}

impl GenericEvaluator {
    pub(crate) fn into_eager_only(mut self) -> Self {
        self.backend_policy = EvaluatorBackendPolicy::EagerOnly;
        self.activate_eager_only();
        self
    }

    fn activate_eager_only(&mut self) {
        self.f64_compiled = None;
        self.activate_eager();
    }

    fn is_eager_only(&self) -> bool {
        matches!(self.backend_policy, EvaluatorBackendPolicy::EagerOnly)
    }

    pub(crate) fn compute_out_size(&self) -> usize {
        let number_type_size = if let Some(dual_shape) = &self.dual_shape {
            dual_shape.len()
        } else {
            1
        };

        number_type_size * self.exprs_len
    }

    pub(crate) fn compile_external(
        &mut self,
        cpp_path: impl AsRef<Path>,
        function_name: impl AsRef<str>,
        lib_path: impl AsRef<Path>,
        frozen_mode: &FrozenCompilationMode,
    ) -> Result<()> {
        if self.is_eager_only() {
            self.activate_eager_only();
            return Ok(());
        }

        let compile_options = frozen_mode
            .to_symbolica_compile_options()
            .ok_or_else(|| eyre!("Frozen mode {frozen_mode} is not externally compiled"))?;
        let compiled = self
            .f64_eager
            .export_cpp::<Complex<f64>>(
                cpp_path.as_ref(),
                function_name.as_ref(),
                frozen_mode.export_settings(),
            )
            .map_err(|err| eyre!(err))?
            .compile(lib_path.as_ref(), compile_options)
            .map_err(|err| eyre!(err))?;
        let loaded = compiled.load().map_err(|err| eyre!(err))?;

        self.f64_compiled = Some(compiled);
        self.loaded_f64_compiled.set(loaded);
        self.symjit_f64.invalidate();
        self.active_f64_backend
            .set(ActiveF64Backend::from_frozen_mode(frozen_mode));
        Ok(())
    }

    pub(crate) fn activate_eager(&mut self) {
        self.loaded_f64_compiled.invalidate();
        self.symjit_f64.invalidate();
        self.active_f64_backend.set(ActiveF64Backend::Eager);
    }

    pub(crate) fn activate_symjit(&mut self) -> Result<()> {
        if self.is_eager_only() {
            self.activate_eager_only();
            return Ok(());
        }

        let rational = self
            .rational
            .as_ref()
            .ok_or_else(|| eyre!("Cannot build symjit backend without the rational evaluator"))?;
        let evaluator = rational
            .jit_compile::<SymComplex<f64>>(JITCompilationSettings::default())
            .map_err(|err| eyre!(err))?;
        self.loaded_f64_compiled.invalidate();
        self.symjit_f64.set(SymjitComplexEvaluatorGL(evaluator));
        self.active_f64_backend.set(ActiveF64Backend::Symjit);
        Ok(())
    }

    pub(crate) fn activate_external_from_artifact(
        &mut self,
        backend: ActiveF64Backend,
    ) -> Result<()> {
        if self.is_eager_only() {
            self.activate_eager_only();
            return Ok(());
        }

        let compiled = self
            .f64_compiled
            .as_ref()
            .ok_or_else(|| eyre!("No external compiled artifact is stored for this evaluator"))?;
        let loaded = compiled.load().map_err(|err| eyre!(err))?;
        self.symjit_f64.invalidate();
        self.loaded_f64_compiled.set(loaded);
        self.active_f64_backend.set(backend);
        Ok(())
    }

    pub(crate) fn has_external_compiled_artifact(&self) -> bool {
        self.is_eager_only() || self.f64_compiled.is_some()
    }

    pub(crate) fn active_f64_backend(&self) -> ActiveF64Backend {
        if self.is_eager_only() {
            return ActiveF64Backend::Eager;
        }

        self.active_f64_backend
            .as_ref()
            .copied()
            .unwrap_or(ActiveF64Backend::Eager)
    }

    pub(crate) fn new_from_builder<I: IntoIterator<Item = Atom>>(
        atoms: I,
        builder: &ParamBuilder<f64>,
        dual_shape: Option<Vec<Vec<usize>>>,
        optimization_settings: OptimizationSettings,
        settings: &EvaluatorSettings,
    ) -> Result<Self> {
        let params: Vec<Atom> = (&builder.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();

        Self::new_from_raw_params(
            atoms,
            &params,
            &builder.fn_map,
            builder.reps.clone(),
            optimization_settings,
            dual_shape,
            settings,
        )
    }

    pub(crate) fn new_from_raw_params<I: IntoIterator<Item = Atom>>(
        atoms: I,
        params: &[Atom],
        fn_map: &FunctionMap,
        fn_map_entries: Vec<FnMapEntry>,
        optimization_settings: OptimizationSettings,
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<Self> {
        let reps = if settings.do_fn_map_replacements {
            fn_map_entries
                .iter()
                .map(|r| r.replacement())
                .collect::<Vec<_>>()
        } else {
            vec![]
        };

        for r in &reps {
            println!("Reps!!{:#}", r)
        }

        let exprs: Vec<Atom> = atoms
            .into_iter()
            .map(|a| a.replace_multiple(&reps).replace_multiple(&reps))
            .collect();

        let mut tree: Option<ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>>> = None;
        for n in exprs.iter() {
            let eval: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = n
                .evaluator(params)
                .function_map(fn_map.clone())
                .optimization_settings(optimization_settings.clone())
                .build()
                .map_err(|e| {
                    let mut settings = SpensoPrintSettings::compact().nice_symbolica();
                    settings.max_line_length = Some(120);
                    settings.hide_all_namespaces = false;
                    eyre!(
                        "Failed to create evaluator for atom: {:120}\n: {}, with params: \n {:120}, and fn_map_entries: \n {}",
                        n.printer(settings),
                        e,
                        params
                            .iter()
                            .map(|a| a.to_string())
                            .collect::<Vec<_>>()
                            .join(", "),
                        fn_map_entries
                            .iter()
                            .map(|a| format!("{}->{}\n",a.lhs,a.rhs))
                            .collect::<Vec<_>>()
                            .join(", "),
                    )
                })?;

            tree = Some(if let Some(mut tree) = tree {
                tree.merge(eval, settings.cpe_iterations)
                    .map_err(|e| eyre!("Failed to merge evaluators: {}", e))?;
                tree
            } else {
                eval
            });
        }

        let mut tree = tree.ok_or_else(|| eyre!("No expressions to evaluate"))?;

        if let Some(dual_shape) = &dual_shape {
            let dual = HyperDual::<SymComplex<Rational>>::new(dual_shape.clone());
            let dualizer = Dualizer::new(dual, vec![]);
            tree = tree.vectorize(&dualizer).unwrap();
        }

        let rational = tree.clone();
        let f64_eager = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let f128 = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));
        let arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
            tree.map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let evaluator = GenericEvaluator {
            exprs_len: exprs.len(),
            fn_map_entries,
            exprs: if settings.store_atom {
                Some(exprs)
            } else {
                None
            },
            backend_policy: EvaluatorBackendPolicy::FollowIntegrand,
            rational: Some(rational),
            f64_compiled: None,
            f64_eager,
            f128,
            dual_shape,
            arb,
            loaded_f64_compiled: RuntimeCache::default(),
            symjit_f64: RuntimeCache::default(),
            active_f64_backend: RuntimeCache::default(),
        };

        let mut evaluator = evaluator;
        evaluator.activate_eager();
        Ok(evaluator)
    }
}

pub enum SliceMut<'a, T: FloatLike> {
    Borrowed(&'a mut [Complex<F<T>>]),
    Owned(Vec<Complex<F<T>>>),
}

pub struct InputParams<'a, T: FloatLike> {
    pub values: SliceMut<'a, T>,
    pub orientations_start: usize,
    pub override_pos: usize,
    pub multiplicative_offset: usize,
}

impl<'a, T: FloatLike> InputParams<'a, T> {
    pub(crate) fn set_orientation_values_impl<A: Clone + Neg<Output = A>, O: GraphOrientation>(
        values: &mut [A],
        one: A,
        zero: A,
        mult_offset: usize,
        start: usize,
        orientation: &O,
    ) {
        let minusone = -(one.clone());
        let mut o_start = start * mult_offset;

        for (_eid, i) in orientation.orientation() {
            // debug!("Setting orientation input for edge {}: {:?}", eid, i);
            match i {
                Orientation::Default => {
                    values[o_start] = one.clone();
                    o_start += mult_offset;
                }
                Orientation::Reversed => {
                    values[o_start] = minusone.clone();
                    o_start += mult_offset;
                }
                Orientation::Undirected => {
                    values[o_start] = zero.clone();
                    o_start += mult_offset;
                }
            }
        }
    }

    pub(crate) fn set_orientation_values<O: GraphOrientation>(&mut self, orientation: &O) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let mult_offset = self.multiplicative_offset;
        let start = self.orientations_start;
        Self::set_orientation_values_impl(
            self.as_mut_slice(),
            one,
            zero,
            mult_offset,
            start,
            orientation,
        );
    }

    pub(crate) fn override_if(&mut self, over_ride: bool) {
        let zero: Complex<F<T>> = Complex::new_re(F(T::from_f64(0.)));
        let one = zero.ref_one();
        let multiplicative_offset = self.multiplicative_offset;
        let start = self.override_pos;
        set_override_if(
            self.as_mut_slice(),
            one,
            zero,
            over_ride,
            start,
            multiplicative_offset,
        );
    }
    pub fn as_mut_slice(&mut self) -> &mut [Complex<F<T>>] {
        match &mut self.values {
            SliceMut::Borrowed(s) => s,
            SliceMut::Owned(v) => v,
        }
    }

    pub fn as_slice(&self) -> &[Complex<F<T>>] {
        match &self.values {
            SliceMut::Borrowed(s) => s,
            SliceMut::Owned(v) => v,
        }
    }
}

impl<T: FloatLike> AsMut<[Complex<F<T>>]> for InputParams<'_, T> {
    fn as_mut(&mut self) -> &mut [Complex<F<T>>] {
        self.as_mut_slice()
    }
}

pub trait GenericEvaluatorFloat<T: FloatLike = Self> {
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Complex<F<T>>;

    #[allow(clippy::type_complexity)]
    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Vec<DualOrNot<Complex<F<T>>>>;

    #[allow(clippy::too_many_arguments)]
    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<T>,
        helicities: &[Helicity],
        additional_params: &[F<T>],
        left_threshold_params: Option<&ThresholdParams<T>>,
        right_threshold_params: Option<&ThresholdParams<T>>,
        lu_params: Option<&LUParams<T>>,
    ) -> InputParams<'a, T>;
}

impl GenericEvaluatorFloat for f64 {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f64>>]) -> Complex<F<f64>> {
        #[inline(always)]
        |params: &[Complex<F<f64>>]| match generic_evaluator.active_f64_backend() {
            ActiveF64Backend::Eager => generic_evaluator.f64_eager.evaluate_single(params),
            ActiveF64Backend::Cpp | ActiveF64Backend::Assembly => {
                let compiled = generic_evaluator
                    .loaded_f64_compiled
                    .as_mut()
                    .expect("compiled f64 backend should be activated before evaluation");
                let mut out = [Complex::default()];

                unsafe {
                    compiled.evaluate(
                        transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                        transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(&mut out),
                    );
                }
                out[0]
            }
            ActiveF64Backend::Symjit => {
                let compiled = generic_evaluator
                    .symjit_f64
                    .as_mut()
                    .expect("symjit f64 backend should be activated before evaluation");
                let mut out = [Complex::default()];
                compiled.evaluate(params, &mut out);
                out[0]
            }
        }
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<Self>>]) -> Vec<DualOrNot<Complex<F<Self>>>> {
        |params: &[Complex<F<f64>>]| {
            let mut out = vec![Complex::default(); generic_evaluator.compute_out_size()];
            match generic_evaluator.active_f64_backend() {
                ActiveF64Backend::Eager => {
                    generic_evaluator.f64_eager.evaluate(params, &mut out);
                }
                ActiveF64Backend::Cpp | ActiveF64Backend::Assembly => {
                    let compiled = generic_evaluator
                        .loaded_f64_compiled
                        .as_mut()
                        .expect("compiled f64 backend should be activated before evaluation");
                    unsafe {
                        compiled.evaluate(
                            transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                            transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(&mut out),
                        );
                    }
                }
                ActiveF64Backend::Symjit => {
                    let compiled = generic_evaluator
                        .symjit_f64
                        .as_mut()
                        .expect("symjit f64 backend should be activated before evaluation");
                    compiled.evaluate(params, &mut out);
                }
            }

            if let Some(dual_shape) = &generic_evaluator.dual_shape {
                let dual_builder = HyperDual::<Complex<F<f64>>>::new(dual_shape.clone());
                let dual_size = dual_builder.values.len();

                out.chunks(dual_size)
                    .map(|chunk| DualOrNot::Dual(new_from_values(&dual_builder, chunk)))
                    .collect()
            } else {
                out.into_iter().map(DualOrNot::NonDual).collect()
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        additional_params: &[F<f64>],
        left_threshold_params: Option<&ThresholdParams<f64>>,
        right_threshold_params: Option<&ThresholdParams<f64>>,
        lu_params: Option<&LUParams<f64>>,
    ) -> InputParams<'a, f64> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            additional_params,
            left_threshold_params,
            right_threshold_params,
            lu_params,
        )
    }

    // fn get_debug_evaluator(
    //     generic_evaluator: &GenericEvaluatorDebug,
    // ) -> impl Fn(&[Complex<F<Self>>]) -> Complex<F<Self>> {
    //     #[inline(always)]
    //     |params: &[Complex<F<f64>>]| {
    //         // generic_evaluator
    //         //     .builder
    //         //     .borrow_mut()
    //         //     .fill_in_values(Vec::from_iter(params.iter().cloned()));

    //         // let a = generic_evaluator
    //         //     .builder
    //         //     .borrow()
    //         //     .replace(&generic_evaluator.expr);

    //         // debug!("Replaced atom:{:+>}", a);
    //         // generic_evaluator
    //         //     .expr
    //         //     .evaluate(
    //         //         |c| Complex::new_re(F::<f64>::from(c)),
    //         //         const_map,
    //         //         function_map,
    //         //     )
    //         //     .unwrap()

    //         // generic_evaluator.expr.evaluate(coeff_map, const_map, function_map)

    //         if let Some(compiled) = &generic_evaluator.f64_compiled {
    //             let mut out = [Complex::default()];
    //             compiled.borrow_mut().evaluate(params, &mut out);
    //             out[0]
    //         } else {
    //             generic_evaluator
    //                 .f64_eager
    //                 .borrow_mut()
    //                 .evaluate_single(params)
    //         }
    //     }
    // }
}

impl GenericEvaluatorFloat for f128 {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f128>>]) -> Complex<F<f128>> {
        // info!("USING COMPLEX F128 SINGLE");
        #[inline(always)]
        |params: &[Complex<F<f128>>]| generic_evaluator.f128.evaluate_single(params)
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<f128>>]) -> Vec<DualOrNot<Complex<F<f128>>>> {
        |params: &[Complex<F<f128>>]| {
            // info!("USING COMPLEX F128 MULTIPLE");
            let mut out = vec![Complex::default(); generic_evaluator.compute_out_size()];
            generic_evaluator.f128.evaluate(params, &mut out);

            if let Some(dual_shape) = &generic_evaluator.dual_shape {
                let dual_builder = HyperDual::<Complex<F<f128>>>::new(dual_shape.clone());
                let dual_size = dual_builder.values.len();

                out.chunks(dual_size)
                    .map(|chunk| DualOrNot::Dual(new_from_values(&dual_builder, chunk)))
                    .collect()
            } else {
                out.into_iter().map(DualOrNot::NonDual).collect()
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        additional_params: &[F<f128>],
        left_threshold_params: Option<&ThresholdParams<f128>>,
        right_threshold_params: Option<&ThresholdParams<f128>>,
        lu_params: Option<&LUParams<f128>>,
    ) -> InputParams<'a, Self> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            additional_params,
            left_threshold_params,
            right_threshold_params,
            lu_params,
        )
    }
}

impl GenericEvaluatorFloat for ArbPrec {
    #[inline(always)]
    fn get_evaluator_single(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<ArbPrec>>]) -> Complex<F<ArbPrec>> {
        #[inline(always)]
        |params: &[Complex<F<ArbPrec>>]| generic_evaluator.arb.evaluate_single(params)
    }

    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<ArbPrec>>]) -> Vec<DualOrNot<Complex<F<ArbPrec>>>> {
        |params: &[Complex<F<ArbPrec>>]| {
            let mut out = vec![Complex::default(); generic_evaluator.compute_out_size()];
            generic_evaluator.arb.evaluate(params, &mut out);

            if let Some(dual_shape) = &generic_evaluator.dual_shape {
                let dual_builder = HyperDual::<Complex<F<ArbPrec>>>::new(dual_shape.clone());
                let dual_size = dual_builder.values.len();

                out.chunks(dual_size)
                    .map(|chunk| DualOrNot::Dual(new_from_values(&dual_builder, chunk)))
                    .collect()
            } else {
                out.into_iter().map(DualOrNot::NonDual).collect()
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        cache: (bool, bool),
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
        additional_params: &[F<ArbPrec>],
        left_threshold_params: Option<&ThresholdParams<ArbPrec>>,
        right_threshold_params: Option<&ThresholdParams<ArbPrec>>,
        lu_params: Option<&LUParams<ArbPrec>>,
    ) -> InputParams<'a, Self> {
        param_builder.update_emr_and_get_params(
            cache,
            sample,
            graph,
            helicities,
            additional_params,
            left_threshold_params,
            right_threshold_params,
            lu_params,
        )
    }
}

#[cfg(test)]
mod tests {
    use idenso::color::CS;
    use symbolica::atom::Symbol;

    use crate::utils::W_;

    use super::*;

    #[test]
    fn test_function_map_summed() {
        fn integer_to_orientation(i: isize) -> Orientation {
            match i {
                1 => Orientation::Default,
                -1 => Orientation::Reversed,
                0 => Orientation::Undirected,
                _ => panic!("Invalid orientation index"),
            }
        }

        let a = [
            EdgeVec::from_iter([1, -1, 1, 1].into_iter().map(integer_to_orientation)),
            EdgeVec::from_iter([1, -1, -1, 1].into_iter().map(integer_to_orientation)),
            EdgeVec::from_iter([-1, -1, 1, 1].into_iter().map(integer_to_orientation)),
        ]
        .iter()
        .map(|a| {
            GS.collect_orientation_if(a.orientation_thetas() * GS.integrand(1, a), true)

            // GS.integrand(a)
        })
        .fold(Atom::Zero, |acc, n| acc + n)
        .replace(
            Symbol::IF.call_args([GS.override_if, W_.b_, W_.c_])
                + Symbol::IF.call_args([GS.override_if, W_.d_, W_.e_]),
        )
        .repeat()
        .with(Symbol::IF.call_args([
            Atom::var(GS.override_if),
            Atom::var(W_.b_) + Atom::var(W_.d_),
            Atom::var(W_.c_) + Atom::var(W_.e_),
        ]));

        println!("{:100}", a);
    }

    #[test]
    fn iterative_subset_output_sums_only_selected_orientation_blocks() {
        let output = vec![10, 20, 30, 40, 50, 60];

        let summed = sum_iterative_outputs_for_selected_orientations(output, 3, &[0, 2]);

        assert_eq!(summed, vec![40, 100]);
    }

    #[test]
    fn iterative_subset_output_keeps_single_selected_orientation() {
        let output = vec![10, 20, 30, 40, 50, 60];

        let summed = sum_iterative_outputs_for_selected_orientations(output, 3, &[1]);

        assert_eq!(summed, vec![20, 50]);
    }

    #[test]
    fn pi_eval() {
        let params = vec![parse!("x"), parse!("y")];
        let _evaluator = (parse!("x + y") + Symbol::PI / CS.cf)
            .evaluator(&params)
            .build()
            .unwrap();
        // assert_eq!(evaluator.evaluate_single(&[1.0, 2.0]), 3.0);
    }
}
