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
    cff::expression::OrientationID,
    cff::orientations::GraphOrientation,
    graph::Graph,
    integrands::{
        evaluation::EvaluationMetaData,
        process::param_builder::{FnMapEntry, LUParams},
    },
    momentum::{Helicity, sample::MomentumSample},
    numerator::symbolica_ext::NumeratorAtomExt,
    processes::{
        ContractionMode, EvaluatorBuildTimings, EvaluatorSettings, ExecutionMode,
        TensorNetworkContractionOrder,
    },
    settings::{
        RuntimeSettings,
        global::{CompilationOptimizationLevel, FrozenCompilationMode},
    },
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
            FrozenCompilationMode::Symjit(_) => ActiveF64Backend::Symjit,
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
    pub(crate) explicit_orientation_sum_only: bool,
    /// Original generalized-3D-representation map key for each dense runtime
    /// orientation channel. Physical edge directions alone do not identify a
    /// raised-energy/contact residue map.
    production_orientation_ids: Vec<OrientationID>,
    pub single_parametric: GenericEvaluator,
    pub iterative: Option<(GenericEvaluator, usize)>,
    // pub iterative_function_map: Option<GenericEvaluator>,
    pub summed_function_map: Option<GenericEvaluator>,
    pub summed: Option<GenericEvaluator>,
}

impl EvaluatorStack {
    fn parametrize_residue_map_selectors<A: AtomCore>(atom: &A, selected_id: Atom) -> Atom {
        atom.as_atom_view()
            .replace(function!(OrientationID::symbol(), W_.a_))
            .with(Symbol::IF.call_args([selected_id - Atom::var(W_.a_), Atom::Zero, Atom::one()]))
    }

    fn sum_residue_map_selectors<A: AtomCore>(atom: &A) -> Atom {
        atom.as_atom_view()
            .replace(function!(OrientationID::symbol(), W_.a_))
            .with(Atom::one())
    }

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

    pub(crate) fn production_orientation_ids(&self) -> &[OrientationID] {
        &self.production_orientation_ids
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
            parametric_atom.iter().map(|atom| {
                GS.collect_orientation_if(Self::parametrize_residue_map_selectors(
                    atom,
                    Atom::var(GS.residue_map_id),
                ))
            }),
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
        production_orientation_ids: &[OrientationID],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(GenericEvaluator, usize)> {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Generating Iterative Evaluator");
        // Each output group contains one entry per generated orientation.

        Ok((
            GenericEvaluator::new_from_builder(
                parametric_atom.iter().flat_map(|atom| {
                    orientations.iter().zip(production_orientation_ids).map(
                        |(orientation, production_id)| {
                            // Select the complete residue-map entry before
                            // resolving its physical-direction metadata. An
                            // inactive entry may contain a selector-local
                            // inverse which becomes `0^-1` in another entry's
                            // physical sector; eliminating it afterwards is too
                            // late because Symbolica has already formed
                            // infinity.
                            let selected = GS.collect_orientation_if(
                                orientation.select(production_id.select(atom.as_atom_view())),
                            );
                            debug!(
                                "Selected iterative residue-map branch {}: {}",
                                production_id.0,
                                selected.log_print(Some(240))
                            );
                            selected
                        },
                    )
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
        production_orientation_ids: &[OrientationID],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard = crate::processes::enter_detailed_progress_span(
            "Generating Summed Function Map Evaluator",
        );
        let first_orientation = orientations.first().ok_or_else(|| {
            eyre!("summed function-map evaluator requires at least one residue-map orientation")
        })?;
        let params: Vec<Atom> = (&param_builder.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();
        let mut fn_map = param_builder.fn_map.clone();

        // The exact residue-map key is an argument independent of physical
        // edge signs: I(map_id, sign(1), sign(2), ...).
        let residue_map_id_arg = symbol!("residue_map_id_arg");

        let entries: Vec<FnMapEntry> = atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let mut args = vec![];
                let mut lhs = FunctionBuilder::new(GS.integrand);
                lhs = lhs.add_arg(i);
                lhs = lhs.add_arg(residue_map_id_arg);
                args.push(residue_map_id_arg.into());
                for (e, _) in first_orientation {
                    lhs = lhs.add_arg(GS.sign(e));
                    args.push(Indeterminate::try_from(GS.sign(e)).unwrap());
                }
                let param_integrand = GS.collect_orientation_if(
                    Self::parametrize_residue_map_selectors(a, Atom::var(residue_map_id_arg)),
                );
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

        // Summed evaluators contain concrete orientation calls; runtime orientation
        // selection stays in the single-parametric evaluator.
        let sum = (0..entries.len())
            .map(|i| {
                orientations
                    .iter()
                    .zip(production_orientation_ids)
                    .map(|(orientation, production_id)| {
                        GS.integrand(i, *production_id, orientation)
                    })
                    .fold(Atom::Zero, |acc, n| acc + n)
            })
            .collect::<Vec<_>>();
        let source_replacements = param_builder
            .reps
            .iter()
            .chain(&entries)
            .map(FnMapEntry::replacement)
            .collect::<Vec<_>>();
        let uses_numerator_sampling_scale =
            GenericEvaluator::source_atoms_use_numerator_sampling_scale(&sum, &source_replacements);
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            sum,
            &params,
            &fn_map,
            entries,
            settings.optimization_settings(),
            dual_shape.clone(),
            settings,
        )?;
        evaluator.uses_numerator_sampling_scale = uses_numerator_sampling_scale;
        Ok(evaluator)
    }

    #[instrument(skip_all)]
    fn new_summed<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        production_orientation_ids: &[OrientationID],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Generating Summed Evaluator");
        let sum = atoms.iter().map(|atom| {
            orientations
                .iter()
                .zip(production_orientation_ids)
                .map(|(orientation, production_id)| {
                    let selected = orientation.select(production_id.select(atom.as_atom_view()));
                    debug!(selected_expr = %selected.log_print(None), "Summed");
                    selected
                })
                .fold(Atom::Zero, |acc, n| acc + n)
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
        let production_orientation_ids = (0..orientations.len())
            .map(OrientationID)
            .collect::<Vec<_>>();
        Ok(Self::new_with_timings(
            atoms,
            param_builder,
            orientations,
            &production_orientation_ids,
            dual_shape,
            settings,
        )?
        .0)
    }

    pub(crate) fn new_explicit_sum_with_timings<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        let mut direct_settings = *settings;
        direct_settings.iterative_orientation_optimization = false;
        direct_settings.summed_function_map = false;
        direct_settings.summed = false;

        let atoms = atoms
            .iter()
            .map(Self::sum_residue_map_selectors)
            .collect::<Vec<_>>();
        let (mut stack, timings) = Self::new_with_timings(
            &atoms,
            param_builder,
            &[],
            &[],
            dual_shape,
            &direct_settings,
        )?;
        stack.explicit_orientation_sum_only = true;
        Ok((stack, timings))
    }

    pub(crate) fn new_deferred_explicit_sum_with_timings(
        compact_atom: &Atom,
        bodies: &[Atom],
        param_builder: &ParamBuilder,
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        let setup_started = std::time::Instant::now();
        let expected_compact = bodies
            .iter()
            .enumerate()
            .map(|(tag, _)| function!(GS.projected_cff_sum, tag))
            .fold(Atom::Zero, |sum, call| sum + call);
        if compact_atom != &expected_compact {
            return Err(eyre!(
                "deferred projected-CFF compact expression does not match its {} function bodies",
                bodies.len()
            ));
        }

        let validation_time = setup_started.elapsed();
        let spenso_started = std::time::Instant::now();
        let bodies = bodies
            .iter()
            .enumerate()
            .map(|(body_index, body)| {
                Self::preprocess_atom(&Self::sum_residue_map_selectors(body), body_index, settings)
            })
            .collect::<Result<Vec<_>>>()?;
        let spenso_time = spenso_started.elapsed();

        let setup_started = std::time::Instant::now();
        let mut param_builder = param_builder.clone();
        for (tag, body) in bodies.iter().enumerate() {
            param_builder
                .add_tagged_function::<Symbol>(
                    GS.projected_cff_sum,
                    vec![Atom::num(tag)],
                    format!("projected_cff_sum_{tag}"),
                    Vec::new(),
                    GS.collect_orientation_if(body.as_view()),
                )
                .map_err(|error| eyre!(error))?;
        }

        let setup_time = validation_time + setup_started.elapsed();
        let (stack, mut timings) = Self::new_explicit_sum_with_timings(
            std::slice::from_ref(compact_atom),
            &param_builder,
            dual_shape,
            settings,
        )?;
        timings.spenso_time += spenso_time;
        timings.symbolica_time += setup_time;
        Ok((stack, timings))
    }

    fn preprocess_atom<A: AtomCore>(
        a: &A,
        atom_index: usize,
        settings: &EvaluatorSettings,
    ) -> Result<Atom> {
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
                format!(
                    "failed to write evaluator pre-network atom to {}",
                    path.display()
                )
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
                    TensorNetworkContractionOrder::SparseAtomAware => net
                        .execute::<$execution_strategy, MinResultRank, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        ),
                    TensorNetworkContractionOrder::AtomAware => {
                        net.execute::<$execution_strategy, MinResultRankWith<
                            { PAIR_SCORE_ATOM_AWARE },
                            { DEFAULT_EXACT_JOIN_LIMIT },
                        >, _, _, _>(
                            TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()
                        )
                    }
                    TensorNetworkContractionOrder::ResultRankOnly => net
                        .execute::<$execution_strategy, MinResultRankWith<
                            { PAIR_SCORE_RESULT_RANK_ONLY },
                            { DEFAULT_EXACT_JOIN_LIMIT },
                        >, _, _, _>(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()),
                    TensorNetworkContractionOrder::EntryAware => {
                        net.execute::<$execution_strategy, MinResultRankWith<
                            { PAIR_SCORE_ENTRY_AWARE },
                            { DEFAULT_EXACT_JOIN_LIMIT },
                        >, _, _, _>(
                            TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()
                        )
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
                Report::from(a).with_note(|| format!("Network looks like: {}", net.dot_pretty()))
            });
        crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
            stage = "evaluator_stack_parse_atom_done",
            atom_index,
            success = result.is_ok(),
            elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );
        result
    }

    #[instrument(skip_all, err)]
    pub fn new_with_timings<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        production_orientation_ids: &[OrientationID],
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        if orientations.len() != production_orientation_ids.len() {
            return Err(eyre!(
                "runtime orientation catalog has {} physical entries but {} exact residue-map IDs",
                orientations.len(),
                production_orientation_ids.len()
            ));
        }
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
            .map(|(atom_index, atom)| Self::preprocess_atom(atom, atom_index, settings))
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
                    production_orientation_ids,
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
                    production_orientation_ids,
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
                    production_orientation_ids,
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
                explicit_orientation_sum_only: false,
                production_orientation_ids: production_orientation_ids.to_vec(),
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
            input.set_residue_map_id(self.production_orientation_ids[usize::from(orientation_id)]);
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
                                "orientation_id={}\norientations_start={}\nmultiplicative_offset={}\nf128_nonfinite_count={}\narb_nonfinite_count={}\n\n# f128_out\n{}\n\n# arb_out\n{}\n\n# params_after_orientation\n{}\n",
                                usize::from(orientation_id),
                                input.orientations_start,
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

    fn evaluate_iterative<'a, T: FloatLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>> {
        let Some((iterative, len)) = &mut self.iterative else {
            return Err(eyre!(
                "Iterative evaluator not available. Regenerate with iterative set to true."
            ));
        };

        let output = evaluate_evaluator(
            iterative,
            input.as_slice(),
            evaluation_metadata,
            record_primary_timing,
        );
        if *len == 0 {
            return Err(eyre!("Iterative evaluator has no generated orientations"));
        }

        let mut values = output.into_iter();
        let mut result = Vec::with_capacity(values.len() / *len);
        while let Some(mut sum) = values.next() {
            for _ in 1..*len {
                sum += values
                    .next()
                    .ok_or_else(|| eyre!("Iterative evaluator returned an incomplete group"))?;
            }
            result.push(sum);
        }
        Ok(result)
    }

    fn evaluate_summed_fnmap<'a, T: FloatLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>> {
        let Some(summed_function_map) = &mut self.summed_function_map else {
            return Err(eyre!(
                "Runtime requested evaluator_method=SummedFunctionMap, but this integrand was generated without a summed function-map evaluator. Regenerate with global.generation.evaluator.summed_function_map=true, or set process runtime general.evaluator_method=SingleParametric."
            ));
        };

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
    }

    fn evaluate_summed<'a, T: FloatLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>> {
        let Some(summed) = &mut self.summed else {
            return Err(eyre!(
                "Summed evaluator not available. Regenerate with summed set to true."
            ));
        };

        Ok(evaluate_evaluator(
            summed,
            input.as_slice(),
            evaluation_metadata,
            record_primary_timing,
        ))
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
        if self.explicit_orientation_sum_only {
            if !orientations.is_all() {
                return Err(eyre!(
                    "An explicit orientation-sum evaluator already contains the complete orientation sum and cannot select an individual orientation"
                ));
            }

            // The atom already contains the complete orientation sum, so
            // applying orientation selection again would double count it.
            return Ok(evaluate_evaluator(
                &mut self.single_parametric,
                input.as_slice(),
                evaluation_metadata,
                record_primary_timing,
            ));
        }

        if !orientations.is_all()
            && !matches!(
                settings.general.evaluator_method,
                EvaluatorMethod::SingleParametric
            )
        {
            return Err(eyre!(
                "Runtime evaluator_method={:?} cannot select individual orientations; use SingleParametric for Monte Carlo sampling or runtime filtering over orientations.",
                settings.general.evaluator_method
            ));
        }

        match settings.general.evaluator_method {
            EvaluatorMethod::SingleParametric => Ok(self.evaluate_parametric(
                input,
                orientations,
                evaluation_metadata,
                record_primary_timing,
            )),
            EvaluatorMethod::Iterative => {
                self.evaluate_iterative(input, evaluation_metadata, record_primary_timing)
            }
            EvaluatorMethod::SummedFunctionMap => {
                self.evaluate_summed_fnmap(input, evaluation_metadata, record_primary_timing)
            }
            EvaluatorMethod::Summed => {
                self.evaluate_summed(input, evaluation_metadata, record_primary_timing)
            }
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
    uses_numerator_sampling_scale: bool,
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
    fn source_atoms_use_numerator_sampling_scale(
        atoms: &[Atom],
        source_replacements: &[Replacement],
    ) -> bool {
        atoms.iter().any(|atom| {
            let mut expanded = atom.clone();
            for _ in 0..=source_replacements.len() {
                if expanded.contains_symbol(GS.numerator_sampling_scale) {
                    return true;
                }
                let next = expanded.replace_multiple(source_replacements);
                if next == expanded {
                    return false;
                }
                expanded = next;
            }
            expanded.contains_symbol(GS.numerator_sampling_scale)
        })
    }

    pub(crate) fn uses_numerator_sampling_scale(&self) -> bool {
        self.uses_numerator_sampling_scale
    }

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

    pub(crate) fn activate_symjit(
        &mut self,
        optimization_level: CompilationOptimizationLevel,
    ) -> Result<()> {
        if self.is_eager_only() {
            self.activate_eager_only();
            return Ok(());
        }

        let rational = self
            .rational
            .as_ref()
            .ok_or_else(|| eyre!("Cannot build symjit backend without the rational evaluator"))?;
        // SymJIT 2.21 supports optimization levels up to O2 and cannot compact some complex
        // temporary layouts.
        let evaluator = rational
            .jit_compile::<SymComplex<f64>>(
                JITCompilationSettings::new()
                    .optimization_level(usize::from(optimization_level).min(2) as u8)
                    .with_option("compact", "false"),
            )
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
        let source_replacements = fn_map_entries
            .iter()
            .map(FnMapEntry::replacement)
            .collect::<Vec<_>>();
        let evaluator_replacements = if settings.do_fn_map_replacements {
            source_replacements.as_slice()
        } else {
            &[]
        };

        // Vakint and older Symbolica states represent the imaginary unit as a
        // symbolic constant. The evaluator domain expects the exact complex
        // coefficient used by current Symbolica instead.
        let mut fn_map = fn_map.clone();
        fn_map
            .add_aliases([
                (Atom::var(vakint::symbols::S.cmplx_i), Atom::i()),
                (Atom::var(symbol!("symbolica::𝑖")), Atom::i()),
            ])
            .map_err(|e| eyre!("Failed to register the imaginary-unit constant: {e}"))?;

        for r in evaluator_replacements {
            println!("Reps!!{:#}", r)
        }

        let exprs: Vec<Atom> = atoms
            .into_iter()
            .map(|a| {
                a.replace_multiple(evaluator_replacements)
                    .replace_multiple(evaluator_replacements)
            })
            .collect();
        // Keep this fact even when `store_atom` is disabled: runtime settings
        // must be validated against the finalized evaluator sources, not the
        // generation mode that may or may not have left M in them. Follow only
        // function bodies reachable from an emitted expression; unrelated
        // entries in its shared function map do not make that evaluator use M.
        let uses_numerator_sampling_scale =
            Self::source_atoms_use_numerator_sampling_scale(&exprs, &source_replacements);

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
            uses_numerator_sampling_scale,
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
    pub residue_map_id_start: usize,
    pub orientations_start: usize,
    pub multiplicative_offset: usize,
}

impl<'a, T: FloatLike> InputParams<'a, T> {
    pub(crate) fn set_residue_map_id(&mut self, id: OrientationID) {
        let zero = F(T::from_f64(0.));
        let value = Complex::new_re(zero.from_usize(id.0));
        let index = self.residue_map_id_start * self.multiplicative_offset;
        self.as_mut_slice()[index] = value;
    }

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
    use std::io::Cursor;

    use idenso::{dirac::AGS, representations::Bispinor};
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{network::tags::SPENSO_TAG, structure::representation::RepName};
    use symbolica::{atom::Symbol, parse_lit, state::State};

    use crate::{
        GammaLoopContextContainer, initialisation::test_initialise,
        integrands::process::param_builder::ParamValuePairs, model::Model,
    };

    use super::*;

    fn scalar_value(result: Vec<DualOrNot<Complex<F<f64>>>>) -> Complex<F<f64>> {
        let [DualOrNot::NonDual(value)] = result.as_slice() else {
            panic!("expected one scalar evaluator output")
        };
        value.clone()
    }

    #[test]
    fn exact_residue_map_keys_distinguish_duplicate_and_undirected_orientations() {
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.residue_map_id = ParamValuePairs::default_from_symbol(GS.residue_map_id);
        builder.pairs.orientations = [GS.sign(EdgeIndex(0)), GS.sign(EdgeIndex(1))]
            .into_iter()
            .collect();
        let parameter_count = builder.pairs.update_ranges();
        builder.values = vec![vec![Complex::new_re(F(0.0)); parameter_count]];

        let orientations = TiVec::<OrientationID, _>::from_iter([
            EdgeVec::from_iter([Orientation::Reversed, Orientation::Undirected]),
            EdgeVec::from_iter([Orientation::Default, Orientation::Undirected]),
            EdgeVec::from_iter([Orientation::Default, Orientation::Undirected]),
            EdgeVec::from_iter([Orientation::Undirected, Orientation::Undirected]),
        ]);
        let production_ids = [
            OrientationID(4),
            OrientationID(9),
            OrientationID(12),
            OrientationID(15),
        ];
        // The first entry owns the inverse selector in its reversed physical
        // sector. Selecting a default-oriented key must discard that complete
        // entry before resolving the physical selector; otherwise the inactive
        // inverse becomes `0^-1` and contaminates the selected expression.
        let reversed_selector_inverse = GS.sign_theta(-GS.sign(EdgeIndex(0))).pow(-1);
        let atom = production_ids[0].atom() * Atom::num(2) * reversed_selector_inverse
            + production_ids[1].atom() * Atom::num(3)
            + production_ids[2].atom() * Atom::num(5)
            + production_ids[3].atom() * Atom::num(7);
        let mut evaluator_settings = EvaluatorSettings::default();
        evaluator_settings.summed = true;
        evaluator_settings.summed_function_map = true;
        let (mut stack, _) = EvaluatorStack::new_with_timings(
            &[atom],
            &builder,
            &orientations.raw,
            &production_ids,
            None,
            &evaluator_settings,
        )
        .unwrap();

        let make_input = || InputParams {
            values: SliceMut::Owned(vec![Complex::new_re(F(0.0)); parameter_count]),
            residue_map_id_start: builder.pairs.residue_map_id.value_range.start,
            orientations_start: builder.pairs.orientations.value_range.start,
            multiplicative_offset: 1,
        };
        let mut metadata = EvaluationMetaData::new_empty();
        for (runtime_id, expected) in [2.0, 3.0, 5.0, 7.0].into_iter().enumerate() {
            let actual = stack.evaluate_parametric(
                make_input(),
                SingleOrAllOrientations::Single {
                    orientation: &orientations[OrientationID(runtime_id)],
                    id: OrientationID(runtime_id),
                },
                &mut metadata,
                false,
            );
            assert_eq!(scalar_value(actual), Complex::new_re(F(expected)));
        }

        let filter = SubSet::full(orientations.len());
        let all = SingleOrAllOrientations::All {
            all: &orientations,
            filter: &filter,
        };
        assert_eq!(
            scalar_value(stack.evaluate_parametric(make_input(), all, &mut metadata, false)),
            Complex::new_re(F(17.0))
        );

        for method in [
            EvaluatorMethod::Iterative,
            EvaluatorMethod::SummedFunctionMap,
            EvaluatorMethod::Summed,
        ] {
            let mut runtime_settings = RuntimeSettings::default();
            runtime_settings.general.evaluator_method = method;
            assert_eq!(
                scalar_value(
                    stack
                        .evaluate(make_input(), all, &runtime_settings, &mut metadata, false,)
                        .unwrap()
                ),
                Complex::new_re(F(17.0))
            );
        }
    }

    #[test]
    fn deferred_explicit_sum_lowering_is_local_and_matches_materialized_value() {
        let builder = ParamBuilder::new_empty();
        let settings = EvaluatorSettings::default();
        let bodies = [Atom::num(2), Atom::num(3)];
        let compact = function!(GS.projected_cff_sum, 0) + function!(GS.projected_cff_sum, 1);
        let (mut deferred, _) = EvaluatorStack::new_deferred_explicit_sum_with_timings(
            &compact, &bodies, &builder, None, &settings,
        )
        .unwrap();
        assert!(builder.reps.is_empty());

        let (mut materialized, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &[Atom::num(5)],
            &builder,
            None,
            &settings,
        )
        .unwrap();
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
            &mut deferred.single_parametric,
        )(&[]);
        let expected = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
            &mut materialized.single_parametric,
        )(&[]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn parameterless_deferred_body_preserves_multiparameter_hyperdual_derivatives() {
        let x = Atom::var(symbol!("evaluator_test::deferred_dual_x"));
        let y = Atom::var(symbol!("evaluator_test::deferred_dual_y"));
        let body = &x * &x + &x * &y + Atom::num(3) * &y;
        let call = function!(GS.projected_cff_sum, 0);
        let entry = FnMapEntry {
            lhs: call.clone(),
            rhs: body.clone(),
            args: Vec::new(),
            tags: vec![Atom::num(0)],
        };
        let mut function_map = FunctionMap::default();
        function_map
            .add_tagged_function(
                GS.projected_cff_sum,
                vec![Atom::num(0)],
                Vec::<Indeterminate>::new(),
                body.clone(),
            )
            .unwrap();
        let dual_shape = Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1));
        let settings = EvaluatorSettings::default();
        let mut deferred = GenericEvaluator::new_from_raw_params(
            [call],
            &[x.clone(), y.clone()],
            &function_map,
            vec![entry],
            OptimizationSettings::default(),
            dual_shape.clone(),
            &settings,
        )
        .unwrap();
        let mut materialized = GenericEvaluator::new_from_raw_params(
            [body],
            &[x, y],
            &FunctionMap::default(),
            vec![],
            OptimizationSettings::default(),
            dual_shape,
            &settings,
        )
        .unwrap();
        let input = [
            Complex::new_re(F(2.0)),
            Complex::new_re(F(7.0)),
            Complex::new_re(F(5.0)),
            Complex::new_re(F(11.0)),
        ];
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut deferred)(&input);
        let expected = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut materialized)(&input);

        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("deferred evaluator did not return the requested dual output")
        };
        let [DualOrNot::Dual(expected)] = expected.as_slice() else {
            panic!("materialized evaluator did not return the requested dual output")
        };
        assert_eq!(actual.values, expected.values);
        assert_eq!(actual.values[0], Complex::new_re(F(29.0)));
        assert_eq!(actual.values[1], Complex::new_re(F(118.0)));
    }

    #[test]
    fn deferred_explicit_sum_preprocesses_tensor_function_bodies() {
        test_initialise().unwrap();
        let builder = ParamBuilder::new_empty();
        let settings = EvaluatorSettings {
            do_algebra: true,
            ..Default::default()
        };
        let body = Bispinor {}.new_rep(4).g(9, 9);
        let compact = function!(GS.projected_cff_sum, 0);
        let (mut deferred, _) = EvaluatorStack::new_deferred_explicit_sum_with_timings(
            &compact,
            &[body],
            &builder,
            None,
            &settings,
        )
        .unwrap();
        let (mut expected, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &[Atom::num(4)],
            &builder,
            None,
            &settings,
        )
        .unwrap();

        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
            &mut deferred.single_parametric,
        )(&[]);
        let expected = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
            &mut expected.single_parametric,
        )(&[]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn evaluator_preprocess_scalarizes_abstract_minkowski_contractions_without_algebra() {
        test_initialise().unwrap();
        let abstract_index = parse_lit!(spenso::mink(4, 1));
        let edge = EdgeIndex(7);
        let mapped_momentum = GS.emr_vec_index(edge, abstract_index.as_view())
            + GS.ose(edge) * GS.energy_delta(abstract_index.as_view());
        let numerator = &mapped_momentum * &mapped_momentum;
        let settings = EvaluatorSettings::default();
        let scalar = EvaluatorStack::preprocess_atom(&numerator, 0, &settings).unwrap();
        let algebraic_scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &EvaluatorSettings {
                do_algebra: true,
                ..settings
            },
        )
        .unwrap();

        assert_eq!(scalar, algebraic_scalar);
    }

    #[test]
    fn evaluator_preprocess_scalarizes_temporal_emr_gamma_slash() {
        test_initialise().unwrap();
        let momentum_index = parse_lit!(spenso::mink(4, 1));
        let compact_minkowski = parse_lit!(spenso::mink(4));
        let left_index = parse_lit!(spenso::bis(4, 2));
        let right_index = parse_lit!(spenso::bis(4, 3));
        let edge = EdgeIndex(7);
        let temporal_momentum = GS.ose(edge) * GS.energy_delta(compact_minkowski.as_view());
        let gamma_factor = FunctionBuilder::new(AGS.gamma)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .add_arg(&temporal_momentum)
            .finish();
        let gamma_slash = FunctionBuilder::new(SPENSO_TAG.chain)
            .add_arg(&left_index)
            .add_arg(&right_index)
            .add_arg(gamma_factor)
            .finish();
        let numerator = function!(GS.vbar, 0, left_index)
            * gamma_slash
            * function!(GS.u, 1, right_index.clone());
        let explicit_numerator = function!(GS.vbar, 0, parse_lit!(spenso::bis(4, 2)))
            * GS.ose(edge)
            * GS.energy_delta(momentum_index.as_view())
            * parse_lit!(spenso::gamma(
                spenso::bis(4, 2),
                spenso::bis(4, 3),
                spenso::mink(4, 1)
            ))
            * function!(GS.u, 1, right_index);
        let settings = EvaluatorSettings::default();
        let scalar = EvaluatorStack::preprocess_atom(&numerator, 0, &settings).unwrap();
        let explicit_scalar =
            EvaluatorStack::preprocess_atom(&explicit_numerator, 0, &settings).unwrap();
        let algebraic_scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &EvaluatorSettings {
                do_algebra: true,
                ..settings
            },
        )
        .unwrap();

        assert!(
            !scalar.contains_symbol(AGS.gamma),
            "gamma slash was componentized before its temporal delta contracted: {scalar}"
        );
        assert_eq!(scalar, explicit_scalar);
        assert!(!algebraic_scalar.contains_symbol(AGS.gamma));
    }

    #[test]
    fn evaluator_archives_sampling_scale_usage_without_stored_atoms() {
        test_initialise().unwrap();
        let source_function = symbol!("evaluator_test::sampled_source");
        let source_call = FunctionBuilder::new(source_function).finish();
        let scale = Atom::var(GS.numerator_sampling_scale);
        let source_body = &scale + Atom::num(1);
        let mut function_map = FunctionMap::default();
        function_map
            .add_function(
                source_function,
                Vec::<Indeterminate>::new(),
                source_body.clone(),
            )
            .unwrap();
        let source_entry = FnMapEntry {
            lhs: source_call.clone(),
            rhs: source_body,
            args: Vec::new(),
            tags: Vec::new(),
        };
        let evaluator = GenericEvaluator::new_from_raw_params(
            [source_call.clone()],
            std::slice::from_ref(&scale),
            &function_map,
            vec![source_entry.clone()],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        assert!(evaluator.exprs.is_none());
        assert!(evaluator.uses_numerator_sampling_scale());

        let encoded = bincode::encode_to_vec(&evaluator, bincode::config::standard()).unwrap();
        let mut state = Vec::new();
        State::export(&mut state).unwrap();
        let state_map = State::import(&mut Cursor::new(state), None).unwrap();
        let model = Model::default();
        let (decoded, _): (GenericEvaluator, _) = bincode::decode_from_slice_with_context(
            &encoded,
            bincode::config::standard(),
            GammaLoopContextContainer {
                state_map: &state_map,
                model: &model,
            },
        )
        .unwrap();

        assert!(decoded.exprs.is_none());
        assert!(decoded.uses_numerator_sampling_scale());

        let independent = GenericEvaluator::new_from_raw_params(
            [Atom::num(1)],
            std::slice::from_ref(&scale),
            &function_map,
            vec![source_entry],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        assert!(!independent.uses_numerator_sampling_scale());
    }

    #[test]
    fn evaluator_treats_vakint_imaginary_symbol_as_exact_constant() {
        test_initialise().unwrap();
        for (index, symbolic_i) in [
            Atom::var(vakint::symbols::S.cmplx_i),
            Atom::var(symbol!("symbolica::𝑖")),
        ]
        .into_iter()
        .enumerate()
        {
            assert_ne!(symbolic_i, Atom::i());
            let source_symbol = symbol!(&format!("evaluator_test::imaginary_source_{index}"));
            let source_call = FunctionBuilder::new(source_symbol).finish();
            let mut function_map = FunctionMap::default();
            function_map
                .add_function(
                    source_symbol,
                    Vec::<Indeterminate>::new(),
                    symbolic_i.clone(),
                )
                .unwrap();

            let mut evaluator = GenericEvaluator::new_from_raw_params(
                [symbolic_i + source_call],
                &[],
                &function_map,
                vec![],
                OptimizationSettings::default(),
                None,
                &EvaluatorSettings::default(),
            )
            .unwrap();
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&[]);

            assert_eq!(actual, Complex::new(F(0.0), F(2.0)));
        }
    }

    #[test]
    fn pi_eval() {
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            [Atom::var(Symbol::PI)],
            &[],
            &FunctionMap::default(),
            vec![],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        // The built-in constant must be created in the active numeric domain;
        // otherwise precision escalation merely pads an f64 approximation.
        let actual = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&[]);
        let zero = F(ArbPrec::default());
        let expected = Complex::new_re(zero.PI());
        assert!((actual - expected).norm().re <= zero.epsilon());
    }
}
