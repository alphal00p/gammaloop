mod parametric_rows;
pub(crate) use parametric_rows::ParametricResidueRows;

use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use eyre::{Context, eyre};
use idenso::{
    IndexTooling,
    color::{ColorSimplifier, ColorSimplifySettings},
    dirac::GammaSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubSetIter, SubSetLike, subset::SubSet},
    typed_vec::IndexLike,
};

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::{
    algebra::{
        algebraic_traits::RefOne,
        complex::{Complex, symbolica_traits::CompiledComplexEvaluatorSpenso},
    },
    iterators::IteratableTensor,
    network::{
        DEFAULT_EXACT_JOIN_LIMIT, ExecutionResult, MAX_EAGER_TENSOR_SUM_BYTES, MinIntermediateCost,
        MinResultRank, MinResultRankWith, PAIR_SCORE_ATOM_AWARE, PAIR_SCORE_ENTRY_AWARE,
        PAIR_SCORE_RESULT_RANK_ONLY, ScalarAliases, Sequential, SequentialExtract, SequentialRef,
        SmallestDegree,
        graph::{NetworkLeaf, NetworkNode, NetworkOp},
        parsing::{AtomStructureExt, ShadowedStructure, StrictTensorFilter},
        store::{TensorScalarStore, TensorScalarStoreMapping},
        tags::SPENSO_TAG,
    },
    shadowing::{
        ANTISYM, CYCLIC, Collectable, SYM,
        symbolica_utils::{LogPrint, SpensoPrintSettings},
    },
    structure::{
        HasStructure, TensorStructure, ToSymbolic,
        abstract_index::AIND_SYMBOLS,
        concrete_index::DualConciousIndex,
        representation::{LibraryRep, RepName},
        slot::{DualSlotTo, IsAbstractSlot, Slot},
    },
    tensors::{
        data::{SparseOrDense, StorageTensor},
        parametric::ParamTensor,
    },
};
use std::{
    collections::{HashMap, HashSet},
    ops::Deref,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};
use std::{mem::transmute, ops::Neg, path::Path};
use symbolica::{
    domains::{dual::HyperDual, float::Complex as SymComplex, rational::Fraction},
    evaluate::{InliningPolicy, JITCompiledEvaluator},
    prelude::*,
};
use tracing::{debug, instrument};
use typed_index_collections::TiSlice;

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
    numerator::{ParsingNet, aind::Aind, symbolica_ext::NumeratorAtomExt},
    processes::{
        ContractionMode, EvaluatorBuildTimings, EvaluatorSettings, ExecutionMode,
        TensorNetworkContractionOrder,
    },
    settings::{
        RuntimeSettings,
        global::{CompilationOptionsSnapshot, FrozenCompilationMode},
    },
    utils::{
        ArbPrec, F, FUN_LIB, FloatLike, GS, Length, RuntimeCache, SamplingFloat, TENSORLIB, W_,
        f128,
        hyperdual_utils::{DualOrNot, new_from_values},
    },
    uv::ParametricIntegrandTerm,
};

type ParsingTensorMap<'a> = dyn Fn(ParamTensor<ShadowedStructure<Aind>>) -> Result<ParamTensor<ShadowedStructure<Aind>>>
    + 'a;

use super::{
    ParamBuilder,
    param_builder::{ThresholdParams, UpdateAndGetParams},
};

const NETWORK_SCALAR_ALIAS_MIN_BYTES: usize = 4096;
static NETWORK_SCALAR_ALIAS_SCOPE: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone, Copy)]
pub enum SingleOrAllOrientations<'a, OID> {
    Single {
        orientation: &'a EdgeVec<Orientation>,
        id: OID,
    },
    All {
        all: &'a TiSlice<OID, EdgeVec<Orientation>>,
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
) -> Complex<F<T>> {
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
) -> Vec<DualOrNot<Complex<F<T>>>> {
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
        all: &'a TiSlice<OID, EdgeVec<Orientation>>,
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
#[cfg_attr(feature = "python", pyo3::pyclass(from_py_object, get_all, set_all))]
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

type SingleParametricAuxiliary = (Vec<Atom>, GenericEvaluator);

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct EvaluatorStack {
    pub(crate) explicit_orientation_sum_only: bool,
    /// Internal coefficient/dispatch rows, independent of physical orientations.
    pub(crate) parametric_rows: Option<ParametricResidueRows>,
    /// Original generalized-3D-representation map key for each dense runtime
    /// orientation channel. Physical edge directions alone do not identify a
    /// raised-energy/contact residue map.
    production_orientation_ids: Vec<OrientationID>,
    pub single_parametric: GenericEvaluator,
    /// Carrier ports and their shared output table, evaluated once per sample.
    /// The row catalogue maps ports to outputs independently of their identities.
    pub(crate) single_parametric_auxiliary: Option<SingleParametricAuxiliary>,
    pub iterative: Option<(GenericEvaluator, usize)>,
    // pub iterative_function_map: Option<GenericEvaluator>,
    pub summed_function_map: Option<GenericEvaluator>,
    pub summed: Option<GenericEvaluator>,
}

impl EvaluatorStack {
    /// Share repeated scalar factors after tensor contraction and factor collection.
    /// Visit the root and retained definitions without unfolding their DAG.
    /// Guards stay at their original use sites inside lazy branches. No subsequent
    /// algebraic collection may move an aliased inverse across them.
    fn share_scalar_factors(atom: AliasedAtom, builder: &ParamBuilder) -> AliasedAtom {
        let started = std::time::Instant::now();
        let input_bytes = atom.get_byte_size();
        let alias_symbol = loop {
            let scope = NETWORK_SCALAR_ALIAS_SCOPE.fetch_add(1, Ordering::Relaxed);
            let name = format!("gammalooprs::evaluator_factor_{scope}");
            if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
            }) {
                break symbol!(name.as_str());
            }
        };
        let (root, aliases) = atom.into_inner_with_aliases();
        // Exact parameter/alias keys and entire functions/powers are opaque.
        // Recursive extraction would also rewrite arguments inside definitions,
        // breaking parameter and function-tag lookup before aliases are resolved.
        let boundaries = (&builder.pairs)
            .into_iter()
            .flat_map(|pair| &pair.params)
            .chain(aliases.keys())
            .map(|key| key.as_view())
            .collect::<HashSet<_>>();
        let max_key_bytes = boundaries
            .iter()
            .map(|key| key.get_data().len())
            .max()
            .unwrap_or(0);
        let is_boundary = |view: AtomView<'_>| {
            view.get_data().len() <= max_key_bytes && boundaries.contains(&view)
        };
        let has_algebra = |view: AtomView<'_>| {
            let mut found = false;
            view.visitor(&mut |part| {
                found |= matches!(part, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_));
                !found
            });
            found
        };
        let sensitive_keys = boundaries
            .iter()
            .copied()
            .filter(|key| has_algebra(*key))
            .collect::<HashSet<_>>();
        // FunctionMap enforces one tag count per head, including registrations
        // absent from reps. Unknown heads conservatively treat all args as tags.
        let tag_counts = builder
            .reps
            .iter()
            .filter_map(|entry| {
                entry
                    .lhs
                    .as_fun_view()
                    .map(|call| (call.get_symbol(), entry.tags.len()))
            })
            .collect::<HashMap<_, _>>();
        let is_factor = |view: AtomView<'_>| {
            matches!(view, AtomView::Fun(_) | AtomView::Pow(_))
                && !matches!(view, AtomView::Fun(call) if call.get_symbol() == Symbol::IF)
        };
        let mut counts = HashMap::<AtomView<'_>, usize>::new();
        for body in std::iter::once(&root).chain(aliases.values()) {
            body.visitor(&mut |view| {
                if is_boundary(view) {
                    return false;
                }
                if is_factor(view) {
                    *counts.entry(view).or_default() += 1;
                    return false;
                }
                true
            });
        }
        // Sort for stable handles, and store every selected body unchanged.
        let mut candidates = counts.into_iter().collect::<Vec<_>>();
        candidates.sort_unstable_by_key(|(view, _)| *view);
        let mut replacements = HashMap::new();
        for (view, count) in candidates {
            if count < 2
                || [
                    OrientationID::symbol(),
                    GS.theta,
                    GS.orientation_delta,
                    Symbol::IF,
                ]
                .into_iter()
                .any(|guard| view.contains_symbol(guard))
            {
                continue;
            }
            // Symbolica also Horners alias bodies. Preserve exact parameter/tag
            // keys, while allowing arithmetic in ordinary arguments, particularly
            // the OnShellEnergy invariant. Atomic momentum/index labels are safe.
            let mut changes_key = false;
            view.visitor(&mut |part| {
                changes_key |=
                    part.get_data().len() <= max_key_bytes && sensitive_keys.contains(&part);
                if let AtomView::Fun(call) = part {
                    let head = call.get_symbol();
                    let tags = tag_counts.get(&head).copied().unwrap_or_else(|| {
                        if head.is_builtin() {
                            0
                        } else {
                            call.get_nargs()
                        }
                    });
                    changes_key |= call.iter().take(tags).any(has_algebra);
                }
                !changes_key
            });
            if changes_key {
                continue;
            }
            let handle = function!(alias_symbol, replacements.len());
            // Include the stored definition and its key, not just root references.
            if view.get_byte_size().saturating_mul(count - 1)
                > handle.as_view().get_byte_size().saturating_mul(count + 1)
            {
                replacements.insert(view, handle);
            }
        }
        let rewrite = |body: &Atom| {
            body.replace_map(|view, _, out| {
                if is_boundary(view) || is_factor(view) {
                    out.set_from_view(
                        &replacements
                            .get(&view)
                            .map_or(view, |handle| handle.as_view()),
                    );
                }
            })
        };
        let mut shared = AliasedAtom::from(rewrite(&root));
        for (alias, body) in &aliases {
            shared.register_alias(alias.clone(), rewrite(body));
        }
        let aliases_created = replacements.len();
        for (body, handle) in replacements {
            shared.register_alias(handle, body.to_owned());
        }
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_stack_share_scalar_factors_done",
            input_bytes,
            result_bytes = shared.get_byte_size(),
            aliases_created,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Shared repeated scalar factors without hiding branch guards"
        );
        shared
    }

    fn parametrize_residue_map_selectors<'a>(
        atom: impl Into<AtomOrView<'a>>,
        selected_id: Atom,
    ) -> Atom {
        let atom = atom.into();
        if !atom.contains_symbol(OrientationID::symbol()) {
            return atom.into_owned();
        }
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
        let mut count = 1 + usize::from(self.single_parametric_auxiliary.is_some());
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
    fn new_single_parametric(
        parametric_atoms: Vec<AliasedAtom>,
        param_builder: &ParamBuilder,
        physical_functions: &[FnMapEntry],
        mut parametric_rows: Option<&mut ParametricResidueRows>,
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(GenericEvaluator, Option<SingleParametricAuxiliary>)> {
        let _progress_guard = crate::processes::enter_detailed_progress_span(
            "Generating Single Parametric Evaluator",
        );
        let opt_settings = settings.optimization_settings();

        let mut atoms = parametric_atoms
            .into_iter()
            .map(|atom| {
                let map = |atom| {
                    GS.collect_orientation_if(Self::parametrize_residue_map_selectors(
                        atom,
                        Atom::var(GS.residue_map_id),
                    ))
                };
                let (root, aliases) = atom.into_inner_with_aliases();
                let mut mapped = AliasedAtom::from(map(root));
                for (alias, body) in aliases {
                    mapped.register_alias(alias, map(body));
                }
                mapped
            })
            .collect::<Vec<_>>();
        let mut params = (&param_builder.pairs)
            .into_iter()
            .flat_map(|pair| pair.params.clone())
            .collect::<Vec<_>>();
        if let Some(rows) = parametric_rows.as_deref_mut() {
            // A catalogue may build several evaluator policies; each starts
            // from the original coefficient interface, including policies
            // which deliberately keep all operations in the row program.
            rows.input_columns = (0..rows.parameters.len()).collect();
            rows.auxiliary_parameters.clear();
            rows.auxiliary_values.clear();
            rows.auxiliary_shared_outputs = 0;
            rows.auxiliary_gathers = vec![Vec::new(); rows.rows.len()];
        }
        // Extract only calls that survive the original symbolic function
        // boundary. Earlier unfolding (including retained-dual lowering) can
        // cancel singular intermediate expressions and must remain unchanged.
        let single_parametric_auxiliary = if let Some(rows) = parametric_rows.as_deref_mut()
            && rows.rows.len() > 1
            && !settings.do_fn_map_replacements
            && !(dual_shape.is_some()
                && param_builder
                    .reps
                    .iter()
                    .any(|entry| entry.inlining == InliningPolicy::Never))
        {
            let port_symbol = loop {
                let scope = NETWORK_SCALAR_ALIAS_SCOPE.fetch_add(1, Ordering::Relaxed);
                let name = format!("gammalooprs::evaluator_physical_{scope}");
                if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                    existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
                }) {
                    break symbol!(name.as_str());
                }
            };
            rows.extract_physical_scalars(
                &mut atoms,
                &params,
                &param_builder.reps,
                physical_functions.len(),
                port_symbol,
            )
            .map(|(ports, roots, definitions)| {
                crate::debug_tags!(#generation, #profile, #compile, #summary;
                    stage = "evaluator_single_parametric_physical_scalars",
                    auxiliary_outputs = roots.len(), carrier_ports = ports.len(),
                    coefficient_rows = rows.rows.len(),
                    "Prepared shared scalar and numerator values for the complete residue sum"
                );
                let mut auxiliary_params = params.clone();
                auxiliary_params.extend(rows.auxiliary_parameters.iter().cloned());
                let dual_config = dual_shape.clone().map(|shape| {
                    let zeros = (params.len()..auxiliary_params.len())
                        .flat_map(|parameter| {
                            (1..shape.len()).map(move |component| (parameter, component))
                        })
                        .collect();
                    (shape, zeros)
                });
                let evaluator = GenericEvaluator::new_from_raw_params(
                    roots,
                    &auxiliary_params,
                    &param_builder.fn_map,
                    definitions,
                    opt_settings.clone(),
                    dual_config,
                    settings,
                )?;
                Ok::<_, eyre::Report>((ports, evaluator))
            })
            .transpose()?
        } else {
            None
        };
        if let Some((ports, _)) = &single_parametric_auxiliary {
            params.extend(ports.iter().cloned());
        }
        // Auxiliary ports carry full physical jets. Only exact coefficient and
        // selector ports have statically zero derivative components.
        let row_start = params.len();
        if let Some(rows) = parametric_rows {
            params.extend(
                rows.input_columns
                    .iter()
                    .map(|&column| rows.parameters[column].clone()),
            );
        }
        let dual_config = dual_shape.clone().map(|shape| {
            let zeros = (row_start..params.len())
                .flat_map(|parameter| (1..shape.len()).map(move |component| (parameter, component)))
                .collect();
            (shape, zeros)
        });
        let evaluator = GenericEvaluator::new_from_raw_params(
            atoms,
            &params,
            &param_builder.fn_map,
            param_builder.reps.clone(),
            opt_settings.clone(),
            dual_config,
            settings,
        )?;
        Ok((evaluator, single_parametric_auxiliary))
    }
    #[instrument(skip_all)]
    fn new_iterative(
        parametric_atom: &[AliasedAtom],
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
                            // Each orientation is a separate output, so handles
                            // keep their names. Only bodies still reached after
                            // selection are specialized and retained.
                            let select = Self::residue_map_selector(orientation, *production_id);
                            let selected = Self::specialize_aliases(
                                atom.get_root(),
                                atom,
                                |atom| GS.collect_orientation_if(select(atom)),
                                None,
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
    fn new_summed_function_map(
        atoms: &[AliasedAtom],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        production_orientation_ids: &[OrientationID],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard = crate::processes::enter_detailed_progress_span(
            "Generating Summed Function Map Evaluator",
        );
        if orientations.is_empty() {
            return Err(eyre!(
                "summed function-map evaluator requires at least one residue-map orientation"
            ));
        }
        // Select each residue-map entry symbolically, as the summed evaluator
        // does, and keep it behind the argument-free alias integrand(i, map_id).
        // Passing the selector and signs as call arguments bound them only
        // during linearization, so no downstream call simplified and the
        // program grew by orders of magnitude before common-subexpression
        // elimination.
        let sum = atoms
            .iter()
            .enumerate()
            .map(|(i, atom)| {
                let mut summed = AliasedAtom::default();
                let mut calls = Vec::with_capacity(orientations.len());
                for (orientation, production_id) in
                    orientations.iter().zip(production_orientation_ids)
                {
                    let (root, aliases) = Self::specialize_aliases(
                        atom.get_root(),
                        atom,
                        Self::residue_map_selector(orientation, *production_id),
                        Some(production_id.0),
                    )
                    .into_inner_with_aliases();
                    let call = FunctionBuilder::new(GS.integrand)
                        .add_arg(i)
                        .add_arg(production_id.0)
                        .finish();
                    summed.register_alias(call.clone(), root);
                    for (alias, body) in aliases {
                        summed.register_alias(alias, body);
                    }
                    calls.push(call);
                }
                summed.map_root(|_| Atom::add_many(calls))
            })
            .collect::<Vec<_>>();
        GenericEvaluator::new_from_builder(
            sum,
            param_builder,
            dual_shape.clone(),
            settings.optimization_settings(),
            settings,
        )
    }

    /// Bind every runtime orientation input after discarding inactive residue maps.
    fn residue_map_selector(
        orientation: &EdgeVec<Orientation>,
        production_id: OrientationID,
    ) -> impl Fn(&Atom) -> Atom + '_ {
        // GraphOrientation::select is a partial selector: undirected edges
        // are wildcards. A compiled evaluator instead binds those inputs to
        // zero, just like InputParams::set_orientation_values_impl.
        let undirected = orientation
            .iter()
            .filter(|(_, direction)| **direction == Orientation::Undirected)
            .map(|(edge, _)| Replacement::new(GS.sign(edge).to_pattern(), Atom::Zero))
            .collect::<Vec<_>>();
        move |atom| {
            let selected = orientation.select(production_id.select(atom));
            if undirected.is_empty() {
                selected
            } else {
                selected.replace_multiple(&undirected)
            }
        }
    }

    /// Specialize the selected alias DAG from its root. Alias handles remain
    /// opaque, and dependencies are discovered only after branch selection, so
    /// inactive scalar bodies never enter orientation or coefficient substitution.
    /// An absent scope preserves the owner's keys when only reachability changes.
    fn specialize_aliases(
        root: &Atom,
        atom: &AliasedAtom,
        specialize: impl Fn(&Atom) -> Atom,
        scope_id: Option<usize>,
    ) -> AliasedAtom {
        let aliases = atom.get_aliases();
        let max_handle_bytes = aliases
            .keys()
            .map(|handle| handle.as_view().get_data().len())
            .max()
            .unwrap_or(0);
        let mut renames = HashMap::<Atom, Atom>::new();
        let mut pending = Vec::new();
        let select =
            |expression: &Atom, renames: &mut HashMap<Atom, Atom>, pending: &mut Vec<Atom>| {
                specialize(expression).replace_map(|view, _, out| {
                    if view.get_data().len() <= max_handle_bytes
                        && let Some((alias, _)) = aliases.get_key_value(view.get_data())
                    {
                        let scoped = renames.entry(alias.clone()).or_insert_with(|| {
                            // Register before traversing the body: shared
                            // dependencies are visited once, including cycles.
                            pending.push(alias.clone());
                            scope_id.map_or_else(
                                || alias.clone(),
                                |scope_id| {
                                    let function = alias.as_fun_view().unwrap();
                                    FunctionBuilder::new(function.get_symbol())
                                        .add_args(function.iter())
                                        .add_arg(scope_id)
                                        .finish()
                                },
                            )
                        });
                        out.set_from_view(&scoped.as_view());
                    }
                })
            };
        let mut selected = AliasedAtom::from(select(root, &mut renames, &mut pending));
        while let Some(alias) = pending.pop() {
            let body = select(&aliases[&alias], &mut renames, &mut pending);
            selected.register_alias(renames[&alias].clone(), body);
        }
        selected
    }

    #[instrument(skip_all)]
    fn new_summed(
        atoms: &[AliasedAtom],
        param_builder: &ParamBuilder,
        orientations: &[EdgeVec<Orientation>],
        production_orientation_ids: &[OrientationID],
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard =
            crate::processes::enter_detailed_progress_span("Generating Summed Evaluator");
        // Atom addition preserves each selected branch's product factors and
        // cancels opposite branches before numerical evaluator construction.
        // Shared numerator bodies remain in the function map.
        let started = std::time::Instant::now();
        let mut selected_alias_bodies = 0;
        let sum = atoms
            .iter()
            .map(|atom| {
                orientations
                    .iter()
                    .zip(production_orientation_ids)
                    .map(|(orientation, production_id)| {
                        let selected = Self::specialize_aliases(
                            atom.get_root(),
                            atom,
                            Self::residue_map_selector(orientation, *production_id),
                            Some(production_id.0),
                        );
                        selected_alias_bodies += selected.get_aliases().len();
                        selected
                    })
                    .fold(AliasedAtom::default(), |acc, atom| {
                        acc.try_add(&atom).unwrap()
                    })
            })
            .collect::<Vec<_>>();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_summed_orientation_selection_done",
            atom_count = atoms.len(), orientation_count = orientations.len(),
            candidate_alias_bodies = atoms.iter().map(|atom| atom.get_aliases().len()).sum::<usize>()
                * orientations.len().min(production_orientation_ids.len()),
            selected_alias_bodies,
            summed_aliases = sum.iter().map(|atom| atom.get_aliases().len()).sum::<usize>(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Selected reachable scalar aliases before evaluator preparation"
        );

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
            &[],
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
        numerator_definitions: &[Arc<FnMapEntry>],
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        let atoms = atoms
            .iter()
            .map(Self::sum_residue_map_selectors)
            .collect::<Vec<_>>();
        let (mut stack, timings) = Self::new_with_timings(
            &atoms,
            param_builder,
            numerator_definitions,
            &[EdgeVec::from_iter(std::iter::empty::<Orientation>())],
            &[OrientationID(0)],
            dual_shape,
            settings,
        )?;
        stack.explicit_orientation_sum_only = true;
        Ok((stack, timings))
    }

    pub(crate) fn from_integrand_with_timings(
        integrand: &Atom,
        param_builder: &ParamBuilder,
        numerator_definitions: &[Arc<FnMapEntry>],
        parametric_terms: &[ParametricIntegrandTerm],
        orientation_catalog: Option<(&[EdgeVec<Orientation>], &[OrientationID])>,
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        if !parametric_terms.is_empty() {
            let (atoms, guards, rows) = ParametricResidueRows::lower(integrand, parametric_terms)?;
            let dummy_orientations = [EdgeVec::from_iter(std::iter::empty::<Orientation>())];
            let dummy_ids = [OrientationID(0)];
            let (orientations, ids) =
                orientation_catalog.unwrap_or((&dummy_orientations, &dummy_ids));
            let atoms = if orientation_catalog.is_none() {
                atoms.iter().map(Self::sum_residue_map_selectors).collect()
            } else {
                atoms
            };
            let definitions = numerator_definitions
                .iter()
                .cloned()
                .chain(
                    parametric_terms
                        .iter()
                        .map(|term| Arc::clone(&term.numerator)),
                )
                .collect::<Vec<_>>();
            let (mut stack, timings) = Self::new_with_rows(
                &atoms,
                param_builder,
                &definitions,
                orientations,
                ids,
                dual_shape,
                settings,
                Some((rows, guards)),
            )?;
            stack.explicit_orientation_sum_only = orientation_catalog.is_none();
            return Ok((stack, timings));
        }
        match orientation_catalog {
            Some((orientations, production_orientation_ids)) => Self::new_with_timings(
                std::slice::from_ref(integrand),
                param_builder,
                numerator_definitions,
                orientations,
                production_orientation_ids,
                dual_shape,
                settings,
            ),
            None => Self::new_explicit_sum_with_timings(
                std::slice::from_ref(integrand),
                param_builder,
                numerator_definitions,
                dual_shape,
                settings,
            ),
        }
    }

    fn preprocess_atom<A: AtomCore>(
        a: &A,
        atom_index: usize,
        settings: &EvaluatorSettings,
        alias_symbol: Symbol,
        tensor_map: Option<&ParsingTensorMap<'_>>,
    ) -> Result<AliasedAtom> {
        Self::preprocess_tensor(
            a,
            atom_index,
            settings,
            tensor_map,
            |_| false,
            |net, scalar_aliases, term_index| {
                let root = match net.result_scalar()? {
                    ExecutionResult::One => Atom::num(1),
                    ExecutionResult::Zero => Atom::Zero,
                    ExecutionResult::Val(value) => value.into_owned(),
                }
                .collect_compact_factors();
                let started = std::time::Instant::now();
                let input_bytes = root.as_view().get_byte_size();
                let (root, aliases) = net
                    .aliased_atom(scalar_aliases, root)
                    .into_inner_with_aliases();
                let mut handles = aliases.keys().cloned().collect::<Vec<_>>();
                handles.sort();
                let renames = handles
                    .into_iter()
                    .enumerate()
                    .map(|(index, handle)| {
                        (
                            handle,
                            function!(alias_symbol, atom_index, term_index, index),
                        )
                    })
                    .collect::<HashMap<_, _>>();
                // Larger subtrees cannot match a handle; avoid hashing
                // their complete serialized contents during renaming.
                let max_handle_bytes = renames
                    .keys()
                    .map(|handle| handle.as_view().get_data().len())
                    .max()
                    .unwrap_or(0);
                let map = |atom: Atom| {
                    if renames.is_empty() {
                        return atom;
                    }
                    atom.replace_map(|view, _, out| {
                        if view.get_data().len() <= max_handle_bytes
                            && let Some(scoped) = renames.get(view.get_data())
                        {
                            out.set_from_view(&scoped.as_view());
                        }
                    })
                };
                let mut retained = AliasedAtom::from(map(root));
                for (alias, body) in aliases {
                    // Spenso may store the entire scalar in a definition, leaving
                    // only a handle in the root. Collect each stored body before
                    // factor sharing, while singular factors are still visible.
                    let body = if body == alias {
                        body
                    } else {
                        map(body.collect_compact_factors())
                    };
                    retained.register_alias(renames[&alias].clone(), body);
                }
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_alias_retention_done",
                    atom_index,
                    term_index,
                    input_bytes,
                    result_bytes = retained.get_byte_size(),
                    elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                Ok(retained)
            },
        )
    }

    fn preprocess_tensor<A: AtomCore>(
        a: &A,
        atom_index: usize,
        settings: &EvaluatorSettings,
        tensor_map: Option<&ParsingTensorMap<'_>>,
        alias_scalar: impl Fn(&Atom) -> bool,
        mut finish: impl FnMut(&ParsingNet, &ScalarAliases, usize) -> Result<AliasedAtom>,
    ) -> Result<AliasedAtom> {
        let atom_started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
            stage = "evaluator_stack_parse_atom_start",
            atom_index,
            do_algebra = settings.do_algebra,
            "Evaluator timing milestone"
        );
        // println!("Parsing {}", a.as_atom_view().log_print(Some(120)));
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
        // CFF sums remain opaque through nested forest operations. Their
        // explicit momentum arguments have already received every deformation.
        let network_input = if network_input.contains_symbol(GS.esurface) {
            GS.lower_esurfaces(&network_input)
        } else {
            network_input
        };
        crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
            stage = "evaluator_stack_parse_atom_normalization_done",
            atom_index,
            elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
            "Normalized evaluator input before independent scalar contractions"
        );
        crate::debug_tags!(#generation, #profile, #compile, #term, #dump;
            stage = "evaluator_stack_parse_atom_before_network_parse",
            atom_index,
            log.atom = network_input,
            "Evaluator atom before network parsing"
        );
        let execute = |net: &mut ParsingNet| -> Result<()> {
            macro_rules! execute_min_result_rank {
                ($execution_strategy:ty) => {
                    match settings.tensor_network_contraction_order {
                        TensorNetworkContractionOrder::IntermediateCost => net
                            .execute::<$execution_strategy, MinIntermediateCost, _, _, _>(
                                TENSORLIB.read().unwrap().deref(),
                                FUN_LIB.deref(),
                            ),
                        TensorNetworkContractionOrder::SparseAtomAware => net
                            .execute::<$execution_strategy, MinResultRank, _, _, _>(
                                TENSORLIB.read().unwrap().deref(),
                                FUN_LIB.deref(),
                            ),
                        TensorNetworkContractionOrder::AtomAware => net
                            .execute::<$execution_strategy, MinResultRankWith<
                                { PAIR_SCORE_ATOM_AWARE },
                                { DEFAULT_EXACT_JOIN_LIMIT },
                            >, _, _, _>(
                                TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()
                            ),
                        TensorNetworkContractionOrder::ResultRankOnly => net
                            .execute::<$execution_strategy, MinResultRankWith<
                                { PAIR_SCORE_RESULT_RANK_ONLY },
                                { DEFAULT_EXACT_JOIN_LIMIT },
                            >, _, _, _>(
                                TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()
                            ),
                        TensorNetworkContractionOrder::EntryAware => net
                            .execute::<$execution_strategy, MinResultRankWith<
                                { PAIR_SCORE_ENTRY_AWARE },
                                { DEFAULT_EXACT_JOIN_LIMIT },
                            >, _, _, _>(
                                TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()
                            ),
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

            Ok(())
        };

        // Materialize only an original additive tensor factor. Scalar spectators,
        // powers and functions stay with the original residual product; this
        // preparation never distributes a graph numerator or changes a strategy.
        // Only exact Gaussian integers qualify; floating coefficients retain
        // their original evaluation order in the ordinary network path.
        let is_gaussian_integer = |atom: AtomView<'_>| {
            let AtomView::Num(number) = atom else {
                return false;
            };
            matches!(number.get_coeff_view().to_owned(),
                symbolica::coefficient::Coefficient::Complex(value)
                    if value.re.is_integer() && value.im.is_integer())
        };
        let constant_factor =
            |factor: AtomView<'_>| -> Result<Option<ParamTensor<ShadowedStructure<Aind>>>> {
                if !matches!(factor, AtomView::Add(_))
                    || factor.get_byte_size() >= MAX_EAGER_TENSOR_SUM_BYTES
                    || spenso::network::profile::lazy_tensor_sums()
                {
                    return Ok(None);
                }
                // Reject scalar parameters and guards before parsing candidate factors.
                // A tensor's index arguments are interpreted by the existing parser.
                let mut pending = vec![factor];
                while let Some(atom) = pending.pop() {
                    match atom {
                        AtomView::Add(sum) => pending.extend(sum.iter()),
                        AtomView::Mul(product) => pending.extend(product.iter()),
                        AtomView::Num(_) if is_gaussian_integer(atom) => {}
                        AtomView::Fun(_) if atom.is_tensorial(StrictTensorFilter::Tagged) => {}
                        _ => return Ok(None),
                    }
                }
                let mut constant = factor.parse_into_net()?;
                let exposed = constant.graph.dangling_indices();
                if exposed.is_empty() || exposed.iter().any(|slot| !slot.matches(slot)) {
                    return Ok(None);
                }
                constant.graph.cache_expr_tree_roots();
                let mut bounds: HashMap<_, (usize, usize)> = HashMap::new();
                for node in constant
                    .graph
                    .cached_expr_preorder_nodes()
                    .into_iter()
                    .rev()
                {
                    let children = constant.graph.cached_expr_children(node);
                    let (entries, bytes) = match &constant.graph.graph[node] {
                        NetworkNode::Op(op @ (NetworkOp::Sum | NetworkOp::Product)) => children
                            .into_iter()
                            .fold((1usize, 0usize), |(entries, bytes), child| {
                                let (child_entries, child_bytes) = bounds[&child];
                                (
                                    if matches!(op, NetworkOp::Sum) {
                                        entries.max(child_entries)
                                    } else {
                                        entries.saturating_mul(child_entries)
                                    },
                                    bytes.saturating_add(child_bytes),
                                )
                            }),
                        NetworkNode::Leaf(NetworkLeaf::Scalar(scalar)) => {
                            let scalar = constant.store.get_scalar_ref(*scalar).as_view();
                            if !is_gaussian_integer(scalar) {
                                return Ok(None);
                            }
                            (1, scalar.get_byte_size() + std::mem::size_of::<Atom>())
                        }
                        NetworkNode::Leaf(
                            leaf @ (NetworkLeaf::LibraryKey { .. } | NetworkLeaf::LocalTensor(_)),
                        ) => {
                            // Bound logical capacity before realizing a library leaf, then
                            // inspect its actual entries rather than its symbol or name.
                            let entries = constant
                                .graph
                                .slots(node)
                                .into_iter()
                                .try_fold(1usize, |size, slot| {
                                    size.checked_mul(usize::try_from(slot.dim()).ok()?)
                                });
                            let Some(entries) = entries.filter(|size| {
                                size.saturating_mul(std::mem::size_of::<Atom>())
                                    < MAX_EAGER_TENSOR_SUM_BYTES
                            }) else {
                                return Ok(None);
                            };
                            let tensor = match leaf {
                                NetworkLeaf::LibraryKey { .. } => constant
                                    .graph
                                    .get_lib_data::<ShadowedStructure<Aind>, _, _>(
                                    TENSORLIB.read().unwrap().deref(),
                                    node,
                                )?,
                                NetworkLeaf::LocalTensor(index) => {
                                    constant.store.get_tensor(*index).clone()
                                }
                                _ => unreachable!(),
                            };
                            let mut bytes = std::mem::size_of::<Atom>();
                            for (_, value) in tensor.iter_flat() {
                                if !is_gaussian_integer(value) {
                                    return Ok(None);
                                }
                                bytes =
                                    bytes.max(value.get_byte_size() + std::mem::size_of::<Atom>());
                            }
                            (entries, bytes)
                        }
                        _ => return Ok(None),
                    };
                    // A product's unreduced Cartesian capacity bounds its partial
                    // contractions; sums preserve the largest child's capacity.
                    // Sum component sizes conservatively and reuse the eager budget.
                    if entries.saturating_mul(bytes) >= MAX_EAGER_TENSOR_SUM_BYTES {
                        return Ok(None);
                    }
                    // Components remain symbolic, so integer arithmetic needs no
                    // floating-point exactness bound.
                    bounds.insert(node, (entries, bytes));
                }
                execute(&mut constant)?;
                let lib = TENSORLIB.read().unwrap();
                let ExecutionResult::Val(tensor) = constant.result_tensor(lib.deref())? else {
                    return Ok(None);
                };
                Ok(Some(tensor.into_owned()))
            };

        // Each existing top-level summand is an independent tensor contraction.
        // Keeping its network local avoids scanning unrelated terms during
        // finite component preparation. Products, powers and nested sums retain
        // their grouping; the resulting expressions are reunited before optimization.
        let terms = if let AtomView::Add(sum) = network_input.as_view() {
            sum.iter().collect::<Vec<_>>()
        } else {
            vec![network_input.as_view()]
        };
        crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
            stage = "evaluator_stack_parse_atom_terms_start",
            atom_index,
            term_count = terms.len(),
            "Contracting independent scalar summands"
        );
        // Network handles are local indices. Give each summand its own scope
        // before combining independently contracted networks.
        let result = terms
            .into_iter()
            .enumerate()
            .map(|(term_index, term)| -> Result<AliasedAtom> {
                let term_started = std::time::Instant::now();
                let mut residual = Vec::new();
                let mut constants = Vec::new();
                if let AtomView::Mul(product) = term {
                    for factor in product.iter() {
                        if let Some(tensor) = constant_factor(factor)? {
                            constants.push(tensor);
                        } else {
                            residual.push(factor);
                        }
                    }
                }
                let mut net = if constants.is_empty() {
                    term.parse_into_net()?
                } else {
                    let mut net = Atom::mul_many(residual).parse_into_net()?;
                    for tensor in constants.into_iter().rev() {
                        net = ParsingNet::from_tensor(tensor) * net;
                    }
                    net
                };
                if let Some(tensor_map) = tensor_map {
                    net = net.map_result(Ok, tensor_map)?;
                }
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_net_done",
                    atom_index,
                    term_index,
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );

                // println!("Net: {}", net.dot_pretty());
                let scalar_aliases = net.alias_scalar_refs(|_, scalar| {
                    (scalar.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES
                        || alias_scalar(scalar))
                        // Keep guards visible until they enclose the complete branch,
                        // including inverses in neighboring scalar factors. These
                        // definitions precede generated handles, so selector structure
                        // cannot be hidden through another registered alias either.
                        && [OrientationID::symbol(), GS.theta, GS.orientation_delta, Symbol::IF]
                            .into_iter().all(|selector| !scalar.contains_symbol(selector))
                });
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_scalar_aliases_done",
                    atom_index,
                    term_index,
                    automatic_scalar_aliasing = true,
                    threshold_bytes = NETWORK_SCALAR_ALIAS_MIN_BYTES,
                    aliases_created = scalar_aliases.aliases_created(),
                    aliased_terms = scalar_aliases.aliased_terms(),
                    aliased_bytes = scalar_aliases.aliased_bytes(),
                    max_aliased_bytes = scalar_aliases.max_aliased_bytes(),
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                // Prepare only the finite component contraction. Raw symbolic networks
                // used by Taylor expansion retain their original product/sum grouping.
                let contraction_preparation_started = std::time::Instant::now();
                let closed_sum_boundaries = net.graph.contract_ready_sum_boundaries();
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_tensor_boundaries_done",
                    atom_index,
                    term_index,
                    closed_sum_boundaries,
                    elapsed_ms = contraction_preparation_started.elapsed().as_secs_f64() * 1000.0,
                    "Prepared finite tensor contractions through pending sums"
                );
                crate::debug_tags!(#generation, #compile, #term, #dump;
                    stage = "evaluator_stack_parse_atom_network_dump",
                    atom_index,
                    term_index,
                    file.atom = %term.to_canonical_string(),
                    file.network = %net.dot_pretty(),
                    "Parsed evaluator network dump"
                );

                let parse_elapsed = term_started.elapsed();
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_parse_elapsed",
                    atom_index,
                    term_index,
                    elapsed_ms = parse_elapsed.as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                let instant = std::time::Instant::now();

                execute(&mut net)?;

                let execute_elapsed = instant.elapsed();
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_execute_elapsed",
                    atom_index,
                    term_index,
                    elapsed_ms = execute_elapsed.as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );
                crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                    stage = "evaluator_stack_parse_atom_execute_done",
                    atom_index,
                    term_index,
                    elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluator timing milestone"
                );

                finish(&net, &scalar_aliases, term_index).map_err(|error| {
                    error.with_note(|| format!("Network looks like: {}", net.dot_pretty()))
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(|terms| {
                terms.into_iter().fold(AliasedAtom::default(), |sum, term| {
                    sum.try_add(&term).unwrap()
                })
            });
        crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
            stage = "evaluator_stack_parse_atom_done",
            atom_index,
            success = result.is_ok(),
            result_bytes = result.as_ref().map_or(0, AliasedAtom::get_byte_size),
            elapsed_ms = atom_started.elapsed().as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );
        result
    }

    /// Close products of shared tensor bodies before enumerating their components.
    /// Each joint body keeps fresh formals per call occurrence and is reused across
    /// physical argument rows. Only residue sums over opaque calls are collected;
    /// scalar weights and factorized numerator sums retain their boundaries.
    fn combine_numerator_families<A: AtomCore>(
        atoms: &[A],
        definitions: &[Arc<FnMapEntry>],
        param_builder: &ParamBuilder,
        families: &HashMap<(Symbol, Vec<Atom>), usize>,
        tag_counts: &HashMap<Symbol, usize>,
        symbols: (Symbol, Symbol, Symbol),
    ) -> (Vec<Atom>, Vec<Arc<FnMapEntry>>) {
        let (joint_symbol, argument_symbol, shadow_symbol) = symbols;
        // Renaming a formal hidden inside a parameter alias would lose its
        // lexical binding. Leave those families at the existing open interface.
        let interfaces = definitions
            .iter()
            .map(|definition| {
                let captured = definition.args.iter().any(|formal| {
                    let formal = Atom::from(formal.clone());
                    param_builder
                        .reps
                        .iter()
                        .any(|entry| entry.rhs.contains(&formal))
                });
                let mut parameterized_slot = false;
                definition.rhs.visitor(&mut |part| {
                    if Slot::<LibraryRep, Aind>::try_from(part).is_ok() {
                        parameterized_slot |= definition
                            .args
                            .iter()
                            .any(|formal| part.contains(Atom::from(formal.clone())));
                    }
                    !parameterized_slot
                });
                // Row-dependent ports must be bound before contracting families.
                // Joint formals describe scalar energies, not tensor interfaces.
                if captured || parameterized_slot {
                    return None;
                }
                // An unknown symbolic interface is not eligible for fusion;
                // the ordinary preparation retains its existing validation.
                definition
                    .rhs
                    .list_dangling::<Aind>()
                    .ok()
                    .map(|slots| slots.into_iter().collect::<HashSet<_>>())
            })
            .collect::<Vec<_>>();
        let family_call = |view: AtomView<'_>| {
            let AtomView::Fun(call) = view else {
                return None;
            };
            let &tag_count = tag_counts.get(&call.get_symbol())?;
            let tags = call.iter().take(tag_count).map(|a| a.to_owned()).collect();
            let &family = families.get(&(call.get_symbol(), tags))?;
            (interfaces[family].is_some()
                && call.get_nargs() == tag_count + definitions[family].args.len())
            .then(|| {
                (
                    family,
                    call.iter()
                        .skip(tag_count)
                        .map(|a| a.to_owned())
                        .collect::<Vec<_>>(),
                )
            })
        };
        let has_family =
            |view: AtomView<'_>| tag_counts.keys().any(|head| view.contains_symbol(*head));
        let structural_function = |view: AtomView<'_>| {
            view.is_tensorial(StrictTensorFilter::Tagged)
                || matches!(view, AtomView::Fun(fun)
                    if [*SYM, *ANTISYM, *CYCLIC].contains(&fun.get_symbol()))
        };
        // Infer closure from small interface shadows, never by copying a body
        // into every physical residue row. Functions keep their argument scope.
        let dangling = |view: AtomView<'_>| {
            let mut valid = true;
            let shadow = view.replace_map(|part, _, out| {
                if let Some((family, _)) = family_call(part) {
                    let mut slots = interfaces[family]
                        .as_ref()
                        .unwrap()
                        .iter()
                        .collect::<Vec<_>>();
                    slots.sort();
                    **out = FunctionBuilder::new(shadow_symbol)
                        .add_arg(family)
                        .add_args(slots)
                        .finish();
                } else if matches!(part, AtomView::Fun(_) | AtomView::Pow(_))
                    && !structural_function(part)
                {
                    valid &= !has_family(part);
                    out.set_from_view(&part);
                }
            });
            valid.then(|| shadow.list_dangling::<Aind>().ok()).flatten()
        };
        let mut occupied_symbols = symbolica::state::State::symbol_iter()
            .flat_map(|(symbol, name)| {
                std::iter::once(name.to_owned()).chain(symbol.get_aliases().to_vec())
            })
            .collect::<HashSet<_>>();
        let mut next_index = 0usize;
        let mut combined = Vec::<Arc<FnMapEntry>>::new();
        let mut cached = HashMap::<Atom, Option<usize>>::new();
        let mut join = |core: AtomView<'_>| -> Atom {
            let mut calls = Vec::new();
            core.visitor(&mut |part| {
                if let Some((family, arguments)) = family_call(part) {
                    calls.push((family, arguments, calls.len()));
                    false
                } else {
                    !matches!(part, AtomView::Fun(_) | AtomView::Pow(_))
                        || structural_function(part)
                }
            });
            if calls.is_empty() {
                return core.to_owned();
            }
            // Physical arguments can reorder a commutative product. Bind in
            // family order, while retaining each occurrence's literal context.
            calls.sort_by_key(|(family, _, _)| *family);
            let mut positions = vec![0; calls.len()];
            for (position, (_, _, original)) in calls.iter().enumerate() {
                positions[*original] = position;
            }
            let mut occurrence = 0;
            let skeleton = core.replace_map(|part, _, out| {
                if let Some((family, _)) = family_call(part) {
                    **out = function!(joint_symbol, family, positions[occurrence]);
                    occurrence += 1;
                } else if matches!(part, AtomView::Fun(_) | AtomView::Pow(_))
                    && !structural_function(part)
                {
                    out.set_from_view(&part);
                }
            });
            let index = *cached.entry(skeleton.clone()).or_insert_with(|| {
                // Do not materialize disconnected open outer products: a joint
                // body must remove the entire finite component interface.
                if !dangling(core).is_some_and(|slots| slots.is_empty()) {
                    return None;
                }
                let index = combined.len();
                let mut formals = Vec::new();
                let mut bodies = Vec::new();
                for (occurrence, (family, _, _)) in calls.iter().enumerate() {
                    let family = *family;
                    let definition = &definitions[family];
                    let arguments = (0..definition.args.len())
                        .map(|position| function!(argument_symbol, index, occurrence, position))
                        .collect::<Vec<_>>();
                    // Each occurrence owns private contractions; only its free
                    // ports can connect it to other factors in this template.
                    let interface = interfaces[family].as_ref().unwrap();
                    let body = definition.rhs.freshen_private_indices(interface, || {
                        loop {
                            let name = format!("{}_index_{next_index}", argument_symbol.get_name());
                            next_index += 1;
                            if occupied_symbols.insert(name.clone()) {
                                break Atom::var(symbol!(&name));
                            }
                        }
                    });

                    bodies.push(
                        body.replace_multiple(
                            definition
                                .args
                                .iter()
                                .cloned()
                                .map(Atom::from)
                                .zip(&arguments)
                                .map(|(formal, argument)| {
                                    Replacement::new(formal.to_pattern(), argument.clone())
                                }),
                        ),
                    );
                    formals.extend(arguments);
                }
                let rhs = skeleton.replace_map(|part, _, out| {
                    if let AtomView::Fun(call) = part
                        && call.get_symbol() == joint_symbol
                    {
                        let position = usize::try_from(call.iter().nth(1).unwrap()).unwrap();
                        **out = bodies[position].clone();
                    }
                });
                // A body may already own a chain/trace binder. Inserting that
                // body inside another binder would change the shared in/out scope;
                // leave such calls at their existing tensor interface instead.
                if rhs.validate_chain_like_nesting().is_err() {
                    return None;
                }
                combined.push(Arc::new(FnMapEntry {
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                    lhs: FunctionBuilder::new(joint_symbol)
                        .add_arg(index)
                        .add_args(&formals)
                        .finish(),
                    rhs,
                    args: formals
                        .into_iter()
                        .map(|arg| Indeterminate::try_from(arg).unwrap())
                        .collect(),
                    tags: vec![Atom::num(index)],
                }));
                Some(index)
            });
            match index {
                Some(index) => FunctionBuilder::new(joint_symbol)
                    .add_arg(index)
                    .add_args(calls.into_iter().flat_map(|(_, arguments, _)| arguments))
                    .finish(),
                None => core.to_owned(),
            }
        };
        let mut prepare_scope = |root: AtomView<'_>| {
            // Collect only the explicit sum over family calls. Complete tensor
            // numerator sums stay opaque, as do independent scalar sums and all
            // function/power scopes. Scalar weights are restored outside the
            // callback by the existing collector, never put in the function map.
            let collected = root.collect_with_map(|part| match part {
                AtomView::Fun(_) | AtomView::Pow(_) => true,
                AtomView::Add(_) => {
                    !has_family(part)
                        || (part != root && dangling(part).is_none_or(|slots| slots.is_empty()))
                }
                _ => false,
            });
            collected
                .map_collects(|wrapped, _, out| {
                    let core = wrapped.as_fun_view().unwrap().iter().next().unwrap();
                    let factors = match core {
                        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                        _ => vec![core],
                    };
                    let mut open = Vec::new();
                    let mut remaining = Vec::new();
                    for factor in factors {
                        let mut eligible = !matches!(factor, AtomView::Pow(_))
                            && (family_call(factor).is_some()
                                || factor.is_tensorial(StrictTensorFilter::Tagged));
                        factor.visitor(&mut |part| {
                            eligible &= match part {
                                AtomView::Pow(power) => {
                                    i64::try_from(power.get_base_exp().1)
                                        .is_ok_and(|exponent| exponent >= 0)
                                        && !has_family(part)
                                }
                                AtomView::Fun(fun) => {
                                    ![
                                        OrientationID::symbol(),
                                        GS.theta,
                                        GS.orientation_delta,
                                        Symbol::IF,
                                    ]
                                    .contains(&fun.get_symbol())
                                        && (!tag_counts.contains_key(&fun.get_symbol())
                                            || family_call(part).is_some())
                                }
                                _ => true,
                            };
                            eligible && family_call(part).is_none()
                        });
                        // Keep an already closed family sum separate from
                        // the surrounding tensor contractions.
                        eligible &= !matches!(factor, AtomView::Add(_)) || !has_family(factor);
                        if !eligible {
                            remaining.push(factor.to_owned());
                        } else if dangling(factor).is_some_and(|slots| slots.is_empty()) {
                            remaining.push(if family_call(factor).is_some() {
                                factor.to_owned()
                            } else {
                                join(factor)
                            });
                        } else {
                            open.push(factor);
                        }
                    }
                    if !open.is_empty() {
                        remaining.push(join(Atom::mul_many(open).as_view()));
                    }
                    **out = Atom::mul_many(remaining);
                })
                .unwrap_collect()
        };
        let atoms = atoms
            .iter()
            .map(|atom| {
                let root = atom.as_atom_view();
                let mut scopes = vec![root.to_owned()];
                root.visitor(&mut |part| {
                    if matches!(part, AtomView::Fun(_) | AtomView::Pow(_)) {
                        return false;
                    }
                    if part != root
                        && matches!(part, AtomView::Add(_))
                        && has_family(part)
                        && dangling(part).is_some_and(|slots| slots.is_empty())
                    {
                        scopes.push(part.to_owned());
                    }
                    true
                });
                // Prepare closed sums independently, from the inside out. They stay
                // opaque in their parent product, but their own residue rows still
                // reuse closed bodies. A scoped worklist preserves function/power
                // barriers that an unrestricted bottom-up traversal would cross.
                let mut prepared = HashMap::<Atom, Atom>::new();
                for scope in scopes.into_iter().rev() {
                    if prepared.contains_key(&scope) {
                        continue;
                    }
                    let inner = scope.replace_map(|part, _, out| {
                        if matches!(part, AtomView::Fun(_) | AtomView::Pow(_)) {
                            out.set_from_view(&part);
                        } else if let Some(replacement) = prepared.get(part.get_data()) {
                            **out = replacement.clone();
                        }
                    });
                    prepared.insert(scope, prepare_scope(inner.as_view()));
                }
                prepared.remove(root.get_data()).unwrap()
            })
            .collect();
        (atoms, combined)
    }

    /// Contract each shared numerator body once, then expose only its scalar
    /// components to the evaluator. Formal branch coefficients remain arguments
    /// of both component functions and the scalar aliases they retain.
    fn preprocess_numerator_families<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        definitions: &[Arc<FnMapEntry>],
        settings: &EvaluatorSettings,
        alias_symbol: Symbol,
        inlining: InliningPolicy,
    ) -> Result<(Vec<AliasedAtom>, ParamBuilder)> {
        let mut builder = param_builder.clone();
        if definitions.is_empty() {
            let atoms = atoms
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    Self::preprocess_atom(atom, index, settings, alias_symbol, None)
                })
                .collect::<Result<_>>()?;
            return Ok((atoms, builder));
        }
        let family_preparation_started = std::time::Instant::now();
        let (shadow_symbol, component_symbol, coefficient_symbol, joint_symbol, argument_symbol) = loop {
            let scope = NETWORK_SCALAR_ALIAS_SCOPE.fetch_add(1, Ordering::Relaxed);
            let names = [
                format!("gammalooprs::numerator_tensor_{scope}"),
                format!("gammalooprs::numerator_component_{scope}"),
                format!("gammalooprs::numerator_coefficient_{scope}"),
                format!("gammalooprs::numerator_joint_{scope}"),
                format!("gammalooprs::numerator_joint_argument_{scope}"),
            ];
            if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                names.iter().any(|name| {
                    existing == name || symbol.get_aliases().iter().any(|alias| alias == name)
                })
            }) {
                break (
                    SPENSO_TAG.tensor_symbol(&names[0]),
                    symbol!(&names[1]),
                    symbol!(&names[2]),
                    symbol!(&names[3]),
                    symbol!(&names[4]; Scalar),
                );
            }
        };
        // Validate the flat physical family catalog before selecting bodies.
        // Unused Taylor coefficients need no tensor contraction, but malformed
        // signatures and conflicting definitions must remain errors.
        let mut families: HashMap<(Symbol, Vec<Atom>), usize> = HashMap::new();
        let mut tag_counts = HashMap::new();
        for (family_index, definition) in definitions.iter().enumerate() {
            let head = definition
                .lhs
                .as_fun_view()
                .ok_or_else(|| eyre!("shared numerator definition requires a function call"))?
                .get_symbol();
            if let Some(tag_count) = tag_counts.insert(head, definition.tags.len())
                && tag_count != definition.tags.len()
            {
                return Err(eyre!("inconsistent shared numerator tag count"));
            }
            if let Some(previous) = families.get(&(head, definition.tags.clone())) {
                if definitions[*previous] != *definition {
                    return Err(eyre!("conflicting shared numerator definitions"));
                }
                continue;
            }
            families.insert((head, definition.tags.clone()), family_index);
            let formal_args = definition
                .args
                .iter()
                .cloned()
                .map(Atom::from)
                .collect::<Vec<_>>();
            if definition.lhs
                != FunctionBuilder::new(head)
                    .add_args(&definition.tags)
                    .add_args(&formal_args)
                    .finish()
            {
                return Err(eyre!(
                    "shared numerator tags and formal arguments do not match its call"
                ));
            }
            if formal_args.iter().collect::<HashSet<_>>().len() != formal_args.len() {
                return Err(eyre!("shared numerator formal arguments must be distinct"));
            }
        }
        if definitions.iter().any(|definition| {
            tag_counts
                .keys()
                .any(|head| definition.rhs.contains_symbol(*head))
        }) {
            return Err(eyre!("nested shared numerator families are unsupported"));
        }
        let combination_started = std::time::Instant::now();
        let (atoms, combined) = Self::combine_numerator_families(
            atoms,
            definitions,
            param_builder,
            &families,
            &tag_counts,
            (joint_symbol, argument_symbol, shadow_symbol),
        );
        crate::debug_tags!(#generation, #profile, #summary;
            stage = "evaluator_stack_numerator_families_combined",
            elapsed_ms = combination_started.elapsed().as_secs_f64() * 1000.0,
            supplied_families = definitions.len(), combined_families = combined.len(),
            "Combined closed numerator tensor families"
        );
        let mut definitions = definitions.to_vec();
        for entry in combined {
            families.insert((joint_symbol, entry.tags.clone()), definitions.len());
            tag_counts.insert(joint_symbol, entry.tags.len());
            definitions.push(entry);
        }
        let mut referenced = HashSet::new();
        let mut reference_error = None;
        for atom in &atoms {
            atom.as_atom_view().visitor(&mut |view| {
                let AtomView::Fun(call) = view else {
                    return true;
                };
                let Some(&tag_count) = tag_counts.get(&call.get_symbol()) else {
                    return true;
                };
                let tags = call
                    .iter()
                    .take(tag_count)
                    .map(|arg| arg.to_owned())
                    .collect::<Vec<_>>();
                let Some(&index) = families.get(&(call.get_symbol(), tags)) else {
                    reference_error = Some("unknown shared numerator family tags");
                    return false;
                };
                if call.get_nargs() != tag_count + definitions[index].args.len() {
                    reference_error = Some("invalid shared numerator family argument count");
                    return false;
                }
                referenced.insert(index);
                true
            });
        }
        if let Some(error) = reference_error {
            return Err(eyre!(error));
        }
        // Family definitions are flat. Calls hidden inside existing parameter
        // functions cannot be lowered through the outer tensor interface.
        if param_builder.reps.iter().any(|entry| {
            tag_counts
                .keys()
                .any(|head| entry.rhs.contains_symbol(*head))
        }) {
            return Err(eyre!(
                "shared numerator calls inside parameter functions are unsupported"
            ));
        }
        // This catalog belongs to this lowering operation. An omitted sparse
        // component is zero only for a known family and a valid component index.
        // Each populated component keeps only its required backend arguments;
        // incoming physical family calls retain their complete signature.
        type FamilyComponents = (Vec<(usize, bool)>, HashMap<Atom, Vec<usize>>, usize);
        let mut components = HashMap::<Vec<Atom>, FamilyComponents>::new();
        let mut shadows = HashMap::new();
        let mut generated = Vec::new();
        let physical_inputs = (&param_builder.pairs)
            .into_iter()
            .flat_map(|pair| pair.params.iter().cloned())
            .collect::<HashSet<_>>();
        let max_input_bytes = physical_inputs
            .iter()
            .map(|input| input.as_view().get_data().len())
            .max()
            .unwrap_or(0);
        for (family_index, definition) in definitions.iter().enumerate() {
            if !referenced.contains(&family_index) {
                continue;
            }
            let formal_args = definition
                .args
                .iter()
                .cloned()
                .map(Atom::from)
                .collect::<Vec<_>>();
            let formal_positions = formal_args
                .iter()
                .cloned()
                .enumerate()
                .map(|(index, atom)| (atom, index))
                .collect::<HashMap<_, _>>();
            let max_formal_bytes = formal_positions
                .keys()
                .map(|formal| formal.as_view().get_data().len())
                .max()
                .unwrap_or(0);
            // An existing zero-argument parameter alias can capture a formal
            // without mentioning it in the family body. Keep such formals
            // conservatively; generated scalar aliases are analyzed exactly below.
            let captured_arguments = formal_args
                .iter()
                .enumerate()
                .filter_map(|(index, formal)| {
                    param_builder
                        .reps
                        .iter()
                        .any(|entry| entry.rhs.contains(formal))
                        .then_some(index)
                })
                .collect::<HashSet<_>>();
            let shadow = Self::preprocess_tensor(
                &definition.rhs,
                family_index,
                settings,
                None,
                |scalar| {
                    // Keep short sums visible during contraction when the final
                    // component is inlined. Aliasing and restoring them afterwards
                    // changes the algebra presented to Horner and loses sharing.
                    // Non-inlined components retain these scalar functions.
                    if inlining == InliningPolicy::Always {
                        return false;
                    }
                    // Functions and physical input keys stay opaque; exact formals
                    // shadow either in this body.
                    let mut has_formal = false;
                    let mut has_sum = false;
                    scalar.visitor(&mut |view| {
                        if has_formal && has_sum {
                            return false;
                        }
                        if view.get_data().len() <= max_formal_bytes
                            && formal_positions.contains_key(view.get_data())
                        {
                            has_formal = true;
                            return false;
                        }
                        if (view.get_data().len() <= max_input_bytes
                            && physical_inputs.contains(view.get_data()))
                            || matches!(view, AtomView::Fun(_))
                        {
                            return false;
                        }
                        has_sum |= matches!(view, AtomView::Add(_));
                        true
                    });
                    has_formal && has_sum
                },
                |net, scalar_aliases, term_index| {
                    let tensor = match net.result_tensor(TENSORLIB.read().unwrap().deref())? {
                        ExecutionResult::One => return Ok(Atom::num(1).into()),
                        ExecutionResult::Zero => return Ok(Atom::Zero.into()),
                        ExecutionResult::Val(tensor) => tensor.into_owned(),
                    };
                    let dimensions = tensor
                        .external_reps_iter()
                        .map(|rep| Ok((usize::try_from(rep.dim)?, !rep.rep.is_base())))
                        .collect::<Result<Vec<_>>>()?;
                    let tags = vec![Atom::num(family_index), Atom::num(term_index)];
                    let (_, mut aliases) = net
                        .aliased_atom(scalar_aliases, Atom::Zero)
                        .into_inner_with_aliases();
                    // Match ordinary scalar preparation before assigning the
                    // retained definitions their minimal formal dependencies.
                    for (alias, body) in &mut aliases {
                        if *body != *alias {
                            *body = body.collect_compact_factors();
                        }
                    }
                    let mut handles = aliases.keys().cloned().collect::<Vec<_>>();
                    handles.sort();
                    let handle_positions = handles
                        .iter()
                        .cloned()
                        .enumerate()
                        .map(|(index, handle)| (handle, index))
                        .collect::<HashMap<_, _>>();
                    let max_handle_bytes = handle_positions
                        .keys()
                        .map(|handle| handle.as_view().get_data().len())
                        .max()
                        .unwrap_or(0);
                    let dependencies = |body: &Atom| {
                        let mut arguments = captured_arguments.clone();
                        let mut referenced_aliases = HashSet::new();
                        body.visitor(&mut |view| {
                            if view.get_data().len() <= max_formal_bytes
                                && let Some(index) = formal_positions.get(view.get_data())
                            {
                                arguments.insert(*index);
                                false
                            } else if view.get_data().len() <= max_handle_bytes
                                && let Some(index) = handle_positions.get(view.get_data())
                            {
                                referenced_aliases.insert(*index);
                                false
                            } else {
                                true
                            }
                        });
                        (arguments, referenced_aliases)
                    };
                    let (mut alias_arguments, alias_dependencies): (Vec<_>, Vec<_>) = handles
                        .iter()
                        .map(|handle| dependencies(&aliases[handle]))
                        .unzip();
                    // Propagate only small argument sets, never inline alias
                    // bodies. At most one pass per alias reaches the fixed point.
                    for _ in 0..handles.len() {
                        let mut changed = false;
                        for index in 0..handles.len() {
                            let mut required = alias_arguments[index].clone();
                            for dependency in &alias_dependencies[index] {
                                required.extend(alias_arguments[*dependency].iter().copied());
                            }
                            changed |= required != alias_arguments[index];
                            alias_arguments[index] = required;
                        }
                        if !changed {
                            break;
                        }
                    }
                    let argument_positions = |body: &Atom| {
                        let (mut required, dependencies) = dependencies(body);
                        for index in dependencies {
                            required.extend(alias_arguments[index].iter().copied());
                        }
                        let mut required = required.into_iter().collect::<Vec<_>>();
                        required.sort_unstable();
                        required
                    };
                    let renames = handles
                        .iter()
                        .cloned()
                        .enumerate()
                        .map(|(index, handle)| {
                            let alias_tags = [tags.clone(), vec![Atom::num(index)]].concat();
                            let mut positions =
                                alias_arguments[index].iter().copied().collect::<Vec<_>>();
                            positions.sort_unstable();
                            let call = FunctionBuilder::new(coefficient_symbol)
                                .add_args(&alias_tags)
                                .add_args(positions.iter().map(|index| &formal_args[*index]))
                                .finish();
                            (handle, (alias_tags, call, positions))
                        })
                        .collect::<HashMap<_, _>>();
                    let rename = |atom: Atom| {
                        if renames.is_empty() {
                            return atom;
                        }
                        atom.replace_map(|view, _, out| {
                            if view.get_data().len() <= max_handle_bytes
                                && let Some((_, call, _)) = renames.get(view.get_data())
                            {
                                out.set_from_view(&call.as_view());
                            }
                        })
                    };
                    for (alias, body) in aliases {
                        let body = if body == alias { body } else { rename(body) };
                        generated.push(FnMapEntry {
                            inlining,
                            is_alias: false,
                            lhs: renames[&alias].1.clone(),
                            rhs: body,
                            tags: renames[&alias].0.clone(),
                            args: renames[&alias]
                                .2
                                .iter()
                                .map(|index| definition.args[*index].clone())
                                .collect(),
                        });
                    }
                    let mut populated = HashMap::new();
                    for (index, body) in tensor.iter_flat() {
                        let body = body.collect_compact_factors();
                        if body.is_zero() {
                            continue;
                        }
                        let index = Atom::from(tensor.co_expanded_index(index)?);
                        let positions = argument_positions(&body);
                        populated.insert(index.clone(), positions.clone());
                        generated.push(FnMapEntry {
                            inlining,
                            is_alias: false,
                            lhs: FunctionBuilder::new(component_symbol)
                                .add_args(&tags)
                                .add_arg(&index)
                                .add_args(positions.iter().map(|index| &formal_args[*index]))
                                .finish(),
                            rhs: rename(body),
                            tags: [tags.clone(), vec![index]].concat(),
                            args: positions
                                .iter()
                                .map(|index| definition.args[*index].clone())
                                .collect(),
                        });
                    }
                    if populated.is_empty() {
                        return Ok(Atom::Zero.into());
                    }
                    if dimensions.is_empty() {
                        return Ok(FunctionBuilder::new(component_symbol)
                            .add_args(&tags)
                            .add_arg(function!(AIND_SYMBOLS.cind))
                            .add_args(
                                populated[&function!(AIND_SYMBOLS.cind)]
                                    .iter()
                                    .map(|index| &formal_args[*index]),
                            )
                            .finish()
                            .into());
                    }
                    let shadow = tensor.structure().to_symbolic_with(
                        shadow_symbol,
                        &[tags.clone(), formal_args.clone()].concat(),
                        None,
                    );
                    components.insert(tags, (dimensions, populated, formal_args.len()));
                    Ok(shadow.into())
                },
            )?;
            shadows.insert(
                family_index,
                (shadow.into_inner(), formal_positions, max_formal_bytes),
            );
        }
        crate::debug_tags!(#generation, #profile, #summary;
            stage = "evaluator_stack_numerator_families_prepared",
            elapsed_ms = family_preparation_started.elapsed().as_secs_f64() * 1000.0,
            referenced_families = referenced.len(), generated_definitions = generated.len(),
            "Prepared shared numerator scalar components"
        );
        let down_symbol = Atom::from(DualConciousIndex::Down(0))
            .as_fun_view()
            .unwrap()
            .get_symbol();
        let lower_components = |atom: Atom| -> Result<Atom> {
            let mut error = None;
            let atom = atom.replace_map(|view, _, out| {
                let AtomView::Fun(function) = view else {
                    return;
                };
                if function.get_symbol() != shadow_symbol {
                    return;
                }
                let args = function
                    .iter()
                    .map(|arg| arg.to_owned())
                    .collect::<Vec<_>>();
                let Some((dimensions, populated, argument_count)) =
                    args.get(..2).and_then(|tags| components.get(tags))
                else {
                    error = Some("unknown numerator tensor component family");
                    return;
                };
                let valid_index = args
                    .last()
                    .and_then(|index| index.as_fun_view())
                    .filter(|index| index.get_symbol() == AIND_SYMBOLS.cind)
                    .is_some_and(|index| {
                        index.get_nargs() == dimensions.len()
                            && index
                                .iter()
                                .zip(dimensions)
                                .all(|(index, (dimension, down))| {
                                    let index = if *down {
                                        let AtomView::Fun(wrapper) = index else {
                                            return false;
                                        };
                                        if wrapper.get_symbol() != down_symbol
                                            || wrapper.get_nargs() != 1
                                        {
                                            return false;
                                        }
                                        wrapper.iter().next().unwrap()
                                    } else {
                                        index
                                    };
                                    matches!(index, AtomView::Num(number)
                                    if matches!(number.get_coeff_view(),
                                        symbolica::coefficient::CoefficientView::Natural(i, 1, 0, 1)
                                            if usize::try_from(i).is_ok_and(|i| i < *dimension)))
                                })
                    });
                if args.len() != argument_count + 3 || !valid_index {
                    error = Some("invalid numerator tensor component index or arguments");
                    return;
                }
                let index = args.last().unwrap();
                let call = if let Some(positions) = populated.get(index) {
                    FunctionBuilder::new(component_symbol)
                        .add_args(&args[..2])
                        .add_arg(index)
                        .add_args(positions.iter().map(|index| &args[2 + index]))
                        .finish()
                } else {
                    Atom::Zero
                };
                out.set_from_view(&call.as_view());
            });
            if let Some(error) = error {
                return Err(eyre!(error));
            }
            Ok(atom)
        };
        let lower_shadow_components = |mut param: ParamTensor<ShadowedStructure<Aind>>| {
            let tensor = param
                .tensor
                .map_data_ref_mut_result(|value| lower_components(std::mem::take(value)))?;
            param.tensor = tensor.to_sparse();
            Ok(param)
        };
        let atoms = atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| {
                let lowering_started = std::time::Instant::now();
                let mut lowered_calls = 0usize;
                let shadowed = atom.replace_map(|view, _, out| {
                    let AtomView::Fun(call) = view else {
                        return;
                    };
                    let Some(&tag_count) = tag_counts.get(&call.get_symbol()) else {
                        return;
                    };
                    let tags = call
                        .iter()
                        .take(tag_count)
                        .map(|arg| arg.to_owned())
                        .collect();
                    let family = families[&(call.get_symbol(), tags)];
                    let (shadow, formal_positions, max_formal_bytes) = &shadows[&family];
                    let arguments = call.iter().skip(tag_count).collect::<Vec<_>>();
                    // Calls and signatures were validated above. Dispatch by their
                    // catalog key, then bind exact formals simultaneously: inserted
                    // arguments must not be visited or substituted a second time.
                    let bound = shadow.replace_map(|part, _, output| {
                        if part.get_data().len() <= *max_formal_bytes
                            && let Some(position) = formal_positions.get(part.get_data())
                        {
                            output.set_from_view(&arguments[*position]);
                        }
                    });
                    out.set_from_view(&bound.as_view());
                    lowered_calls += 1;
                });
                crate::debug_tags!(#generation, #profile, #summary;
                    stage = "evaluator_stack_numerator_calls_lowered",
                    atom_index = index, lowered_calls,
                    input_bytes = atom.as_view().get_byte_size(),
                    output_bytes = shadowed.as_view().get_byte_size(),
                    elapsed_ms = lowering_started.elapsed().as_secs_f64() * 1000.0,
                    "Lowered shared numerator calls by their catalog keys"
                );
                let (root, aliases) = Self::preprocess_atom(
                    &shadowed,
                    index,
                    settings,
                    alias_symbol,
                    Some(&lower_shadow_components),
                )?
                .into_inner_with_aliases();
                let mut result = AliasedAtom::from(lower_components(root)?);
                for (alias, body) in aliases {
                    result.register_alias(alias, lower_components(body)?);
                }
                result.prune();
                Ok(result)
            })
            .collect::<Result<Vec<_>>>()?;
        // Outer projectors can discard populated tensor components. Keep only
        // generated functions reached by the final roots and retained aliases;
        // every original parameter function remains in the cloned builder.
        let generated_keys = generated
            .iter()
            .enumerate()
            .map(|(index, entry)| {
                (
                    (
                        entry.lhs.as_fun_view().unwrap().get_symbol(),
                        entry.tags.clone(),
                    ),
                    index,
                )
            })
            .collect::<HashMap<_, _>>();
        let mut pending = atoms
            .iter()
            .flat_map(|atom| std::iter::once(atom.get_root()).chain(atom.get_aliases().values()))
            .collect::<Vec<_>>();
        let mut reachable = HashSet::new();
        while let Some(body) = pending.pop() {
            let mut error = None;
            body.visitor(&mut |view| {
                let AtomView::Fun(call) = view else {
                    return true;
                };
                if ![component_symbol, coefficient_symbol].contains(&call.get_symbol()) {
                    return true;
                }
                // Generated functions use family, term and component/alias tags.
                let tags = call
                    .iter()
                    .take(3)
                    .map(|arg| arg.to_owned())
                    .collect::<Vec<_>>();
                let Some(&index) = generated_keys.get(&(call.get_symbol(), tags)) else {
                    error = Some("unknown generated numerator function");
                    return false;
                };
                if call.get_nargs() != 3 + generated[index].args.len() {
                    error = Some("invalid generated numerator function arguments");
                    return false;
                }
                if reachable.insert(index) {
                    pending.push(&generated[index].rhs);
                }
                true
            });
            if let Some(error) = error {
                return Err(eyre!(error));
            }
        }
        crate::debug_tags!(#generation, #profile, #summary;
            stage = "evaluator_stack_numerator_components_reachable",
            supplied_families = definitions.len(), referenced_families = referenced.len(),
            generated_definitions = generated.len(), retained_generated_definitions = reachable.len(),
            "Retained numerator functions after outer tensor contraction"
        );
        // Keep each scalar numerator body as one reusable program. The
        // coefficient rows at call sites remain source-local, including in a
        // complete LTD sum and in independently projected UV children.
        for (index, entry) in generated.into_iter().enumerate() {
            if reachable.contains(&index) {
                builder
                    .add_function_entry(entry)
                    .map_err(|error| eyre!(error))?;
            }
        }
        Ok((atoms, builder))
    }

    #[instrument(skip_all, err)]
    pub fn new_with_timings<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        numerator_definitions: &[Arc<FnMapEntry>],
        orientations: &[EdgeVec<Orientation>],
        production_orientation_ids: &[OrientationID],
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        Self::new_with_rows(
            atoms,
            param_builder,
            numerator_definitions,
            orientations,
            production_orientation_ids,
            dual_shape,
            settings,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_with_rows<A: AtomCore>(
        atoms: &[A],
        param_builder: &ParamBuilder,
        numerator_definitions: &[Arc<FnMapEntry>],
        orientations: &[EdgeVec<Orientation>],
        production_orientation_ids: &[OrientationID],
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
        parametric_rows: Option<(ParametricResidueRows, Vec<Atom>)>,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
        let (mut parametric_rows, parametric_guards) = parametric_rows.unzip();
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
            inline_numerator_functions = settings.inline_numerator_functions,
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
        let alias_symbol = loop {
            let scope = NETWORK_SCALAR_ALIAS_SCOPE.fetch_add(1, Ordering::Relaxed);
            let name = format!("gammalooprs::evaluator_scalar_{scope}");
            if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
            }) {
                break symbol!(name.as_str());
            }
        };
        // Borrow the original physical definitions, before preprocessing appends
        // generated numerator components to its cloned builder.
        let physical_functions = &param_builder.reps;
        let (parsed_atoms, prepared_builder) = Self::preprocess_numerator_families(
            atoms,
            param_builder,
            numerator_definitions,
            settings,
            alias_symbol,
            // Symbolica's hyperdual vectorizer requires inlined function bodies.
            if dual_shape.is_some() || settings.inline_numerator_functions {
                InliningPolicy::Always
            } else {
                InliningPolicy::Never
            },
        )?;
        let parsed_atoms = parsed_atoms
            .into_iter()
            .map(|atom| Self::share_scalar_factors(atom, &prepared_builder))
            .collect::<Vec<_>>();
        let parsed_atoms = match &parametric_rows {
            Some(rows) => vec![
                rows.join_preprocessed(
                    parsed_atoms,
                    parametric_guards
                        .as_deref()
                        .expect("parametric row guards are present"),
                ),
            ],
            None => parsed_atoms,
        };
        let param_builder = &prepared_builder;
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
        // Specialization is deliberately below the shared tensor preprocessing
        // boundary. Explicit evaluators consume this derived view, including the
        // argument-free summed function map; SingleParametric retains the
        // original generic products.
        let specialized_atoms = if (settings.iterative_orientation_optimization
            || settings.summed
            || settings.summed_function_map)
            && let Some(rows) = &parametric_rows
        {
            Some(
                parsed_atoms
                    .iter()
                    .map(|atom| {
                        rows.rows
                            .iter()
                            .enumerate()
                            .map(|(id, row)| rows.specialize(atom, row, id))
                            .try_fold(AliasedAtom::default(), |sum, atom| {
                                sum.try_add(&atom).map_err(|error| eyre!(error))
                            })
                    })
                    .collect::<Result<Vec<_>>>()?,
            )
        } else {
            None
        };
        let explicit_atoms = specialized_atoms.as_deref().unwrap_or(&parsed_atoms);
        let iterative_started = std::time::Instant::now();
        let iterative = if settings.iterative_orientation_optimization {
            Some(
                Self::new_iterative(
                    explicit_atoms,
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
                    explicit_atoms,
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
                    explicit_atoms,
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
        // The optional variants have finished borrowing these scalars. Transfer
        // them to the final evaluator instead of retaining a second full copy.
        let (single_parametric, single_parametric_auxiliary) = Self::new_single_parametric(
            parsed_atoms,
            param_builder,
            physical_functions,
            parametric_rows.as_mut(),
            &dual_shape,
            settings,
        )
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
            atom_count = atoms.len(),
            orientation_count = orientations.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            spenso_ms = timings.spenso_time.as_secs_f64() * 1000.0,
            symbolica_ms = timings.symbolica_time.as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );

        Ok((
            EvaluatorStack {
                explicit_orientation_sum_only: false,
                parametric_rows,
                production_orientation_ids: production_orientation_ids.to_vec(),
                single_parametric,
                single_parametric_auxiliary,
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
    ) -> Vec<DualOrNot<Complex<F<T>>>>
    where
        usize: From<OID>,
    {
        let mut result: Option<Vec<DualOrNot<Complex<F<T>>>>> = None;
        for (orientation_id, e) in orientations.iter() {
            input.set_residue_map_id(self.production_orientation_ids[usize::from(orientation_id)]);
            input.set_orientation_values(e);
            let output = self.evaluate_single_parametric(input.as_slice(), evaluation_metadata);
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

    fn evaluate_single_parametric<T: FloatLike>(
        &mut self,
        input: &[Complex<F<T>>],
        metadata: &mut EvaluationMetaData,
    ) -> Vec<DualOrNot<Complex<F<T>>>> {
        let prepared_input;
        let mut auxiliary_output = Vec::new();
        let input = if let Some((ports, evaluator)) = &mut self.single_parametric_auxiliary {
            let rows = self
                .parametric_rows
                .as_ref()
                .expect("single-parametric auxiliary requires coefficient rows");
            let shape_len = evaluator.dual_shape.as_ref().map_or(1, Vec::len);
            let output_size = evaluator.compute_out_size();
            let prototype = input
                .first()
                .map_or_else(|| F(T::new_zero()), |value| value.re.clone());
            let zero = Complex::new_re(prototype.zero());
            let mut auxiliary_input = input.to_vec();
            auxiliary_input.resize(
                input.len() + rows.auxiliary_values.len() * shape_len,
                zero.clone(),
            );
            for (index, value) in rows.auxiliary_values.iter().enumerate() {
                auxiliary_input[input.len() + index * shape_len] =
                    Complex::new_re(prototype.from_rational(value));
            }
            auxiliary_output.resize(output_size, zero.clone());
            let start = std::time::Instant::now();
            <T as GenericEvaluatorFloat>::evaluate_into(
                evaluator,
                &auxiliary_input,
                &mut auxiliary_output,
            );
            metadata.evaluator_evaluation_time = metadata
                .evaluator_evaluation_time
                .saturating_add(start.elapsed());
            let mut prepared = Vec::with_capacity(input.len() + ports.len() * shape_len);
            prepared.extend_from_slice(input);
            // Ports have their own identity in the generic carrier, even when
            // several ports gather the same native auxiliary output. Copy the
            // row-invariant prefix once and let the original row loop gather
            // the remaining full jets in its original accumulation order.
            prepared
                .extend_from_slice(&auxiliary_output[..rows.auxiliary_shared_outputs * shape_len]);
            prepared.resize(input.len() + ports.len() * shape_len, zero);
            prepared_input = prepared;
            prepared_input.as_slice()
        } else {
            input
        };
        match &self.parametric_rows {
            Some(rows) => rows.evaluate(
                &mut self.single_parametric,
                input,
                &auxiliary_output,
                metadata,
            ),
            None => evaluate_evaluator(&mut self.single_parametric, input, metadata),
        }
    }

    fn evaluate_iterative<'a, T: FloatLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        evaluation_metadata: &mut EvaluationMetaData,
    ) -> Result<Vec<DualOrNot<Complex<F<T>>>>> {
        let Some((iterative, len)) = &mut self.iterative else {
            return Err(eyre!(
                "Iterative evaluator not available. Regenerate with iterative set to true."
            ));
        };

        let output = evaluate_evaluator(iterative, input.as_slice(), evaluation_metadata);
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
        ))
    }

    fn evaluate_summed<'a, T: FloatLike>(
        &'a mut self,
        input: InputParams<'a, T>,
        evaluation_metadata: &mut EvaluationMetaData,
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

            // The complete expression owns each independent source sum. Its
            // internal numerator calls retain their small coefficient rows;
            // no production-orientation loop may broadcast UV child sums.
            return match settings.general.evaluator_method {
                EvaluatorMethod::SingleParametric => {
                    Ok(self.evaluate_single_parametric(input.as_slice(), evaluation_metadata))
                }
                EvaluatorMethod::Iterative => self.evaluate_iterative(input, evaluation_metadata),
                EvaluatorMethod::SummedFunctionMap => {
                    self.evaluate_summed_fnmap(input, evaluation_metadata)
                }
                EvaluatorMethod::Summed => self.evaluate_summed(input, evaluation_metadata),
            };
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
            EvaluatorMethod::SingleParametric => {
                Ok(self.evaluate_parametric(input, orientations, evaluation_metadata))
            }
            EvaluatorMethod::Iterative => self.evaluate_iterative(input, evaluation_metadata),
            EvaluatorMethod::SummedFunctionMap => {
                self.evaluate_summed_fnmap(input, evaluation_metadata)
            }
            EvaluatorMethod::Summed => self.evaluate_summed(input, evaluation_metadata),
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

        if let Some((_, auxiliary)) = &mut self.single_parametric_auxiliary {
            let auxiliary_name = format!("{name}_single_parametric_auxiliary");
            auxiliary.compile_external(
                path.as_ref().join(&auxiliary_name).with_extension("cpp"),
                &auxiliary_name,
                path.as_ref().join(&auxiliary_name).with_extension("so"),
                frozen_mode,
            )?;
        }

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
        if let Some((_, auxiliary)) = &mut self.single_parametric_auxiliary {
            f(auxiliary)?;
        }

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

/// Dual shape and statically zero `(parameter, derivative component)` seeds.
type EvaluatorDualConfig = (Vec<Vec<usize>>, Vec<(usize, usize)>);

#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenericEvaluator {
    /// Stored roots when `store_atom` is enabled; evaluate them with the
    /// associated function map, including the definitions in `fn_map_entries`.
    pub exprs: Option<Vec<Atom>>,
    /// Stored function-map entries, including retained scalar alias definitions.
    pub fn_map_entries: Vec<FnMapEntry>,
    pub exprs_len: usize,
    pub backend_policy: EvaluatorBackendPolicy,
    pub rational: Option<ExpressionEvaluator<symbolica::domains::float::Complex<Rational>>>,
    pub f64_compiled: Option<CompiledCode<Complex<f64>>>,
    pub f64_eager: ExpressionEvaluator<Complex<F<f64>>>,
    pub f128: ExpressionEvaluator<Complex<F<f128>>>,
    pub dual_shape: Option<Vec<Vec<usize>>>,
    /// Derivative inputs known to vanish; these also protect singular unused derivatives.
    pub zero_components: Vec<(usize, usize)>,
    pub arb: ExpressionEvaluator<Complex<F<ArbPrec>>>,
    /// Only sampling programs warm this source lane; physical evaluators keep it empty.
    pub(crate) sampling_fixed256: RuntimeCache<ExpressionEvaluator<Complex<F<SamplingFloat>>>>,
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

    pub(crate) fn activate_symjit(&mut self, options: &CompilationOptionsSnapshot) -> Result<()> {
        if self.is_eager_only() {
            self.activate_eager_only();
            return Ok(());
        }

        // Use the same numeric program as eager and external compilation. Rational
        // constant slots (such as pi) remain placeholders until domain mapping.
        // SymJIT cannot compact some complex temporary layouts. Direct translation
        // emits real branches and forces O2. The configured level (capped at O2)
        // and cse=false only affect indirect translation: its function-result cache
        // crosses inactive branches, and optimization can drop stores that later
        // calls still read. Symbolica already shares expressions in this program.
        let program = self
            .f64_eager
            .clone()
            .map_coeff(&|c| SymComplex::new(c.re.0, c.im.0));
        let evaluator = super::standalone::compile_symjit(
            &program,
            JITCompilationSettings::new()
                .direct_translation(options.jit_direct_translation)
                .optimization_level(usize::from(options.optimization_level).min(2) as u8)
                .with_option("compact", "false")
                .with_option("cse", "false"),
        )?;
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

    pub(crate) fn new_from_builder<I: IntoIterator<Item: Into<AliasedAtom>>>(
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
            dual_shape.map(|shape| (shape, Vec::new())),
            settings,
        )
    }

    pub(crate) fn new_from_raw_params<I: IntoIterator<Item: Into<AliasedAtom>>>(
        atoms: I,
        params: &[Atom],
        fn_map: &FunctionMap,
        mut fn_map_entries: Vec<FnMapEntry>,
        optimization_settings: OptimizationSettings,
        dual_config: Option<EvaluatorDualConfig>,
        settings: &EvaluatorSettings,
    ) -> Result<Self> {
        // Known-zero seed components belong to the compiled program. Supplying
        // them to Symbolica avoids evaluating a zero tangent times a singular
        // derivative of a prepared-only subexpression, such as sqrt(m) at m=0.
        // The serialized shape and input layout stay unchanged; simplification
        // is retained in the rational program and every native specialization.
        let (dual_shape, zero_components) =
            dual_config.map_or((None, Vec::new()), |(shape, zeros)| (Some(shape), zeros));
        if let Some(shape) = &dual_shape
            && zero_components.iter().any(|&(parameter, component)| {
                parameter >= params.len() || component == 0 || component >= shape.len()
            })
        {
            return Err(eyre!(
                "statically zero dual seeds must refer to valid derivative components"
            ));
        }
        let evaluator_replacements = if settings.do_fn_map_replacements {
            fn_map_entries
                .iter()
                .filter(|entry| entry.inlining != InliningPolicy::Never)
                .map(FnMapEntry::replacement)
                .collect::<Vec<_>>()
        } else {
            Vec::new()
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

        for replacement in &evaluator_replacements {
            crate::debug_tags!(#generation, #compile, #term, #dump;
                stage = "evaluator_function_map_replacement",
                file.replacement = %replacement,
                "Evaluator function-map replacement"
            );
        }

        let preparation_started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_expression_preparation_start",
            replacement_count = evaluator_replacements.len(),
            "Evaluator timing milestone"
        );
        let atoms: Vec<AliasedAtom> = atoms.into_iter().map(Into::into).collect();
        // A fresh symbol also avoids collisions with definitions supplied only
        // through the opaque FunctionMap.
        let alias_symbol = loop {
            let scope = NETWORK_SCALAR_ALIAS_SCOPE.fetch_add(1, Ordering::Relaxed);
            let name = format!("gammalooprs::evaluator_retained_scalar_{scope}");
            if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
            }) {
                break symbol!(name.as_str());
            }
        };
        let exprs: Vec<AliasedAtom> = atoms
            .into_iter()
            .enumerate()
            .map(|(atom_index, atom)| {
                let (root, aliases) = atom.into_inner_with_aliases();
                let mut handles = aliases.keys().cloned().collect::<Vec<_>>();
                handles.sort();
                let renames = handles
                    .into_iter()
                    .enumerate()
                    .map(|(index, handle)| (handle, function!(alias_symbol, atom_index, index)))
                    .collect::<HashMap<_, _>>();
                // Alias keys may have any Atom kind. Their maximum byte length
                // excludes impossible matches without hashing large subtrees.
                let max_handle_bytes = renames
                    .keys()
                    .map(|handle| handle.as_view().get_data().len())
                    .max()
                    .unwrap_or(0);
                let map = |atom: Atom| {
                    let atom = if renames.is_empty() {
                        atom
                    } else {
                        atom.replace_map(|view, _, out| {
                            if view.get_data().len() <= max_handle_bytes
                                && let Some(scoped) = renames.get(view.get_data())
                            {
                                out.set_from_view(&scoped.as_view());
                            }
                        })
                    };
                    // An empty replacement pass still copies the entire Atom.
                    // Preserve ownership when function-map substitution is disabled.
                    if evaluator_replacements.is_empty() {
                        atom
                    } else {
                        atom.replace_multiple(&evaluator_replacements)
                            .replace_multiple(&evaluator_replacements)
                    }
                };
                let mut retained = AliasedAtom::from(map(root));
                for (alias, body) in aliases {
                    let rhs = if body == alias { body } else { map(body) };
                    retained.register_alias(renames[&alias].clone(), rhs);
                }
                retained
            })
            .collect();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_expression_preparation_done",
            atom_count = exprs.len(),
            elapsed_ms = preparation_started.elapsed().as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );
        // M remains an ordinary runtime input even when `store_atom` is disabled
        // or the emitted sources do not use it. Runtime validation requires a
        // nonzero value without inspecting or expanding shared function bodies.

        let recorded_function_body_bytes = fn_map_entries
            .iter()
            .map(|entry| entry.rhs.as_view().get_byte_size())
            .sum::<usize>();
        let max_recorded_function_body_bytes = fn_map_entries
            .iter()
            .map(|entry| entry.rhs.as_view().get_byte_size())
            .max()
            .unwrap_or(0);
        let shared_function_count = fn_map_entries
            .iter()
            .filter(|entry| entry.inlining == InliningPolicy::Never && !entry.is_alias)
            .count();
        let retained_dual = dual_shape.is_some()
            && fn_map_entries
                .iter()
                .any(|entry| entry.inlining == InliningPolicy::Never);
        let mut tree: Option<ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>>> =
            if retained_dual {
                use super::retained_dual::{RetainedFunctionDefinition, build_dual_evaluator};
                let mut definitions = fn_map_entries
                    .iter()
                    .map(|entry| RetainedFunctionDefinition {
                        lhs: entry.lhs.clone(),
                        rhs: entry.rhs.clone(),
                        tags: entry.tags.clone(),
                        args: entry.args.clone(),
                        retained: entry.inlining == InliningPolicy::Never,
                        is_alias: entry.is_alias,
                    })
                    .collect::<Vec<_>>();
                let mut dual_fn_map = fn_map.clone();
                for expression in &exprs {
                    for (alias, body) in expression.get_aliases() {
                        dual_fn_map
                            .add_aliases([(alias.clone(), body.clone())])
                            .map_err(|error| eyre!(error))?;
                        definitions.push(RetainedFunctionDefinition {
                            lhs: alias.clone(),
                            rhs: body.clone(),
                            tags: alias
                                .as_fun_view()
                                .unwrap()
                                .iter()
                                .map(|tag| tag.to_owned())
                                .collect(),
                            args: Vec::new(),
                            retained: false,
                            is_alias: true,
                        });
                    }
                }
                let roots = exprs
                    .iter()
                    .map(|expression| expression.get_root().clone())
                    .collect::<Vec<_>>();
                Some(
                    build_dual_evaluator(
                        &roots,
                        params,
                        &dual_fn_map,
                        &definitions,
                        dual_shape.clone().unwrap(),
                        zero_components.clone(),
                        optimization_settings.clone(),
                    )
                    .map_err(|error| {
                        eyre!("Failed to dualize retained numerator programs: {error}")
                    })?,
                )
            } else {
                None
            };
        for (atom_index, n) in exprs.iter().enumerate().filter(|_| !retained_dual) {
            // Optional source capture includes recorded function bodies as well
            // as the aliased outer atom. Keep all numerators factorized.
            crate::debug_tags!(#generation, #compile, #evaluator_input, #dump;
                stage = "evaluator_symbolica_input_dump",
                atom_index,
                file.evaluator_input = %serde_json::to_string(&serde_json::json!({
                    "root": n.get_root().to_canonical_string(),
                    "aliases": n.get_aliases().iter().map(|(key, body)|
                        (key.to_canonical_string(), body.to_canonical_string())).collect::<Vec<_>>(),
                    "functions": fn_map_entries.iter().map(|entry| entry.archive::<String>()).collect::<Result<Vec<_>>>()?,
                    "parameters": params.iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>(),
                    "dual_shape": dual_shape,
                    "optimization_settings": optimization_settings,
                }))?,
                "Factorized evaluator source before construction"
            );
        }
        if !retained_dual && !exprs.is_empty() {
            let build_started = std::time::Instant::now();
            crate::debug_tags!(#generation, #profile, #compile, #summary;
                stage = "evaluator_symbolica_build_start",
                atom_count = exprs.len(),
                atom_bytes = exprs.iter().map(AliasedAtom::get_byte_size).sum::<usize>(),
                root_bytes = exprs.iter().map(|atom| atom.get_root().as_view().get_byte_size()).sum::<usize>(),
                scalar_alias_count = exprs.iter().map(|atom| atom.get_aliases().len()).sum::<usize>(),
                recorded_function_count = fn_map_entries.len(),
                recorded_function_body_bytes,
                max_recorded_function_body_bytes,
                shared_function_count,
                dual = dual_shape.is_some(),
                "Evaluator timing milestone"
            );
            // The aliases have already been scoped per output above. Register
            // them after the function map, whose setter replaces existing entries.
            // Build all ordered outputs together instead of repeatedly optimizing
            // and merging the growing program; numerator bodies remain shared.
            let roots = exprs.iter().map(AliasedAtom::get_root).collect::<Vec<_>>();
            let eval: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> =
                Atom::evaluator_multiple(&roots, params)
                    .function_map(fn_map.clone())
                    .add_aliases(exprs.iter().flat_map(|atom| {
                        atom.get_aliases()
                            .iter()
                            .map(|(alias, body)| (alias.clone(), body.clone()))
                    }))
                    .map_err(|error| eyre!(error))?
                    .optimization_settings(optimization_settings.clone())
                    .build()
                    .map_err(|error| {
                        // Keep output indices in an error, without formatting all
                        // source expressions during successful construction.
                        let outputs = exprs
                            .iter()
                            .enumerate()
                            .map(|(index, atom)| {
                                let mut settings =
                                    SpensoPrintSettings::compact().nice_symbolica();
                                settings.max_line_length = Some(120);
                                settings.hide_all_namespaces = false;
                                format!("output[{index}]: {:120}", atom.printer(settings))
                            })
                            .collect::<Vec<_>>()
                            .join("\n");
                        eyre!(
                            "Failed to create evaluator for {} ordered outputs: {error}\n{outputs}\nwith params: {}, and fn_map_entries: {}",
                            exprs.len(),
                            params
                                .iter()
                                .map(|atom| atom.to_string())
                                .collect::<Vec<_>>()
                                .join(", "),
                            fn_map_entries
                                .iter()
                                .map(|entry| format!("{}->{}\n", entry.lhs, entry.rhs))
                                .collect::<Vec<_>>()
                                .join(", "),
                        )
                    })?;

            let operations = eval.count_operations();
            crate::debug_tags!(#generation, #profile, #compile, #summary;
                stage = "evaluator_symbolica_build_done",
                atom_count = exprs.len(),
                root_additions = operations.additions,
                root_multiplications = operations.multiplications,
                root_inversions = operations.inversions,
                root_function_calls = operations.function_calls,
                elapsed_ms = build_started.elapsed().as_secs_f64() * 1000.0,
                "Evaluator timing milestone"
            );
            tree = Some(eval);
        }

        let mut tree = tree.ok_or_else(|| eyre!("No expressions to evaluate"))?;
        let exprs_len = exprs.len();
        // The program owns its data. Release unstored source expressions before
        // constructing dual and numeric programs, which can themselves be large.
        let exprs = if settings.store_atom {
            Some(
                exprs
                    .into_iter()
                    .map(|atom| {
                        let (root, aliases) = atom.into_inner_with_aliases();
                        for (alias, rhs) in aliases {
                            fn_map_entries.push(FnMapEntry {
                                inlining: InliningPolicy::Always,
                                is_alias: true,
                                tags: alias
                                    .as_fun_view()
                                    .unwrap()
                                    .iter()
                                    .map(|arg| arg.to_owned())
                                    .collect(),
                                lhs: alias,
                                rhs,
                                args: Vec::new(),
                            });
                        }
                        root
                    })
                    .collect(),
            )
        } else {
            drop(exprs);
            None
        };
        drop(fn_map);

        let domains_started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_numeric_programs_start",
            "Evaluator timing milestone"
        );
        if let Some(dual_shape) = &dual_shape
            && !retained_dual
        {
            let dual = HyperDual::<SymComplex<Rational>>::new(dual_shape.clone());
            let dualizer = Dualizer::new(dual, zero_components.clone());
            tree = tree
                .vectorize(&dualizer)
                .map_err(|error| eyre!("Failed to vectorize evaluator: {error}"))?;
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
        crate::debug_tags!(#generation, #profile, #compile, #summary;
            stage = "evaluator_numeric_programs_done",
            elapsed_ms = domains_started.elapsed().as_secs_f64() * 1000.0,
            "Evaluator timing milestone"
        );

        let evaluator = GenericEvaluator {
            exprs_len,
            fn_map_entries,
            exprs,
            backend_policy: EvaluatorBackendPolicy::FollowIntegrand,
            rational: Some(rational),
            f64_compiled: None,
            f64_eager,
            f128,
            dual_shape,
            zero_components,
            arb,
            sampling_fixed256: RuntimeCache::default(),
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
        let index = self.residue_map_id_start * self.multiplicative_offset;
        let value = Complex::new_re(self.as_slice()[index].re.from_usize(id.0));
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

    /// Write all scalar or flattened dual outputs into a caller-owned buffer.
    fn evaluate_into(
        generic_evaluator: &mut GenericEvaluator,
        params: &[Complex<F<T>>],
        out: &mut [Complex<F<T>>],
    );

    #[allow(clippy::type_complexity)]
    fn get_evaluator(
        generic_evaluator: &mut GenericEvaluator,
    ) -> impl FnMut(&[Complex<F<T>>]) -> Vec<DualOrNot<Complex<F<T>>>> {
        let out_size = generic_evaluator.compute_out_size();
        let mut out = vec![Complex::default(); out_size];
        let dual_builder = generic_evaluator
            .dual_shape
            .as_ref()
            .map(|shape| HyperDual::<Complex<F<T>>>::new(shape.clone()));
        move |params: &[Complex<F<T>>]| {
            out.resize_with(out_size, Complex::default);
            Self::evaluate_into(generic_evaluator, params, &mut out);
            if let Some(dual_builder) = &dual_builder {
                out.chunks(dual_builder.values.len())
                    .map(|chunk| DualOrNot::Dual(new_from_values(dual_builder, chunk)))
                    .collect()
            } else {
                std::mem::take(&mut out)
                    .into_iter()
                    .map(DualOrNot::NonDual)
                    .collect()
            }
        }
    }

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

    fn evaluate_into(
        generic_evaluator: &mut GenericEvaluator,
        params: &[Complex<F<f64>>],
        out: &mut [Complex<F<f64>>],
    ) {
        assert_eq!(
            out.len(),
            generic_evaluator.compute_out_size(),
            "Incorrect evaluator output length"
        );
        match generic_evaluator.active_f64_backend() {
            ActiveF64Backend::Eager => generic_evaluator.f64_eager.evaluate(params, out),
            ActiveF64Backend::Cpp | ActiveF64Backend::Assembly => {
                let compiled = generic_evaluator
                    .loaded_f64_compiled
                    .as_mut()
                    .expect("compiled f64 backend should be activated before evaluation");
                unsafe {
                    compiled.evaluate(
                        transmute::<&[Complex<F<f64>>], &[Complex<f64>]>(params),
                        transmute::<&mut [Complex<F<f64>>], &mut [Complex<f64>]>(out),
                    );
                }
            }
            ActiveF64Backend::Symjit => generic_evaluator
                .symjit_f64
                .as_mut()
                .expect("symjit f64 backend should be activated before evaluation")
                .evaluate(params, out),
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

// Eager scalar and dual dispatch share the same output shape at every native lane.
macro_rules! impl_eager_evaluator_float {
    ($scalar:ty, $evaluator:ident => $program:expr) => {
        impl GenericEvaluatorFloat for $scalar {
            #[inline(always)]
            fn get_evaluator_single(
                $evaluator: &mut GenericEvaluator,
            ) -> impl FnMut(&[Complex<F<Self>>]) -> Complex<F<Self>> {
                // info!("USING COMPLEX EAGER SINGLE");
                #[inline(always)]
                |params: &[Complex<F<Self>>]| ($program).evaluate_single(params)
            }

            fn evaluate_into(
                $evaluator: &mut GenericEvaluator,
                params: &[Complex<F<Self>>],
                out: &mut [Complex<F<Self>>],
            ) {
                assert_eq!(
                    out.len(),
                    $evaluator.compute_out_size(),
                    "Incorrect evaluator output length"
                );
                ($program).evaluate(params, out);
            }

            fn get_parameters<'a>(
                param_builder: &'a mut ParamBuilder,
                cache: (bool, bool),
                graph: &'a Graph,
                sample: &'a MomentumSample<Self>,
                helicities: &[Helicity],
                additional_params: &[F<Self>],
                left_threshold_params: Option<&ThresholdParams<Self>>,
                right_threshold_params: Option<&ThresholdParams<Self>>,
                lu_params: Option<&LUParams<Self>>,
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
    };
}

impl_eager_evaluator_float!(f128, evaluator => &mut evaluator.f128);
impl_eager_evaluator_float!(ArbPrec, evaluator => &mut evaluator.arb);
impl_eager_evaluator_float!(SamplingFloat, evaluator => evaluator.sampling_fixed256.as_mut()
    .expect("fixed256 sampling evaluator must be prepared during warmup"));

#[cfg(test)]
mod hyperbolic_tests;

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use typed_index_collections::TiVec;

    use idenso::{dirac::AGS, representations::Bispinor};
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{network::tags::SPENSO_TAG, structure::representation::RepName};
    use symbolica::{atom::Symbol, parse_lit, state::State};

    use crate::{
        GammaLoopContextContainer, initialisation::test_initialise,
        integrands::process::param_builder::ParamValuePairs, model::Model,
        settings::global::CompilationOptimizationLevel,
    };

    use super::*;

    fn scalar_value(result: Vec<DualOrNot<Complex<F<f64>>>>) -> Complex<F<f64>> {
        let [DualOrNot::NonDual(value)] = result.as_slice() else {
            panic!("expected one scalar evaluator output")
        };
        *value
    }

    #[test]
    fn generated_numerator_bodies_remain_shared_through_rebuild_and_serialization() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::retained_body_family");
        let q = symbol!("evaluator_test::retained_body_formal");
        let x = parse_lit!(evaluator_test::retained_body_x);
        let g = parse_lit!(evaluator_test::retained_body_global);
        let definition = Arc::new(FnMapEntry {
            lhs: function!(family, 0, q),
            rhs: Atom::add_many((1..=16).map(|i| (Atom::var(q) + &g + i).pow(-1))),
            args: vec![q.into()],
            tags: vec![Atom::num(0)],
            inlining: InliningPolicy::Always,
            is_alias: false,
        });
        let source = Atom::add_many((1..=32).map(|j| function!(family, 0, &x + j)));
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone(), g.clone()].into_iter().collect();
        builder.pairs.update_ranges();
        let expected: f64 = (1..=32)
            .flat_map(|j| (1..=16).map(move |i| 1.0 / (3 + i + j) as f64))
            .sum();
        for do_fn_map_replacements in [false, true] {
            let settings = EvaluatorSettings {
                store_atom: true,
                inline_numerator_functions: false,
                do_fn_map_replacements,
                ..Default::default()
            };
            let (stack, _) = EvaluatorStack::new_explicit_sum_with_timings(
                std::slice::from_ref(&source),
                &builder,
                std::slice::from_ref(&definition),
                None,
                &settings,
            )
            .unwrap();
            let mut evaluator = stack.single_parametric;
            assert!(
                evaluator
                    .fn_map_entries
                    .iter()
                    .any(|entry| entry.inlining == InliningPolicy::Never)
            );
            let instructions = evaluator.f64_eager.export_instructions();
            assert!(!instructions.sub_evaluators.is_empty());
            let shared_count = instructions.sub_evaluators.len();
            let root_count = instructions.instructions.len();
            drop(instructions);
            let roots = evaluator.exprs.as_ref().unwrap().clone();
            let mut function_map = FunctionMap::default();
            for entry in &evaluator.fn_map_entries {
                entry.register(&mut function_map).unwrap();
            }
            let mut rebuilt = GenericEvaluator::new_from_raw_params(
                roots.clone(),
                &[x.clone(), g.clone()],
                &function_map,
                evaluator.fn_map_entries.clone(),
                settings.optimization_settings(),
                None,
                &settings,
            )
            .unwrap();
            let mut inline_map = FunctionMap::default();
            let mut inline_entries = evaluator.fn_map_entries.clone();
            for entry in &mut inline_entries {
                entry.inlining = InliningPolicy::Always;
                entry.register(&mut inline_map).unwrap();
            }
            let mut inline = GenericEvaluator::new_from_raw_params(
                roots,
                &[x.clone(), g.clone()],
                &inline_map,
                inline_entries,
                settings.optimization_settings(),
                None,
                &EvaluatorSettings {
                    do_fn_map_replacements: false,
                    ..settings
                },
            )
            .unwrap();
            let inline_instructions = inline.f64_eager.export_instructions();
            assert!(inline_instructions.sub_evaluators.is_empty());
            assert!(root_count < inline_instructions.instructions.len());
            let (inline_stack, _) = EvaluatorStack::new_explicit_sum_with_timings(
                std::slice::from_ref(&source),
                &builder,
                std::slice::from_ref(&definition),
                None,
                &EvaluatorSettings {
                    store_atom: true,
                    do_fn_map_replacements: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let mut configured_inline = inline_stack.single_parametric;
            assert!(
                configured_inline
                    .fn_map_entries
                    .iter()
                    .all(|entry| entry.inlining == InliningPolicy::Always)
            );
            assert!(
                configured_inline
                    .f64_eager
                    .export_instructions()
                    .sub_evaluators
                    .is_empty()
            );
            let encoded = bincode::encode_to_vec(&evaluator, bincode::config::standard()).unwrap();
            let mut state = Vec::new();
            State::export(&mut state).unwrap();
            let state_map = State::import(&mut Cursor::new(state), None).unwrap();
            let model = Model::default();
            let (mut decoded, _): (GenericEvaluator, _) = bincode::decode_from_slice_with_context(
                &encoded,
                bincode::config::standard(),
                GammaLoopContextContainer {
                    state_map: &state_map,
                    model: &model,
                },
            )
            .unwrap();
            let values = [Complex::new_re(F(1.0)), Complex::new_re(F(2.0))];
            for shared in [&mut evaluator, &mut rebuilt, &mut decoded] {
                assert_eq!(
                    shared.f64_eager.export_instructions().sub_evaluators.len(),
                    shared_count
                );
                for compiled in [false, true] {
                    if compiled {
                        shared
                            .activate_symjit(&CompilationOptionsSnapshot {
                                optimization_level: CompilationOptimizationLevel::O2,
                                ..Default::default()
                            })
                            .unwrap();
                    }
                    let actual =
                        <f64 as GenericEvaluatorFloat>::get_evaluator_single(shared)(&values);
                    assert!((actual.re.0 - expected).abs() < 1e-12);
                    assert_eq!(actual.im.0, 0.0);
                }
                let precise = <f128 as GenericEvaluatorFloat>::get_evaluator_single(shared)(&[
                    Complex::new_re(F::<f128>::from_f64(1.0)),
                    Complex::new_re(F::<f128>::from_f64(2.0)),
                ]);
                assert!((precise.re.to_f64() - expected).abs() < 1e-12);
                let arbitrary =
                    <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(shared)(&[
                        Complex::new_re(F::<ArbPrec>::from_f64(1.0)),
                        Complex::new_re(F::<ArbPrec>::from_f64(2.0)),
                    ]);
                assert!((arbitrary.re.to_f64() - expected).abs() < 1e-12);
            }
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut inline)(&values);
            assert!((actual.re.0 - expected).abs() < 1e-12);
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                &mut configured_inline,
            )(&values);
            assert!((actual.re.0 - expected).abs() < 1e-12);
            assert_eq!(actual.im.0, 0.0);
            let zero = F(ArbPrec::default());
            let arbitrary_values = [
                Complex::new_re(zero.one()),
                Complex::new_re(zero.from_usize(2)),
            ];
            let retained_arbitrary = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(
                &mut evaluator,
            )(&arbitrary_values);
            let inline_arbitrary = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(
                &mut configured_inline,
            )(&arbitrary_values);
            assert!(
                (inline_arbitrary.re - retained_arbitrary.re).abs()
                    < zero.from_usize(10).powi(-290)
            );
            assert_eq!(inline_arbitrary.im, zero);
        }
        // The same production entry point keeps function bodies inlined for the
        // hyperdual backend, which cannot vectorize retained sub-evaluators.
        let (mut dual, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &[source],
            &builder,
            &[definition],
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            &EvaluatorSettings {
                store_atom: true,
                inline_numerator_functions: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            dual.single_parametric
                .fn_map_entries
                .iter()
                .all(|entry| entry.inlining == InliningPolicy::Always)
        );
        assert!(
            dual.single_parametric
                .f64_eager
                .export_instructions()
                .sub_evaluators
                .is_empty()
        );
        let values = [1.0, 1.0, 2.0, 0.0].map(|v| Complex::new_re(F(v)));
        let actual =
            <f64 as GenericEvaluatorFloat>::get_evaluator(&mut dual.single_parametric)(&values);
        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("expected dual result")
        };
        let derivative: f64 = (1..=32)
            .flat_map(|j| (1..=16).map(move |i| -1.0 / ((3 + i + j) as f64).powi(2)))
            .sum();
        assert!((actual.values[0].re.0 - expected).abs() < 1e-12);
        assert!((actual.values[1].re.0 - derivative).abs() < 1e-12);
    }

    #[test]
    fn shared_numerator_outer_projectors_prune_unreachable_component_dependencies() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::reachable_numerator_family");
        let q = Atom::var(symbol!("evaluator_test::reachable_numerator_q"));
        let x = Atom::var(symbol!("evaluator_test::reachable_numerator_x"));
        let vector = SPENSO_TAG.tensor_symbol("evaluator_test::reachable_numerator_vector");
        let original = symbol!("evaluator_test::reachable_original_function");
        let index = parse_lit!(spenso::mink(4, 1));
        let temporal = GS.energy_delta(index.as_view());
        let spatial = GS.emr_vec(EdgeIndex(7), index.as_view());
        let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(2)));
        assert!(weight.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES);
        let definitions = [
            // All four components exist, but the outer temporal projector uses one.
            &q * &weight * function!(vector, index),
            // The outer projector kills the entire large spatial summand, including
            // its populated components and their otherwise retained scalar aliases.
            &temporal + &q * weight * spatial,
        ]
        .into_iter()
        .enumerate()
        .map(|(i, rhs)| {
            Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, i, q.clone()),
                rhs,
                tags: vec![Atom::num(i)],
                args: vec![Indeterminate::try_from(q.clone()).unwrap()],
            })
        })
        .collect::<Vec<_>>();
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [
            x.clone(),
            function!(vector, function!(AIND_SYMBOLS.cind, 0)),
        ]
        .into_iter()
        .collect();
        builder.pairs.update_ranges();
        builder
            .add_tagged_function(
                original,
                vec![Atom::num(0)],
                Default::default(),
                Vec::<Indeterminate>::new(),
                x + 11,
            )
            .unwrap();
        let settings = EvaluatorSettings {
            store_atom: true,
            ..Default::default()
        };
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[
                function!(family, 0, 3) * &temporal,
                function!(family, 1, 3) * temporal,
                function!(original, 0),
            ],
            &builder,
            &definitions,
            &settings,
            symbol!("evaluator_test::reachable_numerator_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        assert_eq!(
            &prepared.reps[..builder.reps.len()],
            builder.reps.as_slice()
        );
        let generated = &prepared.reps[builder.reps.len()..];
        assert_eq!(
            generated
                .iter()
                .filter(|entry| entry
                    .lhs
                    .as_fun_view()
                    .unwrap()
                    .get_symbol()
                    .get_name()
                    .starts_with("gammalooprs::numerator_component_"))
                .count(),
            2
        );
        assert_eq!(
            generated
                .iter()
                .filter(
                    |entry| entry.rhs.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES
                )
                .count(),
            1
        );
        assert!(
            generated
                .iter()
                .any(|entry| entry.lhs == atoms[1].get_root()
                    && entry.rhs.is_one()
                    && entry.args.is_empty())
        );
        for component in 1..4 {
            assert!(generated.iter().all(|entry| {
                !entry
                    .rhs
                    .contains(&function!(vector, function!(AIND_SYMBOLS.cind, component)))
            }));
        }
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let w = (1..512).map(|i| i * i).sum::<i64>() as f64;
        let dw = (1..512).map(|i| 2 * i).sum::<i64>() as f64;
        let values = [0.0, 1.0, 7.0, 0.0].map(|value| Complex::new_re(F(value)));
        for compile in [false, true] {
            if compile {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O0,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(evaluator.active_f64_backend(), ActiveF64Backend::Symjit);
            }
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
            assert_eq!(actual.len(), 3);
            for (actual, expected) in
                actual
                    .iter()
                    .zip([[21.0 * w, 21.0 * dw], [1.0, 0.0], [11.0, 1.0]])
            {
                let DualOrNot::Dual(actual) = actual else {
                    panic!("expected dual result")
                };
                assert_eq!(
                    actual.values,
                    expected.map(|value| Complex::new_re(F(value))).to_vec()
                );
            }
        }
    }

    #[test]
    fn shared_numerator_families_skip_unused_bodies_and_validate_calls() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::filtered_numerator_family");
        let q = Atom::var(symbol!("evaluator_test::filtered_numerator_q"));
        let x = Atom::var(symbol!("evaluator_test::filtered_numerator_x"));
        let tensor = SPENSO_TAG.tensor_symbol("evaluator_test::unused_numerator_tensor");
        let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(2)));
        let definitions = [
            Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, 0, q.clone()),
                rhs: &q * &x,
                tags: vec![Atom::num(0)],
                args: vec![Indeterminate::try_from(q.clone()).unwrap()],
            }),
            Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, 1, q.clone()),
                rhs: weight
                    * function!(
                        tensor,
                        parse_lit!(spenso::mink(evaluator_test::unused_dimension, 1))
                    ),
                tags: vec![Atom::num(1)],
                args: vec![Indeterminate::try_from(q).unwrap()],
            }),
        ];
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x].into_iter().collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings::default();
        let alias = symbol!("evaluator_test::filtered_numerator_scalar");
        let (zero, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[Atom::Zero],
            &builder,
            &definitions,
            &settings,
            alias,
            InliningPolicy::Never,
        )
        .unwrap();
        assert!(zero[0].get_root().is_zero());
        assert!(prepared.reps.is_empty());
        let (live, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[function!(family, 0, 3)],
            &builder,
            &definitions,
            &settings,
            alias,
            InliningPolicy::Never,
        )
        .unwrap();
        assert!(
            prepared
                .reps
                .iter()
                .all(|entry| entry.tags[0] == Atom::num(0))
        );
        let mut evaluator = GenericEvaluator::new_from_builder(
            live,
            &prepared,
            None,
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        assert_eq!(
            <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&[
                Complex::new_re(F(2.0))
            ]),
            Complex::new_re(F(6.0))
        );
        // The unsupported dimension is harmless only while the body is unreferenced.
        assert!(
            EvaluatorStack::preprocess_numerator_families(
                &[function!(family, 1, 3)],
                &builder,
                &definitions,
                &settings,
                alias,
                InliningPolicy::Never,
            )
            .is_err()
        );
        for invalid in [
            function!(family, 99, 3),
            function!(family, 0),
            function!(family, 0, 3, 4),
        ] {
            assert!(
                EvaluatorStack::preprocess_numerator_families(
                    &[invalid],
                    &builder,
                    &definitions,
                    &settings,
                    alias,
                    InliningPolicy::Never,
                )
                .is_err()
            );
        }
        let mut conflict = (*definitions[1]).clone();
        conflict.rhs = Atom::Zero;
        assert!(
            EvaluatorStack::preprocess_numerator_families(
                &[Atom::Zero],
                &builder,
                &[definitions[1].clone(), Arc::new(conflict)],
                &settings,
                alias,
                InliningPolicy::Never,
            )
            .is_err()
        );
        let mut repeated_formal = (*definitions[0]).clone();
        repeated_formal.args.push(repeated_formal.args[0].clone());
        repeated_formal.lhs = function!(
            family,
            0,
            Atom::from(repeated_formal.args[0].clone()),
            Atom::from(repeated_formal.args[1].clone())
        );
        let error = EvaluatorStack::preprocess_numerator_families(
            &[Atom::Zero],
            &builder,
            &[Arc::new(repeated_formal)],
            &settings,
            alias,
            InliningPolicy::Never,
        )
        .err()
        .expect("duplicate formal arguments must be rejected even for unused families");
        assert!(
            error
                .to_string()
                .contains("formal arguments must be distinct")
        );
        let mut nested = (*definitions[0]).clone();
        nested.rhs = function!(family, 1, 3);
        assert!(
            EvaluatorStack::preprocess_numerator_families(
                &[function!(family, 0, 3)],
                &builder,
                &[Arc::new(nested), definitions[1].clone()],
                &settings,
                alias,
                InliningPolicy::Never,
            )
            .is_err()
        );
    }

    #[test]
    fn shared_numerator_call_lookup_binds_function_arguments_simultaneously() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::lookup_numerator_family");
        let other_family = symbol!("evaluator_test::lookup_other_numerator_family");
        let formal = symbol!("evaluator_test::lookup_numerator_argument"; Scalar);
        let nested = symbol!("evaluator_test::lookup_nested_argument"; Scalar);
        let q = function!(formal, 0);
        let r = function!(formal, 1);
        let nested_argument = function!(nested, &q + 1);
        let definitions = [(family, &q + 10 * &r), (other_family, 2 * &q + 3 * &r)]
            .into_iter()
            .map(|(head, rhs)| {
                Arc::new(FnMapEntry {
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                    lhs: function!(head, 7, &q, &r),
                    rhs,
                    tags: vec![Atom::num(7)],
                    args: [&q, &r]
                        .into_iter()
                        .map(|arg| Indeterminate::try_from(arg.clone()).unwrap())
                        .collect(),
                })
            })
            .collect::<Vec<_>>();
        let roots = [
            function!(family, 7, &r, &q),
            function!(family, 7, &nested_argument, &r),
            function!(other_family, 7, &r, &q),
            function!(family, 7, &nested_argument, &q),
        ];
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &roots,
            &ParamBuilder::new_empty(),
            &definitions,
            &EvaluatorSettings::default(),
            symbol!("evaluator_test::lookup_numerator_scalar"),
            InliningPolicy::Never,
        )
        .unwrap();
        let replacements = prepared
            .reps
            .iter()
            .map(FnMapEntry::replacement)
            .collect::<Vec<_>>();
        let expected = [
            &r + 10 * &q,
            &nested_argument + 10 * &r,
            2 * &r + 3 * &q,
            nested_argument + 10 * &q,
        ];
        for (atom, expected) in atoms.iter().zip(expected) {
            // Swapping formals and inserting an argument that contains another
            // formal must not substitute into an already inserted argument.
            assert_eq!(atom.get_root().replace_multiple(&replacements), expected);
        }
    }

    #[test]
    fn shared_numerator_components_compact_before_pruning() {
        use spenso::network::library::symbolic::ETS;

        test_initialise().unwrap();
        let family = symbol!("evaluator_test::compact_numerator_family");
        let q = Atom::var(symbol!("evaluator_test::compact_numerator_q"));
        let r = Atom::var(symbol!("evaluator_test::compact_numerator_r"));
        let x = Atom::var(symbol!("evaluator_test::compact_numerator_x"));
        let y = Atom::var(symbol!("evaluator_test::compact_numerator_y"));
        let mu = parse_lit!(spenso::mink(4, 1));
        let nu = parse_lit!(spenso::mink(4, 2));
        let contraction = GS.energy_delta(mu.as_view())
            * GS.energy_delta(nu.as_view())
            * function!(ETS.metric, &mu, &nu);
        let contractions = [
            contraction.clone(),
            contraction,
            GS.delta_vec.call_args([GS.cind(1), mu.clone()])
                * GS.energy_delta(nu.as_view())
                * function!(ETS.metric, mu, nu),
        ];
        let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(2)));
        let large = Atom::add_many((1..512).flat_map(|i| {
            let term = (&x + i).pow(2);
            [&q * &term, &y * term]
        }));
        assert!(large.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES);
        let sources = [&q * &x + &q * &y, large, r.clone()];
        let expected_bodies = [&q * (&x + &y), (&q + &y) * weight];
        assert_ne!(sources[0], expected_bodies[0]);
        assert_ne!(sources[1], expected_bodies[1]);
        let definitions = sources
            .iter()
            .enumerate()
            .map(|(tag, scalar)| {
                Arc::new(FnMapEntry {
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                    lhs: function!(family, tag, &q, &r),
                    rhs: scalar * &contractions[tag],
                    tags: vec![Atom::num(tag)],
                    args: [&q, &r]
                        .into_iter()
                        .map(|arg| Indeterminate::try_from(arg.clone()).unwrap())
                        .collect(),
                })
            })
            .collect::<Vec<_>>();
        let roots = [
            function!(family, 0, 3, 11),
            function!(family, 0, 5, 13),
            function!(family, 1, 3, 17),
            function!(family, 1, 5, 19),
            function!(family, 2, 3, 23),
        ];
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x, y].into_iter().collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &roots,
            &builder,
            &definitions,
            &settings,
            symbol!("evaluator_test::compact_numerator_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        // Both a contracted component and a retained large scalar alias are
        // compacted before assigning their minimal formal dependencies. The
        // temporal-spatial contraction separately verifies sparse-zero pruning.
        for (expected, source) in expected_bodies.iter().zip(&sources) {
            assert!(prepared.reps.iter().any(|entry| entry.rhs == expected));
            assert!(expected.as_view().get_byte_size() < source.as_view().get_byte_size());
        }
        assert!(
            prepared
                .reps
                .iter()
                .all(|entry| entry.args.iter().all(|arg| *arg == q))
        );
        assert!(atoms[4].get_root().is_zero());
        assert!(
            prepared
                .reps
                .iter()
                .all(|entry| entry.tags[0] != Atom::num(2))
        );
        for ((definition, source), contraction) in definitions.iter().zip(sources).zip(contractions)
        {
            assert_eq!(definition.rhs, source * contraction);
        }
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let values = [0.0, 1.0, 7.0, 0.0].map(|v| Complex::new_re(F(v)));
        let w = (1..512).map(|i| i * i).sum::<i64>() as f64;
        let dw = (1..512).map(|i| 2 * i).sum::<i64>() as f64;
        for compiled in [false, true] {
            if compiled {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O0,
                        ..Default::default()
                    })
                    .unwrap();
            }
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
            assert_eq!(actual.len(), roots.len());
            for (actual, expected) in actual.iter().zip([
                [21.0, 3.0],
                [35.0, 5.0],
                [10.0 * w, 10.0 * dw],
                [12.0 * w, 12.0 * dw],
                [0.0, 0.0],
            ]) {
                let DualOrNot::Dual(actual) = actual else {
                    panic!("expected dual result")
                };
                assert_eq!(
                    actual.values,
                    expected.map(|v| Complex::new_re(F(v))).to_vec()
                );
            }
        }
    }

    #[test]
    fn shared_numerator_components_prune_only_unused_formals_through_aliases() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::pruned_numerator_family");
        // Function-valued formals must be matched as complete atoms, not by head.
        let q1 = function!(symbol!("evaluator_test::pruned_numerator_q"), 0);
        let q2 = function!(symbol!("evaluator_test::pruned_numerator_q"), 1);
        let x = Atom::var(symbol!("evaluator_test::pruned_numerator_x"));
        let vector = SPENSO_TAG.tensor_symbol("evaluator_test::pruned_numerator_vector");
        let index = parse_lit!(spenso::mink(4, 1));
        let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(2)));
        assert!(weight.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES);
        let formals = [q1.clone(), q2.clone()]
            .into_iter()
            .map(|arg| Indeterminate::try_from(arg).unwrap())
            .collect::<Vec<_>>();
        let definitions = [
            Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, 0, q1.clone(), q2.clone()),
                rhs: (&q1 + &x) * &weight * GS.energy_delta(index.as_view()),
                tags: vec![Atom::num(0)],
                args: formals.clone(),
            }),
            Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, 1, q1.clone(), q2.clone()),
                rhs: weight,
                tags: vec![Atom::num(1)],
                args: formals,
            }),
        ];
        let roots = [
            function!(family, 0, 3, 5) * function!(vector, index.clone()),
            function!(family, 0, 3, 99) * function!(vector, index),
            function!(family, 1, 3, 5),
            function!(family, 1, 17, 99),
        ];
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x, function!(vector, function!(AIND_SYMBOLS.cind, 0))]
            .into_iter()
            .collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings {
            store_atom: true,
            ..Default::default()
        };
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &roots,
            &builder,
            &definitions,
            &settings,
            symbol!("evaluator_test::pruned_numerator_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        assert_eq!(atoms[0].get_root(), atoms[1].get_root());
        assert_eq!(atoms[2].get_root(), atoms[3].get_root());
        assert!(prepared.reps.iter().any(|entry| {
            entry
                .lhs
                .as_fun_view()
                .unwrap()
                .get_symbol()
                .get_name()
                .starts_with("gammalooprs::numerator_coefficient_")
                && entry.rhs.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES
        }));
        let surviving_formals = prepared
            .reps
            .iter()
            .flat_map(|entry| entry.args.iter().cloned())
            .collect::<HashSet<_>>();
        assert_eq!(surviving_formals.len(), 1);
        assert!(prepared.reps.iter().all(|entry| entry.args.len() <= 1));
        // The already scalar family keeps its original boundary. Its additive
        // body can yield several components; none may retain either formal.
        let scalar_entries = prepared
            .reps
            .iter()
            .filter(|entry| entry.tags[0] == Atom::num(1))
            .collect::<Vec<_>>();
        assert!(!scalar_entries.is_empty());
        assert!(scalar_entries.iter().all(|entry| entry.args.is_empty()));
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let values = [0.0, 1.0, 7.0, 0.0].map(|v| Complex::new_re(F(v)));
        let w = (1..512).map(|i| i * i).sum::<i64>() as f64;
        let dw = (1..512).map(|i| 2 * i).sum::<i64>() as f64;
        for compile in [false, true] {
            if compile {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O0,
                        ..Default::default()
                    })
                    .unwrap();
                assert_eq!(evaluator.active_f64_backend(), ActiveF64Backend::Symjit);
            }
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
            for (i, result) in actual.iter().enumerate() {
                let DualOrNot::Dual(result) = result else {
                    panic!("expected dual result")
                };
                let expected = if i < 2 {
                    [21.0 * w, 7.0 * (w + 3.0 * dw)]
                } else {
                    [w, dw]
                };
                assert_eq!(
                    result.values,
                    expected.map(|v| Complex::new_re(F(v))).to_vec()
                );
            }
        }

        // Short affine energy sums stay visible for inlined components and
        // share scalar definitions otherwise, including sums under Mul.
        let x = Atom::var(symbol!("evaluator_test::short_affine_x"));
        let y = Atom::var(symbol!("evaluator_test::short_affine_y"));
        let p = Atom::var(symbol!("evaluator_test::short_affine_p"));
        let r = Atom::var(symbol!("evaluator_test::short_affine_r"));
        let gate = Atom::var(symbol!("evaluator_test::short_affine_gate"));
        let q = function!(symbol!("evaluator_test::short_affine_q"), 0);
        let unused = function!(symbol!("evaluator_test::short_affine_q"), 1);
        let index = parse_lit!(spenso::mink(4, 19));
        let temporal = GS.energy_delta(index.as_view());
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params =
            [x.clone(), y.clone(), p.clone(), r.clone(), gate.clone()]
                .into_iter()
                .collect();
        builder.pairs.update_ranges();
        builder
            .add_function(
                GS.tree_denom_wrapper,
                vec![
                    Indeterminate::try_from(Atom::var(symbol!(
                        "evaluator_test::short_affine_tree_arg"
                    )))
                    .unwrap(),
                ],
                Atom::var(symbol!("evaluator_test::short_affine_tree_arg")),
            )
            .unwrap();
        let tree = function!(GS.tree_denom_wrapper, (&x * &x + 4).pow(-1));
        for spectator in [false, true] {
            let affine = &q * &x + &y;
            let weight = if spectator { affine * (&x + 1) } else { affine };
            assert!(weight.as_view().get_byte_size() < NETWORK_SCALAR_ALIAS_MIN_BYTES);
            let definition = Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, 2, q.clone(), unused.clone()),
                rhs: &weight * &temporal,
                tags: vec![Atom::num(2)],
                args: [q.clone(), unused.clone()]
                    .into_iter()
                    .map(|arg| Indeterminate::try_from(arg).unwrap())
                    .collect(),
            });
            let roots = [
                &tree * function!(family, 2, p.clone(), 17) * &temporal,
                function!(Symbol::IF, gate.clone(), x.pow(-1), 0)
                    * function!(family, 2, r.clone(), 23)
                    * &temporal,
            ];
            for inlining in [InliningPolicy::Always, InliningPolicy::Never] {
                let settings = EvaluatorSettings {
                    store_atom: true,
                    ..Default::default()
                };
                let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
                    &roots,
                    &builder,
                    std::slice::from_ref(&definition),
                    &settings,
                    symbol!("evaluator_test::short_affine_scalar"),
                    inlining,
                )
                .unwrap();
                let coefficients = prepared
                    .reps
                    .iter()
                    .filter(|entry| {
                        entry
                            .lhs
                            .as_fun_view()
                            .unwrap()
                            .get_symbol()
                            .get_name()
                            .starts_with("gammalooprs::numerator_coefficient_")
                    })
                    .collect::<Vec<_>>();
                let components = prepared
                    .reps
                    .iter()
                    .filter(|entry| {
                        entry
                            .lhs
                            .as_fun_view()
                            .unwrap()
                            .get_symbol()
                            .get_name()
                            .starts_with("gammalooprs::numerator_component_")
                    })
                    .collect::<Vec<_>>();
                assert_eq!(components.len(), 1);
                let scalar_body =
                    if inlining == InliningPolicy::Always {
                        assert!(coefficients.is_empty());
                        components[0]
                    } else {
                        assert_eq!(coefficients.len(), 1);
                        assert!(components[0].rhs.contains_symbol(
                            coefficients[0].lhs.as_fun_view().unwrap().get_symbol()
                        ));
                        coefficients[0]
                    };
                assert_eq!(scalar_body.args.len(), 1);
                let restored = scalar_body
                    .rhs
                    .replace(Atom::from(scalar_body.args[0].clone()).to_pattern())
                    .with(q.clone());
                assert_eq!(restored, weight.collect_compact_factors());
                assert!(prepared.reps.iter().all(|entry| entry.args.len() <= 1));
                assert_eq!(scalar_body.inlining, inlining);
                assert!(prepared.reps.iter().all(|entry| {
                    !entry
                        .args
                        .contains(&Indeterminate::try_from(unused.clone()).unwrap())
                }));
                assert_eq!(prepared.reps[0], builder.reps[0]);
                assert!(atoms[0].get_root().contains_symbol(GS.tree_denom_wrapper));
                assert_eq!(definition.rhs, &weight * &temporal);
                for dual_shape in [
                    None,
                    Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(2)),
                ] {
                    let width = dual_shape.as_ref().map_or(1, |shape| shape.len());
                    let mut evaluator = GenericEvaluator::new_from_builder(
                        atoms.clone(),
                        &prepared,
                        dual_shape,
                        settings.optimization_settings(),
                        &settings,
                    )
                    .unwrap();
                    if width == 1 && inlining == InliningPolicy::Always {
                        let operations = evaluator.f64_eager.count_operations();
                        assert!(
                            operations.additions <= 12 && operations.multiplications <= 24,
                            "short affine sums lost their compact native boundary: {operations:?}"
                        );
                    }
                    for symjit in [false, true] {
                        if symjit {
                            evaluator
                                .activate_symjit(&CompilationOptionsSnapshot {
                                    optimization_level: CompilationOptimizationLevel::O0,
                                    ..Default::default()
                                })
                                .unwrap();
                        }
                        for (x_value, r_value, active) in [
                            (0.0, 1.0, false),
                            (0.0, 1.0, true),
                            (2.0, 1.0, true),
                            (2.0, 2.0, true),
                        ] {
                            let mut values =
                                [x_value, 3.0, 1.0, r_value, if active { 1.0 } else { 0.0 }]
                                    .into_iter()
                                    .flat_map(|value| {
                                        std::iter::once(Complex::new_re(F(value)))
                                            .chain((1..width).map(|_| Complex::new_re(F(0.0))))
                                    })
                                    .collect::<Vec<_>>();
                            if width > 1 {
                                values[1] = Complex::new_re(F(1.0));
                            }
                            let mut expected = Vec::new();
                            for (row, slope) in [1.0, r_value].into_iter().enumerate() {
                                let base = slope * x_value + 3.0;
                                let (f, df, ddf) = if spectator {
                                    (
                                        base * (x_value + 1.0),
                                        base + slope * (x_value + 1.0),
                                        2.0 * slope,
                                    )
                                } else {
                                    (base, slope, 0.0)
                                };
                                expected.push(if row == 0 {
                                    let d = x_value * x_value + 4.0;
                                    [
                                        f / d,
                                        df / d - 2.0 * x_value * f / (d * d),
                                        ddf / d - 4.0 * x_value * df / (d * d) - 2.0 * f / (d * d)
                                            + 8.0 * x_value * x_value * f / (d * d * d),
                                    ]
                                } else if active {
                                    [
                                        f / x_value,
                                        df / x_value - f / (x_value * x_value),
                                        ddf / x_value - 2.0 * df / (x_value * x_value)
                                            + 2.0 * f / (x_value * x_value * x_value),
                                    ]
                                } else {
                                    [0.0; 3]
                                });
                            }
                            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(
                                &mut evaluator,
                            )(&values);
                            assert_eq!(actual.len(), expected.len());
                            for (row, (actual, expected)) in actual.iter().zip(expected).enumerate()
                            {
                                if row == 1 && active && x_value == 0.0 {
                                    let values = match actual {
                                        DualOrNot::NonDual(value) => std::slice::from_ref(value),
                                        DualOrNot::Dual(value) => value.values.as_slice(),
                                    };
                                    assert!(values.iter().all(|value| {
                                        !value.re.0.is_finite() || !value.im.0.is_finite()
                                    }));
                                    continue;
                                }
                                let expected = expected[..width]
                                    .iter()
                                    .map(|value| Complex::new_re(F(*value)))
                                    .collect::<Vec<_>>();
                                match actual {
                                    DualOrNot::NonDual(value) => assert_eq!(*value, expected[0]),
                                    // Native jets store Taylor coefficients; the oracle uses derivatives.
                                    DualOrNot::Dual(value) => assert_eq!(
                                        crate::utils::hyperdual_utils::extract_t_derivatives_complex(
                                            value.clone(),
                                        ),
                                        expected
                                    ),
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn nested_retained_numerator_coefficients_compile_with_distinct_bodies() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::compiled_nested_numerator_family");
        let q = Atom::var(symbol!("evaluator_test::compiled_nested_numerator_q"));
        let x = Atom::var(symbol!("evaluator_test::compiled_nested_numerator_x"));
        let vector = SPENSO_TAG.tensor_symbol("evaluator_test::compiled_nested_numerator_vector");
        let index = parse_lit!(spenso::mink(4, 1));
        let definitions = [2, 3]
            .into_iter()
            .enumerate()
            .map(|(tag, power)| {
                let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(power)));
                assert!(weight.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES);
                Arc::new(FnMapEntry {
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                    lhs: function!(family, tag, q.clone()),
                    rhs: &q * weight * GS.energy_delta(index.as_view()),
                    tags: vec![Atom::num(tag)],
                    args: vec![Indeterminate::try_from(q.clone()).unwrap()],
                })
            })
            .collect::<Vec<_>>();
        let roots = [
            function!(family, 0, 3) * function!(vector, index.clone()),
            function!(family, 1, 5) * function!(vector, index),
        ];
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x, function!(vector, function!(AIND_SYMBOLS.cind, 0))]
            .into_iter()
            .collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &roots,
            &builder,
            &definitions,
            &settings,
            symbol!("evaluator_test::compiled_nested_numerator_scalar"),
            InliningPolicy::Never,
        )
        .unwrap();
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            None,
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let program = evaluator.f64_eager.export_instructions();
        let nested_definitions = program
            .sub_evaluators
            .iter()
            .flat_map(|sub| &sub.instructions.sub_evaluators)
            .map(|sub| (sub.symbol, sub.tags.clone()))
            .collect::<HashSet<_>>();
        assert!(
            nested_definitions.len() >= 2,
            "distinct large scalar coefficients must remain separate nested programs"
        );
        let values = [2.0, 7.0].map(|value| Complex::new_re(F(value)));
        let expected = [
            21 * (1..512).map(|i: i64| (i + 2).pow(2)).sum::<i64>(),
            35 * (1..512).map(|i: i64| (i + 2).pow(3)).sum::<i64>(),
        ]
        .map(|value| Complex::new_re(F(value as f64)));
        let directory = std::env::temp_dir().join(format!(
            "gammaloop-nested-numerator-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let options = crate::settings::global::CompilationOptionsSnapshot {
            optimization_level: CompilationOptimizationLevel::O0,
            compiler: "g++".into(),
            ..Default::default()
        };
        for mode in [
            None,
            Some(FrozenCompilationMode::Cpp(options.clone())),
            Some(FrozenCompilationMode::Assembly(options)),
        ] {
            let name = mode
                .as_ref()
                .map_or("eager", |mode| mode.active_backend_name());
            if let Some(mode) = &mode {
                evaluator
                    .compile_external(
                        directory.join(format!("{name}.cpp")),
                        "nested_numerator",
                        directory.join(format!("{name}.so")),
                        mode,
                    )
                    .unwrap();
            }
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
            assert_eq!(actual.len(), expected.len());
            for (actual, expected) in actual.iter().zip(&expected) {
                let DualOrNot::NonDual(actual) = actual else {
                    panic!("expected scalar numerator components")
                };
                assert_eq!(actual, expected, "{name} nested numerator coefficient");
            }
        }
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn shared_numerator_pruning_preserves_formals_captured_by_parameter_aliases() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::captured_numerator_family");
        let alias = function!(symbol!("evaluator_test::captured_numerator_alias"), 0);
        let q = Atom::var(symbol!("evaluator_test::captured_numerator_q"));
        let x = Atom::var(symbol!("evaluator_test::captured_numerator_x"));
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        builder.pairs.update_ranges();
        builder
            .fn_map
            .add_aliases([(alias.clone(), &q * &x)])
            .unwrap();
        builder.reps.push(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: alias.clone(),
            rhs: &q * x,
            tags: vec![Atom::num(0)],
            args: vec![],
        });
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 0, q.clone()),
            rhs: alias,
            tags: vec![Atom::num(0)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[function!(family, 0, 3) + function!(family, 0, 5)],
            &builder,
            &[definition],
            &settings,
            symbol!("evaluator_test::captured_numerator_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        assert_eq!(prepared.reps.last().unwrap().args.len(), 1);
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
            Complex::new_re(F(2.0)),
            Complex::new_re(F(1.0)),
        ]);
        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("expected dual result")
        };
        assert_eq!(
            actual.values,
            vec![Complex::new_re(F(16.0)), Complex::new_re(F(8.0))]
        );
    }

    #[test]
    fn shared_numerator_products_close_before_componentization() {
        use spenso::network::library::symbolic::ETS;

        test_initialise().unwrap();
        let family = symbol!("evaluator_test::joint_numerator_family");
        let x = Atom::var(symbol!("evaluator_test::joint_numerator_x"));
        let q = Atom::var(symbol!("evaluator_test::joint_numerator_q"));
        let r = Atom::var(symbol!("evaluator_test::joint_numerator_r"));
        let mu = parse_lit!(spenso::mink(4, 1));
        let nu = parse_lit!(spenso::mink(4, 2));
        let projector = function!(ETS.metric, &mu, &nu);
        let definitions = [
            (&q, (&q + &x) * (&x + 1) * GS.energy_delta(mu.as_view())),
            (&r, (&r - 2 * &x) * GS.energy_delta(nu.as_view())),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (formal, rhs))| {
            Arc::new(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs: function!(family, index, formal.clone()),
                rhs,
                tags: vec![Atom::num(index)],
                args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
            })
        })
        .collect::<Vec<_>>();
        let mut roots = vec![
            (&x + 7) * function!(family, 0, &x + 2) * function!(family, 1, 2 * &x + 3) * &projector,
            (&x + 8) * function!(family, 0, 3 * &x + 4) * function!(family, 1, &x + 5) * projector,
        ];
        // This sum is already scalar, but its interior still needs joint
        // preparation before the outer scalar product makes it opaque.
        roots.push((&x + 9) * (&roots[0] + &roots[1]));
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x].into_iter().collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &roots,
            &builder,
            &definitions,
            &settings,
            symbol!("evaluator_test::joint_numerator_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        // Both argument rows reuse one closed scalar body, rather than the
        // separately materialized open components of either numerator.
        let components = prepared
            .reps
            .iter()
            .filter(|entry| {
                entry
                    .lhs
                    .as_fun_view()
                    .unwrap()
                    .get_symbol()
                    .get_name()
                    .starts_with("gammalooprs::numerator_component_")
            })
            .collect::<Vec<_>>();
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].args.len(), 2);
        assert_eq!(
            components[0].tags.last(),
            Some(&function!(AIND_SYMBOLS.cind))
        );
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
            Complex::new_re(F(1.0)),
            Complex::new_re(F(1.0)),
        ]);
        for (actual, expected) in
            actual
                .iter()
                .zip([[192.0, 216.0], [576.0, 496.0], [7680.0, 7888.0]])
        {
            let DualOrNot::Dual(actual) = actual else {
                panic!("expected dual result")
            };
            assert_eq!(
                actual.values,
                expected.map(|v| Complex::new_re(F(v))).to_vec()
            );
        }
    }

    #[test]
    fn shared_numerator_residue_sum_closes_through_color_trace_and_outer_gamma() {
        use idenso::{color::CS, dirac::gamma_tensor, shorthands::chain::Chain};

        test_initialise().unwrap();
        let family = symbol!("evaluator_test::trace_joint_family");
        let q = Atom::var(symbol!("evaluator_test::trace_joint_q"));
        let x = Atom::var(symbol!("evaluator_test::trace_joint_x"));
        let mu = parse_lit!(spenso::mink(4, 1));
        let left = parse_lit!(spenso::bis(4, 3));
        let right = parse_lit!(spenso::bis(4, 4));
        let generator = CS.chain_t(parse_lit!(spenso::coad(8, 7)));
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 0, &q),
            rhs: (&q + &x) * gamma_tensor(left.clone(), right.clone(), mu.clone()),
            tags: vec![Atom::num(0)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let traced = |argument: Atom| {
            function!(
                SPENSO_TAG.trace,
                parse_lit!(spenso::cof(3)),
                spenso::shadowing::cyclic([
                    generator.clone(),
                    generator.clone() * function!(family, 0, argument),
                ])
            )
        };
        let root =
            ((&x + 3) * traced(&x + 2) + 2 * traced(3 * &x + 4)) * gamma_tensor(right, left, mu);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x].into_iter().collect();
        builder.pairs.update_ranges();
        // Existing compact bodies already own in/out. They must keep their
        // ordinary interface rather than become a nested color/spin binder.
        let compact = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            rhs: definition.rhs.chainify(Bispinor {}.into()),
            ..(*definition).clone()
        });
        let (unfused, joints) = EvaluatorStack::combine_numerator_families(
            std::slice::from_ref(&root),
            &[compact],
            &builder,
            &HashMap::from([((family, vec![Atom::num(0)]), 0)]),
            &HashMap::from([(family, 1)]),
            (
                symbol!("evaluator_test::trace_binder_joint"),
                symbol!("evaluator_test::trace_binder_argument"),
                SPENSO_TAG.tensor_symbol("evaluator_test::trace_binder_shadow"),
            ),
        );
        assert!(joints.is_empty());
        assert!(unfused[0].contains_symbol(family));
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[root],
            &builder,
            &[definition],
            &settings,
            symbol!("evaluator_test::trace_joint_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        let components = prepared
            .reps
            .iter()
            .filter(|entry| {
                entry
                    .lhs
                    .as_fun_view()
                    .unwrap()
                    .get_symbol()
                    .get_name()
                    .starts_with("gammalooprs::numerator_component_")
            })
            .collect::<Vec<_>>();
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].args.len(), 1);
        assert_eq!(
            components[0].tags.last(),
            Some(&function!(AIND_SYMBOLS.cind))
        );
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
            Complex::new_re(F(1.0)),
            Complex::new_re(F(1.0)),
        ]);
        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("expected dual result")
        };
        // Tr_color(T^a T^a)=4 and Tr_spin(gamma_mu gamma^mu)=16.
        // The weighted affine rows give 32 at x=1, with derivative 20.
        assert_eq!(
            actual.values,
            vec![Complex::new_re(F(2048.0)), Complex::new_re(F(1280.0))]
        );
    }

    #[test]
    fn shared_numerator_core_collection_preserves_tensor_sums_and_scalar_boundaries() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::core_boundary_family");
        let joint = symbol!("evaluator_test::core_boundary_joint");
        let argument = symbol!("evaluator_test::core_boundary_argument");
        let shadow = SPENSO_TAG.tensor_symbol("evaluator_test::core_boundary_shadow");
        let q = Atom::var(symbol!("evaluator_test::core_boundary_q"));
        let x = Atom::var(symbol!("evaluator_test::core_boundary_x"));
        let vector = SPENSO_TAG.tensor_symbol("evaluator_test::core_boundary_vector");
        let mu = parse_lit!(spenso::mink(4, 1));
        let temporal = GS.energy_delta(mu.as_view());
        let companion = &temporal + &x * function!(vector, mu);
        let definitions = [&q * &temporal, &q + &x]
            .into_iter()
            .enumerate()
            .map(|(index, rhs)| {
                Arc::new(FnMapEntry {
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                    lhs: function!(family, index, &q),
                    rhs,
                    tags: vec![Atom::num(index)],
                    args: vec![Indeterminate::try_from(q.clone()).unwrap()],
                })
            })
            .collect::<Vec<_>>();
        let calls = HashMap::from([
            ((family, vec![Atom::num(0)]), 0),
            ((family, vec![Atom::num(1)]), 1),
        ]);
        let tag_counts = HashMap::from([(family, 1)]);
        let scalar_product = (function!(family, 1, 1) + function!(family, 1, 2))
            * (function!(family, 1, 3) + function!(family, 1, 4));
        let guarded = Symbol::IF.call_args([
            x.clone(),
            function!(family, 0, 1) * &temporal / x.clone(),
            Atom::Zero,
        ]);
        let roots = [
            (function!(family, 0, 1) + function!(family, 0, 2)) * &companion,
            scalar_product.clone(),
            guarded.clone(),
        ];
        let (prepared, joints) = EvaluatorStack::combine_numerator_families(
            &roots,
            &definitions,
            &ParamBuilder::new_empty(),
            &calls,
            &tag_counts,
            (joint, argument, shadow),
        );
        assert_eq!(joints.len(), 1);
        assert!(joints[0].rhs.contains(&companion));
        assert_eq!(prepared[1], scalar_product);
        assert_eq!(prepared[2], guarded);
        assert!(!joints[0].rhs.contains_symbol(Symbol::IF));
        assert!(!joints[0].rhs.contains(x.pow(-1)));
    }

    #[test]
    fn shared_numerator_product_preserves_private_contractions_per_call() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::repeated_joint_family");
        let q = Atom::var(symbol!("evaluator_test::repeated_joint_q"));
        let x = Atom::var(symbol!("evaluator_test::repeated_joint_x"));
        let p = SPENSO_TAG.tensor_symbol("evaluator_test::repeated_joint_p");
        let r = SPENSO_TAG.tensor_symbol("evaluator_test::repeated_joint_r");
        let internal = parse_lit!(spenso::mink(4, 17));
        let external = parse_lit!(spenso::mink(4, 23));
        let other_internal = parse_lit!(spenso::mink(4, 19));
        let trace = function!(
            SPENSO_TAG.trace,
            parse_lit!(spenso::bis(4)),
            spenso::shadowing::sym([&internal, &other_internal].map(|slot| {
                function!(
                    idenso::dirac::AGS.gamma,
                    SPENSO_TAG.chain_in,
                    SPENSO_TAG.chain_out,
                    slot
                )
            }))
        );
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 0, q.clone()),
            rhs: (&q + &x)
                * function!(p, &internal)
                * function!(r, other_internal)
                * trace
                * GS.energy_delta(external.as_view()),
            tags: vec![Atom::num(0)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let root = function!(family, 0, 2 * &x + 3) * function!(family, 0, &x + 2);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = std::iter::once(x.clone())
            .chain((0..4).flat_map(|i| {
                [
                    function!(p, function!(AIND_SYMBOLS.cind, i)),
                    function!(r, function!(AIND_SYMBOLS.cind, i)),
                ]
            }))
            .collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings {
            summed_function_map: true,
            summed: true,
            iterative_orientation_optimization: true,
            ..Default::default()
        };
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[root],
            &builder,
            &[Arc::clone(&definition)],
            &settings,
            symbol!("evaluator_test::repeated_joint_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        let components = prepared
            .reps
            .iter()
            .filter(|entry| {
                entry
                    .lhs
                    .as_fun_view()
                    .unwrap()
                    .get_symbol()
                    .get_name()
                    .starts_with("gammalooprs::numerator_component_")
            })
            .collect::<Vec<_>>();
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].args.len(), 2);
        assert_eq!(
            components[0].tags.last(),
            Some(&function!(AIND_SYMBOLS.cind))
        );
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let values = std::iter::once([1.0, 1.0])
            .chain([1.0, 5.0, 2.0, 6.0, 3.0, 7.0, 4.0, 8.0].map(|value| [value, 0.0]))
            .flatten()
            .map(|v| Complex::new_re(F(v)))
            .collect::<Vec<_>>();
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("expected dual result")
        };
        // Tr(sym(gamma_mu, gamma_nu)) = 4 g_mu_nu and p.r = -60.
        // Each body owns its contraction, including slots inside sym.
        // Sharing those indices incorrectly replaces (p.r)^2 by p^2 r^2.
        assert_eq!(actual.values, vec![Complex::new_re(F(1382400.0)); 2]);

        // The typed Cartesian-row product must preserve the same private
        // contractions before its single generic body reaches preprocessing.
        let cut = crate::cff::CutCFFIndex::new_all_none();
        let slope = Atom::from(definition.args[0].clone());
        let shift = Atom::var(symbol!("evaluator_test::repeated_joint_shift"));
        let typed_args = vec![
            definition.args[0].clone(),
            Indeterminate::try_from(shift.clone()).unwrap(),
        ];
        let typed_definition = Arc::new(FnMapEntry {
            lhs: symbol!("gammalooprs::uv::numerator_family").call_args(
                definition
                    .tags
                    .iter()
                    .cloned()
                    .chain(typed_args.iter().cloned().map(Atom::from)),
            ),
            rhs: definition
                .rhs
                .replace(slope.to_pattern())
                .with(&slope * &x + shift),
            args: typed_args,
            ..definition.as_ref().clone()
        });
        let factor = |slope: i64, shift: i64| {
            crate::uv::Integrands::from_iter([(cut, Atom::Zero)])
                .with_parametric_terms([(
                    cut,
                    ParametricIntegrandTerm::new(
                        Arc::clone(&typed_definition),
                        vec![crate::uv::ParametricIntegrandRow {
                            coefficients: vec![Rational::from(slope), Rational::from(shift)],
                            carrier: Atom::one(),
                            source_keys: Vec::new(),
                        }],
                    )
                    .unwrap(),
                )])
                .unwrap()
        };
        let product = factor(2, 3).zip_mul(&factor(1, 2)).unwrap();
        let (mut stack, _) = EvaluatorStack::from_integrand_with_timings(
            product.atom(&cut).unwrap(),
            &builder,
            product.numerators(),
            product.parametric_terms(&cut),
            None,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            &settings,
        )
        .unwrap();
        let orientations =
            TiVec::<OrientationID, _>::from_iter([EdgeVec::from_iter(std::iter::empty::<
                Orientation,
            >())]);
        let filter = SubSet::full(1);
        for symjit in [false, true] {
            if symjit {
                stack
                    .for_each_generic_evaluator_mut(|evaluator| {
                        evaluator.activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: CompilationOptimizationLevel::O0,
                            ..Default::default()
                        })
                    })
                    .unwrap();
            }
            for method in [
                EvaluatorMethod::SingleParametric,
                EvaluatorMethod::SummedFunctionMap,
                EvaluatorMethod::Summed,
                EvaluatorMethod::Iterative,
            ] {
                let mut runtime = RuntimeSettings::default();
                runtime.general.evaluator_method = method.clone();
                let output = stack
                    .evaluate(
                        InputParams {
                            values: SliceMut::Owned(values.clone()),
                            residue_map_id_start: 0,
                            orientations_start: 0,
                            multiplicative_offset: 2,
                        },
                        SingleOrAllOrientations::All {
                            all: &orientations,
                            filter: &filter,
                        },
                        &runtime,
                        &mut EvaluationMetaData::new_empty(),
                    )
                    .unwrap();
                let [DualOrNot::Dual(actual)] = output.as_slice() else {
                    panic!("expected typed product dual output")
                };
                assert_eq!(
                    actual.values,
                    vec![Complex::new_re(F(1382400.0)); 2],
                    "{method:?}, SymJIT={symjit}",
                );
            }
        }
    }

    #[test]
    fn shared_numerator_products_keep_alias_captured_formals_bound() {
        test_initialise().unwrap();
        let family = symbol!("evaluator_test::joint_captured_family");
        let alias = function!(symbol!("evaluator_test::joint_captured_alias"), 0);
        let q = Atom::var(symbol!("evaluator_test::joint_captured_q"));
        let x = Atom::var(symbol!("evaluator_test::joint_captured_x"));
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        builder.pairs.update_ranges();
        builder
            .fn_map
            .add_aliases([(alias.clone(), &q * &x)])
            .unwrap();
        builder.reps.push(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: alias.clone(),
            rhs: &q * x,
            tags: vec![Atom::num(0)],
            args: vec![],
        });
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 0, q.clone()),
            rhs: alias,
            tags: vec![Atom::num(0)],
            args: vec![Indeterminate::try_from(q.clone()).unwrap()],
        });
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[function!(family, 0, 3) * function!(family, 0, 5)],
            &builder,
            &[definition],
            &settings,
            symbol!("evaluator_test::joint_captured_scalar"),
            InliningPolicy::Always,
        )
        .unwrap();
        assert_eq!(
            prepared.reps.last().unwrap().args,
            vec![Indeterminate::try_from(q).unwrap()]
        );
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
            Complex::new_re(F(2.0)),
            Complex::new_re(F(1.0)),
        ]);
        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("expected dual result")
        };
        assert_eq!(actual.values, vec![Complex::new_re(F(60.0)); 2]);
    }

    #[test]
    fn scalar_products_of_independent_sums_keep_linear_instruction_growth() {
        test_initialise().unwrap();
        for terms in [2, 4, 8] {
            let parameters = (0..3 * terms)
                .map(|i| Atom::var(symbol!(format!("evaluator_test::independent_sum_{i}"))))
                .collect::<Vec<_>>();
            let sums = parameters
                .chunks(terms)
                .map(|chunk| Atom::add_many(chunk.iter().cloned()))
                .collect::<Vec<_>>();
            let product = Atom::mul_many(sums.iter().cloned());
            assert!(matches!(product.as_view(), AtomView::Mul(mul) if mul.get_nargs() == 3));
            for alias_kind in 0..3 {
                let mut source = AliasedAtom::from(product.clone());
                let mut function_map = FunctionMap::default();
                if alias_kind != 0 {
                    let aliases = (0..3)
                        .map(|i| function!(symbol!("evaluator_test::independent_sum_alias"), i))
                        .collect::<Vec<_>>();
                    source = Atom::mul_many(aliases.iter().cloned()).into();
                    for (alias, sum) in aliases.into_iter().zip(&sums) {
                        if alias_kind == 1 {
                            source.register_alias(alias, sum.clone());
                        } else {
                            function_map.add_aliases([(alias, sum.clone())]).unwrap();
                        }
                    }
                }
                for direct_translation in [false, true] {
                    for horner_iterations in [0, 1] {
                        let settings = EvaluatorSettings {
                            direct_translation,
                            horner_iterations,
                            ..Default::default()
                        };
                        let mut evaluator = GenericEvaluator::new_from_raw_params(
                            [source.clone()],
                            &parameters,
                            &function_map,
                            vec![],
                            settings.optimization_settings(),
                            None,
                            &settings,
                        )
                        .unwrap();
                        let operations = evaluator.f64_eager.count_operations();
                        assert_eq!(operations.additions, 3 * (terms - 1));
                        assert_eq!(operations.multiplications, 2);
                        let values = vec![Complex::new_re(F(1.0)); parameters.len()];
                        assert_eq!(
                            <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(
                                &values
                            ),
                            Complex::new_re(F((terms * terms * terms) as f64))
                        );
                    }
                }
            }
        }
    }

    #[test]
    #[ignore = "bounded manual preparation comparison; run explicitly with --ignored --nocapture"]
    fn shared_numerator_preparation_comparison() {
        use spenso::network::library::symbolic::ETS;

        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::preparation_comparison_x"));
        let q1 = Atom::var(symbol!("evaluator_test::preparation_comparison_q1"));
        let q2 = Atom::var(symbol!("evaluator_test::preparation_comparison_q2"));
        let family = symbol!("evaluator_test::preparation_comparison_family");
        let weight_alias = function!(symbol!("evaluator_test::preparation_comparison_weight"), 0);
        let weight = Atom::add_many((1..17).map(|i| (&x + i).pow(2)));
        assert!(weight.as_view().get_byte_size() < NETWORK_SCALAR_ALIAS_MIN_BYTES);
        let first_index = parse_lit!(spenso::mink(4, 1));
        let second_index = parse_lit!(spenso::mink(4, 2));
        let metric = function!(ETS.metric, first_index.clone(), second_index.clone());
        let temporal =
            GS.energy_delta(first_index.as_view()) * GS.energy_delta(second_index.as_view());
        let first: Atom = &weight * ((&q1 + &x) * &metric + &q1 * &x * &temporal);
        let second: Atom = &weight * ((&q2 * &x + 1) * &metric + (&q2 + 2 * &x) * &temporal);
        // Both coefficient arguments depend on the same joint (a,b) row. Their
        // residue sums therefore cannot be separated into independent factors.
        let rows = [-1i64, 1]
            .into_iter()
            .flat_map(|a| {
                [-1i64, 0, 1]
                    .into_iter()
                    .map(move |b| (a + 2 * b, 2 * a - b))
            })
            .collect::<Vec<_>>();
        let orientations = TiVec::<OrientationID, _>::from_iter(
            rows.iter()
                .map(|_| EdgeVec::from_iter([Orientation::Default])),
        );
        let production_ids = (0..rows.len())
            .map(|i| OrientationID(10 + 3 * i))
            .collect::<Vec<_>>();
        let mut base_builder = ParamBuilder::new_empty();
        base_builder.pairs.residue_map_id = ParamValuePairs::default_from_symbol(GS.residue_map_id);
        base_builder.pairs.orientations = [GS.sign(EdgeIndex(0))].into_iter().collect();
        base_builder.pairs.additional_params = [x.clone()].into_iter().collect();
        let parameter_count = base_builder.pairs.update_ranges();
        base_builder.values = vec![vec![Complex::new_re(F(0.0)); parameter_count]];
        let weight0 = (1..17).map(|i| i * i).sum::<i64>();
        let weight_derivative0 = (1..17).map(|i| 2 * i).sum::<i64>();
        let expected_sum = rows
            .iter()
            .enumerate()
            .map(|(i, (a, b))| (weight0 * weight0 * a * (4 + b)) as f64 / (10 + i) as f64)
            .sum::<f64>();
        let settings = EvaluatorSettings {
            iterative_orientation_optimization: true,
            summed_function_map: true,
            summed: true,
            store_atom: true,
            do_fn_map_replacements: false,
            ..Default::default()
        };
        // Round zero is warm-up only. Rotate the order of the three variants, and
        // report every raw round without asserting a wall-time threshold.
        for round in 0..4 {
            for offset in 0..3 {
                let variant = (round + offset) % 3;
                let variant_name = [
                    "combined_closed_tensor",
                    "separate_open_tensors",
                    "combined_manual_scalar_alias",
                ][variant];
                let mut builder = base_builder.clone();
                let mut definitions = if variant == 1 {
                    vec![
                        Arc::new(FnMapEntry {
                            inlining: InliningPolicy::Always,
                            is_alias: false,
                            lhs: function!(family, 1, q1.clone()),
                            rhs: first.clone(),
                            tags: vec![Atom::num(1)],
                            args: vec![Indeterminate::try_from(q1.clone()).unwrap()],
                        }),
                        Arc::new(FnMapEntry {
                            inlining: InliningPolicy::Always,
                            is_alias: false,
                            lhs: function!(family, 2, q2.clone()),
                            rhs: second.clone(),
                            tags: vec![Atom::num(2)],
                            args: vec![Indeterminate::try_from(q2.clone()).unwrap()],
                        }),
                    ]
                } else {
                    vec![Arc::new(FnMapEntry {
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                        lhs: function!(family, 0, q1.clone(), q2.clone()),
                        rhs: &first * &second,
                        tags: vec![Atom::num(0)],
                        args: [q1.clone(), q2.clone()]
                            .into_iter()
                            .map(|arg| Indeterminate::try_from(arg).unwrap())
                            .collect(),
                    })]
                };
                if variant == 2 {
                    builder
                        .fn_map
                        .add_aliases([(weight_alias.clone(), weight.clone())])
                        .unwrap();
                    builder.reps.push(FnMapEntry {
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                        lhs: weight_alias.clone(),
                        rhs: weight.clone(),
                        tags: vec![Atom::num(0)],
                        args: vec![],
                    });
                    let entry = Arc::make_mut(&mut definitions[0]);
                    entry.rhs = entry
                        .rhs
                        .replace(weight.to_pattern())
                        .with(weight_alias.clone());
                }
                let root = Atom::add_many(rows.iter().zip(&production_ids).enumerate().map(
                    |(i, ((a, b), id))| {
                        let numerator = if variant == 1 {
                            function!(family, 1, *a) * function!(family, 2, *b)
                        } else {
                            function!(family, 0, *a, *b)
                        };
                        id.atom() * numerator / (&x + (10 + i))
                    },
                ));
                let input_bytes = root.as_view().get_byte_size();
                let definition_bytes = definitions
                    .iter()
                    .map(|entry| entry.rhs.as_view().get_byte_size())
                    .sum::<usize>()
                    + builder
                        .reps
                        .iter()
                        .map(|entry| entry.rhs.as_view().get_byte_size())
                        .sum::<usize>();
                let (mut stack, timings) = EvaluatorStack::new_with_timings(
                    std::slice::from_ref(&root),
                    &builder,
                    &definitions,
                    &orientations.raw,
                    &production_ids,
                    None,
                    &settings,
                )
                .unwrap();
                let make_input = || InputParams {
                    values: SliceMut::Owned(vec![Complex::new_re(F(0.0)); parameter_count]),
                    residue_map_id_start: builder.pairs.residue_map_id.value_range.start,
                    orientations_start: builder.pairs.orientations.value_range.start,
                    multiplicative_offset: 1,
                };
                let mut metadata = EvaluationMetaData::new_empty();
                let filter = SubSet::full(orientations.len());
                let all = SingleOrAllOrientations::All {
                    all: &orientations,
                    filter: &filter,
                };
                for method in [
                    EvaluatorMethod::SingleParametric,
                    EvaluatorMethod::Iterative,
                    EvaluatorMethod::SummedFunctionMap,
                    EvaluatorMethod::Summed,
                ] {
                    let mut runtime = RuntimeSettings::default();
                    runtime.general.evaluator_method = method.clone();
                    let actual = scalar_value(
                        stack
                            .evaluate(make_input(), all, &runtime, &mut metadata)
                            .unwrap(),
                    );
                    assert!(
                        (actual.re.0 - expected_sum).abs() <= 1e-11 * expected_sum.abs().max(1.0)
                    );
                }
                // Exact rational values for every selected row, independently of
                // source specialization and floating point operation order.
                // This real fixture maps exactly to Q; nonreal coefficients are
                // rejected before evaluating its selector branches in that ring.
                let mut exact_evaluator = stack
                    .single_parametric
                    .rational
                    .as_ref()
                    .unwrap()
                    .clone()
                    .map_to_ring(&Q)
                    .unwrap();
                for (i, ((a, b), id)) in rows.iter().zip(&production_ids).enumerate() {
                    let mut values = vec![
                        SymComplex::new(Rational::from(0), Rational::from(0));
                        parameter_count
                    ];
                    values[builder.pairs.residue_map_id.value_range.start].re =
                        Rational::from(id.0 as i64);
                    values[builder.pairs.orientations.value_range.start].re = Rational::from(1);
                    let mut result = [SymComplex::new(Rational::from(0), Rational::from(0))];
                    result[0].re = exact_evaluator.evaluate_single_in_ring(
                        &values
                            .iter()
                            .map(|value| value.re.clone())
                            .collect::<Vec<_>>(),
                        &Q,
                    );
                    assert_eq!(
                        result[0],
                        SymComplex::new(
                            Rational::new(weight0 * weight0 * a * (4 + b), (10 + i) as i64),
                            Rational::from(0),
                        )
                    );
                }
                let mut values = vec![Complex::new_re(F(0.0)); parameter_count];
                values[builder.pairs.residue_map_id.value_range.start] =
                    Complex::new_re(F(production_ids[0].0 as f64));
                values[builder.pairs.orientations.value_range.start] = Complex::new_re(F(1.0));
                for _ in 0..16 {
                    std::hint::black_box(<f64 as GenericEvaluatorFloat>::get_evaluator_single(
                        &mut stack.single_parametric,
                    )(&values));
                }
                let runtime_started = std::time::Instant::now();
                for _ in 0..512 {
                    std::hint::black_box(<f64 as GenericEvaluatorFloat>::get_evaluator_single(
                        &mut stack.single_parametric,
                    )(&values));
                }
                let eager_ns_per_call = runtime_started.elapsed().as_nanos() as f64 / 512.0;
                for (mode, evaluator) in [
                    ("single_parametric", &stack.single_parametric),
                    ("iterative", &stack.iterative.as_ref().unwrap().0),
                    (
                        "summed_function_map",
                        stack.summed_function_map.as_ref().unwrap(),
                    ),
                    ("summed", stack.summed.as_ref().unwrap()),
                ] {
                    let operations = evaluator.f64_eager.count_operations();
                    let program = evaluator.f64_eager.export_instructions();
                    crate::debug_tags!(#generation, #profile, #summary;
                        stage = "shared_numerator_preparation_mode_structure", variant = variant_name,
                        round, warmup = round == 0, mode,
                        stored_root_bytes = evaluator.exprs.as_ref().unwrap().iter()
                            .map(|atom| atom.as_view().get_byte_size()).sum::<usize>(),
                        retained_definitions = evaluator.fn_map_entries.len(),
                        retained_argument_slots = evaluator.fn_map_entries.iter()
                            .map(|entry| entry.args.len()).sum::<usize>(),
                        zero_argument_definitions = evaluator.fn_map_entries.iter()
                            .filter(|entry| entry.args.is_empty()).count(),
                        instructions = program.instructions.len(), temporaries = program.temporary_count,
                        additions = operations.additions, multiplications = operations.multiplications,
                        "Shared numerator preparation mode structure"
                    );
                }
                let operations = stack.single_parametric.f64_eager.count_operations();
                let instructions = stack.single_parametric.f64_eager.export_instructions();
                crate::debug_tags!(#generation, #profile, #summary;
                    stage = "shared_numerator_preparation_comparison", variant = variant_name, round,
                    warmup = round == 0, input_bytes, definition_bytes, family_count = definitions.len(),
                    retained_definitions = stack.single_parametric.fn_map_entries.len(),
                    preprocessing_ms = timings.spenso_time.as_secs_f64() * 1000.0,
                    symbolica_ms = timings.symbolica_time.as_secs_f64() * 1000.0,
                    instructions = instructions.instructions.len(), temporaries = instructions.temporary_count,
                    additions = operations.additions, multiplications = operations.multiplications,
                    eager_ns_per_call, "Shared numerator preparation comparison"
                );
                if round == 3 {
                    // Archive the complete stack so all four modes retain the same
                    // small component calls and the same common definitions.
                    let encoded =
                        bincode::encode_to_vec(&stack, bincode::config::standard()).unwrap();
                    let mut state = Vec::new();
                    State::export(&mut state).unwrap();
                    let state_map = State::import(&mut Cursor::new(state), None).unwrap();
                    let model = Model::default();
                    let (mut decoded, _): (EvaluatorStack, _) =
                        bincode::decode_from_slice_with_context(
                            &encoded,
                            bincode::config::standard(),
                            GammaLoopContextContainer {
                                state_map: &state_map,
                                model: &model,
                            },
                        )
                        .unwrap();
                    for method in [
                        EvaluatorMethod::SingleParametric,
                        EvaluatorMethod::Iterative,
                        EvaluatorMethod::SummedFunctionMap,
                        EvaluatorMethod::Summed,
                    ] {
                        let mut runtime = RuntimeSettings::default();
                        runtime.general.evaluator_method = method.clone();
                        let actual = scalar_value(
                            decoded
                                .evaluate(make_input(), all, &runtime, &mut metadata)
                                .unwrap(),
                        );
                        assert!(
                            (actual.re.0 - expected_sum).abs()
                                <= 1e-11 * expected_sum.abs().max(1.0)
                        );
                    }
                    let expected = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                        &mut stack.single_parametric,
                    )(&values);
                    let compile_started = std::time::Instant::now();
                    stack
                        .single_parametric
                        .activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: CompilationOptimizationLevel::O0,
                            ..Default::default()
                        })
                        .unwrap();
                    let compile_ms = compile_started.elapsed().as_secs_f64() * 1000.0;
                    assert_eq!(
                        stack.single_parametric.active_f64_backend(),
                        ActiveF64Backend::Symjit
                    );
                    let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                        &mut stack.single_parametric,
                    )(&values);
                    assert!(
                        (actual.re.0 - expected.re.0).abs() <= 1e-11 * expected.re.0.abs().max(1.0)
                    );
                    for _ in 0..16 {
                        std::hint::black_box(<f64 as GenericEvaluatorFloat>::get_evaluator_single(
                            &mut stack.single_parametric,
                        )(&values));
                    }
                    let compiled_started = std::time::Instant::now();
                    for _ in 0..512 {
                        std::hint::black_box(<f64 as GenericEvaluatorFloat>::get_evaluator_single(
                            &mut stack.single_parametric,
                        )(&values));
                    }
                    crate::debug_tags!(#generation, #profile, #summary;
                        stage = "shared_numerator_preparation_compiled_comparison", variant = variant_name,
                        compile_ms, compiled_ns_per_call = compiled_started.elapsed().as_nanos() as f64 / 512.0,
                        "Shared numerator compiled sanity and runtime"
                    );
                    let dual_settings = EvaluatorSettings {
                        iterative_orientation_optimization: false,
                        summed_function_map: false,
                        summed: false,
                        ..settings
                    };
                    let (mut dual, _) = EvaluatorStack::new_with_timings(
                        &[root],
                        &builder,
                        &definitions,
                        &orientations.raw,
                        &production_ids,
                        Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
                        &dual_settings,
                    )
                    .unwrap();
                    for (i, ((a, b), id)) in rows.iter().zip(&production_ids).enumerate() {
                        let mut values = vec![Complex::new_re(F(0.0)); parameter_count * 2];
                        values[2 * builder.pairs.residue_map_id.value_range.start] =
                            Complex::new_re(F(id.0 as f64));
                        values[2 * builder.pairs.orientations.value_range.start] =
                            Complex::new_re(F(1.0));
                        values[2 * builder.pairs.additional_params.value_range.start + 1] =
                            Complex::new_re(F(1.0));
                        let result = <f64 as GenericEvaluatorFloat>::get_evaluator(
                            &mut dual.single_parametric,
                        )(&values);
                        let [DualOrNot::Dual(actual)] = result.as_slice() else {
                            panic!("expected dual result")
                        };
                        let k = a * (4 + b);
                        let k_derivative = 4 + b + 3 * a + 5 * a * b;
                        let d = (10 + i) as f64;
                        let expected = (2 * weight0 * weight_derivative0 * k) as f64 / d
                            + (weight0 * weight0) as f64
                                * (k_derivative as f64 / d - k as f64 / (d * d));
                        assert!(
                            (actual.values[1].re.0 - expected).abs()
                                <= 1e-11 * expected.abs().max(1.0)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn shared_numerator_family_preserves_free_kinematic_derivatives_and_alias_arguments() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::numerator_family_x"));
        let q = function!(symbol!("evaluator_test::numerator_coefficient"), 2, 0);
        let family = symbol!("evaluator_test::numerator_family");
        let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(2)));
        assert!(weight.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES);
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 7, q.clone()),
            rhs: &q * &weight + q.pow(2) * &x,
            tags: vec![Atom::num(7)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let atom = function!(family, 7, 3) + 2 * function!(family, 7, -1);
        let expected_body: Atom = weight + 11 * &x;
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        builder.pairs.update_ranges();
        let dual_shape = Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1));
        for do_fn_map_replacements in [false, true] {
            let settings = EvaluatorSettings {
                store_atom: true,
                do_fn_map_replacements,
                ..Default::default()
            };
            let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
                std::slice::from_ref(&atom),
                &builder,
                std::slice::from_ref(&definition),
                &settings,
                symbol!("evaluator_test::numerator_family_scalar"),
                InliningPolicy::Always,
            )
            .unwrap();
            assert!(prepared.reps.iter().any(|entry| {
                entry.rhs.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES
                    && entry.args == definition.args
            }));
            let mut evaluator = GenericEvaluator::new_from_builder(
                atoms,
                &prepared,
                dual_shape.clone(),
                settings.optimization_settings(),
                &settings,
            )
            .unwrap();
            let mut expected = GenericEvaluator::new_from_raw_params(
                [expected_body.clone()],
                std::slice::from_ref(&x),
                &FunctionMap::default(),
                vec![],
                settings.optimization_settings(),
                dual_shape.clone().map(|shape| (shape, Vec::new())),
                &settings,
            )
            .unwrap();
            let values = [Complex::new_re(F(0.0)), Complex::new_re(F(1.0))];
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
            let expected = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut expected)(&values);
            let [DualOrNot::Dual(actual)] = actual.as_slice() else {
                panic!("expected dual result")
            };
            let [DualOrNot::Dual(expected)] = expected.as_slice() else {
                panic!("expected dual result")
            };
            assert_eq!(actual.values, expected.values);
            assert_eq!(
                actual.values[1],
                Complex::new_re(F((1..512).map(|i| 2 * i).sum::<i32>() as f64 + 11.0))
            );
        }
    }

    #[test]
    fn two_shared_color_families_preserve_exact_fundamental_casimir() {
        use idenso::{color::CS, representations::ColorFundamental};

        test_initialise().unwrap();
        let family = symbol!("evaluator_test::color_numerator_family");
        let q = Atom::var(symbol!("evaluator_test::color_numerator_q"));
        let left = SPENSO_TAG.tensor_symbol("evaluator_test::color_numerator_left");
        let right = SPENSO_TAG.tensor_symbol("evaluator_test::color_numerator_right");
        let fundamental = ColorFundamental {}.new_rep(3);
        let first = CS.t_pattern(3, 8, 0, 1, 2);
        let second = CS.t_pattern(3, 8, 0, 2, 3);
        let definitions = [first.clone(), second.clone()]
            .into_iter()
            .enumerate()
            .map(|(i, tensor)| {
                Arc::new(FnMapEntry {
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                    lhs: function!(family, i, q.clone()),
                    rhs: &q * tensor,
                    tags: vec![Atom::num(i)],
                    args: vec![Indeterminate::try_from(q.clone()).unwrap()],
                })
            })
            .collect::<Vec<_>>();
        let outer = function!(
            left,
            fundamental
                .dual()
                .slot::<Aind, _>(Aind::Normal(1))
                .to_atom()
        ) * function!(
            right,
            fundamental.slot::<Aind, _>(Aind::Normal(3)).to_atom()
        );
        let shared = function!(family, 0, Atom::num(3) / Atom::num(2))
            * function!(family, 1, Atom::num(2) / Atom::num(3))
            * &outer;
        let unshared = first * second * outer;
        for do_algebra in [false, true] {
            let settings = EvaluatorSettings {
                do_algebra,
                ..Default::default()
            };
            let (roots, builder) = EvaluatorStack::preprocess_numerator_families(
                std::slice::from_ref(&shared),
                &ParamBuilder::new_empty(),
                &definitions,
                &settings,
                symbol!("evaluator_test::shared_color_scalar"),
                InliningPolicy::Never,
            )
            .unwrap();
            let direct = EvaluatorStack::preprocess_atom(
                &unshared,
                0,
                &settings,
                symbol!("evaluator_test::direct_color_scalar"),
                None,
            )
            .unwrap()
            .into_inner();
            let replacements = builder
                .reps
                .iter()
                .map(FnMapEntry::replacement)
                .collect::<Vec<_>>();
            let mut resolved = roots[0].clone().into_inner();
            // Resolve only this finite SU(3) tensor-component fixture. The bound
            // covers its flat definitions and retained aliases without expanding
            // any momentum numerator or introducing an unbounded rewrite loop.
            for _ in 0..=builder.reps.len() {
                let next = resolved.replace_multiple(&replacements);
                if next == resolved {
                    break;
                }
                resolved = next;
            }
            assert!(builder.reps.iter().all(|entry| {
                !resolved.contains_symbol(entry.lhs.as_fun_view().unwrap().get_symbol())
            }));
            // Every diagonal and off-diagonal matrix element is checked exactly:
            // sum_a T^a T^a = (4/3) identity in the fundamental representation.
            for i in 0..3 {
                for j in 0..3 {
                    let values = (0..3)
                        .flat_map(|k| {
                            [
                                Replacement::new(
                                    function!(
                                        left,
                                        function!(
                                            AIND_SYMBOLS.cind,
                                            Atom::from(DualConciousIndex::Down(k))
                                        )
                                    )
                                    .to_pattern(),
                                    Atom::num(usize::from(k == i)),
                                ),
                                Replacement::new(
                                    function!(right, function!(AIND_SYMBOLS.cind, k)).to_pattern(),
                                    Atom::num(usize::from(k == j)),
                                ),
                            ]
                        })
                        .collect::<Vec<_>>();
                    let expected = if i == j {
                        Atom::num(4) / Atom::num(3)
                    } else {
                        Atom::Zero
                    };
                    assert_eq!(resolved.replace_multiple(&values), expected);
                    assert_eq!(direct.replace_multiple(&values), expected);
                }
            }
        }
    }

    #[test]
    fn shared_numerator_family_preserves_dual_component_indices() {
        use spenso::structure::representation::Lorentz;

        test_initialise().unwrap();
        let q = Atom::var(symbol!("evaluator_test::dual_numerator_q"));
        let family = symbol!("evaluator_test::dual_numerator_family");
        let lower = SPENSO_TAG.tensor_symbol("evaluator_test::dual_numerator_lower");
        let upper = SPENSO_TAG.tensor_symbol("evaluator_test::dual_numerator_upper");
        let representation = Lorentz {}.new_rep(4);
        let lower_slot = representation
            .dual()
            .slot::<Aind, _>(Aind::Normal(1))
            .to_atom();
        let upper_slot = representation.slot::<Aind, _>(Aind::Normal(1)).to_atom();
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 5, q.clone()),
            rhs: &q * function!(lower, lower_slot),
            tags: vec![Atom::num(5)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let atom = function!(family, 5, 3) * function!(upper, upper_slot);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = (0..4)
            .flat_map(|i| {
                [
                    function!(
                        lower,
                        function!(AIND_SYMBOLS.cind, Atom::from(DualConciousIndex::Down(i)))
                    ),
                    function!(upper, function!(AIND_SYMBOLS.cind, i)),
                ]
            })
            .collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[atom],
            &builder,
            &[definition],
            &settings,
            symbol!("evaluator_test::dual_numerator_scalar"),
            InliningPolicy::Never,
        )
        .unwrap();
        let mut evaluator = GenericEvaluator::new_from_builder(
            atoms,
            &prepared,
            None,
            settings.optimization_settings(),
            &settings,
        )
        .unwrap();
        let values =
            [1.0, 5.0, 2.0, 6.0, 3.0, 7.0, 4.0, 8.0].map(|value| Complex::new_re(F(value)));
        assert_eq!(
            <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&values),
            Complex::new_re(F(210.0)),
        );
    }

    #[test]
    fn complete_source_sums_share_numerators_and_count_each_row_once() {
        test_initialise().unwrap();
        let q = Atom::var(symbol!("evaluator_test::complete_sum_q"));
        let r = Atom::var(symbol!("evaluator_test::complete_sum_r"));
        let x = Atom::var(symbol!("evaluator_test::complete_sum_x"));
        let family = symbol!("evaluator_test::complete_sum_family");
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 0, q.clone(), r.clone()),
            rhs: (&q * &x + 1).pow(4) * (&r * &x + 2),
            tags: vec![Atom::num(0)],
            args: vec![
                Indeterminate::try_from(q.clone()).unwrap(),
                Indeterminate::try_from(r).unwrap(),
            ],
        });
        // The narrow family's single row runs beside the wide Cartesian sum,
        // retaining its own numerator and open tensor interface.
        let temporal = GS.energy_delta(parse_lit!(spenso::mink(4, 1)));
        let narrow_definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 1, q.clone()),
            rhs: (&q * &x + 3).pow(2) * &temporal,
            tags: vec![Atom::num(1)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let constant_definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 2),
            rhs: x.pow(2) + 7,
            tags: vec![Atom::num(2)],
            args: Vec::new(),
        });
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        let tree_argument = Atom::var(symbol!("evaluator_test::complete_sum_tree_argument"));
        builder
            .add_function_entry(FnMapEntry {
                lhs: function!(GS.tree_denom_wrapper, tree_argument.clone()),
                rhs: tree_argument.clone(),
                args: vec![Indeterminate::try_from(tree_argument).unwrap()],
                tags: Vec::new(),
                inlining: InliningPolicy::Always,
                is_alias: false,
            })
            .unwrap();
        let parameter_count = builder.pairs.update_ranges();
        let orientations =
            TiVec::<OrientationID, _>::from_iter([EdgeVec::from_iter([Orientation::Default])]);
        let filter = SubSet::full(orientations.len());
        for (inline_numerator_functions, with_tree_factor) in
            [(true, false), (false, false), (true, true), (false, true)]
        {
            let settings = EvaluatorSettings {
                iterative_orientation_optimization: true,
                inline_numerator_functions,
                summed_function_map: true,
                summed: true,
                ..Default::default()
            };
            let tree_factor = if with_tree_factor {
                function!(GS.tree_denom_wrapper, (x.pow(2) + 4).pow(-1))
            } else {
                Atom::one()
            };
            let mut numerator_program_shapes = None;
            for rows in [2, 4, 8] {
                // The UV child owns its sum, while this addback has no production
                // row at all. A global loop over production rows duplicates both.
                let term = ParametricIntegrandTerm {
                    numerator: Arc::clone(&definition),
                    rows: (1..=rows)
                        .flat_map(|i| {
                            let carrier = &tree_factor * i;
                            (0..=3).map(move |j| crate::uv::ParametricIntegrandRow {
                                coefficients: vec![Rational::from(i), Rational::from(j)],
                                carrier: carrier.clone(),
                                source_keys: Vec::new(),
                            })
                        })
                        .collect(),
                };
                let terms = [
                    term,
                    ParametricIntegrandTerm {
                        numerator: Arc::clone(&narrow_definition),
                        rows: vec![crate::uv::ParametricIntegrandRow {
                            coefficients: vec![Rational::one()],
                            carrier: 3 * &temporal,
                            source_keys: Vec::new(),
                        }],
                    },
                    ParametricIntegrandTerm {
                        numerator: Arc::clone(&constant_definition),
                        // Duplicate rows and cancellation must keep this
                        // independent source sum once beside the scalar 5.
                        rows: [2, 3, 7, -7]
                            .into_iter()
                            .map(|carrier| crate::uv::ParametricIntegrandRow {
                                coefficients: Vec::new(),
                                carrier: Atom::num(carrier),
                                source_keys: Vec::new(),
                            })
                            .collect(),
                    },
                ];
                // At x=2 the nontrivial tree denominator is exactly 1/8.
                // It multiplies only the wide source, leaving both independent
                // sources and the scalar addback unchanged.
                let wide = (1..=rows).map(|i: i64| i * (2 * i + 1).pow(4)).sum::<i64>()
                    * (0..=3).map(|i| 2 * i + 2).sum::<i64>();
                let expected = wide as f64 / if with_tree_factor { 8.0 } else { 1.0 }
                    + (3 * 5i64.pow(2) + 5 * (4 + 7) + 5) as f64;
                let (mut stack, _) = EvaluatorStack::from_integrand_with_timings(
                    &Atom::num(5),
                    &builder,
                    &[],
                    &terms,
                    None,
                    None,
                    &settings,
                )
                .unwrap();
                assert_eq!(
                    stack.parametric_rows.as_ref().unwrap().rows.len(),
                    rows as usize * 4
                );
                let program = stack.single_parametric.f64_eager.export_instructions();
                // Final inlining is independent of the earlier shared tensor
                // lowering and complete source rows. When bodies remain retained,
                // their shapes must not grow with the number of calling rows.
                let mut shapes = program
                    .sub_evaluators
                    .iter()
                    .map(|sub| (sub.tags.clone(), sub.instructions.instructions.len()))
                    .collect::<Vec<_>>();
                shapes.sort();
                if inline_numerator_functions {
                    assert!(shapes.is_empty());
                } else {
                    assert!(!shapes.is_empty());
                    assert_eq!(
                        numerator_program_shapes.get_or_insert_with(|| shapes.clone()),
                        &shapes
                    );
                }
                let program_count = shapes.len();
                if rows == 8 {
                    let encoded =
                        bincode::encode_to_vec(&stack, bincode::config::standard()).unwrap();
                    let mut state = Vec::new();
                    State::export(&mut state).unwrap();
                    let state_map = State::import(&mut Cursor::new(state), None).unwrap();
                    let model = Model::default();
                    (stack, _) = bincode::decode_from_slice_with_context(
                        &encoded,
                        bincode::config::standard(),
                        GammaLoopContextContainer {
                            state_map: &state_map,
                            model: &model,
                        },
                    )
                    .unwrap();
                    assert_eq!(
                        stack
                            .single_parametric
                            .f64_eager
                            .export_instructions()
                            .sub_evaluators
                            .len(),
                        program_count
                    );
                }
                let input = || {
                    let mut values = vec![Complex::new_re(F(0.0)); parameter_count];
                    values[builder.pairs.additional_params.value_range.start] =
                        Complex::new_re(F(2.0));
                    InputParams {
                        values: SliceMut::Owned(values),
                        residue_map_id_start: builder.pairs.residue_map_id.value_range.start,
                        orientations_start: builder.pairs.orientations.value_range.start,
                        multiplicative_offset: 1,
                    }
                };
                let mut metadata = EvaluationMetaData::new_empty();
                for compiled in [false, true] {
                    if compiled {
                        stack
                            .for_each_generic_evaluator_mut(|evaluator| {
                                evaluator.activate_symjit(&CompilationOptionsSnapshot {
                                    optimization_level: CompilationOptimizationLevel::O0,
                                    ..Default::default()
                                })
                            })
                            .unwrap();
                    }
                    for method in [
                        EvaluatorMethod::SingleParametric,
                        EvaluatorMethod::Iterative,
                        EvaluatorMethod::SummedFunctionMap,
                        EvaluatorMethod::Summed,
                    ] {
                        let mut runtime = RuntimeSettings::default();
                        runtime.general.evaluator_method = method.clone();
                        let result = stack
                            .evaluate(
                                input(),
                                SingleOrAllOrientations::All {
                                    all: &orientations,
                                    filter: &filter,
                                },
                                &runtime,
                                &mut metadata,
                            )
                            .unwrap();
                        assert_eq!(
                            scalar_value(result),
                            Complex::new_re(F(expected)),
                            "{method:?}, {rows} rows, inline={inline_numerator_functions}, tree={with_tree_factor}, compiled={compiled}"
                        );
                    }
                }
                if rows == 8 {
                    let directory = std::env::temp_dir().join(format!(
                        "gammaloop-complete-sum-{}-{}",
                        std::process::id(),
                        std::time::SystemTime::now()
                            .duration_since(std::time::UNIX_EPOCH)
                            .unwrap()
                            .as_nanos()
                    ));
                    std::fs::create_dir_all(&directory).unwrap();
                    let options = crate::settings::global::CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O0,
                        compiler: "g++".into(),
                        ..Default::default()
                    };
                    for mode in [
                        FrozenCompilationMode::Cpp(options.clone()),
                        FrozenCompilationMode::Assembly(options),
                    ] {
                        let name = mode.active_backend_name();
                        let mut program_index = 0;
                        stack
                            .for_each_generic_evaluator_mut(|evaluator| {
                                let label = format!("{name}-{program_index}");
                                program_index += 1;
                                evaluator.compile_external(
                                    directory.join(format!("{label}.cpp")),
                                    "complete_sum",
                                    directory.join(format!("{label}.so")),
                                    &mode,
                                )
                            })
                            .unwrap();
                        assert_eq!(program_index, stack.generic_evaluator_count());
                        assert_eq!(
                            scalar_value(
                                stack
                                    .evaluate(
                                        input(),
                                        SingleOrAllOrientations::All {
                                            all: &orientations,
                                            filter: &filter,
                                        },
                                        &RuntimeSettings::default(),
                                        &mut metadata,
                                    )
                                    .unwrap()
                            ),
                            Complex::new_re(F(expected)),
                            "{name} complete source sum, inline={inline_numerator_functions}, tree={with_tree_factor}"
                        );
                    }
                    std::fs::remove_dir_all(directory).unwrap();
                }
            }
        }
    }

    #[test]
    fn flat_native_outputs_preserve_buffers_precision_and_mixed_jets() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::flat_x"));
        let y = Atom::var(symbol!("evaluator_test::flat_y"));
        let gate = Atom::var(symbol!("evaluator_test::flat_gate"));
        let precise = Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64);
        let roots = [
            x.pow(2) * &y + Atom::num(precise.clone()),
            Symbol::IF.call_args([gate.clone(), x.pow(-1), y.pow(2)]),
        ];
        fn check<T: FloatLike>(evaluator: &mut GenericEvaluator, precise: &Rational) {
            let zero = F(T::new_zero());
            let width = evaluator.dual_shape.as_ref().map_or(1, Vec::len);
            assert!(width == 1 || width == 4);
            let sentinel = Complex::new_re(zero.from_usize(17));
            let mut output = vec![sentinel.clone(); evaluator.compute_out_size() + 2];
            for point in [0, 2, -2, 0] {
                let x = zero.from_i64(point);
                let y = zero.from_usize(3);
                let mut input = vec![Complex::new_re(zero.clone()); 3 * width];
                input[0] = Complex::new_re(x.clone());
                input[width] = Complex::new_re(y.clone());
                input[2 * width] = Complex::new_re(zero.from_usize(usize::from(point != 0)));
                if width > 1 {
                    input[1] = Complex::new_re(zero.one());
                    input[width + 2] = Complex::new_re(zero.one());
                }
                let output_end = output.len() - 1;
                <T as GenericEvaluatorFloat>::evaluate_into(
                    evaluator,
                    &input,
                    &mut output[1..output_end],
                );
                assert_eq!(output[0], sentinel);
                assert_eq!(output[output_end], sentinel);
                let wrapped = <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&input)
                    .into_iter()
                    .flat_map(|value| match value {
                        DualOrNot::Dual(dual) => dual.values,
                        DualOrNot::NonDual(value) => vec![value],
                    })
                    .collect::<Vec<_>>();
                assert_eq!(output[1..output_end], wrapped);
                let first = [
                    x.powi(2) * &y + zero.from_rational(precise),
                    zero.from_usize(2) * &x * &y,
                    x.powi(2),
                    zero.from_usize(2) * &x,
                ];
                // x=0 must never evaluate the pole, including its physical jets.
                let second = if point == 0 {
                    [
                        y.powi(2),
                        zero.clone(),
                        zero.from_usize(2) * &y,
                        zero.clone(),
                    ]
                } else {
                    [x.powi(-1), -x.powi(-2), zero.clone(), zero.clone()]
                };
                for (actual, expected) in wrapped.iter().zip(
                    first
                        .into_iter()
                        .take(width)
                        .chain(second.into_iter().take(width)),
                ) {
                    assert!(
                        (&actual.re - &expected).abs()
                            <= (expected.abs() + zero.one())
                                * zero.epsilon()
                                * zero.from_usize(128),
                        "flat output differs at x={point}, width={width}: {actual:?} != {expected:?}"
                    );
                    assert_eq!(actual.im, zero);
                }
                if point == 0 && zero.epsilon() < zero.from_usize(2).powi(-60) {
                    assert!(wrapped[0].re > zero.one());
                }
            }
            // Caller-owned extents must be checked before a compiled backend
            // receives the pointer, not merely by the eager implementation.
            let input = vec![Complex::new_re(zero); 3 * width];
            assert!(
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    <T as GenericEvaluatorFloat>::evaluate_into(evaluator, &input, &mut []);
                }))
                .is_err()
            );
        }
        let directory = std::env::temp_dir().join(format!(
            "gammaloop-flat-native-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
        ));
        std::fs::create_dir_all(&directory).unwrap();
        for shape in [
            None,
            Some(vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]]),
        ] {
            let settings = EvaluatorSettings::default();
            let mut evaluator = GenericEvaluator::new_from_raw_params(
                roots.clone(),
                &[x.clone(), y.clone(), gate.clone()],
                &FunctionMap::new(),
                vec![],
                settings.optimization_settings(),
                shape.map(|shape| {
                    let zeros = (1..shape.len()).map(|component| (2, component)).collect();
                    (shape, zeros)
                }),
                &settings,
            )
            .unwrap();
            check::<f64>(&mut evaluator, &precise);
            check::<f128>(&mut evaluator, &precise);
            check::<ArbPrec>(&mut evaluator, &precise);
            // Sampling owns this warmup in production. Exercise its existing
            // source-only lane because the shared wrapper dispatches it too.
            let sampling = evaluator
                .rational
                .as_ref()
                .unwrap()
                .clone()
                .map_coeff(&|r| {
                    Complex::new(
                        F::<SamplingFloat>::from(&r.re),
                        F::<SamplingFloat>::from(&r.im),
                    )
                });
            evaluator.sampling_fixed256.set(sampling);
            check::<SamplingFloat>(&mut evaluator, &precise);
            for optimization in [
                CompilationOptimizationLevel::O0,
                CompilationOptimizationLevel::O2,
            ] {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: optimization,
                        ..Default::default()
                    })
                    .unwrap();
                check::<f64>(&mut evaluator, &precise);
            }
            let options = crate::settings::global::CompilationOptionsSnapshot {
                optimization_level: CompilationOptimizationLevel::O0,
                compiler: "g++".into(),
                ..Default::default()
            };
            for mode in [
                FrozenCompilationMode::Cpp(options.clone()),
                FrozenCompilationMode::Assembly(options),
            ] {
                let label = format!(
                    "{}-{}",
                    mode.active_backend_name(),
                    evaluator.compute_out_size()
                );
                evaluator
                    .compile_external(
                        directory.join(format!("{label}.cpp")),
                        "flat_native_outputs",
                        directory.join(format!("{label}.so")),
                        &mode,
                    )
                    .unwrap();
                check::<f64>(&mut evaluator, &precise);
            }
            evaluator.activate_eager();
            check::<f64>(&mut evaluator, &precise);
        }
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn generated_numerator_tables_preserve_argument_roles_and_native_jets() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::batched_role_x"));
        let y = Atom::var(symbol!("evaluator_test::batched_role_y"));
        let columns = (0..4)
            .map(|i| function!(symbol!("evaluator_test::batched_role_column"), i))
            .collect::<Vec<_>>();
        let [q, r, s, t] = columns.as_slice() else {
            unreachable!()
        };
        let formals = (0..3)
            .map(|i| function!(symbol!("evaluator_test::batched_role_formal"), i))
            .collect::<Vec<_>>();
        let [a, b, c] = formals.as_slice() else {
            unreachable!()
        };
        let n = symbol!("evaluator_test::batched_role_numerator");
        let definition = FnMapEntry {
            lhs: function!(n, a.clone(), b.clone(), c.clone()),
            rhs: a * x.pow(2) * &y + b * &x * y.pow(2) + c * &x * &y,
            tags: Vec::new(),
            args: formals
                .iter()
                .cloned()
                .map(|a| Indeterminate::try_from(a).unwrap())
                .collect(),
            inlining: InliningPolicy::Always,
            is_alias: false,
        };
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x, y].into_iter().collect();
        let input_count = builder.pairs.update_ranges();
        let x_index = builder.pairs.additional_params.value_range.start;
        builder.add_function_entry(definition.clone()).unwrap();
        let atoms = [
            function!(n, q.clone(), r.clone(), q.clone()),
            function!(n, s.clone(), t.clone(), s.clone()),
            function!(n, q.clone(), r.clone(), r.clone()),
            function!(n, q.clone(), 1, q.clone()),
        ];
        fn check<T: FloatLike>(stack: &mut EvaluatorStack, input_count: usize, x_index: usize) {
            let zero = F(T::new_zero());
            let width = stack
                .single_parametric
                .dual_shape
                .as_ref()
                .map_or(1, Vec::len);
            for point in [0, 2, -2, 2] {
                let x = zero.from_i64(point);
                let y = zero.from_usize(3);
                let mut input = vec![Complex::new_re(zero.clone()); input_count * width];
                input[x_index * width] = Complex::new_re(x.clone());
                input[(x_index + 1) * width] = Complex::new_re(y.clone());
                if width > 1 {
                    input[x_index * width + 1] = Complex::new_re(zero.one());
                    input[(x_index + 1) * width + 2] = Complex::new_re(zero.one());
                }
                let outputs =
                    stack.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
                assert_eq!(outputs.len(), 4);
                // Independent coefficient sums for the four original calls.
                for (output, [a, b, c]) in
                    outputs
                        .iter()
                        .zip([[3, 5, 3], [3, 5, 3], [3, 5, 5], [3, 3, 3]])
                {
                    let a = zero.from_usize(a);
                    let b = zero.from_usize(b);
                    let c = zero.from_usize(c);
                    let expected = [
                        &a * x.powi(2) * &y + &b * &x * y.powi(2) + &c * &x * &y,
                        zero.from_usize(2) * &a * &x * &y + &b * y.powi(2) + &c * &y,
                        &a * x.powi(2) + zero.from_usize(2) * &b * &x * &y + &c * &x,
                        zero.from_usize(2) * (&a * &x + &b * &y) + &c,
                    ];
                    let values = match output {
                        DualOrNot::Dual(value) => value.values.as_slice(),
                        DualOrNot::NonDual(value) => std::slice::from_ref(value),
                    };
                    assert_eq!(values.len(), width);
                    for (actual, expected) in values.iter().zip(expected) {
                        assert_eq!(*actual, Complex::new_re(expected));
                    }
                }
            }
        }
        for shape in [
            None,
            Some(vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]]),
        ] {
            let mut rows = ParametricResidueRows {
                parameters: columns.clone(),
                input_columns: (0..4).collect(),
                rows: [[1, 1, 1, 1], [2, 3, 2, 3], [0, 1, 0, 1]]
                    .map(|row| row.map(Rational::from).to_vec())
                    .to_vec(),
                ..Default::default()
            };
            let (main, auxiliary) = EvaluatorStack::new_single_parametric(
                atoms.iter().cloned().map(AliasedAtom::from).collect(),
                &builder,
                &[],
                Some(&mut rows),
                &shape,
                &EvaluatorSettings {
                    store_atom: true,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(builder.reps, vec![definition.clone()]);
            let (ports, program) = auxiliary
                .as_ref()
                .expect("repeated native numerator calls are batched");
            assert_eq!(
                ports.len(),
                4,
                "original carrier call identities stay distinct"
            );
            assert_eq!(rows.auxiliary_shared_outputs, 0);
            assert!(rows.input_columns.is_empty());
            assert_eq!(program.exprs_len, 9);
            assert_eq!(rows.auxiliary_gathers.len(), 3);
            let carrier_columns = main
                .exprs
                .as_ref()
                .unwrap()
                .iter()
                .map(|root| {
                    ports
                        .iter()
                        .position(|port| port == root)
                        .expect("each output keeps its original caller port")
                })
                .collect::<Vec<_>>();
            for gather in &rows.auxiliary_gathers {
                let gather = carrier_columns
                    .iter()
                    .map(|&column| gather[column])
                    .collect::<Vec<_>>();
                assert_eq!(gather.len(), 4);
                assert!(gather.iter().all(Option::is_some));
                assert_eq!(
                    gather[0], gather[1],
                    "equal argument-role patterns share outputs"
                );
                assert_ne!(
                    gather[0], gather[2],
                    "distinct argument columns do not become equal"
                );
                assert_ne!(
                    gather[0], gather[3],
                    "a dynamic value of one remains distinct from literal one"
                );
            }
            let mut stack = EvaluatorStack {
                explicit_orientation_sum_only: false,
                parametric_rows: Some(rows),
                production_orientation_ids: Vec::new(),
                single_parametric: main,
                single_parametric_auxiliary: auxiliary,
                iterative: None,
                summed_function_map: None,
                summed: None,
            };
            check::<f64>(&mut stack, input_count, x_index);
            check::<f128>(&mut stack, input_count, x_index);
            check::<ArbPrec>(&mut stack, input_count, x_index);
            for optimization in [
                CompilationOptimizationLevel::O0,
                CompilationOptimizationLevel::O2,
            ] {
                stack
                    .for_each_generic_evaluator_mut(|evaluator| {
                        evaluator.activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: optimization,
                            ..Default::default()
                        })
                    })
                    .unwrap();
                check::<f64>(&mut stack, input_count, x_index);
            }
        }
    }

    #[test]
    fn row_scalar_leaves_preserve_alias_scopes_precision_and_mixed_jets() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::row_scalar_x"));
        let y = Atom::var(symbol!("evaluator_test::row_scalar_y"));
        let q = function!(symbol!("evaluator_test::row_scalar_q"));
        let gate = Atom::var(symbol!("evaluator_test::row_scalar_gate"));
        let alias = function!(symbol!("evaluator_test::row_scalar_alias"));
        let physical_alias = function!(symbol!("evaluator_test::row_scalar_physical_alias"));
        let formal = Atom::var(symbol!("evaluator_test::row_scalar_formal"));
        let other_formal = Atom::var(symbol!("evaluator_test::row_scalar_other_formal"));
        let one = symbol!("evaluator_test::row_scalar_one");
        let two = symbol!("evaluator_test::row_scalar_two");
        let precise = Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone(), y.clone()].into_iter().collect();
        let input_count = builder.pairs.update_ranges();
        let x_index = builder.pairs.additional_params.value_range.start;
        // Exact native inputs take precedence over ordinary definitions with
        // the same key, including during scalar dependency rebinding.
        builder
            .add_function_entry(FnMapEntry {
                lhs: q.clone(),
                rhs: x.clone(),
                tags: Vec::new(),
                args: Vec::new(),
                inlining: InliningPolicy::Always,
                is_alias: false,
            })
            .unwrap();
        let scalar_functions = [
            FnMapEntry {
                lhs: function!(GS.tree_denom_wrapper, formal.clone()),
                rhs: formal.clone(),
                tags: Vec::new(),
                args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            },
            FnMapEntry {
                lhs: function!(one, formal.clone()),
                rhs: formal.pow(2) * &y,
                tags: Vec::new(),
                args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            },
            FnMapEntry {
                lhs: function!(two, 7, formal.clone(), other_formal.clone()),
                rhs: &formal * &other_formal + formal.pow(2),
                tags: vec![Atom::num(7)],
                args: vec![
                    Indeterminate::try_from(formal.clone()).unwrap(),
                    Indeterminate::try_from(other_formal).unwrap(),
                ],
                inlining: InliningPolicy::Always,
                is_alias: false,
            },
        ];
        for definition in &scalar_functions {
            builder.add_function_entry(definition.clone()).unwrap();
        }
        let physical_functions = builder.reps.clone();
        let n1 = symbol!("evaluator_test::row_scalar_generated_one");
        let n2 = symbol!("evaluator_test::row_scalar_generated_two");
        let numerators = [
            FnMapEntry {
                lhs: function!(n1, formal.clone()),
                rhs: &formal * &x,
                tags: Vec::new(),
                args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            },
            FnMapEntry {
                lhs: function!(n2, formal.clone()),
                rhs: &formal + &y,
                tags: Vec::new(),
                args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            },
        ];
        for definition in &numerators {
            builder.add_function_entry(definition.clone()).unwrap();
        }
        let row_tree_factor = function!(
            GS.tree_denom_wrapper,
            (q.pow(2) * x.pow(2) + 1).pow(-2) * (y.pow(2) + 2).pow(-1)
        );
        let mut selected =
            AliasedAtom::from(Symbol::IF.call_args([gate.clone(), alias.clone(), -&q]));
        selected.register_alias(alias.clone(), &q * &physical_alias);
        selected.register_alias(physical_alias, x.pow(2) * &y);
        let mut other_owner = AliasedAtom::from(alias.clone());
        other_owner.register_alias(alias.clone(), &q * &x * y.pow(2));
        let tree_factor = function!(
            GS.tree_denom_wrapper,
            (x.pow(2) + 1).pow(-2) * (y.pow(2) + 2).pow(-1)
        );
        let atoms = vec![
            q.clone().into(),
            (-&q).into(),
            selected,
            other_owner,
            // This unchanged identity wrapper blocked the complete scalar
            // branch in the physical massless soft-counterterm fixture.
            Symbol::IF
                .call_args([
                    gate.clone(),
                    function!(GS.tree_denom_wrapper, 1) * function!(one, &q * &x),
                    -&q * &y,
                ])
                .into(),
            (function!(two, 7, &q * &x, y.clone()) * function!(GS.tree_denom_wrapper, 1)).into(),
            // Keep the nontrivial product of tree propagators outside the
            // complete source sum, including both ordinary scalar calls.
            (&tree_factor * (function!(one, &q * &x) + function!(two, 7, &q * &x, y.clone())))
                .into(),
            // The wrapper factor multiplies the complete selected scalar sum
            // and two generated numerator calls; none is distributed or inlined.
            (&row_tree_factor
                * function!(n1, q.clone())
                * function!(n2, q.clone())
                * (Symbol::IF.call_args([gate.clone(), &x + &y, &x - &y]) + 1))
                .into(),
        ];
        fn check<T: FloatLike>(
            stack: &mut EvaluatorStack,
            baseline: &mut GenericEvaluator,
            original_rows: &ParametricResidueRows,
            input_count: usize,
            x_index: usize,
            precise: &Rational,
        ) {
            let zero = F(T::new_zero());
            let width = stack
                .single_parametric
                .dual_shape
                .as_ref()
                .map_or(1, Vec::len);
            let p = zero.from_rational(precise);
            let sum = &p + zero.from_usize(2);
            let squares = p.powi(2) + zero.from_usize(4);
            for point in [0, 2, -2, 2] {
                let x = zero.from_i64(point);
                let y = zero.from_usize(3);
                let two = zero.from_usize(2);
                let mut input = vec![Complex::new_re(zero.clone()); input_count * width];
                input[x_index * width] = Complex::new_re(x.clone());
                input[(x_index + 1) * width] = Complex::new_re(y.clone());
                if width > 1 {
                    input[x_index * width + 1] = Complex::new_re(zero.one());
                    input[(x_index + 1) * width + 2] = Complex::new_re(zero.one());
                }
                let output =
                    stack.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
                let unbatched = original_rows.evaluate(
                    baseline,
                    &input,
                    &[],
                    &mut EvaluationMetaData::new_empty(),
                );
                // The seventh source is F(x,y) P(x,y). Compute both jets
                // independently, then use the mixed product rule.
                let xden = x.powi(2) + zero.one();
                let yden = y.powi(2) + &two;
                let factor = [
                    xden.powi(-2) * yden.powi(-1),
                    -zero.from_usize(4) * &x * xden.powi(-3) * yden.powi(-1),
                    -&two * &y * xden.powi(-2) * yden.powi(-2),
                    zero.from_usize(8) * &x * &y * xden.powi(-3) * yden.powi(-2),
                ];
                let polynomial = [
                    &sum * &x * &y + &squares * x.powi(2) * (&y + zero.one()),
                    &sum * &y + &two * &squares * &x * (&y + zero.one()),
                    &sum * &x + &squares * x.powi(2),
                    &sum + &two * &squares * &x,
                ];
                // Sum the two complete graph products independently. Here
                // N1(q)=q*x, N2(q)=q+y and the scalar branch is x+sign*y+1.
                let mut graph_product: [_; 4] = std::array::from_fn(|_| zero.clone());
                for (coefficient, sign) in [(p.clone(), -1), (two.clone(), 1)] {
                    let sign = zero.from_i64(sign);
                    let q2 = coefficient.powi(2);
                    let u = &q2 * x.powi(2) + zero.one();
                    let v = y.powi(2) + &two;
                    let f = [
                        u.powi(-2) * v.powi(-1),
                        -zero.from_usize(4) * &q2 * &x * u.powi(-3) * v.powi(-1),
                        -&two * &y * u.powi(-2) * v.powi(-2),
                        zero.from_usize(8) * &q2 * &x * &y * u.powi(-3) * v.powi(-2),
                    ];
                    let h = &coefficient + &y;
                    let k = &x + &sign * &y + zero.one();
                    let p = [
                        &coefficient * &x * &h * &k,
                        &coefficient * &h * (&k + &x),
                        &coefficient * &x * (&k + &sign * &h),
                        &coefficient * (&k + &x + &sign * &h),
                    ];
                    for (sum, value) in graph_product.iter_mut().zip([
                        &f[0] * &p[0],
                        &f[1] * &p[0] + &f[0] * &p[1],
                        &f[2] * &p[0] + &f[0] * &p[2],
                        &f[3] * &p[0] + &f[1] * &p[2] + &f[2] * &p[1] + &f[0] * &p[3],
                    ]) {
                        *sum += value;
                    }
                }
                let expected = [
                    [sum.clone(), zero.clone(), zero.clone(), zero.clone()],
                    [-&sum, zero.clone(), zero.clone(), zero.clone()],
                    [
                        &two * x.powi(2) * &y - &p,
                        zero.from_usize(4) * &x * &y,
                        &two * x.powi(2),
                        zero.from_usize(4) * &x,
                    ],
                    [
                        &sum * &x * y.powi(2),
                        &sum * y.powi(2),
                        &two * &sum * &x * &y,
                        &two * &sum * &y,
                    ],
                    [
                        zero.from_usize(4) * x.powi(2) * &y - &p * &y,
                        zero.from_usize(8) * &x * &y,
                        zero.from_usize(4) * x.powi(2) - &p,
                        zero.from_usize(8) * &x,
                    ],
                    [
                        &sum * &x * &y + &squares * x.powi(2),
                        &sum * &y + &two * &squares * &x,
                        &sum * &x,
                        sum.clone(),
                    ],
                    [
                        &factor[0] * &polynomial[0],
                        &factor[1] * &polynomial[0] + &factor[0] * &polynomial[1],
                        &factor[2] * &polynomial[0] + &factor[0] * &polynomial[2],
                        &factor[3] * &polynomial[0]
                            + &factor[1] * &polynomial[2]
                            + &factor[2] * &polynomial[1]
                            + &factor[0] * &polynomial[3],
                    ],
                    graph_product,
                ];
                assert_eq!(output.len(), expected.len());
                assert_eq!(output.len(), unbatched.len());
                for ((output, unbatched), expected) in output.iter().zip(&unbatched).zip(expected) {
                    let values = match output {
                        DualOrNot::Dual(value) => value.values.as_slice(),
                        DualOrNot::NonDual(value) => std::slice::from_ref(value),
                    };
                    let reference = match unbatched {
                        DualOrNot::Dual(value) => value.values.as_slice(),
                        DualOrNot::NonDual(value) => std::slice::from_ref(value),
                    };
                    assert_eq!(values.len(), width);
                    assert_eq!(reference.len(), width);
                    for ((actual, unbatched), expected) in
                        values.iter().zip(reference).zip(expected)
                    {
                        assert!(
                            (actual.re.clone() - &unbatched.re).abs()
                                <= (expected.abs() + zero.one())
                                    * zero.epsilon()
                                    * zero.from_usize(128)
                        );
                        assert_eq!(actual.im, unbatched.im);
                        assert!(
                            (actual.re.clone() - &expected).abs()
                                <= (expected.abs() + zero.one())
                                    * zero.epsilon()
                                    * zero.from_usize(128)
                        );
                        assert_eq!(actual.im, zero);
                    }
                }
            }
        }
        let directory = std::env::temp_dir().join(format!(
            "gammaloop-row-scalars-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        for shape in [
            None,
            Some(vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]]),
        ] {
            let mut rows = ParametricResidueRows {
                parameters: vec![q.clone(), gate.clone()],
                input_columns: vec![0, 1],
                rows: vec![
                    vec![precise.clone(), Rational::zero()],
                    vec![Rational::from(2), Rational::one()],
                ],
                ..Default::default()
            };
            let original_rows = rows.clone();
            let parameters = (&builder.pairs)
                .into_iter()
                .flat_map(|pair| pair.params.clone())
                .chain(rows.parameters.iter().cloned())
                .collect::<Vec<_>>();
            let settings = EvaluatorSettings {
                store_atom: true,
                ..Default::default()
            };
            let mut baseline = GenericEvaluator::new_from_raw_params(
                atoms.clone(),
                &parameters,
                &builder.fn_map,
                builder.reps.clone(),
                settings.optimization_settings(),
                shape.clone().map(|shape| {
                    let zeros = (input_count..parameters.len())
                        .flat_map(|p| (1..shape.len()).map(move |component| (p, component)))
                        .collect();
                    (shape, zeros)
                }),
                &settings,
            )
            .unwrap();
            let (main, auxiliary) = EvaluatorStack::new_single_parametric(
                atoms.clone(),
                &builder,
                &physical_functions,
                Some(&mut rows),
                &shape,
                &settings,
            )
            .unwrap();
            let (ports, program) = auxiliary.as_ref().expect("row scalar leaves are gathered");
            assert_eq!(ports.len(), 11);
            assert_eq!(program.exprs_len, 22);
            assert!(
                main.exprs
                    .as_ref()
                    .unwrap()
                    .iter()
                    .take(7)
                    .all(|root| ports.contains(root)),
                "each complete output moves to one port, including the identity-wrapped scalar branch"
            );
            let whole_source_port = main.exprs.as_ref().unwrap()[6].clone();
            let whole_source_column = ports
                .iter()
                .position(|port| port == &whole_source_port)
                .unwrap();
            for gather in &rows.auxiliary_gathers {
                let root = &program.exprs.as_ref().unwrap()[gather[whole_source_column].unwrap()];
                let AtomView::Mul(product) = root.as_view() else {
                    panic!("the tree denominator must factorize the complete source");
                };
                assert!(product.iter().any(|part| part == tree_factor.as_view()));
                assert!(product.iter().any(|part| matches!(part, AtomView::Add(_))));
                let mut wrappers = 0;
                root.visitor(&mut |view| {
                    if let AtomView::Fun(call) = view
                        && call.get_symbol() == GS.tree_denom_wrapper
                    {
                        wrappers += 1;
                    }
                    true
                });
                assert_eq!(
                    wrappers, 1,
                    "the common factor is not distributed into its source sum"
                );
            }
            let AtomView::Mul(graph_factors) = main.exprs.as_ref().unwrap()[7].as_view() else {
                panic!("the graph numerator product must remain factored in the carrier");
            };
            assert_eq!(graph_factors.iter().count(), 4);
            let graph_columns = graph_factors
                .iter()
                .filter_map(|factor| ports.iter().position(|port| port.as_view() == factor))
                .collect::<Vec<_>>();
            assert_eq!(graph_columns.len(), 3);
            let scalar_factor = graph_factors
                .iter()
                .find(|factor| matches!(factor, AtomView::Add(_)))
                .expect("the selected scalar sum remains one carrier factor");
            let scalar_column = ports
                .iter()
                .position(|port| scalar_factor == (port + 1).as_view())
                .expect("the sum retains its distinct selector port and constant");
            for (row, gather) in rows.auxiliary_gathers.iter().enumerate() {
                let roots = graph_columns
                    .iter()
                    .map(|&column| &program.exprs.as_ref().unwrap()[gather[column].unwrap()])
                    .collect::<Vec<_>>();
                let calls = [n1, n2].map(|head| {
                    let found = roots
                        .iter()
                        .filter_map(|root| root.as_fun_view())
                        .filter(|call| call.get_symbol() == head)
                        .collect::<Vec<_>>();
                    assert_eq!(
                        found.len(),
                        1,
                        "each numerator retains one distinct call factor"
                    );
                    found[0]
                });
                let role = calls[0].iter().next().unwrap().to_owned();
                assert_eq!(calls[1].iter().next().unwrap(), role.as_view());
                assert!(rows.auxiliary_parameters.contains(&role));
                let expected_tree = row_tree_factor
                    .replace(q.to_pattern())
                    .with(role.to_pattern());
                assert_eq!(
                    roots.iter().filter(|root| ***root == expected_tree).count(),
                    1,
                    "the exact coefficient stays a native input in the single outer denominator factor"
                );
                assert_eq!(
                    program.exprs.as_ref().unwrap()[gather[scalar_column].unwrap()],
                    if row == 0 { &x - &y } else { &x + &y },
                    "only the selector leaf is gathered; the tree factor and numerator product are not distributed"
                );
            }
            for definition in scalar_functions.iter().chain(&numerators) {
                assert!(
                    program.fn_map_entries.contains(definition),
                    "ordinary function bodies and signatures remain shared and unchanged"
                );
                assert!(
                    program.exprs.as_ref().unwrap().iter().any(|root| root
                        .as_view()
                        .contains_symbol(definition.lhs.as_fun_view().unwrap().get_symbol())),
                    "auxiliary roots keep the ordinary call instead of expanding its body"
                );
            }
            assert!(
                program.fn_map_entries.iter().any(|entry| entry.is_alias)
                    && program.fn_map_entries.iter().all(|entry| entry.is_alias
                        || scalar_functions.contains(entry)
                        || numerators.contains(entry)),
                "owned aliases are retained; the colliding ordinary q() definition is shadowed by its input"
            );
            assert_eq!(rows.auxiliary_shared_outputs, 0);
            assert!(rows.input_columns.is_empty());
            assert!(rows.auxiliary_values.contains(&precise));
            for gather in &rows.auxiliary_gathers {
                assert_eq!(gather.len(), 11);
                assert!(gather.iter().all(Option::is_some));
            }
            for column in 0..11 {
                assert_ne!(
                    rows.auxiliary_gathers[0][column], rows.auxiliary_gathers[1][column],
                    "an unchanged alias handle still has row-dependent RHS bindings"
                );
            }
            let mut stack = EvaluatorStack {
                explicit_orientation_sum_only: false,
                parametric_rows: Some(rows),
                production_orientation_ids: Vec::new(),
                single_parametric: main,
                single_parametric_auxiliary: auxiliary,
                iterative: None,
                summed_function_map: None,
                summed: None,
            };
            check::<f64>(
                &mut stack,
                &mut baseline,
                &original_rows,
                input_count,
                x_index,
                &precise,
            );
            check::<f128>(
                &mut stack,
                &mut baseline,
                &original_rows,
                input_count,
                x_index,
                &precise,
            );
            check::<ArbPrec>(
                &mut stack,
                &mut baseline,
                &original_rows,
                input_count,
                x_index,
                &precise,
            );
            for optimization in [
                CompilationOptimizationLevel::O0,
                CompilationOptimizationLevel::O2,
            ] {
                stack
                    .for_each_generic_evaluator_mut(|evaluator| {
                        evaluator.activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: optimization,
                            ..Default::default()
                        })
                    })
                    .unwrap();
                check::<f64>(
                    &mut stack,
                    &mut baseline,
                    &original_rows,
                    input_count,
                    x_index,
                    &precise,
                );
            }
            let options = crate::settings::global::CompilationOptionsSnapshot {
                optimization_level: CompilationOptimizationLevel::O0,
                compiler: "g++".into(),
                ..Default::default()
            };
            for (backend, mode) in [
                ("cpp", FrozenCompilationMode::Cpp(options.clone())),
                ("asm", FrozenCompilationMode::Assembly(options)),
            ] {
                let name = format!(
                    "row_scalars_{backend}_{}",
                    shape.as_ref().map_or(1, Vec::len)
                );
                stack.compile(&name, &directory, &mode).unwrap();
                check::<f64>(
                    &mut stack,
                    &mut baseline,
                    &original_rows,
                    input_count,
                    x_index,
                    &precise,
                );
            }
        }
        for case in [
            "captured_row",
            "unknown_call",
            "cycle",
            "caller_alias",
            "generated_dependency",
            "row_tag",
            "alias_tag",
        ] {
            let opaque = symbol!("evaluator_test::row_scalar_opaque");
            let generated = symbol!("evaluator_test::row_scalar_generated");
            let tags = match case {
                "row_tag" => vec![function!(opaque, q.clone())],
                "alias_tag" => vec![alias.clone()],
                _ => Vec::new(),
            };
            let mut definitions = vec![FnMapEntry {
                lhs: one.call_args(tags.iter().chain(std::iter::once(&formal))),
                rhs: match case {
                    "captured_row" => &formal * &q,
                    "unknown_call" => function!(opaque, formal.clone()),
                    "cycle" => function!(one, formal.clone()),
                    "caller_alias" => &formal * &alias,
                    "generated_dependency" => function!(generated, formal.clone()),
                    "row_tag" | "alias_tag" => &formal * &x,
                    _ => unreachable!(),
                },
                tags: tags.clone(),
                args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            }];
            if case == "generated_dependency" {
                definitions.push(FnMapEntry {
                    lhs: function!(generated, formal.clone()),
                    rhs: &formal * &x,
                    tags: Vec::new(),
                    args: vec![Indeterminate::try_from(formal.clone()).unwrap()],
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                });
            }
            let call = one.call_args(tags.iter().chain(std::iter::once(&q)));
            let mut atom = AliasedAtom::from(&call * &x + x.pow(2));
            if matches!(case, "caller_alias" | "alias_tag") {
                atom.register_alias(alias.clone(), y.clone());
            }
            let mut rows = ParametricResidueRows {
                parameters: vec![q.clone()],
                input_columns: vec![0],
                rows: vec![vec![Rational::zero()], vec![Rational::one()]],
                ..Default::default()
            };
            let original = atom.get_root().clone();
            let extracted = rows.extract_physical_scalars(
                std::slice::from_mut(&mut atom),
                &[x.clone(), y.clone()],
                &definitions,
                1,
                symbol!("evaluator_test::row_scalar_negative_port"),
            );
            if case == "alias_tag" {
                // A definition key containing a caller alias makes the whole
                // owner ambiguous before individual calls are classified.
                assert!(extracted.is_none());
                assert_eq!(atom.get_root(), &original);
                continue;
            }
            let (_, roots, retained) =
                extracted.expect("the independent physical scalar still moves");
            assert_eq!(roots.len(), 1, "{case}: ordinary call must remain opaque");
            assert_eq!(roots[0].get_root(), &x.pow(2), "{case}");
            assert!(
                retained.is_empty(),
                "{case}: unsafe ordinary body is not moved"
            );
            assert!(rows.auxiliary_values.is_empty(), "{case}");
            assert!(atom.get_root().as_view().contains_symbol(one), "{case}");
        }
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn generated_numerator_tables_preserve_dynamic_domains_and_fallbacks() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::batched_domain_x"));
        let q = Atom::var(symbol!("evaluator_test::batched_domain_q"));
        let r = Atom::var(symbol!("evaluator_test::batched_domain_r"));
        let dispatch = Atom::var(symbol!("evaluator_test::batched_domain_dispatch"));
        let a = Atom::var(symbol!("evaluator_test::batched_domain_formal"));
        let n = symbol!("evaluator_test::batched_domain_numerator");
        let opaque = symbol!("evaluator_test::batched_domain_opaque");
        let metadata = function!(symbol!("evaluator_test::batched_domain_metadata"), 7);
        let capture = function!(symbol!("evaluator_test::batched_domain_capture"));
        let alias = function!(symbol!("evaluator_test::batched_domain_alias"));
        let terminal = Symbol::IF.call_args([x.clone(), Atom::one(), Atom::Zero]);
        for case in [
            "zero_factor",
            "near_one",
            "compound",
            "captured",
            "opaque_tagged_capture",
            "recorded_tagged_capture",
            "physical_guard",
            "shared_physical_guard",
            "distinct_calls",
            "alias_call",
            "nested_guard",
            "physical_if_terminal",
            "alias_if_terminal",
            "row_if_terminal",
            "row_if_boolean_terminal",
            "scalar_ratio",
            "scalar_repeated_column",
            "scalar_captured",
            "scalar_physical_guard",
            "alias_selector",
            "wrapped_inactive_scalar",
            "compound_distinct_guards",
            "never",
            "replacement",
        ] {
            let terminal = if case == "row_if_boolean_terminal" {
                Symbol::IF.call_args([q.clone(), Atom::one(), Atom::Zero])
            } else {
                terminal.clone()
            };
            let mut builder = ParamBuilder::new_empty();
            builder.pairs.additional_params = if case == "physical_if_terminal" {
                [x.clone(), terminal.clone()].into_iter().collect()
            } else {
                [x.clone()].into_iter().collect()
            };
            let input_count = builder.pairs.update_ranges();
            let x_index = builder.pairs.additional_params.value_range.start;
            if matches!(case, "captured" | "scalar_captured") {
                builder
                    .add_function_entry(FnMapEntry {
                        lhs: capture.clone(),
                        rhs: r.clone(),
                        tags: Vec::new(),
                        args: Vec::new(),
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
            }
            if case == "opaque_tagged_capture" {
                builder
                    .add_function_entry(FnMapEntry {
                        lhs: function!(opaque, 0),
                        rhs: x.pow(2),
                        tags: vec![Atom::num(0)],
                        args: Vec::new(),
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
                // This supported direct registration is intentionally absent
                // from reps. Sharing a head with a recorded overload does not
                // certify its signature or its captured global row parameter.
                FnMapEntry {
                    lhs: function!(opaque, 1),
                    rhs: q.clone(),
                    tags: vec![Atom::num(1)],
                    args: Vec::new(),
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                }
                .register(&mut builder.fn_map)
                .unwrap();
            }
            if case == "recorded_tagged_capture" {
                builder
                    .add_function_entry(FnMapEntry {
                        lhs: function!(opaque, metadata.clone()),
                        rhs: r.clone(),
                        tags: vec![metadata.clone()],
                        args: Vec::new(),
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
            }
            if case == "wrapped_inactive_scalar" {
                builder
                    .add_function_entry(FnMapEntry {
                        lhs: function!(GS.tree_denom_wrapper, a.clone()),
                        rhs: a.clone(),
                        tags: Vec::new(),
                        args: vec![Indeterminate::try_from(a.clone()).unwrap()],
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
            }
            let physical_functions = builder.reps.clone();
            let body = match case {
                "zero_factor" => &a / &x,
                "nested_guard" => &x / &a,
                "near_one" => Symbol::IF.call_args([&a - 1, x.pow(-1), &x + 7]),
                "captured" => &a * &x + &capture,
                "physical_guard" | "distinct_calls" => x.pow(-1),
                _ => &a * x.pow(2),
            };
            let definition = FnMapEntry {
                lhs: function!(n, a.clone()),
                rhs: body,
                tags: Vec::new(),
                args: vec![Indeterminate::try_from(a.clone()).unwrap()],
                inlining: if case == "never" {
                    InliningPolicy::Never
                } else {
                    InliningPolicy::Always
                },
                is_alias: false,
            };
            builder.add_function_entry(definition.clone()).unwrap();
            let mut atom = AliasedAtom::from(match case {
                "compound" => function!(n, &q - &r),
                "opaque_tagged_capture" => function!(opaque, 1) + x.pow(2),
                "recorded_tagged_capture" => function!(opaque, metadata.clone()) + x.pow(2),
                "physical_guard" => {
                    Symbol::IF.call_args([x.clone(), function!(n, q.clone()), Atom::Zero])
                }
                "shared_physical_guard" => {
                    function!(n, q.clone())
                        + Symbol::IF.call_args([x.clone(), function!(n, q.clone()), Atom::Zero])
                }
                "nested_guard" => Symbol::IF.call_args([
                    Symbol::IF.call_args([r.clone(), q.clone(), Atom::Zero]),
                    function!(n, q.clone()),
                    Atom::Zero,
                ]),
                "physical_if_terminal"
                | "alias_if_terminal"
                | "row_if_terminal"
                | "row_if_boolean_terminal" => {
                    Symbol::IF.call_args([q.clone(), terminal.clone(), Atom::Zero])
                }
                "compound_distinct_guards" => {
                    Symbol::IF.call_args([q.clone(), x.pow(-1), Atom::Zero])
                        - Symbol::IF.call_args([r.clone(), x.pow(-1), Atom::Zero])
                        + 1
                }
                "distinct_calls" => alias.clone(),
                "scalar_ratio" | "scalar_repeated_column" => alias.clone(),
                "scalar_captured" => &capture + x.pow(2),
                "scalar_physical_guard" => {
                    Symbol::IF.call_args([x.clone(), &q / &x, Atom::Zero]) + x.pow(2)
                }
                "alias_selector" | "wrapped_inactive_scalar" => {
                    Symbol::IF.call_args([dispatch.clone(), alias.clone(), Atom::Zero])
                }
                _ => function!(n, q.clone()),
            });
            if matches!(case, "scalar_ratio" | "scalar_repeated_column") {
                // Distinct columns must remain distinct native inputs even
                // when every sampled value agrees, including the 0/0 row.
                atom.register_alias(
                    alias.clone(),
                    &q / if case == "scalar_ratio" { &r } else { &q },
                );
            }
            if matches!(case, "alias_selector" | "wrapped_inactive_scalar") {
                // Follow only the pure Boolean selector behind the handle;
                // its singular arm is inactive even at x=0 and q=0.
                atom.register_alias(
                    alias.clone(),
                    Symbol::IF.call_args([
                        r.clone(),
                        if case == "wrapped_inactive_scalar" {
                            function!(GS.tree_denom_wrapper, &q / &x)
                        } else {
                            &q / &x
                        },
                        &q * x.pow(2) + 1,
                    ]),
                );
            }
            if case == "distinct_calls" {
                // Do not merge these ORIGINAL call sites, even when their
                // gathered numeric output indices coincide: Inf-Inf must not
                // be replaced by a normalized literal zero.
                atom.register_alias(
                    alias.clone(),
                    function!(n, q.clone()) - function!(n, r.clone()),
                );
            }
            if case == "alias_if_terminal" {
                // An IF-shaped alias is an opaque terminal of the outer
                // selector. Its body and lexical owner must remain intact.
                atom.register_alias(terminal.clone(), &x + 7);
            }
            if case == "alias_call" {
                atom.register_alias(function!(n, q.clone()), &q * &x + 4);
            }
            let rows = ParametricResidueRows {
                parameters: vec![
                    q.clone(),
                    r.clone(),
                    if matches!(case, "row_if_terminal" | "row_if_boolean_terminal") {
                        terminal.clone()
                    } else {
                        dispatch.clone()
                    },
                ],
                input_columns: vec![0, 1, 2],
                rows: if matches!(case, "row_if_terminal" | "row_if_boolean_terminal") {
                    // The exact IF-shaped row input is not its apparent syntax.
                    // Both rows activate the selector, whose values sum to five.
                    [[1, 0, 2], [1, 1, 3]]
                        .map(|row| row.map(Rational::from).to_vec())
                        .to_vec()
                } else if case == "near_one" {
                    vec![
                        vec![
                            Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64),
                            Rational::zero(),
                            Rational::zero(),
                        ],
                        vec![
                            Rational::new(18_014_398_509_481_985i64, 18_014_398_509_481_984i64),
                            Rational::one(),
                            Rational::one(),
                        ],
                    ]
                } else if matches!(case, "alias_selector" | "wrapped_inactive_scalar") {
                    [[0, 0, 1], [1, 0, 1]]
                        .map(|row| row.map(Rational::from).to_vec())
                        .to_vec()
                } else if case == "nested_guard" {
                    [[1, 1, 0], [1, 1, 1], [0, 1, 0], [0, 0, 1]]
                        .map(|row| row.map(Rational::from).to_vec())
                        .to_vec()
                } else {
                    [[0, 0, 0], [1, 1, 1]]
                        .map(|row| row.map(Rational::from).to_vec())
                        .to_vec()
                },
                ..Default::default()
            };
            let parameters = (&builder.pairs)
                .into_iter()
                .flat_map(|pair| pair.params.clone())
                .chain(rows.parameters.iter().cloned())
                .collect::<Vec<_>>();
            for shape in [None, Some(vec![vec![0], vec![1]])] {
                let width = shape.as_ref().map_or(1, Vec::len);
                let settings = EvaluatorSettings {
                    do_fn_map_replacements: case == "replacement",
                    store_atom: true,
                    ..Default::default()
                };
                let mut baseline = GenericEvaluator::new_from_raw_params(
                    [atom.clone()],
                    &parameters,
                    &builder.fn_map,
                    builder.reps.clone(),
                    settings.optimization_settings(),
                    shape.clone().map(|shape| {
                        let zeros = (input_count..parameters.len())
                            .flat_map(|p| (1..shape.len()).map(move |component| (p, component)))
                            .collect();
                        (shape, zeros)
                    }),
                    &settings,
                )
                .unwrap();
                let mut prepared_rows = rows.clone();
                let (main, auxiliary) = EvaluatorStack::new_single_parametric(
                    vec![atom.clone()],
                    &builder,
                    &physical_functions,
                    Some(&mut prepared_rows),
                    &shape,
                    &settings,
                )
                .unwrap();
                let numerator_batched = matches!(
                    case,
                    "zero_factor"
                        | "near_one"
                        | "distinct_calls"
                        | "nested_guard"
                        | "shared_physical_guard"
                );
                let scalar_batched = matches!(
                    case,
                    "physical_if_terminal"
                        | "alias_if_terminal"
                        | "opaque_tagged_capture"
                        | "recorded_tagged_capture"
                        | "scalar_captured"
                        | "scalar_physical_guard"
                );
                let row_scalar_batched = matches!(
                    case,
                    "scalar_ratio"
                        | "alias_call"
                        | "alias_selector"
                        | "wrapped_inactive_scalar"
                        | "row_if_terminal"
                        | "row_if_boolean_terminal"
                );
                assert_eq!(
                    auxiliary.is_some(),
                    numerator_batched
                        || scalar_batched
                        || row_scalar_batched
                        || case == "compound_distinct_guards",
                    "{case}, shape={shape:?}"
                );
                if scalar_batched {
                    assert!(prepared_rows.auxiliary_values.is_empty());
                    assert_eq!(auxiliary.as_ref().unwrap().0.len(), 1);
                    assert_eq!(auxiliary.as_ref().unwrap().1.exprs_len, 1);
                }
                if numerator_batched {
                    assert!(!prepared_rows.auxiliary_values.is_empty());
                    assert!(
                        auxiliary
                            .as_ref()
                            .unwrap()
                            .1
                            .fn_map_entries
                            .contains(&definition)
                    );
                }
                if row_scalar_batched {
                    assert!(!prepared_rows.auxiliary_values.is_empty());
                    assert!(prepared_rows.input_columns.is_empty());
                    assert_eq!(auxiliary.as_ref().unwrap().0.len(), 1);
                    assert_eq!(auxiliary.as_ref().unwrap().1.exprs_len, 2);
                }
                if case == "compound_distinct_guards" {
                    let (ports, program) = auxiliary.as_ref().unwrap();
                    assert_eq!(ports.len(), 2);
                    assert_ne!(ports[0], ports[1]);
                    assert_eq!(program.exprs_len, 1);
                    assert_eq!(program.exprs.as_ref().unwrap(), &vec![x.pow(-1)]);
                    assert_eq!(prepared_rows.auxiliary_shared_outputs, 0);
                    assert!(prepared_rows.auxiliary_values.is_empty());
                    assert!(prepared_rows.input_columns.is_empty());
                    assert_eq!(
                        prepared_rows.auxiliary_gathers,
                        vec![vec![None, None], vec![Some(0), Some(0)]]
                    );
                    let root = &main.exprs.as_ref().unwrap()[0];
                    assert!(
                        *root == &ports[0] - &ports[1] + 1 || *root == &ports[1] - &ports[0] + 1,
                        "distinct guards keep distinct caller ports even when their pole output is shared"
                    );
                }
                if case == "scalar_captured" {
                    assert!(prepared_rows.input_columns.contains(&1));
                }
                if matches!(case, "alias_selector" | "wrapped_inactive_scalar") {
                    assert!(
                        !auxiliary
                            .as_ref()
                            .unwrap()
                            .1
                            .exprs
                            .as_ref()
                            .unwrap()
                            .iter()
                            .any(|root| root.as_view().contains_symbol(Symbol::IF))
                    );
                }
                if case == "opaque_tagged_capture" {
                    assert_eq!(
                        prepared_rows.input_columns, rows.input_columns,
                        "an unrecorded overload can capture any original row input"
                    );
                    assert_eq!(
                        prepared_rows.auxiliary_shared_outputs, 1,
                        "the independent x² scalar forces signature pruning to run"
                    );
                }
                if case == "recorded_tagged_capture" {
                    assert_eq!(
                        prepared_rows.input_columns,
                        [1],
                        "metadata tags do not capture columns; the definition's free r does"
                    );
                }
                if case == "shared_physical_guard" {
                    assert!(
                        prepared_rows.input_columns.contains(&0),
                        "the hidden physical-guard call still needs its original coefficient"
                    );
                }
                if case == "nested_guard" {
                    assert!(
                        prepared_rows
                            .auxiliary_values
                            .iter()
                            .all(|value| value == &Rational::one())
                    );
                    assert!(
                        prepared_rows.auxiliary_gathers[2..]
                            .iter()
                            .all(|row| row.iter().any(Option::is_none))
                    );
                }
                if case == "distinct_calls" {
                    assert_eq!(auxiliary.as_ref().unwrap().0.len(), 2);
                    for gather in &prepared_rows.auxiliary_gathers {
                        assert_eq!(gather.len(), 2);
                        assert!(gather[0].is_some());
                        assert_eq!(gather[0], gather[1]);
                    }
                }
                let mut stack = EvaluatorStack {
                    explicit_orientation_sum_only: false,
                    parametric_rows: Some(prepared_rows),
                    production_orientation_ids: Vec::new(),
                    single_parametric: main,
                    single_parametric_auxiliary: auxiliary,
                    iterative: None,
                    summed_function_map: None,
                    summed: None,
                };
                for compiled in [false, true] {
                    if compiled {
                        baseline
                            .activate_symjit(&CompilationOptionsSnapshot {
                                optimization_level: CompilationOptimizationLevel::O2,
                                ..Default::default()
                            })
                            .unwrap();
                        stack
                            .for_each_generic_evaluator_mut(|evaluator| {
                                evaluator.activate_symjit(&CompilationOptionsSnapshot {
                                    optimization_level: CompilationOptimizationLevel::O2,
                                    ..Default::default()
                                })
                            })
                            .unwrap();
                    }
                    for point in [0., 0.5, -0.5] {
                        let mut input = vec![Complex::new_re(F(0.0)); input_count * width];
                        input[x_index * width] = Complex::new_re(F(point));
                        if case == "physical_if_terminal" {
                            // Deliberately independent of x: at x=0 the raw
                            // IF syntax would be zero, while this input is 7.
                            input[(x_index + 1) * width] = Complex::new_re(F(7.0));
                            if width > 1 {
                                input[(x_index + 1) * width + 1] = Complex::new_re(F(2.0));
                            }
                        }
                        if width > 1 {
                            input[x_index * width + 1] = Complex::new_re(F(1.0));
                        }
                        let old = rows.evaluate(
                            &mut baseline,
                            &input,
                            &[],
                            &mut EvaluationMetaData::new_empty(),
                        );
                        let new = stack.evaluate_single_parametric(
                            &input,
                            &mut EvaluationMetaData::new_empty(),
                        );
                        assert_eq!(old.len(), new.len());
                        if case == "wrapped_inactive_scalar" {
                            let values = match &new[0] {
                                DualOrNot::Dual(value) => value.values.as_slice(),
                                DualOrNot::NonDual(value) => std::slice::from_ref(value),
                            };
                            assert_eq!(
                                values,
                                [2. + point * point, 2. * point]
                                    .into_iter()
                                    .take(width)
                                    .map(|v| Complex::new_re(F(v)))
                                    .collect::<Vec<_>>(),
                                "the inactive tree_denoms(q/x) argument is never evaluated, including 0/0 at x=0"
                            );
                        }
                        if case == "compound_distinct_guards" && point == 0. {
                            let values = match &old[0] {
                                DualOrNot::Dual(value) => value.values.as_slice(),
                                DualOrNot::NonDual(value) => std::slice::from_ref(value),
                            };
                            assert!(
                                values.iter().any(|value| !value.re.0.is_finite() || !value.im.0.is_finite()),
                                "the original distinct guards must actually exercise active Inf-Inf"
                            );
                        }
                        if matches!(case, "row_if_terminal" | "row_if_boolean_terminal") {
                            let values = match &new[0] {
                                DualOrNot::Dual(value) => value.values.as_slice(),
                                DualOrNot::NonDual(value) => std::slice::from_ref(value),
                            };
                            assert_eq!(
                                values,
                                [5., 0.]
                                    .into_iter()
                                    .take(width)
                                    .map(|value| Complex::new_re(F(value)))
                                    .collect::<Vec<_>>()
                            );
                        }
                        if matches!(case, "opaque_tagged_capture" | "recorded_tagged_capture") {
                            let values = match &new[0] {
                                DualOrNot::Dual(value) => value.values.as_slice(),
                                DualOrNot::NonDual(value) => std::slice::from_ref(value),
                            };
                            assert_eq!(
                                values,
                                [1. + 2. * point * point, 4. * point]
                                    .into_iter()
                                    .take(width)
                                    .map(|v| Complex::new_re(F(v)))
                                    .collect::<Vec<_>>()
                            );
                        }
                        for (old, new) in old.iter().zip(&new) {
                            let old = match old {
                                DualOrNot::Dual(v) => v.values.as_slice(),
                                DualOrNot::NonDual(v) => std::slice::from_ref(v),
                            };
                            let new = match new {
                                DualOrNot::Dual(v) => v.values.as_slice(),
                                DualOrNot::NonDual(v) => std::slice::from_ref(v),
                            };
                            assert_eq!(old.len(), new.len());
                            for (old, new) in old.iter().zip(new) {
                                for (old, new) in [(old.re.0, new.re.0), (old.im.0, new.im.0)] {
                                    assert_eq!(
                                        old.is_finite(),
                                        new.is_finite(),
                                        "{case}, x={point}, shape={shape:?}, compiled={compiled}"
                                    );
                                    assert_eq!(old.is_nan(), new.is_nan(), "{case}");
                                    if old.is_infinite() {
                                        assert_eq!(old, new, "{case}");
                                    }
                                    if old.is_finite() {
                                        assert!(
                                            (old - new).abs()
                                                <= (old.abs() + 1.) * f64::EPSILON * 64.,
                                            "{case}: {old} != {new}"
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn physical_functions_preserve_complete_rows_and_native_higher_jets() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::physical_function_x"));
        let q = Atom::var(symbol!("evaluator_test::physical_function_q"));
        let energy = function!(symbol!("evaluator_test::physical_energy"));
        let family = symbol!("evaluator_test::physical_family");
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        let input_count = builder.pairs.update_ranges();
        let x_index = builder.pairs.additional_params.value_range.start;
        let energy_definition = FnMapEntry {
            lhs: energy.clone(),
            rhs: (x.pow(2) + 1).pow((1, 2)),
            tags: Vec::new(),
            args: Vec::new(),
            inlining: InliningPolicy::Always,
            is_alias: false,
        };
        builder
            .add_function_entry(energy_definition.clone())
            .unwrap();
        let terms = [
            (q.clone(), energy.clone()),
            (Atom::one(), (&energy + &q).pow(-1)),
        ]
        .into_iter()
        .enumerate()
        .map(|(tag, (rhs, carrier))| ParametricIntegrandTerm {
            numerator: Arc::new(FnMapEntry {
                lhs: function!(family, tag, q.clone()),
                rhs,
                tags: vec![Atom::num(tag)],
                args: vec![Indeterminate::try_from(q.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            }),
            rows: (1..=3)
                .map(|coefficient| crate::uv::ParametricIntegrandRow {
                    coefficients: vec![Rational::from(coefficient)],
                    carrier: carrier.clone(),
                    source_keys: Vec::new(),
                })
                .collect(),
        })
        .collect::<Vec<_>>();
        fn check<T: FloatLike>(stack: &mut EvaluatorStack, input_count: usize, x_index: usize) {
            let zero = F(T::new_zero());
            let orientations =
                TiVec::<OrientationID, _>::from_iter([EdgeVec::from_iter(std::iter::empty::<
                    Orientation,
                >())]);
            let filter = SubSet::full(1);
            for method in [
                EvaluatorMethod::SingleParametric,
                EvaluatorMethod::SummedFunctionMap,
                EvaluatorMethod::Summed,
                EvaluatorMethod::Iterative,
            ] {
                for (value, first, second) in
                    [(0, 1, 0), (1, -2, 1), (-2, 3, -1), (0, 0, 2), (1, -2, 1)]
                {
                    let x = zero.from_i64(value);
                    let dx = zero.from_i64(first);
                    let d2x = zero.from_i64(second);
                    let energy = (x.powi(2) + zero.one()).sqrt();
                    let de = &x * &dx / &energy;
                    let d2e = (dx.powi(2) + zero.from_usize(2) * &x * &d2x)
                        / (zero.from_usize(2) * &energy)
                        - x.powi(2) * dx.powi(2) / (zero.from_usize(2) * energy.powi(3));
                    let mut expected_value = zero.from_usize(6) * &energy;
                    let mut slope = zero.from_usize(6);
                    let mut curvature = zero.clone();
                    for coefficient in 1..=3 {
                        let denominator = &energy + zero.from_usize(coefficient);
                        expected_value += denominator.powi(-1);
                        slope -= denominator.powi(-2);
                        curvature += denominator.powi(-3);
                    }
                    let expected = [
                        expected_value,
                        &slope * &de,
                        slope * d2e + curvature * de.powi(2),
                    ];
                    let mut values = vec![Complex::new_re(zero.clone()); input_count * 3];
                    values[x_index * 3] = Complex::new_re(x);
                    values[x_index * 3 + 1] = Complex::new_re(dx);
                    values[x_index * 3 + 2] = Complex::new_re(d2x);
                    let input = InputParams {
                        values: SliceMut::Owned(values),
                        residue_map_id_start: 0,
                        orientations_start: 0,
                        multiplicative_offset: 3,
                    };
                    let mut runtime = RuntimeSettings::default();
                    runtime.general.evaluator_method = method.clone();
                    let output = stack
                        .evaluate(
                            input,
                            SingleOrAllOrientations::All {
                                all: &orientations,
                                filter: &filter,
                            },
                            &runtime,
                            &mut EvaluationMetaData::new_empty(),
                        )
                        .unwrap();
                    let [DualOrNot::Dual(output)] = output.as_slice() else {
                        panic!("expected one complete higher-jet output")
                    };
                    for (actual, expected) in output.values.iter().zip(expected) {
                        assert!(
                            (actual.re.clone() - &expected).abs()
                                <= (expected.abs() + zero.one())
                                    * zero.epsilon()
                                    * zero.from_usize(128),
                            "{method:?}, x={value}, seed=({first},{second}): {actual:?} != {expected:?}"
                        );
                        assert_eq!(actual.im, zero);
                    }
                }
            }
        }
        for inline_numerator_functions in [false, true] {
            for store_atom in [false, true] {
                let settings = EvaluatorSettings {
                    inline_numerator_functions,
                    store_atom,
                    summed: true,
                    summed_function_map: true,
                    iterative_orientation_optimization: true,
                    ..Default::default()
                };
                let (mut stack, _) = EvaluatorStack::from_integrand_with_timings(
                    &Atom::Zero,
                    &builder,
                    &[],
                    &terms,
                    None,
                    Some(vec![vec![0], vec![1], vec![2]]),
                    &settings,
                )
                .unwrap();
                let (ports, auxiliary) = stack
                    .single_parametric_auxiliary
                    .as_ref()
                    .expect("energy is reached on all three rows");
                // The guarded energy, numerator, and complete reciprocal
                // scalar keep distinct ports. The reciprocal now binds q as
                // a native input inside each gathered leaf.
                assert_eq!(ports.len(), 3);
                assert_eq!(auxiliary.fn_map_entries.len(), 2);
                assert_eq!(auxiliary.fn_map_entries[0], energy_definition);
                let numerator = &auxiliary.fn_map_entries[1];
                assert_eq!(numerator.args.len(), 1);
                assert_eq!(numerator.rhs, Atom::from(numerator.args[0].clone()));
                assert_eq!(auxiliary.exprs_len, 7);
                assert_eq!(
                    stack
                        .parametric_rows
                        .as_ref()
                        .unwrap()
                        .auxiliary_shared_outputs,
                    0
                );
                assert_eq!(stack.generic_evaluator_count(), 5);
                let rows = stack.parametric_rows.as_ref().unwrap();
                assert_eq!(rows.rows.len(), 3);
                assert!(rows.auxiliary_gathers.iter().all(|row| row.len() == 3));
                let mut gathered_outputs = (0..3)
                    .map(|port| {
                        rows.auxiliary_gathers
                            .iter()
                            .map(|row| {
                                row[port].expect("all three source ports are active on every row")
                            })
                            .collect::<HashSet<_>>()
                            .len()
                    })
                    .collect::<Vec<_>>();
                gathered_outputs.sort_unstable();
                assert_eq!(gathered_outputs, [1, 3, 3]);
                assert_eq!(
                    stack.single_parametric.zero_components,
                    (input_count + ports.len()
                        ..input_count + ports.len() + rows.input_columns.len())
                        .flat_map(|p| (1..3).map(move |component| (p, component)))
                        .collect::<Vec<_>>()
                );
                check::<f64>(&mut stack, input_count, x_index);
                check::<f128>(&mut stack, input_count, x_index);
                check::<ArbPrec>(&mut stack, input_count, x_index);
                let encoded = bincode::encode_to_vec(&stack, bincode::config::standard()).unwrap();
                let mut state = Vec::new();
                State::export(&mut state).unwrap();
                let state_map = State::import(&mut Cursor::new(state), None).unwrap();
                let model = Model::default();
                (stack, _) = bincode::decode_from_slice_with_context(
                    &encoded,
                    bincode::config::standard(),
                    GammaLoopContextContainer {
                        state_map: &state_map,
                        model: &model,
                    },
                )
                .unwrap();
                for optimization in [
                    CompilationOptimizationLevel::O0,
                    CompilationOptimizationLevel::O2,
                ] {
                    stack
                        .for_each_generic_evaluator_mut(|evaluator| {
                            evaluator.activate_symjit(&CompilationOptionsSnapshot {
                                optimization_level: optimization,
                                ..Default::default()
                            })
                        })
                        .unwrap();
                    assert_eq!(
                        stack
                            .single_parametric_auxiliary
                            .as_ref()
                            .unwrap()
                            .1
                            .active_f64_backend(),
                        ActiveF64Backend::Symjit
                    );
                    check::<f64>(&mut stack, input_count, x_index);
                }
            }
        }
        // Literal normalization removed sqrt(x²) before evaluator construction.
        // It must remain finite at x=0; extraction never recovers that root.
        let literal = x.pow(2).pow((1, 2)).pow(2);
        assert_eq!(literal, x.pow(2));
        let rows = ParametricResidueRows {
            input_columns: vec![0],
            parameters: vec![q],
            rows: vec![vec![Rational::one()], vec![Rational::from(2)]],
            ..Default::default()
        };
        for do_fn_map_replacements in [false, true] {
            let settings = EvaluatorSettings {
                do_fn_map_replacements,
                ..Default::default()
            };
            let mut prepared_rows = rows.clone();
            let (evaluator, auxiliary) = EvaluatorStack::new_single_parametric(
                vec![literal.clone().into()],
                &builder,
                &builder.reps,
                Some(&mut prepared_rows),
                &Some(vec![vec![0], vec![1], vec![2]]),
                &settings,
            )
            .unwrap();
            let mut stack = EvaluatorStack {
                explicit_orientation_sum_only: false,
                parametric_rows: Some(prepared_rows),
                production_orientation_ids: Vec::new(),
                single_parametric: evaluator,
                single_parametric_auxiliary: auxiliary,
                iterative: None,
                summed_function_map: None,
                summed: None,
            };
            let mut values = vec![Complex::new_re(F(0.0)); input_count * 3];
            values[x_index * 3 + 1] = Complex::new_re(F(1.0));
            let output =
                stack.evaluate_single_parametric(&values, &mut EvaluationMetaData::new_empty());
            let [DualOrNot::Dual(output)] = output.as_slice() else {
                panic!("expected literal jet")
            };
            assert_eq!(
                output.values,
                [0., 0., 2.].map(|value| Complex::new_re(F(value)))
            );
        }
    }

    #[test]
    fn physical_functions_preserve_opaque_numerators_and_lazy_fallbacks() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::physical_boundary_x"));
        let q = Atom::var(symbol!("evaluator_test::physical_boundary_q"));
        let energy = function!(symbol!("evaluator_test::physical_boundary_energy"));
        let pole = function!(symbol!("evaluator_test::physical_boundary_pole"));
        let numerator = function!(symbol!("evaluator_test::physical_boundary_numerator"));
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        let input_count = builder.pairs.update_ranges();
        let x_index = builder.pairs.additional_params.value_range.start;
        for (lhs, rhs) in [
            (energy.clone(), (x.pow(2) + 1).pow((1, 2))),
            (pole.clone(), x.pow(-1)),
        ] {
            builder
                .add_function_entry(FnMapEntry {
                    lhs,
                    rhs,
                    tags: Vec::new(),
                    args: Vec::new(),
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                })
                .unwrap();
        }
        let mut prepared_builder = builder.clone();
        prepared_builder
            .add_function_entry(FnMapEntry {
                lhs: numerator.clone(),
                rhs: x.pow(2) + 2,
                tags: Vec::new(),
                args: Vec::new(),
                inlining: InliningPolicy::Always,
                is_alias: false,
            })
            .unwrap();
        let rows = ParametricResidueRows {
            input_columns: vec![0],
            parameters: vec![q.clone()],
            rows: vec![vec![Rational::one()], vec![Rational::from(2)]],
            ..Default::default()
        };
        // A closed numerator can have the same physical-only dependency pattern
        // as an OSE. Its owner, not its function name or formal arity, excludes it.
        // The pole is active only under a physical predicate and must stay lazy.
        let atom =
            &energy + &numerator + &q + Symbol::IF.call_args([x.clone(), pole.clone(), Atom::Zero]);
        for dual in [false, true] {
            let shape = dual.then(|| vec![vec![0], vec![1], vec![2]]);
            let size = if dual { 3 } else { 1 };
            let mut prepared_rows = rows.clone();
            let (evaluator, auxiliary) = EvaluatorStack::new_single_parametric(
                vec![atom.clone().into()],
                &prepared_builder,
                &builder.reps,
                Some(&mut prepared_rows),
                &shape,
                &EvaluatorSettings::default(),
            )
            .unwrap();
            let (ports, program) = auxiliary.as_ref().unwrap();
            assert_eq!(ports.len(), 2);
            assert_eq!(program.exprs_len, 3);
            assert_eq!(prepared_rows.auxiliary_shared_outputs, 1);
            assert_eq!(prepared_rows.auxiliary_gathers.len(), 2);
            assert_ne!(
                prepared_rows.auxiliary_gathers[0],
                prepared_rows.auxiliary_gathers[1]
            );
            assert_eq!(program.fn_map_entries, vec![builder.reps[0].clone()]);
            assert!(
                evaluator
                    .fn_map_entries
                    .iter()
                    .any(|entry| entry.lhs == numerator)
            );
            prepared_rows.rows.truncate(1);
            prepared_rows.auxiliary_gathers.truncate(1);
            let mut stack = EvaluatorStack {
                explicit_orientation_sum_only: false,
                parametric_rows: Some(prepared_rows),
                production_orientation_ids: Vec::new(),
                single_parametric: evaluator,
                single_parametric_auxiliary: auxiliary,
                iterative: None,
                summed_function_map: None,
                summed: None,
            };
            for compiled in [false, true] {
                if compiled {
                    stack
                        .for_each_generic_evaluator_mut(|evaluator| {
                            evaluator.activate_symjit(&CompilationOptionsSnapshot {
                                optimization_level: CompilationOptimizationLevel::O2,
                                ..Default::default()
                            })
                        })
                        .unwrap();
                }
                for (point, expected) in [
                    (0., [4., 0., 1.5]),
                    (
                        1.,
                        [
                            2.0f64.sqrt() + 5.,
                            1. + 1. / 2.0f64.sqrt(),
                            2. + 1. / (4. * 2.0f64.sqrt()),
                        ],
                    ),
                ] {
                    let mut input = vec![Complex::new_re(F(0.0)); input_count * size];
                    input[x_index * size] = Complex::new_re(F(point));
                    if dual {
                        input[x_index * size + 1] = Complex::new_re(F(1.0));
                    }
                    let output = stack
                        .evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
                    let values = match &output[0] {
                        DualOrNot::Dual(value) => value.values.as_slice(),
                        DualOrNot::NonDual(value) => std::slice::from_ref(value),
                    };
                    for (actual, expected) in values.iter().zip(expected) {
                        assert!(
                            (actual.re.0 - expected).abs() < 1.0e-13,
                            "{actual:?} != {expected}; x={point}, dual={dual}, compiled={compiled}"
                        );
                        assert_eq!(actual.im, F(0.0));
                    }
                }
            }
        }
        // Explicit unfolding and retained-dual lowering have different symbolic
        // cancellation boundaries. They retain the existing unsplit evaluator.
        for retained in [false, true] {
            let mut fallback_builder = builder.clone();
            if retained {
                fallback_builder.reps[0].inlining = InliningPolicy::Never;
                fallback_builder.fn_map = FunctionMap::new();
                for entry in &fallback_builder.reps {
                    entry.register(&mut fallback_builder.fn_map).unwrap();
                }
            }
            let settings = EvaluatorSettings {
                do_fn_map_replacements: !retained,
                ..Default::default()
            };
            let mut prepared_rows = rows.clone();
            let (mut evaluator, auxiliary) = EvaluatorStack::new_single_parametric(
                vec![(&energy + &q).into()],
                &fallback_builder,
                &fallback_builder.reps,
                Some(&mut prepared_rows),
                &Some(vec![vec![0], vec![1], vec![2]]),
                &settings,
            )
            .unwrap();
            assert!(auxiliary.is_none());
            if retained {
                assert!(
                    evaluator
                        .fn_map_entries
                        .iter()
                        .any(|entry| entry.inlining == InliningPolicy::Never)
                );
            }
            let mut input = vec![Complex::new_re(F(0.0)); (input_count + 1) * 3];
            input[x_index * 3 + 1] = Complex::new_re(F(1.0));
            input[input_count * 3] = Complex::new_re(F(1.0));
            let output = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&input);
            let [DualOrNot::Dual(output)] = output.as_slice() else {
                panic!("expected fallback jet")
            };
            assert_eq!(
                output.values,
                [2., 0., 0.5].map(|value| Complex::new_re(F(value)))
            );
        }
    }

    #[test]
    fn physical_function_activity_requires_native_boolean_selectors() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::activity_physical_x"));
        let q = Atom::var(symbol!("evaluator_test::activity_row_q"));
        let dispatch = Atom::var(symbol!("evaluator_test::activity_dispatch"));
        let energy = function!(symbol!("evaluator_test::activity_energy"));
        let pole = function!(symbol!("evaluator_test::activity_pole"));
        let composite = Symbol::IF.call_args([q.clone(), Atom::one(), Atom::Zero]);
        for case in [
            "local_alias",
            "physical_input",
            "alias_inside_selector",
            "rounded_rational",
        ] {
            let mut builder = ParamBuilder::new_empty();
            builder.pairs.additional_params = if case == "physical_input" {
                [x.clone(), composite.clone()].into_iter().collect()
            } else {
                [x.clone()].into_iter().collect()
            };
            let input_count = builder.pairs.update_ranges();
            let parameters = (&builder.pairs)
                .into_iter()
                .flat_map(|pair| pair.params.clone())
                .collect::<Vec<_>>();
            let x_index = parameters
                .iter()
                .position(|parameter| parameter == &x)
                .unwrap();
            for (lhs, rhs) in [(energy.clone(), x.pow(2) + 1), (pole.clone(), x.pow(-1))] {
                builder
                    .add_function_entry(FnMapEntry {
                        lhs,
                        rhs,
                        tags: Vec::new(),
                        args: Vec::new(),
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
            }
            let wrapper = symbol!("evaluator_test::activity_selector_wrapper");
            let local = Atom::var(symbol!("evaluator_test::activity_selector_local"));
            let formal = Atom::var(symbol!("evaluator_test::activity_selector_formal"));
            let selector = if case == "alias_inside_selector" {
                function!(wrapper, local.clone())
            } else {
                q.clone()
            };
            if case == "alias_inside_selector" {
                builder
                    .add_function_entry(FnMapEntry {
                        lhs: function!(wrapper, formal.clone()),
                        rhs: formal.clone(),
                        tags: Vec::new(),
                        args: vec![Indeterminate::try_from(formal).unwrap()],
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
            }
            let rows = ParametricResidueRows {
                input_columns: vec![0, 1],
                parameters: vec![selector.clone(), dispatch.clone()],
                rows: if case == "rounded_rational" {
                    // Distinct exact coefficients both round to one in Double.
                    vec![
                        vec![
                            Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64),
                            Rational::zero(),
                        ],
                        vec![
                            Rational::new(18_014_398_509_481_985i64, 18_014_398_509_481_984i64),
                            Rational::one(),
                        ],
                    ]
                } else {
                    vec![
                        vec![Rational::one(), Rational::zero()],
                        vec![Rational::one(), Rational::one()],
                    ]
                },
                ..Default::default()
            };
            let condition = if case == "rounded_rational" {
                &q - 1
            } else if case == "alias_inside_selector" {
                selector
            } else {
                composite.clone()
            };
            let mut atom = AliasedAtom::from(
                &energy + Symbol::IF.call_args([condition, pole.clone(), Atom::Zero]),
            );
            if case == "local_alias" {
                atom.register_alias(composite.clone(), x.clone());
            }
            if case == "alias_inside_selector" {
                atom.register_alias(local, x.clone());
            }
            for dual in [false, true] {
                let shape = dual.then(|| vec![vec![0], vec![1]]);
                let width = if dual { 2 } else { 1 };
                let mut prepared_rows = rows.clone();
                let (main, auxiliary) = EvaluatorStack::new_single_parametric(
                    vec![atom.clone()],
                    &builder,
                    &builder.reps,
                    Some(&mut prepared_rows),
                    &shape,
                    &EvaluatorSettings::default(),
                )
                .unwrap();
                // This owner shadows a subatom of a row key. The scalar
                // extractor conservatively keeps the entire owner unsplit.
                if case == "alias_inside_selector" {
                    assert!(auxiliary.is_none());
                } else {
                    assert_eq!(auxiliary.as_ref().unwrap().0.len(), 1, "{case}");
                    assert_eq!(
                        auxiliary.as_ref().unwrap().1.fn_map_entries,
                        vec![builder.reps[0].clone()]
                    );
                }
                let mut stack = EvaluatorStack {
                    explicit_orientation_sum_only: false,
                    parametric_rows: Some(prepared_rows),
                    production_orientation_ids: Vec::new(),
                    single_parametric: main,
                    single_parametric_auxiliary: auxiliary,
                    iterative: None,
                    summed_function_map: None,
                    summed: None,
                };
                for compiled in [false, true] {
                    if compiled {
                        stack
                            .for_each_generic_evaluator_mut(|evaluator| {
                                evaluator.activate_symjit(&CompilationOptionsSnapshot {
                                    optimization_level: CompilationOptimizationLevel::O2,
                                    ..Default::default()
                                })
                            })
                            .unwrap();
                    }
                    let mut input = vec![Complex::new_re(F(0.0)); input_count * width];
                    if dual {
                        input[x_index * width + 1] = Complex::new_re(F(1.0));
                    }
                    let result = stack
                        .evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
                    let values = match &result[0] {
                        DualOrNot::Dual(value) => value.values.clone(),
                        DualOrNot::NonDual(value) => vec![*value],
                    };
                    assert_eq!(
                        values,
                        [2., 0.]
                            .into_iter()
                            .take(width)
                            .map(|value| Complex::new_re(F(value)))
                            .collect::<Vec<_>>(),
                        "{case}, dual={dual}, compiled={compiled}"
                    );
                }
            }
        }
    }

    #[test]
    fn physical_aliases_preserve_native_domains_and_mixed_jets() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::physical_alias_x"));
        let y = Atom::var(symbol!("evaluator_test::physical_alias_y"));
        let q = Atom::var(symbol!("evaluator_test::physical_alias_q"));
        let a = function!(symbol!("evaluator_test::physical_alias"), 0);
        let b = function!(symbol!("evaluator_test::physical_alias"), 1);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone(), y.clone()].into_iter().collect();
        let input_count = builder.pairs.update_ranges();
        let x_index = builder.pairs.additional_params.value_range.start;
        let y_index = x_index + 1;
        let parameters = (&builder.pairs)
            .into_iter()
            .flat_map(|pair| pair.params.clone())
            .chain([q.clone()])
            .collect::<Vec<_>>();
        let rows = ParametricResidueRows {
            input_columns: vec![0],
            parameters: vec![q.clone()],
            rows: vec![vec![Rational::zero()], vec![Rational::one()]],
            ..Default::default()
        };
        // These are caller-local aliases, not literal algebraic subexpressions.
        // Compare against the existing unsplit GenericEvaluator boundary, with
        // the same exact row table and zero derivative seeds for its coefficients.
        let cases = [
            ("reciprocal", vec![(a.clone(), x.pow(-1))], a.pow(-1), 1),
            (
                "root_squared",
                vec![(a.clone(), x.pow(2).pow((1, 2)))],
                a.pow(2),
                1,
            ),
            ("zero_row_pole", vec![(a.clone(), x.pow(-1))], &q * &a, 1),
            (
                "cancel_distinct",
                vec![(a.clone(), x.pow(-1)), (b.clone(), x.pow(-1))],
                &a - &b + 1,
                1,
            ),
            (
                "cancel_nested",
                vec![(a.clone(), x.pow(-1)), (b.clone(), &a + 1)],
                &b - &a,
                1,
            ),
            (
                "parent_only",
                vec![
                    (a.clone(), x.pow(2) + y.pow(2) + 1),
                    (b.clone(), (&a + 1).pow(-1) * &x * &y),
                ],
                &q * &b,
                1,
            ),
            (
                "parent_and_child",
                vec![
                    (a.clone(), x.pow(2) + y.pow(2) + 1),
                    (b.clone(), (&a + 1).pow(-1) * &x * &y),
                ],
                &q * &b + &a,
                1,
            ),
            (
                "compound_mul",
                vec![(
                    a.clone(),
                    &x * &y * Symbol::IF.call_args([q.clone(), x.pow(2) + 1, y.pow(2) + 2]),
                )],
                a.clone(),
                1,
            ),
            (
                "compound_add",
                vec![(
                    a.clone(),
                    Symbol::IF.call_args([q.clone(), x.pow(2) * &y, &x * y.pow(2)]) + &x * &y,
                )],
                a.clone(),
                2,
            ),
            (
                "compound_inactive_pole",
                vec![(
                    a.clone(),
                    &x * &y
                        * Symbol::IF.call_args([
                            q.clone(),
                            Symbol::IF.call_args([q.clone(), x.pow(2) + 1, x.pow(-1)]),
                            &y + 1,
                        ]),
                )],
                a.clone(),
                1,
            ),
            // Keep both native operations: the zero branch must not turn
            // physical Inf * 0 into a symbolically normalized zero.
            (
                "compound_zero_inf",
                vec![(a.clone(), x.pow(-1))],
                &a * Symbol::IF.call_args([q.clone(), Atom::one(), Atom::Zero]),
                2,
            ),
            (
                "compound_overflow",
                vec![(a.clone(), x.pow(2)), (b.clone(), x.pow(2))],
                Symbol::IF.call_args([q.clone(), a.clone(), Atom::Zero])
                    - Symbol::IF.call_args([q.clone(), b.clone(), Atom::Zero]),
                2,
            ),
            // An alias eliminated by normalization has no reachable use.
            ("dead_alias", vec![(a.clone(), x.pow(-1))], &a / &a, 0),
            ("input_copy", vec![(a.clone(), x.clone())], &a + &q, 1),
            (
                "physical_guard",
                vec![(a.clone(), x.pow(-1))],
                Symbol::IF.call_args([x.clone(), a.clone(), Atom::Zero]),
                0,
            ),
        ];
        #[allow(clippy::too_many_arguments)]
        fn compare<T: FloatLike>(
            baseline: &mut GenericEvaluator,
            stack: &mut EvaluatorStack,
            rows: &ParametricResidueRows,
            input_count: usize,
            x_index: usize,
            y_index: usize,
            shape: &Option<Vec<Vec<usize>>>,
            case: &str,
        ) {
            let zero = F(T::new_zero());
            let width = shape.as_ref().map_or(1, Vec::len);
            let mixed = shape.as_ref().is_some_and(|shape| shape[0].len() == 2);
            for point in [0, 1, -1, 2]
                .into_iter()
                .filter(|point| *point != 2 || case == "compound_overflow")
            {
                let mut input = vec![Complex::new_re(zero.clone()); input_count * width];
                input[x_index * width] = Complex::new_re(if point == 2 {
                    zero.from_usize(10).powi(200) // Finite binary64 input; its square overflows.
                } else {
                    zero.from_i64(point) / zero.from_usize(2)
                });
                input[y_index * width] = Complex::new_re(zero.from_usize(2));
                if width > 1 {
                    input[x_index * width + 1] = Complex::new_re(zero.one());
                }
                if mixed {
                    input[y_index * width + 2] = Complex::new_re(zero.one());
                }
                let old =
                    rows.evaluate(baseline, &input, &[], &mut EvaluationMetaData::new_empty());
                let new =
                    stack.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
                assert_eq!(old.len(), new.len());
                if matches!(
                    case,
                    "compound_mul" | "compound_add" | "compound_inactive_pole"
                ) && (width == 1 || mixed)
                {
                    let x = input[x_index * width].re.clone();
                    let y = input[y_index * width].re.clone();
                    let two = zero.from_usize(2);
                    let three = zero.from_usize(3);
                    // Independent complete two-row sums and mixed derivatives;
                    // all chosen probes and coefficients are exactly dyadic.
                    let expected = match case {
                        "compound_mul" => [
                            &x * &y * (x.powi(2) + y.powi(2) + &three),
                            &three * x.powi(2) * &y + y.powi(3) + &three * &y,
                            x.powi(3) + &three * &x * y.powi(2) + &three * &x,
                            &three * (x.powi(2) + y.powi(2) + zero.one()),
                        ],
                        "compound_add" => [
                            &x * &y * (&x + &y + &two),
                            &two * &x * &y + y.powi(2) + &two * &y,
                            x.powi(2) + &two * &x * &y + &two * &x,
                            &two * (&x + &y + zero.one()),
                        ],
                        _ => [
                            &x * &y * (x.powi(2) + &y + &two),
                            &three * x.powi(2) * &y + y.powi(2) + &two * &y,
                            x.powi(3) + &two * &x * &y + &two * &x,
                            &three * x.powi(2) + &two * &y + &two,
                        ],
                    };
                    let values = match &new[0] {
                        DualOrNot::Dual(value) => value.values.as_slice(),
                        DualOrNot::NonDual(value) => std::slice::from_ref(value),
                    };
                    for (actual, expected) in values.iter().zip(expected) {
                        assert_eq!(*actual, Complex::new_re(expected), "{case}");
                    }
                }
                for (old, new) in old.iter().zip(&new) {
                    let values = |value: &DualOrNot<Complex<F<T>>>| match value {
                        DualOrNot::Dual(value) => value.values.clone(),
                        DualOrNot::NonDual(value) => vec![value.clone()],
                    };
                    for (old, new) in values(old).into_iter().zip(values(new)) {
                        for (old, new) in [(old.re, new.re), (old.im, new.im)] {
                            let finite = old.to_f64().is_finite();
                            assert_eq!(
                                finite,
                                new.to_f64().is_finite(),
                                "domain changed at x={point}/2, shape={shape:?}"
                            );
                            if case.starts_with("compound_") {
                                assert_eq!(old.to_f64().is_nan(), new.to_f64().is_nan(), "{case}");
                                if old.to_f64().is_infinite() {
                                    assert_eq!(old.to_f64(), new.to_f64(), "{case}");
                                }
                            }
                            if finite {
                                assert!(
                                    (&old - &new).abs()
                                        <= (old.abs() + zero.one())
                                            * zero.epsilon()
                                            * zero.from_usize(128),
                                    "alias jet differs at x={point}/2, shape={shape:?}: {old:?} != {new:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
        for (name, aliases, root, expected_ports) in cases {
            let mut atom = AliasedAtom::from(root);
            for (key, body) in aliases {
                atom.register_alias(key, body);
            }
            for shape in [
                None,
                Some(vec![vec![0], vec![1]]),
                Some(vec![vec![0], vec![1], vec![2]]),
                Some(vec![vec![0], vec![1], vec![2], vec![3]]),
                Some(vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]]),
            ] {
                let settings = EvaluatorSettings {
                    store_atom: name.starts_with("compound_"),
                    ..Default::default()
                };
                let mut baseline = GenericEvaluator::new_from_raw_params(
                    [atom.clone()],
                    &parameters,
                    &builder.fn_map,
                    builder.reps.clone(),
                    settings.optimization_settings(),
                    shape.clone().map(|shape| {
                        let zeros = (1..shape.len())
                            .map(|component| (input_count, component))
                            .collect();
                        (shape, zeros)
                    }),
                    &settings,
                )
                .unwrap();
                let mut prepared_rows = rows.clone();
                let (main, auxiliary) = EvaluatorStack::new_single_parametric(
                    vec![atom.clone()],
                    &builder,
                    &builder.reps,
                    Some(&mut prepared_rows),
                    &shape,
                    &settings,
                )
                .unwrap();
                assert_eq!(
                    auxiliary.as_ref().map_or(0, |(ports, _)| ports.len()),
                    expected_ports,
                    "{name}"
                );
                if name.starts_with("compound_") {
                    let (ports, auxiliary) = auxiliary.as_ref().unwrap();
                    let (output_count, shared) = match name {
                        "compound_add" => (3, 1),
                        "compound_zero_inf" => (2, 1),
                        _ => (2, 0),
                    };
                    assert_eq!(auxiliary.exprs_len, output_count);
                    assert_eq!(prepared_rows.auxiliary_shared_outputs, shared);
                    assert!(prepared_rows.input_columns.is_empty());
                    assert!(prepared_rows.auxiliary_values.is_empty());
                    assert_eq!(prepared_rows.auxiliary_gathers.len(), 2);
                    assert!(
                        prepared_rows
                            .auxiliary_gathers
                            .iter()
                            .all(|row| row.len() == ports.len() - shared)
                    );
                    let root = &main.exprs.as_ref().unwrap()[0];
                    let aliases = main.fn_map_entries.iter().filter(|entry| entry.is_alias);
                    assert_eq!(
                        aliases.clone().count(),
                        usize::from(matches!(
                            name,
                            "compound_mul" | "compound_add" | "compound_inactive_pole"
                        ))
                    );
                    let carrier = aliases
                        .clone()
                        .find(|entry| entry.lhs == *root)
                        .map_or(root, |entry| &entry.rhs);
                    let roots = auxiliary
                        .exprs
                        .as_ref()
                        .unwrap()
                        .iter()
                        .map(|root| {
                            auxiliary
                                .fn_map_entries
                                .iter()
                                .find(|entry| entry.is_alias && entry.lhs == *root)
                                .map_or_else(|| root.clone(), |entry| entry.rhs.clone())
                        })
                        .collect::<Vec<_>>();
                    let expected_carrier = match name {
                        "compound_mul" | "compound_inactive_pole" => &x * &y * &ports[0],
                        "compound_add" => &ports[0] + &ports[1],
                        "compound_zero_inf" => &ports[0] * &ports[1],
                        _ => {
                            assert_ne!(ports[0], ports[1]);
                            assert!(
                                *carrier == &ports[0] - &ports[1]
                                    || *carrier == &ports[1] - &ports[0]
                            );
                            carrier.clone()
                        }
                    };
                    assert_eq!(carrier, &expected_carrier);
                    if shared == 1 {
                        assert_eq!(
                            roots[0],
                            if name == "compound_add" {
                                &x * &y
                            } else {
                                x.pow(-1)
                            }
                        );
                    }
                    if matches!(name, "compound_zero_inf" | "compound_overflow") {
                        assert!(
                            prepared_rows.auxiliary_gathers[0]
                                .iter()
                                .all(Option::is_none)
                        );
                        for output in &prepared_rows.auxiliary_gathers[1] {
                            assert_eq!(
                                roots[output.unwrap()],
                                if name == "compound_zero_inf" {
                                    Atom::one()
                                } else {
                                    x.pow(2)
                                }
                            );
                        }
                        if name == "compound_overflow" {
                            assert_ne!(
                                auxiliary.exprs.as_ref().unwrap()[0],
                                auxiliary.exprs.as_ref().unwrap()[1],
                                "distinct owned alias handles are retained"
                            );
                            assert_ne!(
                                prepared_rows.auxiliary_gathers[1][0],
                                prepared_rows.auxiliary_gathers[1][1]
                            );
                        }
                    } else {
                        let leaves = match name {
                            "compound_mul" => [y.pow(2) + 2, x.pow(2) + 1],
                            "compound_add" => [&x * y.pow(2), x.pow(2) * &y],
                            _ => [&y + 1, x.pow(2) + 1],
                        };
                        for (gather, leaf) in prepared_rows.auxiliary_gathers.iter().zip(leaves) {
                            assert_eq!(gather.len(), 1);
                            assert_eq!(roots[gather[0].unwrap()], leaf);
                        }
                    }
                }
                let mut stack = EvaluatorStack {
                    explicit_orientation_sum_only: false,
                    parametric_rows: Some(prepared_rows),
                    production_orientation_ids: Vec::new(),
                    single_parametric: main,
                    single_parametric_auxiliary: auxiliary,
                    iterative: None,
                    summed_function_map: None,
                    summed: None,
                };
                compare::<f64>(
                    &mut baseline,
                    &mut stack,
                    &rows,
                    input_count,
                    x_index,
                    y_index,
                    &shape,
                    name,
                );
                // Singular-point status is exercised in binary64 above; these
                // regular mixed-axis aliases also certify full native precision.
                if matches!(
                    name,
                    "parent_and_child" | "compound_mul" | "compound_add" | "compound_inactive_pole"
                ) && shape.as_ref().is_some_and(|shape| shape[0].len() == 2)
                {
                    compare::<f128>(
                        &mut baseline,
                        &mut stack,
                        &rows,
                        input_count,
                        x_index,
                        y_index,
                        &shape,
                        name,
                    );
                    compare::<ArbPrec>(
                        &mut baseline,
                        &mut stack,
                        &rows,
                        input_count,
                        x_index,
                        y_index,
                        &shape,
                        name,
                    );
                }
                baseline
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O2,
                        ..Default::default()
                    })
                    .unwrap();
                stack
                    .for_each_generic_evaluator_mut(|evaluator| {
                        evaluator.activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: CompilationOptimizationLevel::O2,
                            ..Default::default()
                        })
                    })
                    .unwrap();
                compare::<f64>(
                    &mut baseline,
                    &mut stack,
                    &rows,
                    input_count,
                    x_index,
                    y_index,
                    &shape,
                    name,
                );
            }
        }
        // Native preparation renames caller aliases before exact input lookup,
        // including an alias key occurring inside a function-shaped input key.
        // Register P(t)=t so the renamed P(q) remains a valid original program.
        let p = symbol!("evaluator_test::physical_alias_input_wrapper");
        let t = Atom::var(symbol!("evaluator_test::physical_alias_formal"));
        for composite in [false, true] {
            let input = if composite {
                function!(p, x.clone())
            } else {
                x.clone()
            };
            let mut scoped_builder = ParamBuilder::new_empty();
            scoped_builder.pairs.additional_params = [input.clone()].into_iter().collect();
            let count = scoped_builder.pairs.update_ranges();
            if composite {
                scoped_builder
                    .add_function_entry(FnMapEntry {
                        lhs: function!(p, t.clone()),
                        rhs: t.clone(),
                        tags: Vec::new(),
                        args: vec![Indeterminate::try_from(t.clone()).unwrap()],
                        inlining: InliningPolicy::Always,
                        is_alias: false,
                    })
                    .unwrap();
            }
            let mut atom = AliasedAtom::from(a.clone());
            atom.register_alias(a.clone(), input);
            atom.register_alias(x.clone(), q.clone());
            let mut prepared_rows = rows.clone();
            let (evaluator, auxiliary) = EvaluatorStack::new_single_parametric(
                vec![atom],
                &scoped_builder,
                &scoped_builder.reps,
                Some(&mut prepared_rows),
                &Some(vec![vec![0], vec![1]]),
                &EvaluatorSettings::default(),
            )
            .unwrap();
            assert!(
                auxiliary.is_none(),
                "a caller alias shadows a physical input key"
            );
            let mut stack = EvaluatorStack {
                explicit_orientation_sum_only: false,
                parametric_rows: Some(prepared_rows),
                production_orientation_ids: Vec::new(),
                single_parametric: evaluator,
                single_parametric_auxiliary: auxiliary,
                iterative: None,
                summed_function_map: None,
                summed: None,
            };
            let input = vec![Complex::new_re(F(0.0)); count * 2];
            let result =
                stack.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
            let [DualOrNot::Dual(result)] = result.as_slice() else {
                panic!("expected dual output")
            };
            assert_eq!(
                result.values,
                [1., 0.].map(|value| Complex::new_re(F(value)))
            );
        }
        // Equal alias handles in different roots belong to different owners.
        // They must remain separate outputs with their own captured definitions.
        let mut first = AliasedAtom::from(a.clone());
        first.register_alias(a.clone(), x.pow(2) + 1);
        let mut second = AliasedAtom::from(a.clone());
        second.register_alias(a.clone(), x.pow(2) + 2);
        let mut prepared_rows = rows.clone();
        let (main, auxiliary) = EvaluatorStack::new_single_parametric(
            vec![first, second],
            &builder,
            &builder.reps,
            Some(&mut prepared_rows),
            &None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        assert_eq!(auxiliary.as_ref().unwrap().0.len(), 2);
        let mut stack = EvaluatorStack {
            explicit_orientation_sum_only: false,
            parametric_rows: Some(prepared_rows),
            production_orientation_ids: Vec::new(),
            single_parametric: main,
            single_parametric_auxiliary: auxiliary,
            iterative: None,
            summed_function_map: None,
            summed: None,
        };
        let input = vec![Complex::new_re(F(0.0)); input_count];
        let result = stack.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
        let result = result
            .into_iter()
            .map(|value| {
                let DualOrNot::NonDual(value) = value else {
                    panic!("expected scalar")
                };
                value
            })
            .collect::<Vec<_>>();
        assert_eq!(result, [2., 4.].map(|value| Complex::new_re(F(value))));

        // Ordinary function bodies resolve global aliases in their own scope.
        // An identically named caller alias cannot certify a captured row port.
        let global = function!(symbol!("evaluator_test::physical_alias_global"));
        let ordinary = function!(symbol!("evaluator_test::physical_alias_global_user"));
        let mut scoped_builder = builder.clone();
        // Keep the captured global alias in the opaque map, not the recorded
        // candidate definitions, so this specifically checks scope reset rather
        // than the conservative owner-key collision rejection.
        scoped_builder
            .fn_map
            .add_aliases([(global.clone(), q.clone())])
            .unwrap();
        scoped_builder
            .add_function_entry(FnMapEntry {
                lhs: ordinary.clone(),
                rhs: global.clone(),
                tags: Vec::new(),
                args: Vec::new(),
                inlining: InliningPolicy::Always,
                is_alias: false,
            })
            .unwrap();
        let mut atom = AliasedAtom::from(a.clone());
        atom.register_alias(a.clone(), &ordinary + 1);
        atom.register_alias(global, x.clone());
        let mut prepared_rows = rows.clone();
        let (evaluator, auxiliary) = EvaluatorStack::new_single_parametric(
            vec![atom],
            &scoped_builder,
            &scoped_builder.reps,
            Some(&mut prepared_rows),
            &None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        assert!(auxiliary.is_none());
        let mut stack = EvaluatorStack {
            explicit_orientation_sum_only: false,
            parametric_rows: Some(prepared_rows),
            production_orientation_ids: Vec::new(),
            single_parametric: evaluator,
            single_parametric_auxiliary: auxiliary,
            iterative: None,
            summed_function_map: None,
            summed: None,
        };
        let input = vec![Complex::new_re(F(0.0)); input_count];
        let result = stack.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
        let [DualOrNot::NonDual(result)] = result.as_slice() else {
            panic!("expected scalar")
        };
        assert_eq!(*result, Complex::new_re(F(3.0)));
    }

    #[test]
    fn joint_rows_preserve_native_precision_derivatives_and_inactive_terms() {
        test_initialise().unwrap();
        let q = Atom::var(symbol!("evaluator_test::joint_precision_q"));
        let x = Atom::var(symbol!("evaluator_test::joint_precision_x"));
        let family = symbol!("evaluator_test::joint_precision_family");
        let precise = Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64);
        let terms = [
            (
                &q * x.pow(2),
                vec![precise.clone(), Rational::one(), Rational::from(4)],
            ),
            // This shorter source must contribute twice, not once per row of
            // the longer source. Its carriers share a singular scalar factor
            // but have distinct rational weights; the inactive q=1 probe below
            // must still skip the singular carrier in the inactive tail.
            (q.pow(-1) * &x, vec![Rational::from(2), Rational::from(3)]),
            // This independent source shares the previous activity guard.
            // Both complete products and their poles must remain lazy after
            // the preprocessed scalar bodies are grouped into one branch.
            (
                (&q - 1).pow(-1) * &x * 2,
                vec![Rational::from(2), Rational::from(3)],
            ),
            // Its coefficient is a constant zero. An unpruned tangent of
            // sqrt(q) would introduce 0 / sqrt(0) into physical derivatives.
            (
                q.pow(Atom::num(Rational::new(1, 2))) * &x,
                vec![Rational::zero()],
            ),
            // This column disappears into an exact literal call argument. Its
            // coefficient must still survive beyond binary64 precision.
            (q.clone(), vec![precise.clone()]),
        ]
        .into_iter()
        .enumerate()
        .map(|(tag, (rhs, rows))| ParametricIntegrandTerm {
            numerator: Arc::new(FnMapEntry {
                lhs: function!(family, tag, q.clone()),
                rhs,
                tags: vec![Atom::num(tag)],
                args: vec![Indeterminate::try_from(q.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            }),
            rows: rows
                .into_iter()
                .map(|coefficient| {
                    let carrier = if tag == 1 {
                        // q equals this exact row coefficient on active rows,
                        // preserving x/(q-1) while exposing shared carrier work.
                        Atom::num(coefficient.clone()) * (&q - 1).pow(-1)
                    } else {
                        Atom::one()
                    };
                    crate::uv::ParametricIntegrandRow {
                        coefficients: vec![coefficient],
                        carrier,
                        source_keys: Vec::new(),
                    }
                })
                .collect(),
        })
        .collect::<Vec<_>>();
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x].into_iter().collect();
        let input_count = builder.pairs.update_ranges();
        let x_index = builder.pairs.additional_params.value_range.start;
        fn check<T: FloatLike>(
            stack: &mut EvaluatorStack,
            input_count: usize,
            x_index: usize,
            precise: &Rational,
        ) {
            let zero = F(T::new_zero());
            let sum_q = zero.from_rational(precise) + zero.from_usize(5);
            let linear = zero.from_rational(&Rational::new(9, 2));
            let orientations =
                TiVec::<OrientationID, _>::from_iter([EdgeVec::from_iter(std::iter::empty::<
                    Orientation,
                >())]);
            let filter = SubSet::full(1);
            for method in [
                EvaluatorMethod::SingleParametric,
                EvaluatorMethod::SummedFunctionMap,
                EvaluatorMethod::Summed,
                EvaluatorMethod::Iterative,
            ] {
                // Repeated and zero probes detect stale scratch values, including tangents.
                for (x, tangent) in [(3, 1), (1, -2), (0, 0), (3, 1)] {
                    let x = zero.from_i64(x);
                    let tangent = zero.from_i64(tangent);
                    let expected = [
                        x.powi(2) * &sum_q
                            + &linear * &x
                            + zero.from_usize(5)
                            + zero.from_rational(precise),
                        (zero.from_usize(2) * &x * &sum_q + &linear) * &tangent,
                    ];
                    let mut values = vec![Complex::new_re(zero.clone()); input_count * 2];
                    values[2 * x_index] = Complex::new_re(x.clone());
                    values[2 * x_index + 1] = Complex::new_re(tangent.clone());
                    let input = InputParams {
                        values: SliceMut::Owned(values),
                        residue_map_id_start: 0,
                        orientations_start: 0,
                        multiplicative_offset: 2,
                    };
                    let mut runtime = RuntimeSettings::default();
                    runtime.general.evaluator_method = method.clone();
                    let output = stack
                        .evaluate(
                            input,
                            SingleOrAllOrientations::All {
                                all: &orientations,
                                filter: &filter,
                            },
                            &runtime,
                            &mut EvaluationMetaData::new_empty(),
                        )
                        .unwrap();
                    let [DualOrNot::Dual(value)] = output.as_slice() else {
                        panic!("expected one dual result")
                    };
                    for (actual, expected) in value.values.iter().zip(&expected) {
                        assert!(
                            (actual.re.clone() - expected).abs()
                                <= expected.abs() * zero.epsilon() * zero.from_usize(64),
                            "{method:?}: {actual:?} != {expected:?}"
                        );
                        assert_eq!(actual.im, zero);
                    }
                    if x == zero.from_usize(3) && zero.epsilon() < zero.from_usize(2).powi(-60) {
                        assert!(value.values[0].re > zero.from_rational(&Rational::new(147, 2)));
                    }
                }
            }
            // At the final ordinal both short families are absent. Force their
            // shared inactive coefficient port to a singular value; the common
            // guard must skip the factored scalar pole, the other numerator's
            // pole, and their physical tangents.
            let mut inactive = stack.clone();
            let rows = inactive.parametric_rows.as_mut().unwrap();
            let inactive_port = (0..rows.parameters.len())
                .find(|&column| {
                    rows.rows
                        .iter()
                        .map(|row| row[column].clone())
                        .collect::<Vec<_>>()
                        == [2, 3, 3].map(Rational::from)
                })
                .unwrap();
            let mut last = rows.rows.last().unwrap().clone();
            last[inactive_port] = Rational::one();
            rows.rows = vec![last];
            if let Some(last) = rows.auxiliary_gathers.last().cloned() {
                rows.auxiliary_gathers = vec![last];
            }
            let mut input = vec![Complex::new_re(zero.clone()); input_count * 2];
            input[x_index * 2] = Complex::new_re(zero.from_usize(3));
            input[x_index * 2 + 1] = Complex::new_re(zero.one());
            let output =
                inactive.evaluate_single_parametric(&input, &mut EvaluationMetaData::new_empty());
            let [DualOrNot::Dual(value)] = output.as_slice() else {
                panic!("expected guarded dual output")
            };
            assert_eq!(
                value.values,
                [36, 24].map(|v| Complex::new_re(zero.from_usize(v)))
            );
        }
        for store_atom in [false, true] {
            for (inline_numerator_functions, do_fn_map_replacements) in
                [(false, false), (true, true), (true, false)]
            {
                let settings = EvaluatorSettings {
                    store_atom,
                    do_fn_map_replacements,
                    inline_numerator_functions,
                    summed_function_map: true,
                    summed: true,
                    iterative_orientation_optimization: true,
                    ..Default::default()
                };
                let (mut stack, _) = EvaluatorStack::from_integrand_with_timings(
                    &Atom::num(5),
                    &builder,
                    &[],
                    &terms,
                    None,
                    Some(vec![vec![0], vec![1]]),
                    &settings,
                )
                .unwrap();
                let rows = stack.parametric_rows.as_ref().unwrap();
                assert_eq!(rows.rows.len(), 3);
                assert_eq!(
                    stack.single_parametric.zero_components,
                    (input_count
                        + stack
                            .single_parametric_auxiliary
                            .as_ref()
                            .map_or(0, |(ports, _)| ports.len())
                        ..input_count
                            + stack
                                .single_parametric_auxiliary
                                .as_ref()
                                .map_or(0, |(ports, _)| ports.len())
                            + rows.input_columns.len())
                        .map(|p| (p, 1))
                        .collect::<Vec<_>>()
                );
                check::<f64>(&mut stack, input_count, x_index, &precise);
                check::<f128>(&mut stack, input_count, x_index, &precise);
                check::<ArbPrec>(&mut stack, input_count, x_index, &precise);
                for optimization in [
                    CompilationOptimizationLevel::O0,
                    CompilationOptimizationLevel::O2,
                ] {
                    stack
                        .for_each_generic_evaluator_mut(|evaluator| {
                            evaluator.activate_symjit(&CompilationOptionsSnapshot {
                                optimization_level: optimization,
                                ..Default::default()
                            })
                        })
                        .unwrap();
                    check::<f64>(&mut stack, input_count, x_index, &precise);
                }
            }
        }
    }

    #[test]
    fn joint_rows_support_an_empty_physical_signature_without_losing_precision() {
        test_initialise().unwrap();
        let q = Atom::var(symbol!("evaluator_test::constant_joint_row"));
        let family = symbol!("evaluator_test::constant_joint_family");
        let term = ParametricIntegrandTerm {
            numerator: Arc::new(FnMapEntry {
                lhs: function!(family, 0, q.clone()),
                rhs: q.clone(),
                tags: vec![Atom::num(0)],
                args: vec![Indeterminate::try_from(q).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            }),
            rows: [
                Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64),
                Rational::from(-1),
            ]
            .into_iter()
            .map(|value| crate::uv::ParametricIntegrandRow {
                coefficients: vec![value],
                carrier: Atom::one(),
                source_keys: Vec::new(),
            })
            .collect(),
        };
        let (mut stack, _) = EvaluatorStack::from_integrand_with_timings(
            &Atom::Zero,
            &ParamBuilder::new_empty(),
            &[],
            &[term],
            None,
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        let actual =
            stack.evaluate_single_parametric::<ArbPrec>(&[], &mut EvaluationMetaData::new_empty());
        let expected =
            F::<ArbPrec>::default().from_rational(&Rational::new(1, 9_007_199_254_740_992i64));
        assert!(expected > expected.zero());
        assert_eq!(
            actual.into_iter().next().unwrap().unwrap_real(),
            Complex::new_re(expected)
        );
    }

    #[test]
    fn shared_sparse_numerator_family_preserves_open_indices_and_residue_map_modes() {
        test_initialise().unwrap();
        let q = Atom::var(symbol!("evaluator_test::sparse_numerator_q"));
        let family = symbol!("evaluator_test::sparse_numerator_family");
        let vector = SPENSO_TAG.tensor_symbol("evaluator_test::sparse_numerator_vector");
        let index = parse_lit!(spenso::mink(4, 1));
        let temporal = GS.energy_delta(index.as_view());
        let definition = Arc::new(FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: function!(family, 3, q.clone()),
            rhs: (&q + 1) * temporal,
            tags: vec![Atom::num(3)],
            args: vec![Indeterminate::try_from(q).unwrap()],
        });
        let orientations = TiVec::<OrientationID, _>::from_iter([
            EdgeVec::from_iter([Orientation::Default]),
            EdgeVec::from_iter([Orientation::Default]),
        ]);
        let production_ids = [OrientationID(7), OrientationID(11)];
        let atom = (production_ids[0].atom() * function!(family, 3, 2)
            + production_ids[1].atom() * function!(family, 3, 5))
            * function!(vector, index);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.residue_map_id = ParamValuePairs::default_from_symbol(GS.residue_map_id);
        builder.pairs.orientations = [GS.sign(EdgeIndex(0))].into_iter().collect();
        builder.pairs.additional_params = (0..4)
            .map(|i| function!(vector, function!(AIND_SYMBOLS.cind, i)))
            .collect();
        let parameter_count = builder.pairs.update_ranges();
        builder.values = vec![vec![Complex::new_re(F(0.0)); parameter_count]];
        for do_fn_map_replacements in [false, true] {
            for store_atom in [false, true] {
                let settings = EvaluatorSettings {
                    iterative_orientation_optimization: true,
                    summed_function_map: true,
                    summed: true,
                    store_atom,
                    do_fn_map_replacements,
                    ..Default::default()
                };
                let (mut stack, _) = EvaluatorStack::new_with_timings(
                    std::slice::from_ref(&atom),
                    &builder,
                    std::slice::from_ref(&definition),
                    &orientations.raw,
                    &production_ids,
                    None,
                    &settings,
                )
                .unwrap();
                let make_input = || {
                    let mut values = vec![Complex::new_re(F(0.0)); parameter_count];
                    for (i, value) in [7.0, 11.0, 13.0, 17.0].into_iter().enumerate() {
                        values[builder.pairs.additional_params.value_range.start + i] =
                            Complex::new_re(F(value));
                    }
                    InputParams {
                        values: SliceMut::Owned(values),
                        residue_map_id_start: builder.pairs.residue_map_id.value_range.start,
                        orientations_start: builder.pairs.orientations.value_range.start,
                        multiplicative_offset: 1,
                    }
                };
                let mut metadata = EvaluationMetaData::new_empty();
                for (id, expected) in [21.0, 42.0].into_iter().enumerate() {
                    assert_eq!(
                        scalar_value(stack.evaluate_parametric(
                            make_input(),
                            SingleOrAllOrientations::Single {
                                orientation: &orientations[OrientationID(id)],
                                id: OrientationID(id),
                            },
                            &mut metadata,
                        )),
                        Complex::new_re(F(expected))
                    );
                }
                let filter = SubSet::full(orientations.len());
                let all = SingleOrAllOrientations::All {
                    all: &orientations,
                    filter: &filter,
                };
                for method in [
                    EvaluatorMethod::SingleParametric,
                    EvaluatorMethod::Iterative,
                    EvaluatorMethod::SummedFunctionMap,
                    EvaluatorMethod::Summed,
                ] {
                    let mut runtime = RuntimeSettings::default();
                    runtime.general.evaluator_method = method.clone();
                    assert_eq!(
                        scalar_value(
                            stack
                                .evaluate(make_input(), all, &runtime, &mut metadata)
                                .unwrap()
                        ),
                        Complex::new_re(F(63.0))
                    );
                }
                if store_atom {
                    // The shared component body must survive archives in every
                    // evaluator mode, including the summed function map.
                    for evaluator in [
                        &stack.single_parametric,
                        &stack.iterative.as_ref().unwrap().0,
                        stack.summed_function_map.as_ref().unwrap(),
                        stack.summed.as_ref().unwrap(),
                    ] {
                        assert!(evaluator.fn_map_entries.iter().any(|entry| {
                            entry
                                .lhs
                                .as_fun_view()
                                .unwrap()
                                .get_symbol()
                                .get_name()
                                .starts_with("gammalooprs::numerator_component_")
                        }));
                    }
                }
            }
        }
    }

    #[test]
    fn selector_free_parametrization_preserves_factorized_scalar() {
        let (a, b, c, d) = symbol!(
            "unselected_scalar_a",
            "unselected_scalar_b",
            "unselected_scalar_c",
            "unselected_scalar_d"
        );
        let expression = (Atom::var(a) + b).pow(7) * (Atom::var(c) + d).pow(5);
        assert_eq!(
            EvaluatorStack::parametrize_residue_map_selectors(
                &expression,
                Atom::var(GS.residue_map_id)
            ),
            expression,
        );
        assert_eq!(
            EvaluatorStack::parametrize_residue_map_selectors(
                expression.clone(),
                Atom::var(GS.residue_map_id)
            ),
            expression,
        );
        let guarded = Symbol::IF.call_args([Atom::var(a), Atom::Zero, expression]);
        assert_eq!(
            EvaluatorStack::parametrize_residue_map_selectors(
                &guarded,
                Atom::var(GS.residue_map_id)
            ),
            guarded,
        );
    }

    #[test]
    fn summed_orientation_selects_reachable_alias_dependencies() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::reachable_alias_x"));
        let y = Atom::var(symbol!("evaluator_test::reachable_alias_y"));
        let head = symbol!("evaluator_test::reachable_alias");
        let outer = function!(head, 0);
        let leaf = function!(head, 1);
        let dead = function!(head, 2);
        // An alias key is opaque: its argument is not a dependency on `dead`.
        let opaque = function!(head, &dead);
        let selected_id = OrientationID(7);
        let other_id = OrientationID(11);
        let sign = GS.sign(EdgeIndex(0));
        let mut source =
            AliasedAtom::from(selected_id.atom() * (&outer + 1).pow(2) + other_id.atom() * &dead);
        source.register_alias(outer, (&leaf + 1) * (&opaque + 2) + other_id.atom() * &dead);
        source.register_alias(leaf, &x * &sign);
        source.register_alias(opaque, y.clone());
        source.register_alias(dead, (&sign - 1).pow(-1));

        let selected = EvaluatorStack::specialize_aliases(
            source.get_root(),
            &source,
            |expression| {
                EdgeVec::from_iter([Orientation::Default]).select(selected_id.select(expression))
            },
            Some(selected_id.0),
        );
        assert_eq!(selected.get_aliases().len(), 3);
        assert!(!selected.get_aliases().contains_key(&function!(head, 2, 7)));
        assert_eq!(selected.into_inner(), ((&x + 1) * (&y + 2) + 1).pow(2));
        // Selection owns a separate DAG; the original inactive body survives.
        assert_eq!(source.get_aliases().len(), 4);

        // Iterative outputs keep their handles. Unlike AliasedAtom::prune, which
        // would keep `dead` for its occurrence inside the opaque key, the
        // inactive body selected to 0^-1 is not retained.
        let orientation = EdgeVec::from_iter([Orientation::Default]);
        let iterative = EvaluatorStack::specialize_aliases(
            source.get_root(),
            &source,
            |atom| GS.collect_orientation_if(orientation.select(selected_id.select(atom))),
            None,
        );
        assert_eq!(iterative.get_aliases().len(), 3);
        assert!(!iterative.get_aliases().contains_key(&function!(head, 2)));
        assert!(
            iterative
                .get_aliases()
                .contains_key(&function!(head, function!(head, 2)))
        );
        assert_eq!(iterative.into_inner(), ((&x + 1) * (&y + 2) + 1).pow(2));
    }

    #[test]
    fn summed_orientation_aliases_preserve_lazy_branches_and_numeric_values() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::selected_alias_x"));
        let head = symbol!("evaluator_test::selected_alias");
        let outer = function!(head, 0);
        let leaf = function!(head, 1);
        let dead = function!(head, 2);
        let production_ids = [OrientationID(7), OrientationID(11)];
        let orientations = [
            EdgeVec::from_iter([Orientation::Default]),
            EdgeVec::from_iter([Orientation::Reversed]),
        ];
        let mut source =
            AliasedAtom::from((production_ids[0].atom() + production_ids[1].atom()) * &outer);
        source.register_alias(
            outer,
            Symbol::IF.call_args([
                production_ids[0].atom(),
                (&leaf + 1).pow(2),
                (&leaf + 3).pow(2),
            ]),
        );
        source.register_alias(leaf, &x * GS.sign(EdgeIndex(0)));
        source.register_alias(dead, (GS.sign(EdgeIndex(0)) - 1).pow(-1));
        for (index, expected) in [(&x + 1).pow(2), (Atom::num(3) - &x).pow(2)]
            .into_iter()
            .enumerate()
        {
            let selected = EvaluatorStack::specialize_aliases(
                source.get_root(),
                &source,
                |expression| orientations[index].select(production_ids[index].select(expression)),
                Some(production_ids[index].0),
            );
            assert_eq!(selected.get_aliases().len(), 2);
            assert_eq!(selected.into_inner(), expected);
        }

        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x].into_iter().collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings {
            direct_translation: true,
            horner_iterations: 0,
            ..Default::default()
        };
        let mut evaluator = EvaluatorStack::new_summed(
            &[source],
            &builder,
            &orientations,
            &production_ids,
            &None,
            &settings,
        )
        .unwrap();
        for (input, expected) in [(0.0, 10.0), (2.0, 10.0), (-3.0, 40.0)] {
            assert_eq!(
                <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&[
                    Complex::new_re(F(input))
                ]),
                Complex::new_re(F(expected))
            );
        }
    }

    #[test]
    fn summed_orientation_preserves_products_and_exact_cancellation() {
        test_initialise().unwrap();
        let left = Atom::var(symbol!("evaluator_test::summed_branch_left"));
        let right = Atom::var(symbol!("evaluator_test::summed_branch_right"));
        let root = (&left + 1).pow(7) * (&right + 2).pow(5);
        let orientations = TiVec::<OrientationID, _>::from_iter([
            EdgeVec::from_iter([Orientation::Default]),
            EdgeVec::from_iter([Orientation::Default]),
        ]);
        let production_ids = [OrientationID(7), OrientationID(11)];
        let cancelling = (production_ids[0].atom() - production_ids[1].atom()) * left.pow(-1);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [left, right].into_iter().collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings {
            store_atom: true,
            direct_translation: true,
            horner_iterations: 0,
            ..Default::default()
        };
        let mut evaluator = EvaluatorStack::new_summed(
            &[
                AliasedAtom::from(root.clone()),
                AliasedAtom::from(cancelling),
            ],
            &builder,
            &orientations.raw,
            &production_ids,
            &None,
            &settings,
        )
        .unwrap();
        assert_eq!(
            evaluator.exprs.as_ref().unwrap(),
            &vec![2 * root, Atom::Zero]
        );
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(
            &[Complex::new_re(F(0.0)); 2],
        );
        let [DualOrNot::NonDual(product), DualOrNot::NonDual(cancelled)] = actual.as_slice() else {
            panic!("expected two scalar outputs");
        };
        assert_eq!(*product, Complex::new_re(F(64.0)));
        assert_eq!(*cancelled, Complex::new_re(F(0.0)));
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
        for ((production_id, orientation), expected) in production_ids
            .iter()
            .zip(orientations.iter())
            .zip([2, 3, 5, 7])
        {
            let selected = production_id.select(&atom);
            assert_eq!(orientation.select(&selected), Atom::num(expected));
        }
        // The public selector keeps an undirected sign symbolic. Compiled
        // orientation sums must instead bind it to zero, independently of a
        // previous evaluation's values in the reusable orientation slots.
        let undirected_sign = GS.sign(EdgeIndex(1));
        assert_eq!(
            orientations[OrientationID(0)].select(&undirected_sign),
            undirected_sign
        );
        let atom =
            atom + Atom::add_many(production_ids.iter().map(|id| id.atom())) * undirected_sign;
        for store_atom in [false, true] {
            let evaluator_settings = EvaluatorSettings {
                iterative_orientation_optimization: true,
                store_atom,
                summed: true,
                summed_function_map: true,
                ..Default::default()
            };
            let (mut stack, _) = EvaluatorStack::new_with_timings(
                std::slice::from_ref(&atom),
                &builder,
                &[],
                &orientations.raw,
                &production_ids,
                None,
                &evaluator_settings,
            )
            .unwrap();

            let make_input = || InputParams {
                values: SliceMut::Owned(vec![Complex::new_re(F(11.0)); parameter_count]),
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
                );
                assert_eq!(scalar_value(actual), Complex::new_re(F(expected)));
            }

            let filter = SubSet::full(orientations.len());
            let all = SingleOrAllOrientations::All {
                all: &orientations,
                filter: &filter,
            };
            assert_eq!(
                scalar_value(stack.evaluate_parametric(make_input(), all, &mut metadata)),
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
                            .evaluate(make_input(), all, &runtime_settings, &mut metadata)
                            .unwrap()
                    ),
                    Complex::new_re(F(17.0))
                );
            }
        }
    }

    #[test]
    fn summed_function_map_dualizes_retained_calls_inside_aliases_without_source_storage() {
        let x = parse_lit!(evaluator_test::summed_retained_x);
        let z = parse_lit!(evaluator_test::summed_retained_z);
        let function = symbol!("evaluator_test::summed_retained_function");
        let alias = function!(symbol!("evaluator_test::summed_retained_alias"), 0);
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        builder.pairs.update_ranges();
        builder
            .add_function_entry(FnMapEntry {
                lhs: function!(function, 0, &z),
                rhs: z.pow(2),
                tags: vec![Atom::num(0)],
                args: vec![Indeterminate::try_from(z.clone()).unwrap()],
                inlining: InliningPolicy::Never,
                is_alias: false,
            })
            .unwrap();
        let ids = [OrientationID(4), OrientationID(9)];
        let mut source = AliasedAtom::from(alias.clone());
        source.register_alias(
            alias,
            (ids[0].atom() + 2 * ids[1].atom()) * function!(function, 0, x),
        );
        let orientations = [
            EdgeVec::from_iter(std::iter::empty::<Orientation>()),
            EdgeVec::from_iter(std::iter::empty::<Orientation>()),
        ];
        // Row substitution must precede function-map lowering as well. The
        // first row discards the singular alias before specializing its body;
        // the second retains it with weight 1/2. Both reuse the numerator.
        let coefficient = parse_lit!(evaluator_test::summed_retained_coefficient);
        let pole = function!(symbol!("evaluator_test::summed_retained_pole"), 0);
        let mut row_source = source
            .clone()
            .map_root(|root| Symbol::IF.call_args([&coefficient - 1, &pole * &root, root]));
        row_source.register_alias(pole, (&coefficient - 1).pow(-1));
        let rows = ParametricResidueRows {
            parameters: vec![coefficient],
            rows: vec![vec![Rational::from(1)], vec![Rational::from(3)]],
            ..Default::default()
        };
        let specialized = rows
            .rows
            .iter()
            .enumerate()
            .map(|(id, row)| rows.specialize(&row_source, row, id))
            .fold(AliasedAtom::default(), |sum, row| {
                sum.try_add(&row).unwrap()
            });
        for (source, expected) in [(source, [27.0, 18.0]), (specialized, [40.5, 27.0])] {
            for store_atom in [false, true] {
                let settings = EvaluatorSettings {
                    store_atom,
                    ..Default::default()
                };
                let mut evaluator = EvaluatorStack::new_summed_function_map(
                    std::slice::from_ref(&source),
                    &builder,
                    &orientations,
                    &ids,
                    &Some(vec![vec![0], vec![1]]),
                    &settings,
                )
                .unwrap();
                let mut output = [Complex::new_re(F(0.0)); 2];
                evaluator.f64_eager.evaluate(
                    &[Complex::new_re(F(3.0)), Complex::new_re(F(1.0))],
                    &mut output,
                );
                assert_eq!(output, expected.map(|value| Complex::new_re(F(value))));
                assert_eq!(evaluator.exprs.is_some(), store_atom);
            }
        }
    }

    #[test]
    fn parameterless_function_body_preserves_multiparameter_hyperdual_derivatives() {
        let x = Atom::var(symbol!("evaluator_test::function_dual_x"));
        let y = Atom::var(symbol!("evaluator_test::function_dual_y"));
        let body = &x * &x + &x * &y + Atom::num(3) * &y;
        let function_symbol = symbol!("evaluator_test::parameterless_dual_function");
        let call = function!(function_symbol, 0);
        let entry = FnMapEntry {
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: call.clone(),
            rhs: body.clone(),
            args: Vec::new(),
            tags: vec![Atom::num(0)],
        };
        let mut function_map = FunctionMap::default();
        function_map
            .add_tagged_function(
                function_symbol,
                vec![Atom::num(0)],
                Vec::<Indeterminate>::new(),
                body.clone(),
            )
            .unwrap();
        let dual_shape = Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1));
        let settings = EvaluatorSettings::default();
        let mut function_evaluator = GenericEvaluator::new_from_raw_params(
            [call],
            &[x.clone(), y.clone()],
            &function_map,
            vec![entry],
            OptimizationSettings::default(),
            dual_shape.clone().map(|shape| (shape, Vec::new())),
            &settings,
        )
        .unwrap();
        let mut materialized = GenericEvaluator::new_from_raw_params(
            [body],
            &[x, y],
            &FunctionMap::default(),
            vec![],
            OptimizationSettings::default(),
            dual_shape.map(|shape| (shape, Vec::new())),
            &settings,
        )
        .unwrap();
        let input = [
            Complex::new_re(F(2.0)),
            Complex::new_re(F(7.0)),
            Complex::new_re(F(5.0)),
            Complex::new_re(F(11.0)),
        ];
        let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut function_evaluator)(&input);
        let expected = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut materialized)(&input);

        let [DualOrNot::Dual(actual)] = actual.as_slice() else {
            panic!("function evaluator did not return the requested dual output")
        };
        let [DualOrNot::Dual(expected)] = expected.as_slice() else {
            panic!("materialized evaluator did not return the requested dual output")
        };
        assert_eq!(actual.values, expected.values);
        assert_eq!(actual.values[0], Complex::new_re(F(29.0)));
        assert_eq!(actual.values[1], Complex::new_re(F(118.0)));
    }

    #[test]
    fn evaluator_source_retention_preserves_nested_functions_and_multiple_outputs() {
        let x = Atom::var(symbol!("evaluator_test::retained_source_x"));
        let inner = symbol!("evaluator_test::retained_inner");
        let outer = symbol!("evaluator_test::retained_outer");
        let inner_call = function!(inner);
        let outer_call = function!(outer);
        let inner_body = (&x + 2).pow(3);
        let outer_body = &inner_call + 1;
        let mut fn_map = FunctionMap::default();
        let mut entries = Vec::new();
        for (symbol, lhs, rhs) in [
            (inner, inner_call.clone(), inner_body),
            (outer, outer_call.clone(), outer_body),
        ] {
            fn_map
                .add_function(symbol, Vec::<Indeterminate>::new(), rhs.clone())
                .unwrap();
            entries.push(FnMapEntry {
                inlining: InliningPolicy::Always,
                is_alias: false,
                lhs,
                rhs,
                args: Vec::new(),
                tags: Vec::new(),
            });
        }
        let source = [outer_call.clone(), inner_call * outer_call];
        for store_atom in [false, true] {
            for do_fn_map_replacements in [false, true] {
                let settings = EvaluatorSettings {
                    store_atom,
                    do_fn_map_replacements,
                    ..Default::default()
                };
                let mut evaluator = GenericEvaluator::new_from_raw_params(
                    source.clone(),
                    std::slice::from_ref(&x),
                    &fn_map,
                    entries.clone(),
                    settings.optimization_settings(),
                    None,
                    &settings,
                )
                .unwrap();
                assert_eq!(evaluator.exprs.is_some(), store_atom);
                if let Some(stored) = &evaluator.exprs
                    && !do_fn_map_replacements
                {
                    assert_eq!(stored, &source);
                }
                for (input, expected) in [(2.0, [65.0, 4160.0]), (-1.0, [2.0, 2.0])] {
                    let values = [Complex::new_re(F(input))];
                    let actual =
                        <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&values);
                    assert_eq!(
                        actual
                            .into_iter()
                            .map(DualOrNot::unwrap_real)
                            .collect::<Vec<_>>(),
                        expected.map(|value| Complex::new_re(F(value)))
                    );
                }
            }
        }
    }

    #[test]
    fn evaluator_retained_aliases_preserve_collisions_substitutions_and_archives() {
        use spenso::network::{
            Network,
            store::{NetworkStore, TensorScalarStore},
            tags::scalar_store_alias,
        };

        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::retained_alias_x"));
        let unregistered = scalar_store_alias(99);
        let original_handle = function!(symbol!("evaluator_test::retained_original"), 0, 0);
        let model_function = symbol!("evaluator_test::retained_alias_model");
        let model_call = function!(model_function);
        let mut builder = ParamBuilder::<f64>::new_empty();
        builder
            .add_function(model_function, Vec::<Indeterminate>::new(), &x + 1)
            .unwrap();
        // A function may be registered without a persistence entry. It still
        // owns its name and must never capture a newly scoped network handle.
        builder
            .fn_map
            .add_tagged_function(
                symbol!("evaluator_retained_scalar_0"),
                vec![Atom::num(0), Atom::num(0)],
                Vec::<Indeterminate>::new(),
                Atom::num(12345),
            )
            .unwrap();

        builder
            .fn_map
            .add_tagged_function(
                symbol!(
                    "gammalooprs::retained_alternate",
                    aliases = ["gammalooprs::evaluator_retained_scalar_1"]
                ),
                vec![Atom::num(0), Atom::num(0)],
                Vec::<Indeterminate>::new(),
                Atom::num(67890),
            )
            .unwrap();

        let mut first: Network<NetworkStore<(), Atom>, i8, i8> =
            Network::from_scalar(model_call.clone());
        first.store.add_scalar(scalar_store_alias(3) * 2);
        first.store.add_scalar(model_call);
        first.store.add_scalar(scalar_store_alias(0) + 3);
        let aliases = first.alias_scalar_refs(|_, _| true);
        let first = first.aliased_atom(
            &aliases,
            scalar_store_alias(1) + scalar_store_alias(2).pow(2) + &unregistered + &original_handle,
        );
        let mut second: Network<NetworkStore<(), Atom>, i8, i8> = Network::from_scalar(&x + 7);
        let aliases = second.alias_scalar_refs(|_, _| true);
        let second = second.aliased_atom(&aliases, scalar_store_alias(0));
        let source = [first, second];
        let original = source.clone().map(AliasedAtom::into_inner);
        let params = [x, unregistered, original_handle];

        for store_atom in [false, true] {
            for do_fn_map_replacements in [false, true] {
                for dual_shape in [
                    None,
                    Some(crate::utils::hyperdual_utils::simple_n_deriv_shape(1)),
                ] {
                    let settings = EvaluatorSettings {
                        store_atom,
                        do_fn_map_replacements,
                        ..Default::default()
                    };
                    let mut retained = GenericEvaluator::new_from_raw_params(
                        source.clone(),
                        &params,
                        &builder.fn_map,
                        builder.reps.clone(),
                        settings.optimization_settings(),
                        dual_shape.clone().map(|shape| (shape, Vec::new())),
                        &settings,
                    )
                    .unwrap();
                    let mut direct = GenericEvaluator::new_from_raw_params(
                        original.clone(),
                        &params,
                        &builder.fn_map,
                        builder.reps.clone(),
                        settings.optimization_settings(),
                        dual_shape.clone().map(|shape| (shape, Vec::new())),
                        &settings,
                    )
                    .unwrap();
                    let mut rebuilt = retained.exprs.as_ref().map(|roots| {
                        let mut function_map = FunctionMap::default();
                        for entry in &retained.fn_map_entries {
                            entry.register(&mut function_map).unwrap();
                        }
                        GenericEvaluator::new_from_raw_params(
                            roots.clone(),
                            &params,
                            &function_map,
                            retained.fn_map_entries.clone(),
                            settings.optimization_settings(),
                            dual_shape.clone().map(|shape| (shape, Vec::new())),
                            &EvaluatorSettings {
                                do_fn_map_replacements: false,
                                ..settings
                            },
                        )
                        .unwrap()
                    });
                    let encoded =
                        bincode::encode_to_vec(&retained, bincode::config::standard()).unwrap();
                    let mut state = Vec::new();
                    State::export(&mut state).unwrap();
                    let state_map = State::import(&mut Cursor::new(state), None).unwrap();
                    let model = Model::default();
                    let (mut decoded, _): (GenericEvaluator, _) =
                        bincode::decode_from_slice_with_context(
                            &encoded,
                            bincode::config::standard(),
                            GammaLoopContextContainer {
                                state_map: &state_map,
                                model: &model,
                            },
                        )
                        .unwrap();
                    let values = if dual_shape.is_some() {
                        vec![2.0, 1.0, 5.0, 0.0, 11.0, 0.0]
                    } else {
                        vec![2.0, 5.0, 11.0]
                    }
                    .into_iter()
                    .map(|value| Complex::new_re(F(value)))
                    .collect::<Vec<_>>();
                    let expected =
                        <f64 as GenericEvaluatorFloat>::get_evaluator(&mut direct)(&values);
                    for evaluator in [&mut retained, &mut decoded]
                        .into_iter()
                        .chain(rebuilt.iter_mut())
                    {
                        let actual =
                            <f64 as GenericEvaluatorFloat>::get_evaluator(evaluator)(&values);
                        for (actual, expected) in actual.iter().zip(&expected) {
                            match (actual, expected) {
                                (DualOrNot::NonDual(actual), DualOrNot::NonDual(expected)) => {
                                    assert_eq!(actual, expected)
                                }
                                (DualOrNot::Dual(actual), DualOrNot::Dual(expected)) => {
                                    assert_eq!(actual.values, expected.values)
                                }
                                _ => panic!("retained aliases changed the output domain"),
                            }
                        }
                    }
                    if !store_atom && !do_fn_map_replacements && dual_shape.is_none() {
                        retained
                            .activate_symjit(&CompilationOptionsSnapshot {
                                optimization_level: CompilationOptimizationLevel::O0,
                                ..Default::default()
                            })
                            .unwrap();
                        let actual =
                            <f64 as GenericEvaluatorFloat>::get_evaluator(&mut retained)(&values);
                        assert_eq!(
                            actual
                                .into_iter()
                                .map(DualOrNot::unwrap_real)
                                .collect::<Vec<_>>(),
                            expected
                                .into_iter()
                                .map(DualOrNot::unwrap_real)
                                .collect::<Vec<_>>()
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn evaluator_retained_aliases_skip_inactive_singular_definitions() {
        let x = Atom::var(symbol!("evaluator_test::retained_singular_x"));
        let y = Atom::var(symbol!("evaluator_test::retained_singular_y"));
        let mass = Atom::var(symbol!("evaluator_test::retained_singular_mass"));
        let active = function!(symbol!("evaluator_test::retained_singular"), 0);
        let finite = function!(symbol!("evaluator_test::retained_singular"), 1);
        let dead = function!(symbol!("evaluator_test::retained_singular"), 2);
        let mut first =
            AliasedAtom::from(Symbol::IF.call_args([x.clone(), active.clone(), finite.clone()]));
        first.register_alias(active.clone(), x.pow(-1));
        first.register_alias(finite.clone(), &x + 7 + mass.pow((1, 2)));
        first.register_alias(dead.clone(), (&x - &y).pow(-1));
        // The same handles have different bodies in the second output owner.
        // Neither their registration nor another output's active branch may
        // make an inactive or entirely unreachable pole execute.
        let mut second =
            AliasedAtom::from(Symbol::IF.call_args([y.clone(), active.clone(), finite.clone()]));
        second.register_alias(active, y.pow(-1));
        second.register_alias(finite, &y + 11);
        second.register_alias(dead, (&x + &y).pow(-1));
        let roots = [second.clone(), first, second];
        for store_atom in [false, true] {
            for do_fn_map_replacements in [false, true] {
                for shape in [
                    None,
                    Some(vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]]),
                ] {
                    let settings = EvaluatorSettings {
                        store_atom,
                        do_fn_map_replacements,
                        ..Default::default()
                    };
                    let width = shape.as_ref().map_or(1, Vec::len);
                    let zeros = (1..width)
                        .map(|component| (2, component))
                        .collect::<Vec<_>>();
                    let mut evaluator = GenericEvaluator::new_from_raw_params(
                        roots.clone(),
                        &[x.clone(), y.clone(), mass.clone()],
                        &FunctionMap::new(),
                        vec![],
                        settings.optimization_settings(),
                        shape.map(|shape| (shape, zeros.clone())),
                        &settings,
                    )
                    .unwrap();
                    assert_eq!(evaluator.exprs_len, roots.len());
                    assert_eq!(evaluator.zero_components, zeros);
                    assert_eq!(evaluator.exprs.is_some(), store_atom);
                    for compiled in [false, true] {
                        if compiled {
                            evaluator
                                .activate_symjit(&CompilationOptionsSnapshot {
                                    optimization_level: CompilationOptimizationLevel::O0,
                                    ..Default::default()
                                })
                                .unwrap();
                        }
                        for (x, y) in [(0.0_f64, 0.0_f64), (2.0, 0.0), (0.0, 2.0), (2.0, 4.0)] {
                            let mut input = vec![Complex::new_re(F(0.0)); 3 * width];
                            input[0] = Complex::new_re(F(x));
                            input[width] = Complex::new_re(F(y));
                            if width > 1 {
                                input[1] = Complex::new_re(F(1.0));
                                input[width + 2] = Complex::new_re(F(1.0));
                            }
                            let first = if x == 0.0 {
                                [7.0, 1.0, 0.0, 0.0]
                            } else {
                                [x.recip(), -x.powi(-2), 0.0, 0.0]
                            };
                            let second = if y == 0.0 {
                                [11.0, 0.0, 1.0, 0.0]
                            } else {
                                [y.recip(), 0.0, -y.powi(-2), 0.0]
                            };
                            let expected = [second, first, second]
                                .into_iter()
                                .flat_map(|jet| jet.into_iter().take(width))
                                .map(|value| Complex::new_re(F(value)))
                                .collect::<Vec<_>>();
                            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(
                                &mut evaluator,
                            )(&input)
                            .into_iter()
                            .flat_map(|value| match value {
                                DualOrNot::NonDual(value) => vec![value],
                                DualOrNot::Dual(dual) => dual.values,
                            })
                            .collect::<Vec<_>>();
                            assert_eq!(actual, expected);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn weighted_orientation_alias_preserves_inactive_pole_laziness() {
        test_initialise().unwrap();
        let key = Atom::var(symbol!("evaluator_test::weighted_guard_key"));
        let x = Atom::var(symbol!("evaluator_test::weighted_guard_x"));
        let weight = function!(symbol!("evaluator_test::weighted_guard_alias"), 7);
        let collected = GS.collect_orientation_if(
            Symbol::IF.call_args([key.clone(), Atom::Zero, Atom::num(3)]) * &weight,
        );
        assert_eq!(
            collected,
            Symbol::IF.call_args([key.clone(), Atom::Zero, Atom::num(3) * &weight])
        );
        let mut retained = AliasedAtom::from(collected);
        retained.register_alias(weight, x.pow(-1));
        let protected = Symbol::IF.call_args([key.clone(), Atom::Zero, Atom::num(3) / &x]);
        let params = [key, x];
        let settings = EvaluatorSettings::default();
        let mut evaluators = [retained, AliasedAtom::from(protected)].map(|atom| {
            GenericEvaluator::new_from_raw_params(
                [atom],
                &params,
                &FunctionMap::new(),
                vec![],
                settings.optimization_settings(),
                None,
                &settings,
            )
            .unwrap()
        });
        for optimization in [
            None,
            Some(CompilationOptimizationLevel::O0),
            Some(CompilationOptimizationLevel::O2),
        ] {
            for evaluator in &mut evaluators {
                if let Some(optimization) = optimization {
                    evaluator
                        .activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: optimization,
                            jit_direct_translation: false,
                            ..Default::default()
                        })
                        .unwrap();
                } else {
                    evaluator.activate_eager();
                }
                // The inverse is singular only in the inactive branch. Compare
                // with an explicitly protected expression, not the original 0/0.
                for (key, x, expected) in [(1.0, 0.0, 0.0), (0.0, 2.0, 1.5)] {
                    assert_eq!(
                        <f64 as GenericEvaluatorFloat>::get_evaluator_single(evaluator)(&[
                            Complex::new_re(F(key)),
                            Complex::new_re(F(x)),
                        ]),
                        Complex::new_re(F(expected))
                    );
                }
            }
        }
    }

    #[test]
    fn evaluator_retained_aliases_preserve_residue_map_modes() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::retained_orientation_x"));
        let weight = Atom::add_many((1..512).map(|i| (&x + i).pow(2)));
        assert!(weight.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES);
        let orientations = TiVec::<OrientationID, _>::from_iter([
            EdgeVec::from_iter([Orientation::Reversed]),
            EdgeVec::from_iter([Orientation::Default]),
        ]);
        let production_ids = [OrientationID(4), OrientationID(9)];
        let temporal = GS.energy_delta(parse_lit!(spenso::mink(4, 1)));
        // Open tensor sums keep the large scalar weights separate from the
        // outer guards during parsing. The second temporal vector closes each
        // sum with T.T = 1 before the retained scalars reach the evaluator.
        let reversed = &weight * &temporal + (&weight + 1) * &temporal;
        let default = (&weight + 5) * &temporal + (&weight + 6) * &temporal;
        let atom = production_ids[0].atom()
            * reversed
            * &temporal
            * GS.sign_theta(-GS.sign(EdgeIndex(0))).pow(-1)
            + production_ids[1].atom() * default * &temporal;
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.residue_map_id = ParamValuePairs::default_from_symbol(GS.residue_map_id);
        builder.pairs.orientations = [GS.sign(EdgeIndex(0))].into_iter().collect();
        builder.pairs.additional_params = [x].into_iter().collect();
        let parameter_count = builder.pairs.update_ranges();
        builder.values = vec![vec![Complex::new_re(F(0.0)); parameter_count]];
        let weight_at_zero = (1..512).map(|i| (i * i) as f64).sum::<f64>();
        for do_fn_map_replacements in [false, true] {
            for store_atom in [false, true] {
                let settings = EvaluatorSettings {
                    summed: true,
                    summed_function_map: true,
                    iterative_orientation_optimization: true,
                    do_fn_map_replacements,
                    store_atom,
                    ..Default::default()
                };
                let (mut stack, _) = EvaluatorStack::new_with_timings(
                    std::slice::from_ref(&atom),
                    &builder,
                    &[],
                    &orientations.raw,
                    &production_ids,
                    None,
                    &settings,
                )
                .unwrap();
                let make_input = || InputParams {
                    values: SliceMut::Owned(vec![Complex::new_re(F(0.0)); parameter_count]),
                    residue_map_id_start: builder.pairs.residue_map_id.value_range.start,
                    orientations_start: builder.pairs.orientations.value_range.start,
                    multiplicative_offset: 1,
                };
                let mut metadata = EvaluationMetaData::new_empty();
                for (id, expected) in [2.0 * weight_at_zero + 1.0, 2.0 * weight_at_zero + 11.0]
                    .into_iter()
                    .enumerate()
                {
                    let selected = SingleOrAllOrientations::Single {
                        orientation: &orientations[OrientationID(id)],
                        id: OrientationID(id),
                    };
                    assert_eq!(
                        scalar_value(stack.evaluate_parametric(
                            make_input(),
                            selected,
                            &mut metadata,
                        )),
                        Complex::new_re(F(expected))
                    );
                }
                let filter = SubSet::full(orientations.len());
                let all = SingleOrAllOrientations::All {
                    all: &orientations,
                    filter: &filter,
                };
                for method in [
                    EvaluatorMethod::SingleParametric,
                    EvaluatorMethod::Iterative,
                    EvaluatorMethod::SummedFunctionMap,
                    EvaluatorMethod::Summed,
                ] {
                    let mut runtime_settings = RuntimeSettings::default();
                    runtime_settings.general.evaluator_method = method;
                    assert_eq!(
                        scalar_value(
                            stack
                                .evaluate(make_input(), all, &runtime_settings, &mut metadata,)
                                .unwrap()
                        ),
                        Complex::new_re(F(4.0 * weight_at_zero + 12.0))
                    );
                }
            }
        }
    }

    #[test]
    fn explicit_sum_preprocesses_tensor_integrands() {
        test_initialise().unwrap();
        let builder = ParamBuilder::new_empty();
        let settings = EvaluatorSettings {
            do_algebra: true,
            ..Default::default()
        };
        let body = Bispinor {}.new_rep(4).g(9, 9);
        let (mut stack, _) =
            EvaluatorStack::new_explicit_sum_with_timings(&[body], &builder, &[], None, &settings)
                .unwrap();
        let (mut expected, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &[Atom::num(4)],
            &builder,
            &[],
            None,
            &settings,
        )
        .unwrap();

        let actual =
            <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut stack.single_parametric)(&[]);
        let expected = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
            &mut expected.single_parametric,
        )(&[]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn evaluator_preprocess_constant_tensor_factors_preserve_numeric_domains() {
        use spenso::network::library::symbolic::ETS;

        test_initialise().unwrap();
        let x = Atom::var(symbol!("evaluator_test::constant_tensor_factor_x"));
        // Deliberately reverse the index labels. The metric and temporal pair
        // are both indexed tensors; their sum must retain its slot ownership.
        let left = parse_lit!(spenso::mink(4, 9));
        let right = parse_lit!(spenso::mink(4, 2));
        let metric = function!(ETS.metric, left.clone(), right.clone());
        let temporal = GS.energy_delta(left.as_view()) * GS.energy_delta(right.as_view());
        let branches = &x * &temporal + (&x + 1) * &temporal;
        for (a, b) in [
            (Atom::num(1) + 2 * Atom::i(), Atom::num(2) - Atom::i()),
            (Atom::num(1) / Atom::num(3), Atom::num(2) / Atom::num(5)),
            (Atom::num(1u64 << 54), Atom::num(1)),
            (Atom::num(1u64 << 51), Atom::num(1)),
            (x.clone(), Atom::num(2)),
        ] {
            let numerator = (a * &metric + b * &temporal) * &branches;
            // The unchanged complete-input route is the numeric oracle for
            // admitted integers and declined rational and symbolic data.
            let mut whole = numerator.parse_into_net().unwrap();
            whole.graph.contract_ready_sum_boundaries();
            whole
                .execute::<SequentialRef, SmallestDegree, _, _, _>(
                    TENSORLIB.read().unwrap().deref(),
                    FUN_LIB.deref(),
                )
                .unwrap();
            let ExecutionResult::Val(expected) = whole.result_scalar().unwrap() else {
                panic!("expected an indexed contraction to produce a scalar");
            };
            for mode in [
                ExecutionMode::Sequential,
                ExecutionMode::SequentialRef,
                ExecutionMode::SequentialExtract,
            ] {
                for contraction in [
                    ContractionMode::SmallestDegree,
                    ContractionMode::MinResultRank,
                ] {
                    let settings = EvaluatorSettings {
                        spenso_execution_mode: (mode, contraction),
                        tensor_network_contraction_order:
                            TensorNetworkContractionOrder::IntermediateCost,
                        ..Default::default()
                    };
                    let actual = EvaluatorStack::preprocess_atom(
                        &numerator,
                        0,
                        &settings,
                        symbol!("evaluator_test::constant_tensor_factor_alias"),
                        None,
                    )
                    .unwrap();
                    let mut evaluator = GenericEvaluator::new_from_raw_params(
                        [actual, expected.as_ref().clone().into()],
                        std::slice::from_ref(&x),
                        &FunctionMap::default(),
                        vec![],
                        settings.optimization_settings(),
                        None,
                        &settings,
                    )
                    .unwrap();
                    let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
                        Complex::new_re(F(2.0)),
                    ])
                    .into_iter()
                    .map(DualOrNot::unwrap_real)
                    .collect::<Vec<_>>();
                    assert_eq!(actual[0], actual[1], "{mode:?}, {contraction:?}");
                    let actual =
                        <ArbPrec as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
                            Complex::new_re(F::<ArbPrec>::from_f64(2.0)),
                        ])
                        .into_iter()
                        .map(DualOrNot::unwrap_real)
                        .collect::<Vec<_>>();
                    assert_eq!(actual[0], actual[1], "Arb: {mode:?}, {contraction:?}");
                }
            }
        }
    }
    #[test]
    fn evaluator_preprocess_scalarizes_abstract_minkowski_contractions_without_algebra() {
        test_initialise().unwrap();
        let abstract_index = parse_lit!(spenso::mink(4, 1));
        let edge = EdgeIndex(7);
        let mapped_momentum = GS.emr_vec(edge, abstract_index.as_view())
            + GS.ose(edge) * GS.energy_delta(abstract_index.as_view());
        let numerator = &mapped_momentum * &mapped_momentum;
        let settings = EvaluatorSettings::default();
        let scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &settings,
            symbol!("evaluator_test::preprocess_scalar"),
            None,
        )
        .unwrap()
        .into_inner();
        let algebraic_scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &EvaluatorSettings {
                do_algebra: true,
                ..settings
            },
            symbol!("evaluator_test::preprocess_scalar"),
            None,
        )
        .unwrap()
        .into_inner();

        assert!((scalar - algebraic_scalar).expand().is_zero());
    }

    #[test]
    fn evaluator_preprocess_additive_tensor_terms_matches_whole_network() {
        test_initialise().unwrap();
        let index = parse_lit!(spenso::mink(4, 1));
        let temporal = GS.energy_delta(index.as_view());
        let momentum = GS.emr_vec(EdgeIndex(7), index.as_view()) + GS.ose(EdgeIndex(7)) * &temporal;
        let numerator = &momentum * &momentum + Atom::num(3) * &temporal * &momentum + Atom::num(5);
        let mut whole = numerator.parse_into_net().unwrap();
        whole
            .execute::<SequentialRef, SmallestDegree, _, _, _>(
                TENSORLIB.read().unwrap().deref(),
                FUN_LIB.deref(),
            )
            .unwrap();
        let ExecutionResult::Val(expected) = whole.result_scalar().unwrap() else {
            panic!("expected a nonconstant scalar");
        };
        for mode in [
            ExecutionMode::Sequential,
            ExecutionMode::SequentialRef,
            ExecutionMode::SequentialExtract,
            ExecutionMode::Parallel,
        ] {
            for (contraction, order) in [
                (
                    ContractionMode::SmallestDegree,
                    TensorNetworkContractionOrder::IntermediateCost,
                ),
                (
                    ContractionMode::MinResultRank,
                    TensorNetworkContractionOrder::IntermediateCost,
                ),
                (
                    ContractionMode::MinResultRank,
                    TensorNetworkContractionOrder::SparseAtomAware,
                ),
                (
                    ContractionMode::MinResultRank,
                    TensorNetworkContractionOrder::AtomAware,
                ),
                (
                    ContractionMode::MinResultRank,
                    TensorNetworkContractionOrder::ResultRankOnly,
                ),
                (
                    ContractionMode::MinResultRank,
                    TensorNetworkContractionOrder::EntryAware,
                ),
            ] {
                let settings = EvaluatorSettings {
                    spenso_execution_mode: (mode, contraction),
                    tensor_network_contraction_order: order,
                    ..Default::default()
                };
                assert_eq!(
                    EvaluatorStack::preprocess_atom(
                        &numerator,
                        0,
                        &settings,
                        symbol!("evaluator_test::preprocess_scalar"),
                        None
                    )
                    .unwrap()
                    .into_inner(),
                    *expected,
                    "{mode:?}, {contraction:?}, {order:?}",
                );
            }
        }
    }

    #[test]
    fn evaluator_preprocess_preserves_large_factorized_coefficients() {
        test_initialise().unwrap();
        let coefficients: Vec<_> = [
            symbol!("evaluator_test::sum_coefficient_left"),
            symbol!("evaluator_test::sum_coefficient_right"),
        ]
        .into_iter()
        .map(|symbol| {
            Atom::add_many((0..512).map(|i| function!(symbol, i)).collect::<Vec<_>>()).pow(7)
        })
        .collect();
        let numerator = Bispinor {}.new_rep(4).g(9, 9) * &coefficients[0]
            + Bispinor {}.new_rep(4).g(8, 8) * &coefficients[1];
        let expected = Atom::num(4) * &coefficients[0] + Atom::num(4) * &coefficients[1];
        let scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &EvaluatorSettings::default(),
            symbol!("evaluator_test::preprocess_scalar"),
            None,
        )
        .unwrap()
        .into_inner();
        // Exact Atom equality checks both aliases without expanding the
        // factorized scalar blocks.
        assert_eq!(scalar, expected);
    }

    #[test]
    fn evaluator_preprocess_collects_retained_scalar_definitions() {
        test_initialise().unwrap();
        let x = parse_lit!(evaluator_test::retained_collection_x);
        let y = parse_lit!(evaluator_test::retained_collection_y);
        let z = parse_lit!(evaluator_test::retained_collection_z);
        let factors = (0..256)
            .map(|i| function!(symbol!("evaluator_test::retained_collection_factor"), i))
            .collect::<Vec<_>>();
        // A product keeps the scalar in one network instead of separate
        // top-level summands, matching the projected-4D input shape.
        let source = &z * Atom::add_many(factors.iter().map(|factor| &x * factor + &y * factor));
        let expected = z * (&x + &y) * Atom::add_many(factors);
        let scalar = EvaluatorStack::preprocess_atom(
            &source,
            0,
            &EvaluatorSettings::default(),
            symbol!("evaluator_test::retained_collection_alias"),
            None,
        )
        .unwrap();
        assert!(!scalar.get_aliases().is_empty());
        assert!(scalar.get_byte_size() < source.as_view().get_byte_size());
        assert_eq!(scalar.into_inner(), expected);
    }

    #[test]
    fn evaluator_preprocess_preserves_symbolic_empty_trace() {
        test_initialise().unwrap();
        let dimension = Atom::var(symbol!("evaluator_test::symbolic_empty_trace_dimension"));
        let trace = function!(
            SPENSO_TAG.trace,
            parse_lit!(spenso::mink(evaluator_test::symbolic_empty_trace_dimension))
        );
        let scalar = EvaluatorStack::preprocess_atom(
            &trace,
            0,
            &EvaluatorSettings {
                do_algebra: false,
                ..Default::default()
            },
            symbol!("evaluator_test::symbolic_empty_trace_scalar"),
            None,
        )
        .unwrap();
        assert_eq!(scalar.get_root(), &dimension);
        assert!(scalar.get_aliases().is_empty());
    }

    #[test]
    fn evaluator_preprocess_preserves_closed_zeros_and_rejects_open_sums() {
        test_initialise().unwrap();
        let settings = EvaluatorSettings::default();
        let index = parse_lit!(spenso::mink(4, 1));
        let spatial = GS.emr_vec(EdgeIndex(7), index.as_view());
        let temporal = GS.energy_delta(index.as_view());
        let closed_zero = &spatial * &temporal;
        let cancellation =
            Bispinor {}.new_rep(4).g(9, 9) - Bispinor {}.new_rep(4).g(8, 8) + closed_zero;
        assert!(
            EvaluatorStack::preprocess_atom(
                &cancellation,
                0,
                &settings,
                symbol!("evaluator_test::preprocess_scalar"),
                None
            )
            .unwrap()
            .into_inner()
            .is_zero()
        );
        assert!(
            EvaluatorStack::preprocess_atom(
                &(spatial + temporal),
                0,
                &settings,
                symbol!("evaluator_test::preprocess_scalar"),
                None
            )
            .is_err()
        );
        assert!(
            EvaluatorStack::preprocess_atom(
                &Atom::Zero,
                0,
                &settings,
                symbol!("evaluator_test::preprocess_scalar"),
                None
            )
            .unwrap()
            .into_inner()
            .is_zero()
        );
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
        let scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &settings,
            symbol!("evaluator_test::preprocess_scalar"),
            None,
        )
        .unwrap()
        .into_inner();
        let explicit_scalar = EvaluatorStack::preprocess_atom(
            &explicit_numerator,
            0,
            &settings,
            symbol!("evaluator_test::preprocess_scalar"),
            None,
        )
        .unwrap()
        .into_inner();
        let algebraic_scalar = EvaluatorStack::preprocess_atom(
            &numerator,
            0,
            &EvaluatorSettings {
                do_algebra: true,
                ..settings
            },
            symbol!("evaluator_test::preprocess_scalar"),
            None,
        )
        .unwrap()
        .into_inner();

        assert!((&scalar - explicit_scalar).expand().is_zero());
        assert!((scalar - algebraic_scalar).expand().is_zero());
    }

    #[test]
    fn evaluator_archives_numerator_sampling_scale_input_without_stored_atoms() {
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
            inlining: InliningPolicy::Always,
            is_alias: false,
            lhs: source_call.clone(),
            rhs: source_body,
            args: Vec::new(),
            tags: Vec::new(),
        };
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            [source_call.clone()],
            std::slice::from_ref(&scale),
            &function_map,
            vec![source_entry.clone()],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        for value in [1.0, -2.0] {
            let input = [Complex::new_re(F(value))];
            assert_eq!(
                <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&input),
                Complex::new_re(F(value + 1.0)),
            );
        }

        let encoded = bincode::encode_to_vec(&evaluator, bincode::config::standard()).unwrap();
        let mut state = Vec::new();
        State::export(&mut state).unwrap();
        let state_map = State::import(&mut Cursor::new(state), None).unwrap();
        let model = Model::default();
        let (mut decoded, _): (GenericEvaluator, _) = bincode::decode_from_slice_with_context(
            &encoded,
            bincode::config::standard(),
            GammaLoopContextContainer {
                state_map: &state_map,
                model: &model,
            },
        )
        .unwrap();

        for value in [1.0, -2.0] {
            let input = [Complex::new_re(F(value))];
            assert_eq!(
                <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut decoded)(&input),
                Complex::new_re(F(value + 1.0)),
            );
        }

        // Argument substitution cannot introduce an absent symbol, but a custom
        // normalizer can introduce M when a formal argument becomes concrete.
        // Reachable bodies can also cancel M during substitution. Neither case
        // changes the registered M input, including for an evaluator with only
        // an unrelated M-bearing shared function body.
        let mut independent = GenericEvaluator::new_from_raw_params(
            [Atom::num(1)],
            std::slice::from_ref(&scale),
            &function_map,
            vec![source_entry],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        for value in [1.0, -2.0] {
            let input = [Complex::new_re(F(value))];
            assert_eq!(
                <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut independent)(&input),
                Complex::new_re(F(1.0)),
            );
        }
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
    fn shared_scalar_factors_preserve_exact_expression_and_existing_aliases() {
        test_initialise().unwrap();
        let x = parse_lit!(evaluator_test::shared_factor_x);
        let factor = (&x + 2).pow(-3);
        let old = function!(symbol!("evaluator_test::existing_factor_alias"));
        let source =
            Atom::add_many((0..32).map(|i| OrientationID(i).atom() * &factor * (&x + i + &old)));
        let mut input = AliasedAtom::from(source.clone());
        input.register_alias(old.clone(), &x + 7);
        let expected = input.clone().into_inner();
        let input_bytes = input.get_byte_size();
        let shared = EvaluatorStack::share_scalar_factors(input, &ParamBuilder::new_empty());
        assert!(shared.get_aliases().len() > 1);
        assert!(shared.get_byte_size() < input_bytes);
        assert_eq!(shared.get_aliases()[&old], &x + 7);
        assert!(shared.get_root().contains_symbol(OrientationID::symbol()));
        assert!(
            shared
                .get_aliases()
                .values()
                .all(|body| { !body.contains_symbol(OrientationID::symbol()) })
        );
        assert_eq!(shared.into_inner(), expected);
    }

    #[test]
    fn shared_scalar_factors_preserve_parameter_tag_and_alias_keys() {
        test_initialise().unwrap();
        let x = parse_lit!(evaluator_test::shared_key_x);
        let y = parse_lit!(evaluator_test::shared_key_y);
        let key = Atom::add_many((1..=8).map(|power| x.pow(power))).pow(-1);
        let free = Atom::add_many((1..=8).map(|power| y.pow(power))).pow(-1);
        let head = symbol!("evaluator_test::shared_key_function");
        let call = function!(head, &key);
        let source = Atom::add_many(
            (1..=4).map(|i| (&x + i) * call.pow(2) + (&y + i) * &key + (&x + &y + i) * &free),
        );
        for mode in 0..5 {
            let mut builder = ParamBuilder::new_empty();
            let mut input = AliasedAtom::from(source.clone());
            let mut params = vec![x.clone(), y.clone()];
            match mode {
                0 => {
                    // Atomic parameters take precedence over functions with the
                    // same head. Its catalog declares an ordinary argument.
                    params.push(call.clone());
                    builder
                        .add_function(
                            head,
                            vec![symbol!("evaluator_test::shared_key_formal")],
                            Atom::num(99),
                        )
                        .unwrap();
                }
                1 => builder
                    .add_tagged_function(
                        head,
                        vec![key.clone()],
                        String::new(),
                        Vec::<Indeterminate>::new(),
                        Atom::num(7),
                    )
                    .unwrap(),
                2 => builder
                    .fn_map
                    .add_tagged_function(
                        head,
                        vec![key.clone()],
                        Vec::<Indeterminate>::new(),
                        Atom::num(7),
                    )
                    .unwrap(),
                3 => input.register_alias(call.clone(), Atom::num(7)),
                4 => builder
                    .add_function(
                        head,
                        vec![symbol!("evaluator_test::shared_key_ordinary_argument")],
                        Atom::num(7),
                    )
                    .unwrap(),
                _ => unreachable!(),
            }
            builder.pairs.additional_params = params.iter().cloned().collect();
            let expected = input.clone().into_inner();
            let shared = EvaluatorStack::share_scalar_factors(input, &builder);
            assert!(shared.get_aliases().values().any(|body| body == &free));
            assert_eq!(shared.clone().into_inner(), expected);
            if mode == 3 {
                assert_eq!(shared.get_aliases()[&call], Atom::num(7));
            }
            if mode == 4 {
                assert!(
                    shared
                        .get_aliases()
                        .values()
                        .any(|body| body == &call.pow(2))
                );
            }
            let settings = EvaluatorSettings {
                do_fn_map_replacements: false,
                ..Default::default()
            };
            let mut evaluator = GenericEvaluator::new_from_raw_params(
                [shared],
                &params,
                &builder.fn_map,
                builder.reps.clone(),
                settings.optimization_settings(),
                None,
                &settings,
            )
            .unwrap();
            let mut values = vec![Complex::new_re(F(1.0)); params.len()];
            if mode == 0 {
                values[2] = Complex::new_re(F(7.0));
            }
            assert_eq!(
                <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&values),
                Complex::new_re(F(690.0)),
                "mode={mode}",
            );
        }
    }

    #[test]
    fn shared_scalar_factors_remain_lazy_in_independent_branches() {
        test_initialise().unwrap();
        let q = parse_lit!(evaluator_test::shared_singular_q);
        let r = parse_lit!(evaluator_test::shared_singular_r);
        let key = Atom::var(GS.residue_map_id);
        let inverse = Atom::add_many((1..=8).map(|power| q.pow(power))).pow(-1);
        let source = Symbol::IF.call_args([&key - 4, Atom::Zero, inverse.clone()])
            + Symbol::IF.call_args([&key - 9, Atom::Zero, (&r + 1) * &inverse])
            + Symbol::IF.call_args([&key - 17, Atom::Zero, Atom::num(5)]);
        let shared =
            EvaluatorStack::share_scalar_factors(source.into(), &ParamBuilder::new_empty());
        assert!(shared.get_aliases().values().any(|body| body == &inverse));
        assert!(shared.get_root().contains_symbol(Symbol::IF));
        let collected = shared.map_root(|root| GS.collect_orientation_if(root));
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            [collected],
            &[key, q, r],
            &FunctionMap::default(),
            vec![],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        for jit_direct_translation in [None, Some(false), Some(true)] {
            if let Some(jit_direct_translation) = jit_direct_translation {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O2,
                        jit_direct_translation,
                        ..Default::default()
                    })
                    .unwrap();
            }
            for (inputs, expected) in [
                ([4.0, 1.0, 0.0], 0.125),
                ([9.0, 1.0, 3.0], 0.5),
                ([17.0, 0.0, 0.0], 5.0),
                ([23.0, 0.0, 0.0], 0.0),
            ] {
                let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(
                    &inputs.map(|value| Complex::new_re(F(value))),
                );
                assert_eq!(
                    actual,
                    Complex::new_re(F(expected)),
                    "jit_direct_translation={jit_direct_translation:?}, inputs={inputs:?}"
                );
            }
        }
    }

    #[test]
    fn shared_scalar_factors_inside_alias_dag_remain_lazy_and_idempotent() {
        test_initialise().unwrap();
        let q = parse_lit!(evaluator_test::shared_dag_q);
        let r = parse_lit!(evaluator_test::shared_dag_r);
        let key = Atom::var(GS.residue_map_id);
        let inverse = Atom::add_many((1..=8).map(|power| q.pow(power))).pow(-1);
        let old_a = function!(symbol!("evaluator_test::shared_dag_a"));
        let old_b = function!(symbol!("evaluator_test::shared_dag_b"));
        let outer = function!(symbol!("evaluator_test::shared_dag_outer"));
        let source = Symbol::IF.call_args([&key - 4, Atom::Zero, old_a.clone()])
            + Symbol::IF.call_args([&key - 9, Atom::Zero, old_b.clone()])
            + Symbol::IF.call_args([&key - 17, Atom::Zero, Atom::num(5)]);
        let mut input = AliasedAtom::from(outer.clone());
        input.register_alias(outer.clone(), source.clone());
        input.register_alias(old_a.clone(), inverse.clone());
        input.register_alias(old_b.clone(), &old_a + &r * &inverse);
        let expected = input.clone().into_inner();
        let input_bytes = input.get_byte_size();
        let builder = ParamBuilder::new_empty();
        let shared = EvaluatorStack::share_scalar_factors(input, &builder);
        assert_eq!(shared.get_aliases().len(), 4);
        assert!(shared.get_byte_size() < input_bytes);
        assert_eq!(shared.get_root(), &outer);
        assert_eq!(shared.get_aliases()[&outer], source);
        let factor_alias = &shared.get_aliases()[&old_a];
        assert_eq!(shared.get_aliases()[factor_alias], inverse);
        assert_eq!(shared.get_aliases()[&old_b], &old_a + &r * factor_alias);
        assert_eq!(shared.clone().into_inner(), expected);
        let repeated = EvaluatorStack::share_scalar_factors(shared.clone(), &builder);
        assert_eq!(repeated.get_root(), shared.get_root());
        assert_eq!(repeated.get_aliases(), shared.get_aliases());
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            [shared],
            &[key, q, r],
            &FunctionMap::default(),
            vec![],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        for jit_direct_translation in [None, Some(false), Some(true)] {
            if let Some(jit_direct_translation) = jit_direct_translation {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level: CompilationOptimizationLevel::O2,
                        jit_direct_translation,
                        ..Default::default()
                    })
                    .unwrap();
            }
            for (inputs, expected) in [
                ([4.0, 1.0, 0.0], 0.125),
                ([9.0, 1.0, 3.0], 0.5),
                ([17.0, 0.0, 0.0], 5.0),
                ([23.0, 0.0, 0.0], 0.0),
            ] {
                let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(
                    &inputs.map(|value| Complex::new_re(F(value))),
                );
                assert_eq!(
                    actual,
                    Complex::new_re(F(expected)),
                    "jit_direct_translation={jit_direct_translation:?}, inputs={inputs:?}"
                );
            }
        }
    }

    #[test]
    fn collected_residue_map_guards_preserve_lazy_values() {
        test_initialise().unwrap();
        let q = parse_lit!(evaluator_test::collected_guard_q);
        let r = parse_lit!(evaluator_test::collected_guard_r);
        let key = Atom::var(GS.residue_map_id);
        let alias = function!(symbol!("evaluator_test::collected_guard_alias"));
        for retained in [false, true] {
            let inverse = if retained { alias.clone() } else { q.pow(-1) };
            let source = Symbol::IF.call_args([&key - 4, Atom::Zero, inverse.clone()])
                + Symbol::IF.call_args([&key - 4, Atom::Zero, inverse.pow(2)])
                + Symbol::IF.call_args([&key - 9, Atom::Zero, r.pow(-1)]);
            let mut source = AliasedAtom::from(source);
            if retained {
                source.register_alias(alias.clone(), q.pow(-1));
            }
            let collected =
                EvaluatorStack::share_scalar_factors(source, &ParamBuilder::new_empty())
                    .map_root(|root| GS.collect_orientation_if(root));
            let mut evaluator = GenericEvaluator::new_from_raw_params(
                [collected],
                &[key.clone(), q.clone(), r.clone()],
                &FunctionMap::default(),
                vec![],
                OptimizationSettings::default(),
                None,
                &EvaluatorSettings::default(),
            )
            .unwrap();
            for jit_direct_translation in [None, Some(false), Some(true)] {
                if let Some(jit_direct_translation) = jit_direct_translation {
                    evaluator
                        .activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: CompilationOptimizationLevel::O2,
                            jit_direct_translation,
                            ..Default::default()
                        })
                        .unwrap();
                } else {
                    evaluator.activate_eager();
                }
                for (inputs, expected) in [
                    ([4.0, 2.0, 0.0], 0.75),
                    ([9.0, 0.0, 2.0], 0.5),
                    ([17.0, 0.0, 0.0], 0.0),
                ] {
                    let actual =
                        <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(
                            &inputs.map(|value| Complex::new_re(F(value))),
                        );
                    assert_eq!(
                        actual,
                        Complex::new_re(F(expected)),
                        "retained={retained}, jit_direct_translation={jit_direct_translation:?}, inputs={inputs:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn symjit_preserves_independently_guarded_function_results() {
        test_initialise().unwrap();
        let q = parse_lit!(evaluator_test::cached_power_q);
        let s = parse_lit!(evaluator_test::cached_power_s);
        let t = parse_lit!(evaluator_test::cached_power_t);
        let power = q.pow(parse_lit!(-1 / 2));
        let sources = [s.clone(), t.clone()]
            .map(|guard| Symbol::IF.call_args([guard, power.clone(), Atom::Zero]));
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            sources,
            &[q, s, t],
            &FunctionMap::default(),
            vec![],
            OptimizationSettings::default(),
            None,
            &EvaluatorSettings::default(),
        )
        .unwrap();
        for compilation in [
            None,
            Some((CompilationOptimizationLevel::O0, false)),
            Some((CompilationOptimizationLevel::O0, true)),
            Some((CompilationOptimizationLevel::O2, false)),
            Some((CompilationOptimizationLevel::O2, true)),
        ] {
            if let Some((optimization_level, jit_direct_translation)) = compilation {
                evaluator
                    .activate_symjit(&CompilationOptionsSnapshot {
                        optimization_level,
                        jit_direct_translation,
                        ..Default::default()
                    })
                    .unwrap();
            } else {
                evaluator.activate_eager();
            }
            for (q, expected_power) in [
                (Complex::new_re(F(4.0)), (0.5, 0.0)),
                (Complex::new(F(0.0), F(2.0)), (0.5, -0.5)),
            ] {
                // The second guard must work when the first call never ran,
                // and when both calls ran but the shared result was reused.
                for guards in [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]] {
                    let actual = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut evaluator)(&[
                        q,
                        Complex::new_re(F(guards[0])),
                        Complex::new_re(F(guards[1])),
                    ]);
                    assert_eq!(actual.len(), 2);
                    for (actual, guard) in actual.into_iter().zip(guards) {
                        let actual = actual.unwrap_real();
                        let expected = if guard == 0.0 {
                            (0.0, 0.0)
                        } else {
                            expected_power
                        };
                        let difference = (actual.re.0 - expected.0).hypot(actual.im.0 - expected.1);
                        let scale = actual
                            .re
                            .0
                            .hypot(actual.im.0)
                            .max(expected.0.hypot(expected.1));
                        assert!(actual.re.0.is_finite() && actual.im.0.is_finite());
                        assert!(
                            difference <= 1e-12 * scale,
                            "compilation={compilation:?}, q={q:?}, guards={guards:?}, actual={actual:?}, expected={expected:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn symjit_direct_translation_survives_reload() {
        use crate::settings::global::GammaloopCompileOptions;

        test_initialise().unwrap();
        let parameter = Atom::var(symbol!("jit_direct_translation_parameter"));
        let leaf = symbol!("jit_capacity_leaf");
        let outer = symbol!("jit_capacity_outer");
        let formals = (0..64)
            .map(|i| Atom::var(symbol!(format!("jit_capacity_argument_{i}"))))
            .collect::<Vec<_>>();
        let mut builder = ParamBuilder::<f64>::new_empty();
        for entry in [
            FnMapEntry {
                lhs: leaf.call_args(&formals),
                rhs: Atom::add_many(
                    formals
                        .iter()
                        .enumerate()
                        .map(|(i, argument)| argument.pow(2) * (i + 2)),
                ),
                tags: vec![],
                args: formals
                    .iter()
                    .cloned()
                    .map(Indeterminate::try_from)
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap(),
                inlining: InliningPolicy::Never,
                is_alias: false,
            },
            FnMapEntry {
                lhs: outer.call(&parameter),
                rhs: leaf.call_args((0..64).map(|i| &parameter + i)) + &parameter,
                tags: vec![],
                args: vec![Indeterminate::try_from(parameter.clone()).unwrap()],
                inlining: InliningPolicy::Never,
                is_alias: false,
            },
        ] {
            builder.add_function_entry(entry).unwrap();
        }
        for jit_direct_translation in [true, false] {
            let settings = EvaluatorSettings {
                compile: true,
                ..Default::default()
            };
            let compilation = GammaloopCompileOptions {
                jit_direct_translation,
                ..Default::default()
            }
            .frozen_mode(&settings);
            let evaluator = GenericEvaluator::new_from_raw_params(
                [outer.call(&parameter)],
                std::slice::from_ref(&parameter),
                &builder.fn_map,
                builder.reps.clone(),
                settings.optimization_settings(),
                None,
                &settings,
            )
            .unwrap();
            // Original signatures have at most64 formals. Hoisted numeric
            // constants widen the nested leaf call beyond that default budget.
            let exported = evaluator.f64_eager.export_instructions();
            let nested_arity = exported
                .sub_evaluators
                .iter()
                .flat_map(|callee| &callee.instructions.instructions)
                .filter_map(|instruction| match instruction {
                    symbolica::evaluate::Instruction::Fun(_, function, _) => Some(function.2.len()),
                    _ => None,
                })
                .max()
                .unwrap_or(0);
            assert!(nested_arity > 64);
            let encoded =
                bincode::encode_to_vec((&evaluator, &compilation), bincode::config::standard())
                    .unwrap();
            let mut state = Vec::new();
            State::export(&mut state).unwrap();
            let state_map = State::import(&mut Cursor::new(state), None).unwrap();
            let model = Model::default();
            let ((mut decoded, decoded_compilation), _): (
                (GenericEvaluator, FrozenCompilationMode),
                _,
            ) = bincode::decode_from_slice_with_context(
                &encoded,
                bincode::config::standard(),
                GammaLoopContextContainer {
                    state_map: &state_map,
                    model: &model,
                },
            )
            .unwrap();
            assert_eq!(decoded_compilation, compilation);
            let FrozenCompilationMode::Symjit(options) = decoded_compilation else {
                panic!("expected frozen JIT compilation settings");
            };
            decoded.activate_symjit(&options).unwrap();
            assert_eq!(decoded.active_f64_backend(), ActiveF64Backend::Symjit);
            // Symbolica serializes the compiled instance as (external functions,
            // IR, settings); inspect those settings to test actual forwarding.
            let compiled = serde_json::to_value(&decoded.symjit_f64.as_ref().unwrap().0).unwrap();
            assert_eq!(
                compiled[2]["direct_translation"].as_bool(),
                Some(jit_direct_translation)
            );
            let actual = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut decoded)(&[
                Complex::new(F(2.0), F(1.0)),
            ]);
            let expected = (0..64).fold(Complex::new(F(2.0), F(1.0)), |sum, i| {
                let real = (i + 2) as f64;
                sum + Complex::new(F(real * real - 1.0), F(2.0 * real)) * F((i + 2) as f64)
            });
            assert_eq!(actual, expected);
            let capacity = compiled[2]["options"]["num_args"]
                .as_str()
                .unwrap()
                .parse::<usize>()
                .unwrap();
            assert!(capacity >= nested_arity && capacity <= 2048);
        }
    }

    #[test]
    fn symjit_rejects_retained_calls_above_supported_capacity() {
        test_initialise().unwrap();
        let head = symbol!("jit_capacity_unsupported");
        let parameters = (0..2049)
            .map(|i| Atom::var(symbol!(format!("jit_capacity_unsupported_{i}"))))
            .collect::<Vec<_>>();
        let mut builder = ParamBuilder::<f64>::new_empty();
        builder
            .add_function_entry(FnMapEntry {
                lhs: head.call_args(&parameters),
                rhs: Atom::add_many(&parameters),
                tags: vec![],
                args: parameters
                    .iter()
                    .cloned()
                    .map(Indeterminate::try_from)
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap(),
                inlining: InliningPolicy::Never,
                is_alias: false,
            })
            .unwrap();
        let settings = EvaluatorSettings::default();
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            [head.call_args(&parameters)],
            &parameters,
            &builder.fn_map,
            builder.reps.clone(),
            settings.optimization_settings(),
            None,
            &settings,
        )
        .unwrap();
        for jit_direct_translation in [false, true] {
            let error = evaluator
                .activate_symjit(&CompilationOptionsSnapshot {
                    jit_direct_translation,
                    ..Default::default()
                })
                .unwrap_err();
            assert!(error.to_string().contains("requires 2049"), "{error}");
            assert!(error.to_string().contains("2048 call arguments"), "{error}");
            assert_eq!(evaluator.active_f64_backend(), ActiveF64Backend::Eager);
        }
    }

    #[test]
    fn symjit_resolves_builtin_constant_slots() {
        test_initialise().unwrap();
        let pi = Atom::var(Symbol::PI);
        for (source, expected) in [
            (pi.clone(), std::f64::consts::PI),
            (pi.pow(-3), std::f64::consts::PI.powi(-3)),
        ] {
            let mut evaluator = GenericEvaluator::new_from_raw_params(
                [source],
                &[],
                &FunctionMap::default(),
                vec![],
                OptimizationSettings::default(),
                None,
                &EvaluatorSettings::default(),
            )
            .unwrap();
            evaluator.activate_eager();
            let eager = scalar_value(<f64 as GenericEvaluatorFloat>::get_evaluator(
                &mut evaluator,
            )(&[]));
            for jit_direct_translation in [false, true] {
                for level in [
                    CompilationOptimizationLevel::O0,
                    CompilationOptimizationLevel::O3,
                ] {
                    evaluator
                        .activate_symjit(&CompilationOptionsSnapshot {
                            optimization_level: level,
                            jit_direct_translation,
                            ..Default::default()
                        })
                        .unwrap();
                    let jit = scalar_value(<f64 as GenericEvaluatorFloat>::get_evaluator(
                        &mut evaluator,
                    )(&[]));
                    for actual in [eager, jit] {
                        assert!(actual.re.0.is_finite() && actual.im.0.is_finite());
                        assert!(
                            (actual.re.0 - expected).abs() <= 4.0 * f64::EPSILON * expected.abs(),
                            "level={level}, jit_direct_translation={jit_direct_translation}, actual={actual:?}, expected={expected}"
                        );
                        assert_eq!(actual.im.0, 0.0);
                    }
                }
            }
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
