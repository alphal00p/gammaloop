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
    numerator::{ParsingNet, aind::Aind, symbolica_ext::NumeratorAtomExt},
    processes::{
        ContractionMode, EvaluatorBuildTimings, EvaluatorSettings, ExecutionMode,
        TensorNetworkContractionOrder,
    },
    settings::{
        RuntimeSettings,
        global::{CompilationOptimizationLevel, FrozenCompilationMode},
    },
    utils::{
        ArbPrec, F, FUN_LIB, FloatLike, GS, Length, RuntimeCache, SamplingFloat, TENSORLIB, W_,
        f128,
        hyperdual_utils::{DualOrNot, new_from_values},
    },
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
    fn new_single_parametric(
        parametric_atoms: Vec<AliasedAtom>,
        param_builder: &ParamBuilder,
        dual_shape: &Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<GenericEvaluator> {
        let _progress_guard = crate::processes::enter_detailed_progress_span(
            "Generating Single Parametric Evaluator",
        );
        let opt_settings = settings.optimization_settings();

        GenericEvaluator::new_from_builder(
            parametric_atoms.into_iter().map(|atom| {
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
            }),
            param_builder,
            dual_shape.clone(),
            opt_settings.clone(),
            settings,
        )
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
                            let map = |atom: &Atom| {
                                GS.collect_orientation_if(
                                    orientation.select(production_id.select(atom)),
                                )
                            };
                            let mut selected = AliasedAtom::from(map(atom.get_root()));
                            for (alias, body) in atom.get_aliases() {
                                selected.register_alias(alias.clone(), map(body));
                            }
                            selected.prune();
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

        let mut alias_entries = Vec::new();
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
                // Explicit arguments keep alias bodies in scope when the
                // integrand call is replaced by its body before compilation.
                let alias_calls = a
                    .get_aliases()
                    .keys()
                    .map(|alias| {
                        let call = FunctionBuilder::from_atom(alias)
                            .add_args(args.iter().cloned().map(Atom::from))
                            .finish();
                        Replacement::new(alias.to_pattern(), call)
                    })
                    .collect::<Vec<_>>();
                let map = |atom: &Atom| {
                    GS.collect_orientation_if(Self::parametrize_residue_map_selectors(
                        atom,
                        Atom::var(residue_map_id_arg),
                    ))
                    .replace_multiple(&alias_calls)
                };
                let param_integrand = map(a.get_root());
                for (alias, body) in a.get_aliases() {
                    let rhs = map(body);
                    let function = alias.as_fun_view().unwrap();
                    let tags = function
                        .iter()
                        .map(|arg| arg.to_owned())
                        .collect::<Vec<_>>();
                    fn_map
                        .add_tagged_function(
                            function.get_symbol(),
                            tags.clone(),
                            args.clone(),
                            rhs.clone(),
                        )
                        .map_err(|e| eyre!(e))?;
                    if settings.store_atom {
                        alias_entries.push(FnMapEntry {
                            lhs: alias.replace_multiple(&alias_calls),
                            rhs,
                            tags,
                            args: args.clone(),
                        });
                    }
                }
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
                    rhs: param_integrand,
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
        let entries = param_builder.reps.iter().cloned().chain(entries).collect();
        let mut evaluator = GenericEvaluator::new_from_raw_params(
            sum,
            &params,
            &fn_map,
            entries,
            settings.optimization_settings(),
            dual_shape.clone().map(|shape| (shape, Vec::new())),
            settings,
        )?;
        evaluator.fn_map_entries.extend(alias_entries);
        Ok(evaluator)
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
        let sum = atoms.iter().map(|atom| {
            orientations
                .iter()
                .zip(production_orientation_ids)
                .map(|(orientation, production_id)| {
                    // Concrete orientations own distinct scalar definitions.
                    // Retain each network scope and append its residue-map key.
                    let renames = atom
                        .get_aliases()
                        .keys()
                        .map(|alias| {
                            let function = alias.as_fun_view().unwrap();
                            let scoped = FunctionBuilder::new(function.get_symbol())
                                .add_args(function.iter())
                                .add_arg(production_id.0)
                                .finish();
                            (alias.clone(), scoped)
                        })
                        .collect::<HashMap<_, _>>();
                    let max_handle_bytes = renames
                        .keys()
                        .map(|handle| handle.as_view().get_data().len())
                        .max()
                        .unwrap_or(0);
                    let map = |atom: &Atom| {
                        orientation.select(production_id.select(atom)).replace_map(
                            |view, _, out| {
                                if view.get_data().len() <= max_handle_bytes
                                    && let Some(scoped) = renames.get(view.get_data())
                                {
                                    out.set_from_view(&scoped.as_view());
                                }
                            },
                        )
                    };
                    let mut selected = AliasedAtom::from(map(atom.get_root()));
                    for (alias, body) in atom.get_aliases() {
                        selected.register_alias(renames[alias].clone(), map(body));
                    }
                    selected.prune();
                    debug!(selected_expr = %selected.log_print(None), "Summed");
                    selected
                })
                .fold(AliasedAtom::default(), |acc, atom| {
                    acc.try_add(&atom).unwrap()
                })
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
            numerator_definitions,
            &[],
            &[],
            dual_shape,
            &direct_settings,
        )?;
        stack.explicit_orientation_sum_only = true;
        Ok((stack, timings))
    }

    pub(crate) fn from_integrand_with_timings(
        integrand: &Atom,
        param_builder: &ParamBuilder,
        numerator_definitions: &[Arc<FnMapEntry>],
        orientation_catalog: Option<(&[EdgeVec<Orientation>], &[OrientationID])>,
        dual_shape: Option<Vec<Vec<usize>>>,
        settings: &EvaluatorSettings,
    ) -> Result<(Self, EvaluatorBuildTimings)> {
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
            |net, scalar_aliases, term_index| {
                let root = match net.result_scalar()? {
                    ExecutionResult::One => Atom::num(1),
                    ExecutionResult::Zero => Atom::Zero,
                    ExecutionResult::Val(value) => value.into_owned(),
                };
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
                    let body = if body == alias { body } else { map(body) };
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
                    scalar.as_view().get_byte_size() >= NETWORK_SCALAR_ALIAS_MIN_BYTES
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
                    let mut private_indices = HashMap::new();
                    let body = definition.rhs.replace_map(|part, _, output| {
                        if let Ok(slot) = Slot::<LibraryRep, Aind>::try_from(part) {
                            if !interface.contains(&part.to_owned()) {
                                let fresh =
                                    private_indices.entry(slot.aind()).or_insert_with(|| {
                                        loop {
                                            let name = format!(
                                                "{}_index_{next_index}",
                                                argument_symbol.get_name()
                                            );
                                            next_index += 1;
                                            if occupied_symbols.insert(name.clone()) {
                                                break Atom::var(symbol!(&name));
                                            }
                                        }
                                    });
                                **output = slot.rep().to_symbolic([fresh.clone()]);
                            } else {
                                // A dual slot owns its inner representation;
                                // do not revisit that as a separate base slot.
                                output.set_from_view(&part);
                            }
                        } else if let AtomView::Fun(function) = part
                            && !part.is_tensorial(StrictTensorFilter::Tagged)
                            && ![*SYM, *ANTISYM, *CYCLIC].contains(&function.get_symbol())
                        {
                            // Scalar function payloads are opaque to the
                            // tensor interface, including parameter aliases.
                            output.set_from_view(&part);
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
                        || (part != root && !dangling(part).is_some_and(|slots| !slots.is_empty()))
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
        }
        if definitions.iter().any(|definition| {
            tag_counts
                .keys()
                .any(|head| definition.rhs.contains_symbol(*head))
        }) {
            return Err(eyre!("nested shared numerator families are unsupported"));
        }
        let (atoms, combined) = Self::combine_numerator_families(
            atoms,
            definitions,
            param_builder,
            &families,
            &tag_counts,
            (joint_symbol, argument_symbol, shadow_symbol),
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
        let mut replacements = Vec::new();
        let mut generated = Vec::new();
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
                    let (_, aliases) = net
                        .aliased_atom(scalar_aliases, Atom::Zero)
                        .into_inner_with_aliases();
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
                        if body.is_zero() {
                            continue;
                        }
                        let index = Atom::from(tensor.co_expanded_index(index)?);
                        let body = body.to_owned();
                        let positions = argument_positions(&body);
                        populated.insert(index.clone(), positions.clone());
                        generated.push(FnMapEntry {
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
            replacements.push(
                FnMapEntry {
                    lhs: definition.lhs.clone(),
                    rhs: shadow.into_inner(),
                    args: definition.args.clone(),
                    tags: definition.tags.clone(),
                }
                .replacement(),
            );
        }
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
                let shadowed = atom.as_atom_view().replace_multiple(&replacements);
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
        for (index, entry) in generated.into_iter().enumerate() {
            if reachable.contains(&index) {
                builder
                    .add_tagged_function(
                        entry.lhs.as_fun_view().unwrap().get_symbol(),
                        entry.tags,
                        String::new(),
                        entry.args,
                        entry.rhs,
                    )
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
        let alias_symbol = loop {
            let scope = NETWORK_SCALAR_ALIAS_SCOPE.fetch_add(1, Ordering::Relaxed);
            let name = format!("gammalooprs::evaluator_scalar_{scope}");
            if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
                existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
            }) {
                break symbol!(name.as_str());
            }
        };
        let (parsed_atoms, prepared_builder) = Self::preprocess_numerator_families(
            atoms,
            param_builder,
            numerator_definitions,
            settings,
            alias_symbol,
        )?;
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
        // The optional variants have finished borrowing these scalars. Transfer
        // them to the final evaluator instead of retaining a second full copy.
        let single_parametric =
            Self::new_single_parametric(parsed_atoms, param_builder, &dual_shape, settings)
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
            );
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

            // The atom already contains the complete orientation sum, so
            // applying orientation selection again would double count it.
            return Ok(evaluate_evaluator(
                &mut self.single_parametric,
                input.as_slice(),
                evaluation_metadata,
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

    pub(crate) fn activate_symjit(
        &mut self,
        optimization_level: CompilationOptimizationLevel,
    ) -> Result<()> {
        if self.is_eager_only() {
            self.activate_eager_only();
            return Ok(());
        }

        // Use the same numeric program as eager and external compilation. Rational
        // constant slots (such as pi) remain placeholders until domain mapping.
        // SymJIT 2.21 supports optimization levels up to O2 and cannot compact some complex
        // temporary layouts. Its function-result cache crosses inactive branches,
        // and optimization can drop stores that later calls still read. Disable
        // that cache; Symbolica already shares expressions in this numeric program.
        let evaluator = self
            .f64_eager
            .clone()
            .map_coeff(&|c| SymComplex::new(c.re.0, c.im.0))
            .jit_compile(
                JITCompilationSettings::new()
                    .optimization_level(usize::from(optimization_level).min(2) as u8)
                    .with_option("compact", "false")
                    .with_option("cse", "false"),
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
        // A fresh symbol also avoids definitions supplied only through the
        // opaque FunctionMap; its registration API silently keeps duplicate keys.
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

        let mut tree: Option<ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>>> = None;
        for (atom_index, n) in exprs.iter().enumerate() {
            let build_started = std::time::Instant::now();
            crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                stage = "evaluator_symbolica_build_start",
                atom_index,
                atom_bytes = n.get_byte_size(),
                "Evaluator timing milestone"
            );
            let eval: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = n
                .get_root()
                .evaluator(params)
                .function_map(fn_map.clone())
                .add_aliases(n.get_aliases().iter().map(|(alias, body)| (alias.clone(), body.clone())))
                .map_err(|e| eyre!(e))?
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

            crate::debug_tags!(#generation, #profile, #compile, #term, #summary;
                stage = "evaluator_symbolica_build_done",
                atom_index,
                elapsed_ms = build_started.elapsed().as_secs_f64() * 1000.0,
                "Evaluator timing milestone"
            );
            tree = Some(if let Some(mut tree) = tree {
                tree.merge(eval, settings.cpe_iterations)
                    .map_err(|e| eyre!("Failed to merge evaluators: {}", e))?;
                tree
            } else {
                eval
            });
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
        if let Some(dual_shape) = &dual_shape {
            let dual = HyperDual::<SymComplex<Rational>>::new(dual_shape.clone());
            let dualizer = Dualizer::new(dual, zero_components);
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

            fn get_evaluator(
                $evaluator: &mut GenericEvaluator,
            ) -> impl FnMut(&[Complex<F<Self>>]) -> Vec<DualOrNot<Complex<F<Self>>>> {
                |params: &[Complex<F<Self>>]| {
                    // info!("USING COMPLEX EAGER MULTIPLE");
                    let mut out = vec![Complex::default(); $evaluator.compute_out_size()];
                    ($program).evaluate(params, &mut out);

                    if let Some(dual_shape) = &$evaluator.dual_shape {
                        let dual_builder = HyperDual::<Complex<F<Self>>>::new(dual_shape.clone());
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
        *value
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
        let spatial = GS.emr_vec_index(EdgeIndex(7), index.as_view());
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
                String::new(),
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
                    .activate_symjit(CompilationOptimizationLevel::O0)
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
                lhs: function!(family, 0, q.clone()),
                rhs: &q * &x,
                tags: vec![Atom::num(0)],
                args: vec![Indeterminate::try_from(q.clone()).unwrap()],
            }),
            Arc::new(FnMapEntry {
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
                alias
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
                    alias
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
                alias
            )
            .is_err()
        );
        let mut nested = (*definitions[0]).clone();
        nested.rhs = function!(family, 1, 3);
        assert!(
            EvaluatorStack::preprocess_numerator_families(
                &[function!(family, 0, 3)],
                &builder,
                &[Arc::new(nested), definitions[1].clone()],
                &settings,
                alias
            )
            .is_err()
        );
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
                lhs: function!(family, 0, q1.clone(), q2.clone()),
                rhs: (&q1 + &x) * &weight * GS.energy_delta(index.as_view()),
                tags: vec![Atom::num(0)],
                args: formals.clone(),
            }),
            Arc::new(FnMapEntry {
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
                    .activate_symjit(CompilationOptimizationLevel::O0)
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
            lhs: alias.clone(),
            rhs: &q * x,
            tags: vec![Atom::num(0)],
            args: vec![],
        });
        let definition = Arc::new(FnMapEntry {
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
        )
        .unwrap();
        // Both argument rows reuse one closed scalar body, rather than the
        // separately materialized open components of either numerator.
        assert_eq!(prepared.reps.len(), 1);
        assert_eq!(prepared.reps[0].args.len(), 2);
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
        )
        .unwrap();
        assert_eq!(prepared.reps.len(), 1);
        assert_eq!(prepared.reps[0].args.len(), 1);
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
        builder.pairs.additional_params = std::iter::once(x)
            .chain((0..4).flat_map(|i| {
                [
                    function!(p, function!(AIND_SYMBOLS.cind, i)),
                    function!(r, function!(AIND_SYMBOLS.cind, i)),
                ]
            }))
            .collect();
        builder.pairs.update_ranges();
        let settings = EvaluatorSettings::default();
        let (atoms, prepared) = EvaluatorStack::preprocess_numerator_families(
            &[root],
            &builder,
            &[definition],
            &settings,
            symbol!("evaluator_test::repeated_joint_scalar"),
        )
        .unwrap();
        assert_eq!(prepared.reps.len(), 1);
        assert_eq!(prepared.reps[0].args.len(), 2);
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
            lhs: alias.clone(),
            rhs: &q * x,
            tags: vec![Atom::num(0)],
            args: vec![],
        });
        let definition = Arc::new(FnMapEntry {
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
                            lhs: function!(family, 1, q1.clone()),
                            rhs: first.clone(),
                            tags: vec![Atom::num(1)],
                            args: vec![Indeterminate::try_from(q1.clone()).unwrap()],
                        }),
                        Arc::new(FnMapEntry {
                            lhs: function!(family, 2, q2.clone()),
                            rhs: second.clone(),
                            tags: vec![Atom::num(2)],
                            args: vec![Indeterminate::try_from(q2.clone()).unwrap()],
                        }),
                    ]
                } else {
                    vec![Arc::new(FnMapEntry {
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
                    runtime.general.evaluator_method = method;
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
                        runtime.general.evaluator_method = method;
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
                        .activate_symjit(CompilationOptimizationLevel::O0)
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
    fn shared_sparse_numerator_family_preserves_open_indices_and_residue_map_modes() {
        test_initialise().unwrap();
        let q = Atom::var(symbol!("evaluator_test::sparse_numerator_q"));
        let family = symbol!("evaluator_test::sparse_numerator_family");
        let vector = SPENSO_TAG.tensor_symbol("evaluator_test::sparse_numerator_vector");
        let index = parse_lit!(spenso::mink(4, 1));
        let temporal = GS.energy_delta(index.as_view());
        let definition = Arc::new(FnMapEntry {
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
                    runtime.general.evaluator_method = method;
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
        let evaluator_settings = EvaluatorSettings {
            summed: true,
            summed_function_map: true,
            ..Default::default()
        };
        let (mut stack, _) = EvaluatorStack::new_with_timings(
            &[atom],
            &builder,
            &[],
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

    #[test]
    fn parameterless_function_body_preserves_multiparameter_hyperdual_derivatives() {
        let x = Atom::var(symbol!("evaluator_test::function_dual_x"));
        let y = Atom::var(symbol!("evaluator_test::function_dual_y"));
        let body = &x * &x + &x * &y + Atom::num(3) * &y;
        let function_symbol = symbol!("evaluator_test::parameterless_dual_function");
        let call = function!(function_symbol, 0);
        let entry = FnMapEntry {
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
                            function_map
                                .add_tagged_function(
                                    entry.lhs.as_fun_view().unwrap().get_symbol(),
                                    entry.tags.clone(),
                                    entry.args.clone(),
                                    entry.rhs.clone(),
                                )
                                .unwrap();
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
                            .activate_symjit(CompilationOptimizationLevel::O0)
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
        let active = function!(symbol!("evaluator_test::retained_singular"), 0);
        let finite = function!(symbol!("evaluator_test::retained_singular"), 1);
        let mut atom =
            AliasedAtom::from(Symbol::IF.call_args([x.clone(), active.clone(), finite.clone()]));
        atom.register_alias(active, x.pow(-1));
        atom.register_alias(finite, &x + 7);
        for do_fn_map_replacements in [false, true] {
            let settings = EvaluatorSettings {
                do_fn_map_replacements,
                ..Default::default()
            };
            let mut evaluator = GenericEvaluator::new_from_raw_params(
                [atom.clone()],
                std::slice::from_ref(&x),
                &FunctionMap::new(),
                vec![],
                settings.optimization_settings(),
                None,
                &settings,
            )
            .unwrap();
            for compiled in [false, true] {
                if compiled {
                    evaluator
                        .activate_symjit(CompilationOptimizationLevel::O0)
                        .unwrap();
                }
                for (input, expected) in [(0.0, 7.0), (2.0, 0.5)] {
                    assert_eq!(
                        <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut evaluator)(&[
                            Complex::new_re(F(input))
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
        let mapped_momentum = GS.emr_vec_index(edge, abstract_index.as_view())
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
        let momentum =
            GS.emr_vec_index(EdgeIndex(7), index.as_view()) + GS.ose(EdgeIndex(7)) * &temporal;
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
        let spatial = GS.emr_vec_index(EdgeIndex(7), index.as_view());
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
            let mut collected = AliasedAtom::from(GS.collect_orientation_if(source));
            if retained {
                collected.register_alias(alias.clone(), q.pow(-1));
            }
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
            for compiled in [false, true] {
                if compiled {
                    evaluator
                        .activate_symjit(CompilationOptimizationLevel::O2)
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
                        "retained={retained}, compiled={compiled}, inputs={inputs:?}"
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
        for level in [
            None,
            Some(CompilationOptimizationLevel::O0),
            Some(CompilationOptimizationLevel::O2),
        ] {
            if let Some(level) = level {
                evaluator.activate_symjit(level).unwrap();
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
                            "level={level:?}, q={q:?}, guards={guards:?}, actual={actual:?}, expected={expected:?}"
                        );
                    }
                }
            }
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
            for level in [
                CompilationOptimizationLevel::O0,
                CompilationOptimizationLevel::O3,
            ] {
                evaluator.activate_symjit(level).unwrap();
                let jit = scalar_value(<f64 as GenericEvaluatorFloat>::get_evaluator(
                    &mut evaluator,
                )(&[]));
                for actual in [eager, jit] {
                    assert!(actual.re.0.is_finite() && actual.im.0.is_finite());
                    assert!(
                        (actual.re.0 - expected).abs() <= 4.0 * f64::EPSILON * expected.abs(),
                        "level={level}, actual={actual:?}, expected={expected}"
                    );
                    assert_eq!(actual.im.0, 0.0);
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
