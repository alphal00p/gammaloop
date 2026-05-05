use crate::{
    DependentMomentaConstructor, GammaLoopContext, GammaLoopContextContainer,
    cff::{
        CutCFFIndex,
        esurface::Esurface,
        expression::{GraphOrientation, OrientationID},
        surface::HybridSurfaceID,
    },
    graph::{
        ExternalConnection, FeynmanGraph, Graph, GraphGroup, GroupId, LmbIndex, LoopMomentumBasis,
    },
    integrands::{
        HasIntegrand,
        evaluation::{EvaluationResult, GraphEvaluationResult},
        process::{
            ChannelIndex, GraphTermEvaluationContext, ParamBuilder,
            evaluators::{ActiveF64Backend, EvaluatorStack, evaluate_evaluator_single},
            param_builder::LUParams,
            prepare_buffered_event,
        },
    },
    model::Model,
    momentum::{
        Energy, FourMomentum, Rotation, RotationMethod, ThreeMomentum,
        sample::{ExternalIndex, LoopMomenta, MomentumSample, Subspace},
    },
    observables::{AdditionalWeightKey, EventProcessingRuntime, GenericEvent, GenericEventGroup},
    processes::{
        CrossSectionCut, CrossSectionGraph, CutId, GraphGenerationStats, RaisedCutData, RaisedCutId,
    },
    settings::{
        GlobalSettings, RuntimeSettings, global::FrozenCompilationMode, runtime::IntegralUnit,
    },
    subtraction::{
        generate_rstar_t_dependence_evaluator,
        lu_counterterm::{LUCTKinematicPoint, LUCounterTerm, LUCounterTermEvaluators},
    },
    utils::{
        F, FloatLike, Length, h, h_dual,
        hyperdual_utils::{
            DualOrNot, extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            simple_n_deriv_shape,
        },
        newton_solver::{NewtonIterationResult, newton_iteration_and_derivative},
        serde_utils::SmartSerde,
    },
};
use bincode::Encode;
use bincode_trait_derive::Decode;
use color_eyre::{Result, owo_colors::OwoColorize};
use eyre::Context;
use eyre::eyre;
use std::{
    collections::{BTreeMap, HashSet},
    slice,
    time::{Duration, Instant},
};

use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Flow, HedgePair, Orientation},
    subgraph::{ModifySubSet, SubSetLike, subset::SubSet},
};
use rayon::{
    ThreadPool,
    iter::{IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::algebra::complex::Complex;
use std::{
    fs::{self},
    path::Path,
    vec,
};
use symbolica::{
    domains::{dual::HyperDual, float::SingleFloat},
    numerical_integration::{Grid, Sample},
};
use tracing::debug;
use typed_index_collections::{TiVec, ti_vec};

use super::{
    GraphTerm, LmbMultiChannelingSetup, ProcessIntegrandImpl, RuntimeCache, create_grid,
    evaluate_sample, filtered_orientation_count, format_lmb_channel_label,
    format_orientation_label, histogram_process_info_for_integrand, resolve_visible_orientation_id,
};

pub mod export;
pub mod load;

#[allow(clippy::excessive_precision)]
const PICOBARN_CONVERSION: F<f64> = F(3.89379372171859372125651613062e8);

fn barn_conversion_factor<T: FloatLike>(unit: IntegralUnit, one: F<T>) -> F<T> {
    let Some(relative_to_picobarn) = unit.relative_to_picobarn_factor(&one) else {
        return one.one();
    };

    F::from_ff64(PICOBARN_CONVERSION) * relative_to_picobarn
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionIntegrand {
    pub settings: RuntimeSettings,
    pub data: CrossSectionIntegrandData,
    pub(crate) event_processing_runtime: RuntimeCache<EventProcessingRuntime>,
    pub(crate) active_f64_backend: RuntimeCache<ActiveF64Backend>,
}
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionIntegrandData {
    pub name: String,
    pub compilation: FrozenCompilationMode,
    pub loop_cache_id: usize,
    pub external_cache_id: usize,
    /// Cache ID for the base (unrotated) external momentum configuration
    pub base_external_cache_id: usize,
    // pub polarizations: Vec<Polarizations>,
    pub rotations: Option<Vec<Rotation>>,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
    pub n_incoming: usize,
    pub external_connections: Vec<ExternalConnection>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub graph_to_group_id: Vec<usize>,
    // pub builder_cache: ParamBuilder<f64>,
}

impl CrossSectionIntegrand {
    pub(crate) fn frozen_compilation(&self) -> &FrozenCompilationMode {
        &self.data.compilation
    }

    pub(crate) fn active_f64_backend(&self) -> ActiveF64Backend {
        self.active_f64_backend
            .as_ref()
            .copied()
            .unwrap_or(ActiveF64Backend::Eager)
    }

    fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        for graph_term in &mut self.data.graph_terms {
            graph_term.for_each_generic_evaluator_mut(&mut f)?;
        }
        Ok(())
    }

    fn has_complete_external_artifacts(&mut self) -> Result<bool> {
        let mut has_all = true;
        self.for_each_generic_evaluator_mut(|evaluator| {
            has_all &= evaluator.has_external_compiled_artifact();
            Ok(())
        })?;
        Ok(has_all)
    }

    pub(crate) fn prepare_runtime_backends_after_generation_with_compile_times(
        &mut self,
    ) -> Result<Vec<Duration>> {
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        match self.data.compilation {
            FrozenCompilationMode::Symjit => {
                let mut compile_times = Vec::with_capacity(self.data.graph_terms.len());
                for graph_term in &mut self.data.graph_terms {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let compile_started = Instant::now();
                    graph_term
                        .for_each_generic_evaluator_mut(|evaluator| evaluator.activate_symjit())?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    compile_times.push(compile_started.elapsed());
                }
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(compile_times)
            }
            FrozenCompilationMode::Eager
            | FrozenCompilationMode::Cpp(_)
            | FrozenCompilationMode::Assembly(_) => {
                self.for_each_generic_evaluator_mut(|evaluator| {
                    evaluator.activate_eager();
                    Ok(())
                })?;
                self.active_f64_backend.set(ActiveF64Backend::Eager);
                Ok(vec![Duration::ZERO; self.data.graph_terms.len()])
            }
        }
    }

    pub(crate) fn prepare_runtime_backends_after_generation(&mut self) -> Result<()> {
        let _ = self.prepare_runtime_backends_after_generation_with_compile_times()?;
        Ok(())
    }

    pub(crate) fn activate_runtime_backends_after_load(
        &mut self,
        allow_symjit_fallback: bool,
    ) -> Result<Option<String>> {
        match self.data.compilation.clone() {
            FrozenCompilationMode::Eager => {
                self.prepare_runtime_backends_after_generation()?;
                Ok(None)
            }
            FrozenCompilationMode::Symjit => {
                self.for_each_generic_evaluator_mut(|evaluator| evaluator.activate_symjit())?;
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(None)
            }
            FrozenCompilationMode::Cpp(_) => {
                self.activate_external_after_load(ActiveF64Backend::Cpp, allow_symjit_fallback)
            }
            FrozenCompilationMode::Assembly(_) => {
                self.activate_external_after_load(ActiveF64Backend::Assembly, allow_symjit_fallback)
            }
        }
    }

    fn activate_external_after_load(
        &mut self,
        backend: ActiveF64Backend,
        allow_symjit_fallback: bool,
    ) -> Result<Option<String>> {
        if !self.has_complete_external_artifacts()? {
            self.prepare_runtime_backends_after_generation()?;
            return Ok(None);
        }

        match self.for_each_generic_evaluator_mut(|evaluator| {
            evaluator.activate_external_from_artifact(backend)
        }) {
            Ok(()) => {
                self.active_f64_backend.set(backend);
                Ok(None)
            }
            Err(err) if allow_symjit_fallback => {
                let error_message = err.to_string();
                self.for_each_generic_evaluator_mut(|evaluator| evaluator.activate_symjit())?;
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(Some(error_message))
            }
            Err(err) => Err(err),
        }
    }

    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self.data, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;

        self.settings
            .to_file(path.as_ref().join("settings.toml"), override_existing)
            .with_context(|| "Error saving settings.toml file for amplitude integrand")?;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _): (CrossSectionIntegrandData, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let settings = SmartSerde::from_file(
            path.as_ref().join("settings.toml"),
            "runtime settings for amplitude integrand",
        )?;

        Ok(CrossSectionIntegrand {
            settings,
            data,
            event_processing_runtime: RuntimeCache::default(),
            active_f64_backend: RuntimeCache::default(),
        })
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path> + Sync,
        override_existing: bool,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<(String, Duration)>> {
        let frozen_mode = self.data.compilation.clone();
        let compile_times = thread_pool.install(|| {
            self.data
                .graph_terms
                .par_iter_mut()
                .map(|term| {
                    term.compile(path.as_ref(), override_existing, &frozen_mode)
                        .map(|duration| (term.graph.name.clone(), duration))
                })
                .collect::<Result<Vec<_>>>()
        })?;

        self.active_f64_backend
            .set(ActiveF64Backend::from_frozen_mode(&self.data.compilation));
        Ok(compile_times)
    }

    pub(crate) fn invalidate_event_processing_runtime(&mut self) {
        self.event_processing_runtime.invalidate();
    }
}

impl ProcessIntegrandImpl for CrossSectionIntegrand {
    type G = CrossSectionGraphTerm;

    fn external_cache_id(&self) -> usize {
        self.data.external_cache_id
    }

    fn increment_external_cache_id(&mut self, val: usize) {
        self.data.external_cache_id += val
    }

    fn signal_external_momenta_changed(&mut self) {
        self.increment_external_cache_id(1);
        // Update base cache ID when the fundamental configuration changes
        self.data.base_external_cache_id = self.data.external_cache_id;
    }

    fn get_current_external_cache_id(&self) -> usize {
        self.external_cache_id()
    }

    /// Revert to the base external cache ID for the current configuration
    fn revert_to_base_external_cache_id(&mut self) {
        self.data.external_cache_id = self.data.base_external_cache_id;
    }

    fn get_base_external_cache_id(&self) -> usize {
        self.data.base_external_cache_id
    }

    fn increment_loop_cache_id(&mut self, val: usize) {
        self.data.loop_cache_id += val
    }

    fn loop_cache_id(&self) -> usize {
        self.data.loop_cache_id
    }
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.as_ref().expect("forgot warmup").iter()
    }

    fn get_group_structure(&self) -> &TiVec<GroupId, GraphGroup> {
        &self.data.graph_group_structure
    }

    fn warm_up(&mut self, model: &Model) -> Result<()> {
        self.data.rotations = Some(
            Some(Rotation::new(RotationMethod::Identity))
                .into_iter()
                .chain(
                    self.settings
                        .stability
                        .rotation_axis
                        .iter()
                        .map(|axis| Rotation::new(axis.rotation_method())),
                )
                .collect(),
        );

        for a in self.data.graph_terms.iter_mut() {
            a.warm_up(&self.settings, model)?;
        }
        self.event_processing_runtime.set(
            EventProcessingRuntime::from_settings_with_model_and_process_info(
                &self.settings,
                model,
                &histogram_process_info_for_integrand(self),
            )?,
        );
        Ok(())
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
    }

    fn graph_count(&self) -> usize {
        self.data.graph_terms.len()
    }

    fn get_group_masters(&self) -> impl Iterator<Item = &Self::G> {
        self.data
            .graph_group_structure
            .iter()
            .map(|group| &self.data.graph_terms[group.master()])
    }

    fn get_settings(&self) -> &RuntimeSettings {
        &self.settings
    }

    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G {
        &mut self.data.graph_terms[graph_id]
    }

    fn graph_group_id_for_graph(&self, graph_id: usize) -> Option<usize> {
        self.data.graph_to_group_id.get(graph_id).copied()
    }

    fn get_graph(&self, graph_id: usize) -> &Self::G {
        &self.data.graph_terms[graph_id]
    }

    fn get_master_graph(&self, group_id: GroupId) -> &Self::G {
        let group_master = self.data.graph_group_structure[group_id].master();

        &self.data.graph_terms[group_master]
    }

    fn get_group(&self, group_id: GroupId) -> &crate::graph::GraphGroup {
        &self.data.graph_group_structure[group_id]
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor<'_> {
        DependentMomentaConstructor::CrossSection
    }

    fn take_event_processing_runtime(&mut self) -> Option<EventProcessingRuntime> {
        self.event_processing_runtime.take()
    }

    fn restore_event_processing_runtime(&mut self, runtime: Option<EventProcessingRuntime>) {
        if let Some(runtime) = runtime {
            self.event_processing_runtime.set(runtime);
        }
    }

    fn event_processing_runtime(&self) -> Option<&EventProcessingRuntime> {
        self.event_processing_runtime.as_ref()
    }

    fn event_processing_runtime_mut(&mut self) -> Option<&mut EventProcessingRuntime> {
        self.event_processing_runtime.as_mut()
    }

    // fn get_builder_cache(&self) -> &ParamBuilder<f64> {
    //     &self.data.builder_cache
    // }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionGraphTerm {
    pub integrand: TiVec<RaisedCutId, BTreeMap<CutCFFIndex, EvaluatorStack>>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub reversed_edges: TiVec<RaisedCutId, Vec<EdgeIndex>>,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder<f64>,
    pub orientations: TiVec<OrientationID, EdgeVec<Orientation>>,
    pub orientation_filter: SubSet<OrientationID>,
    #[allow(private_interfaces)]
    pub counterterm: LUCounterTerm,
    pub raised_data: RaisedCutData,
}

impl CrossSectionGraphTerm {
    pub fn threshold_esurface_ids_for_raised_cut(
        &self,
        raised_cut_id: RaisedCutId,
    ) -> (Vec<usize>, Vec<usize>) {
        if self.counterterm.thresholds.is_empty() {
            return (vec![], vec![]);
        }

        let (left_thresholds, right_thresholds) = &self.counterterm.thresholds[raised_cut_id];
        let resolve_ids = |thresholds: &[Esurface]| {
            thresholds
                .iter()
                .map(|threshold| {
                    self.graph
                        .surface_cache
                        .esurface_cache
                        .iter()
                        .position(|candidate| candidate == threshold)
                        .expect("threshold esurface should resolve in graph")
                })
                .collect::<Vec<_>>()
        };

        (
            resolve_ids(left_thresholds.raw.as_slice()),
            resolve_ids(right_thresholds.raw.as_slice()),
        )
    }

    pub fn from_cross_section_graph(
        graph: &CrossSectionGraph,
        settings: &GlobalSettings,
    ) -> Result<(Self, GraphGenerationStats)> {
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let mut stats = GraphGenerationStats::default();
        let selected_generation_orientations = graph
            .derived_data
            .global_cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .filter(|orientation| {
                settings.generation.orientation_pattern.filter(*orientation)
                    && orientation.expression.iter_nodes().any(|tree_node| {
                        graph.cut_esurface_id_map.iter().any(|cut_esurface_id| {
                            tree_node.data == HybridSurfaceID::Esurface(*cut_esurface_id)
                        })
                    })
            })
            .collect_vec();
        let orientations: TiVec<OrientationID, EdgeVec<Orientation>> =
            selected_generation_orientations
                .iter()
                .map(|data| data.orientation().clone())
                .collect();
        if orientations.is_empty() {
            let pattern = settings
                .generation
                .orientation_pattern
                .pat
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "<empty>".to_string());
            return Err(eyre!(
                "Generation orientation pattern {pattern} matched no orientations for graph {}",
                graph.graph.name
            ));
        }

        let selected_generation_esurfaces = selected_generation_orientations
            .iter()
            .flat_map(|orientation| {
                orientation.expression.iter_nodes().filter_map(|tree_node| {
                    if let HybridSurfaceID::Esurface(esurface_id) = tree_node.data {
                        Some(esurface_id)
                    } else {
                        None
                    }
                })
            })
            .collect::<HashSet<_>>();

        let active_cuts: TiVec<RaisedCutId, bool> = graph
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .map(|raised_cut_group| {
                raised_cut_group
                    .related_esurface_group
                    .esurface_ids
                    .iter()
                    .any(|esurface_id| selected_generation_esurfaces.contains(esurface_id))
            })
            .collect();

        let masked_cut_parametric_integrand: TiVec<RaisedCutId, _> = graph
            .derived_data
            .cut_paramatric_integrand
            .iter_enumerated()
            .map(|(raised_cut_id, integrands)| {
                if active_cuts[raised_cut_id] {
                    integrands.clone()
                } else {
                    integrands.zero_like()
                }
            })
            .collect();

        let mut active_left_thresholds: TiVec<_, TiVec<_, bool>> = TiVec::new();
        let mut active_right_thresholds: TiVec<_, TiVec<_, bool>> = TiVec::new();
        let mut active_iterated_thresholds: TiVec<_, _> = TiVec::new();
        let mut masked_threshold_counterterms: TiVec<RaisedCutId, _> = TiVec::new();
        for (raised_cut_id, counterterm_data) in
            graph.derived_data.threshold_counterterms.iter_enumerated()
        {
            let left_active: TiVec<_, bool> = counterterm_data
                .left_thresholds
                .iter()
                .map(|esurface_id| {
                    active_cuts[raised_cut_id]
                        && selected_generation_esurfaces.contains(esurface_id)
                })
                .collect();
            let right_active: TiVec<_, bool> = counterterm_data
                .right_thresholds
                .iter()
                .map(|esurface_id| {
                    active_cuts[raised_cut_id]
                        && selected_generation_esurfaces.contains(esurface_id)
                })
                .collect();
            let mut iterated_active = counterterm_data.iterated.map_ref(|_| false);
            for (left_id, _) in counterterm_data.left_thresholds.iter_enumerated() {
                for (right_id, _) in counterterm_data.right_thresholds.iter_enumerated() {
                    iterated_active[(left_id, right_id)] =
                        left_active[left_id] && right_active[right_id];
                }
            }

            let mut masked_counterterm_data = counterterm_data.clone();
            for (left_id, integrands) in masked_counterterm_data.left_atoms.iter_mut_enumerated() {
                if !left_active[left_id] {
                    *integrands = integrands.zero_like();
                }
            }
            for (right_id, integrands) in masked_counterterm_data.right_atoms.iter_mut_enumerated()
            {
                if !right_active[right_id] {
                    *integrands = integrands.zero_like();
                }
            }
            for (integrands, is_active) in masked_counterterm_data
                .iterated
                .iter_mut()
                .zip(iterated_active.iter())
            {
                if !*is_active {
                    *integrands = integrands.zero_like();
                }
            }

            active_left_thresholds.push(left_active);
            active_right_thresholds.push(right_active);
            active_iterated_thresholds.push(iterated_active);
            masked_threshold_counterterms.push(masked_counterterm_data);
        }

        let mut integrand = TiVec::new();
        for integrand_for_cut in &masked_cut_parametric_integrand {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            let mut cut_integrands = BTreeMap::new();
            for (cut_cff_index, integrand_for_subset) in integrand_for_cut.integrands.iter() {
                let num_derivatives = cut_cff_index.lu_cut_order.unwrap_or(0);
                if crate::is_interrupted() {
                    return Err(eyre!("Generation interrupted by user"));
                }
                let dual_shape = if num_derivatives > 0 {
                    Some(graph.derived_data.raised_data.dual_shapes[num_derivatives - 1].clone())
                } else {
                    None
                };

                let (evaluator_stack, evaluator_timings) = EvaluatorStack::new_with_timings(
                    slice::from_ref(integrand_for_subset),
                    &graph.graph.param_builder,
                    &orientations.raw,
                    dual_shape,
                    &settings.generation.evaluator,
                )
                .with_context(|| {
                    format!(
                        "Failed to create evaluator for graph{}",
                        graph.graph.debug_dot()
                    )
                })?;
                if crate::is_interrupted() {
                    return Err(eyre!("Generation interrupted by user"));
                }
                stats.add_evaluator_build_timings(evaluator_timings);
                stats.evaluator_count += evaluator_stack.generic_evaluator_count();
                cut_integrands.insert(*cut_cff_index, evaluator_stack);
            }
            integrand.push(cut_integrands);
        }

        let mut ct_evaluators = TiVec::new();
        for (raised_cut_id, ct_data) in masked_threshold_counterterms.iter_enumerated() {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            let (evaluators, evaluator_timings) = LUCounterTermEvaluators::from_atoms(
                ct_data,
                graph.derived_data.raised_data.raised_cut_groups[raised_cut_id]
                    .related_esurface_group
                    .max_occurence,
                &graph.graph.param_builder,
                settings,
                &orientations,
            );
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            stats.add_evaluator_build_timings(evaluator_timings);
            stats.evaluator_count += evaluators.generic_compileable_evaluator_count();
            ct_evaluators.push(evaluators);
        }

        let mut thresholds = TiVec::new();
        for ct_data in &graph.derived_data.threshold_counterterms {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            thresholds.push((
                ct_data
                    .left_thresholds
                    .iter()
                    .map(|esurface_id| {
                        graph.graph.surface_cache.esurface_cache[*esurface_id].clone()
                    })
                    .collect(),
                ct_data
                    .right_thresholds
                    .iter()
                    .map(|esurface_id| {
                        graph.graph.surface_cache.esurface_cache[*esurface_id].clone()
                    })
                    .collect(),
            ));
        }

        let rstar_dependence_calculator = graph
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .map(|raised_cut_group| {
                generate_rstar_t_dependence_evaluator(
                    raised_cut_group.related_esurface_group.max_occurence - 1,
                )
            })
            .collect::<Result<TiVec<RaisedCutId, _>>>()?;

        let counterterm = LUCounterTerm {
            evaluators: ct_evaluators,
            thresholds,
            subspaces: graph.derived_data.subspace_data.clone(),
            rstar_dependence_calculator,
            active_cuts,
            active_left_thresholds,
            active_right_thresholds,
            active_iterated_thresholds,
        };

        let reversed_edges = graph
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .map(|cut_group| {
                let mut reversed_edges = HashSet::new();

                cut_group.cuts.iter().for_each(|cut_id| {
                    let cut = &graph.cuts[*cut_id];
                    cut.cut
                        .iter_edges(&graph.graph)
                        .for_each(|(orientation, edge_data)| {
                            if orientation == Orientation::Reversed {
                                reversed_edges.insert(
                                    graph
                                        .graph
                                        .edge_name_to_index(&edge_data.data.name)
                                        .unwrap(),
                                );
                            }
                        });
                });

                reversed_edges.into_iter().sorted().collect()
            })
            .collect();

        Ok((
            Self {
                integrand,
                graph: graph.graph.clone(),
                cut_esurface: graph.cut_esurface.clone(),
                cuts: graph.cuts.clone(),
                multi_channeling_setup: LmbMultiChannelingSetup {
                    channels: TiVec::new(),
                    graph: graph.graph.clone(), // will be overwritten later,
                    all_bases: TiVec::new(),
                },
                lmbs: graph.derived_data.lmbs.as_ref().unwrap().clone(),
                estimated_scale: None,
                param_builder: graph.graph.param_builder.clone(),
                orientation_filter: SubSet::full(orientations.len()),
                orientations,
                counterterm,
                reversed_edges,
                raised_data: graph.derived_data.raised_data.clone(),
            },
            stats,
        ))
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        _override_existing: bool,
        frozen_mode: &FrozenCompilationMode,
    ) -> Result<Duration> {
        let compile_started = Instant::now();
        let graph_path = path.as_ref().join(&self.graph.name);

        fs::create_dir_all(&graph_path).with_context(|| {
            format!(
                "Failed to create directory for cross section graph {} at {}",
                self.graph.name,
                graph_path.display()
            )
        })?;

        for (raised_cut_id, integrands) in self.integrand.iter_mut().enumerate() {
            for (cut_cff_index, integrand) in integrands.iter_mut() {
                let n_derivatives = cut_cff_index.lu_cut_order.unwrap_or(0);
                integrand.compile(
                    format!(
                        "integrand_zen_cut_{}_deriv_{}",
                        raised_cut_id, n_derivatives
                    ),
                    graph_path.clone(),
                    frozen_mode,
                )?;
            }
        }

        self.counterterm.compile(&graph_path, frozen_mode)?;

        for (index, evaluator) in self.raised_data.pass_two_evaluators.iter_mut().enumerate() {
            evaluator.compile_external(
                graph_path
                    .join(format!("pass_two_{index}"))
                    .with_extension("cpp"),
                format!("pass_two_{index}"),
                graph_path
                    .join(format!("pass_two_{index}"))
                    .with_extension("so"),
                frozen_mode,
            )?;
        }

        Ok(compile_started.elapsed())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        for raised_cut_integrands in self.integrand.iter_mut() {
            for (_cut_cff_index, evaluator_stack) in raised_cut_integrands.iter_mut() {
                evaluator_stack.for_each_generic_evaluator_mut(&mut f)?;
            }
        }

        self.counterterm.for_each_generic_evaluator_mut(&mut f)?;

        for evaluator in self.raised_data.pass_two_evaluators.iter_mut() {
            f(evaluator)?;
        }

        Ok(())
    }

    fn generate_event_for_cut<T: FloatLike>(
        &self,
        model: &Model,
        t_scaling_solution: &NewtonIterationResult<T>,
        momentum_sample: &MomentumSample<T>,
        cut_id: CutId,
        cut: &CrossSectionCut,
        channel_id: Option<ChannelIndex>,
    ) -> Result<GenericEvent<T>> {
        let rescaled_momenta =
            momentum_sample.rescaled_loop_momenta(&t_scaling_solution.solution, Subspace::None);

        let mut new_event = GenericEvent::<T>::default();
        new_event.cut_info.cut_id = cut_id.0;
        new_event.cut_info.orientation_id = momentum_sample.sample.orientation;
        new_event.cut_info.lmb_channel_id = channel_id.map(usize::from);
        new_event.cut_info.lmb_channel_edge_ids =
            channel_id.map(|channel_id| self.multi_channeling_setup.channel_edge_ids(channel_id));
        // Set initial momenta and PDGs for the event
        new_event
            .kinematic_configuration
            .0
            .extend(rescaled_momenta.external_moms().clone());

        let mut incoming_pdgs = self
            .graph
            .underlying
            .iter_edges_of(&self.graph.initial_state_cut)
            .map(|(_, edge_index, edge_data)| {
                edge_data
                    .data
                    .particle()
                    .map(|particle| (edge_index, particle.pdg_code))
                    .ok_or_else(|| {
                        eyre!(
                            "Initial-state cut edge {edge_index:?} in graph {} has no particle specifier",
                            self.graph.name
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        incoming_pdgs.sort_by_key(|(edge_index, _)| edge_index.0);
        new_event
            .cut_info
            .particle_pdgs
            .0
            .extend(incoming_pdgs.into_iter().map(|(_, pdg)| pdg));

        let initial_state_cut_edges = self.graph.get_edges_in_initial_state_cut();

        for (p, eid, d) in self.graph.iter_edges_of(&cut.cut) {
            if initial_state_cut_edges.contains(&eid) {
                continue;
            }

            let cut_flow = match p {
                HedgePair::Split {
                    source: _,
                    sink: _,
                    split,
                } => split,
                HedgePair::Unpaired { hedge: _, flow } => flow,
                HedgePair::Paired { .. } => {
                    return Err(eyre!(
                        "Found paired edge {eid:?} while building an event for cut {:?} in graph {}",
                        cut.cut,
                        self.graph.name
                    ));
                }
            };

            let mut edge_spatial_momentum = self.graph.loop_momentum_basis.edge_signatures[eid]
                .compute_three_momentum_from_four(
                    rescaled_momenta.loop_moms(),
                    momentum_sample.external_moms(),
                );

            let edge_pdg = d.data.particle().map(|p| p.pdg_code).ok_or_else(|| {
                eyre!("Cut legs in Local Unitarity must have a particle specifier.")
            })?;

            let cut_pdg = match cut_flow {
                Flow::Source => edge_pdg,
                Flow::Sink => {
                    edge_spatial_momentum = -edge_spatial_momentum;
                    model
                        .get_particle_from_pdg(edge_pdg)
                        .get_anti_particle(model)
                        .pdg_code
                }
            };

            let mass_value = if let Some(mass) = d.data.mass.value(model, &self.param_builder) {
                if !mass.im.is_zero() {
                    return Err(eyre!(
                        "Cut particles should have real-valued masses ({})",
                        edge_pdg
                    ));
                }
                Some(mass.re)
            } else {
                None
            };

            let cut_four_momentum = edge_spatial_momentum.into_on_shell_four_momentum(mass_value);

            debug!(
                "event cut leg: edge={eid:?} name={} pair={} graph_orientation={:?} cut_flow={:?} pdg={} p={}",
                d.data.name, p, d.orientation, cut_flow, cut_pdg, cut_four_momentum,
            );

            new_event.kinematic_configuration.1.push(cut_four_momentum);
            new_event.cut_info.particle_pdgs.1.push(cut_pdg);
        }

        Ok(new_event)
    }
}

impl GraphTerm for CrossSectionGraphTerm {
    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64> {
        &mut self.param_builder
    }

    fn get_num_channels(&self) -> usize {
        self.multi_channeling_setup.channels.len()
    }

    fn orientation_label(&self, orientation_id: usize) -> Option<String> {
        self.orientations
            .get(resolve_visible_orientation_id(
                &self.orientation_filter,
                orientation_id,
            )?)
            .map(format_orientation_label)
    }

    fn lmb_channel_label(&self, channel_id: ChannelIndex) -> Option<String> {
        Some(format_lmb_channel_label(
            &self.multi_channeling_setup.channel_edge_ids(channel_id),
        ))
    }

    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()> {
        self.estimated_scale = Some(
            self.graph
                .expected_scale(F(settings.kinematics.e_cm), model),
        );

        self.orientation_filter = SubSet::empty(self.orientations.len());

        for (i, or) in self.orientations.iter_enumerated() {
            if settings.general.orientation_pat.filter(or) {
                self.orientation_filter.add(i);
            }
        }
        if self.orientation_filter.included_iter().next().is_none() {
            let pattern = settings
                .general
                .orientation_pat
                .pat
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "<empty>".to_string());
            return Err(eyre!(
                "Runtime orientation pattern {pattern} matched no orientations for graph {}",
                self.graph.name
            ));
        }

        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::CrossSection)
            .with_context(|| {
                format!(
                    "Failed to get dependent external momenta for graph {}",
                    self.graph.name
                )
            })?;

        if externals.len() != self.graph.loop_momentum_basis.ext_edges.len() {
            return Err(eyre!(
                "Number Externals supplied in the settings {} do not match number of externals {} in graph {}",
                externals.len(),
                self.graph.loop_momentum_basis.ext_edges.len(),
                self.graph.name
            ));
        }
        self.graph
            .param_builder
            .add_external_four_mom_all_derivatives(&externals);

        let pols = self.graph.param_builder.pairs.polarizations_values(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );
        self.graph.param_builder.pairs.warn_zero_polarizations(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );

        for (value_index, values) in self.graph.param_builder.values.iter_mut().enumerate() {
            let multiplicative_offset = value_index + 1;
            let mut pol_start = self
                .graph
                .param_builder
                .pairs
                .polarizations
                .value_range
                .start
                * multiplicative_offset;

            for pol in &pols {
                values[pol_start] = *pol;
                pol_start += multiplicative_offset;
            }
        }

        self.graph
            .param_builder
            .m_uv_value(Complex::new_re(F(settings.general.m_uv)));
        self.graph
            .param_builder
            .mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq())));
        self.graph.param_builder.update_model_values(model);

        self.param_builder = self.graph.param_builder.clone();

        Ok(())
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        mut context: GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<GraphEvaluationResult<T>> {
        let orientations =
            momentum_sample.orientations(&self.orientation_filter, &self.orientations);

        // let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        let center = LoopMomenta::from_iter(vec![
            ThreeMomentum {
                px: momentum_sample.zero(),
                py: momentum_sample.zero(),
                pz: momentum_sample.zero(),
            };
            momentum_sample.loop_moms().len()
        ]);
        let masses = self.graph.get_real_mass_vector(context.model);
        let hel = context.settings.kinematics.externals.get_helicities();
        let mut cut_results: TiVec<RaisedCutId, Vec<Complex<F<T>>>> =
            ti_vec![Vec::new(); self.raised_data.raised_cut_groups.len()];
        let mut cut_threshold_counterterms = TiVec::<RaisedCutId, Complex<F<T>>>::new();
        let mut differential_result = GraphEvaluationResult::zero(momentum_sample.zero());
        let mut accepted_event_group = GenericEventGroup::default();

        let momentum_sample = if let Some((channel_id, _alpha)) = &context.channel_id {
            MomentumSample {
                sample: self.multi_channeling_setup.reinterpret_loop_momenta_impl(
                    *channel_id,
                    &momentum_sample.sample,
                    momentum_sample.sample.loop_mom_cache_id,
                ),
            }
        } else {
            momentum_sample.clone()
        };

        crate::debug_tags!(#integration, #sample, #inspect;
            "loop moms: {}",
            momentum_sample.loop_moms()
        );

        for (raised_cut, raised_cut_group) in self.raised_data.raised_cut_groups.iter_enumerated() {
            let max_occurance = raised_cut_group.related_esurface_group.max_occurence;
            if !self.counterterm.cut_is_active(raised_cut) {
                let zero = Complex::new_re(momentum_sample.zero());
                for _ in 1..=max_occurance {
                    cut_results[raised_cut].push(zero.clone());
                }
                cut_threshold_counterterms.push(zero);
                continue;
            }
            crate::debug_tags!(#integration, #cut;
                "\n =====START EVALUTAION FOR CUT {}=====",
                raised_cut.0
            );
            let representative_esurface = &self.cut_esurface[raised_cut_group.cuts[0]];

            crate::debug_tags!(#integration, #cut, #inspect;
                "representative esurface: {:#?}",
                representative_esurface
            );

            let function = |t: &F<T>| {
                representative_esurface.compute_self_and_r_derivative(
                    t,
                    momentum_sample.loop_moms(),
                    &center,
                    momentum_sample.external_moms(),
                    &masses,
                    &self.graph.loop_momentum_basis,
                )
            };

            let (guess, _) = representative_esurface.get_radius_guess(
                momentum_sample.loop_moms(),
                momentum_sample.external_moms(),
                &self.graph.loop_momentum_basis,
            );

            let solution = newton_iteration_and_derivative(
                &guess,
                function,
                &F::from_f64(1.0),
                2000,
                &F::from_f64(context.settings.kinematics.e_cm),
            );

            crate::debug_tags!(#integration, #cut, #solver;
                "tolerance for newton solver: {}",
                F::from_f64(context.settings.kinematics.e_cm) * guess.epsilon()
            );

            crate::debug_tags!(#integration, #cut, #solver;
                "solution: {:?}",
                solution
            );

            let prepared_event = prepare_buffered_event(
                context.settings,
                context.rotation,
                context.event_processing_runtime.as_deref_mut(),
                || {
                    let mut generated = self.generate_event_for_cut::<T>(
                        context.model,
                        &solution,
                        &momentum_sample,
                        self.raised_data.raised_cut_groups[raised_cut].cuts[0],
                        &self.cuts[self.raised_data.raised_cut_groups[raised_cut].cuts[0]],
                        context
                            .channel_id
                            .as_ref()
                            .map(|(channel_id, _)| *channel_id),
                    )?;
                    generated.inverse_rotate(context.rotation);
                    Ok(generated)
                },
            )?;
            differential_result.generated_event_count += prepared_event.generated_event_count;
            differential_result.accepted_event_count += prepared_event.accepted_event_count;
            differential_result.event_processing_time += prepared_event.event_processing_time;

            if !prepared_event.selectors_pass {
                let zero = Complex::new_re(momentum_sample.zero());
                for _ in 1..=max_occurance {
                    cut_results[raised_cut].push(zero.clone());
                }
                cut_threshold_counterterms.push(zero);
                continue;
            }

            let accepted_event = prepared_event.buffered_event;
            let mut bare_cut_total = Complex::new_re(momentum_sample.zero());
            let mut threshold_counterterm_weights = Vec::with_capacity(max_occurance);
            let mut kinematic_point = LUCTKinematicPoint::new(momentum_sample.clone());
            for num_esurfaces in 1..=max_occurance {
                let dual_shape = if num_esurfaces > 1 {
                    Some(HyperDual::<F<T>>::new(
                        self.raised_data.dual_shapes[num_esurfaces - 2].clone(),
                    ))
                } else {
                    None
                };

                let (tstar, h_function, esurface_derivatives, rescaled_momenta) =
                    if let Some(dual_shape) = dual_shape {
                        let dual_t_for_integrand =
                            dual_shape.variable(0, solution.solution.clone());
                        let dual_h_function = h_dual(
                            &dual_t_for_integrand,
                            None,
                            None,
                            &context.settings.lu_h_function,
                        );
                        let dual_momenta_for_integrand = momentum_sample
                            .loop_moms()
                            .rescale_with_hyper_dual(&dual_t_for_integrand, None);

                        let dual_shape_for_esurface =
                            HyperDual::<F<T>>::new(simple_n_deriv_shape(num_esurfaces));

                        let dual_t_for_esurface =
                            dual_shape_for_esurface.variable(0, solution.solution.clone());
                        let dual_momenta_for_esurface = momentum_sample
                            .loop_moms()
                            .rescale_with_hyper_dual(&dual_t_for_esurface, None);

                        let dual_externals = momentum_sample
                            .external_moms()
                            .iter()
                            .map(|mom| FourMomentum {
                                temporal: Energy {
                                    value: new_constant(&dual_t_for_esurface, &mom.temporal.value),
                                },
                                spatial: ThreeMomentum {
                                    px: new_constant(&dual_t_for_esurface, &mom.spatial.px),
                                    py: new_constant(&dual_t_for_esurface, &mom.spatial.py),
                                    pz: new_constant(&dual_t_for_esurface, &mom.spatial.pz),
                                },
                            })
                            .collect();

                        let dual_e_surface = representative_esurface.compute_from_dual_momenta(
                            &self.graph.loop_momentum_basis,
                            &masses,
                            &dual_momenta_for_esurface,
                            &dual_externals,
                        );

                        let mut momentum_sample_with_duals = momentum_sample.clone();
                        momentum_sample_with_duals.sample.dual_loop_moms =
                            Some(dual_momenta_for_integrand);

                        (
                            DualOrNot::Dual(dual_t_for_integrand),
                            DualOrNot::Dual(dual_h_function),
                            DualOrNot::Dual(dual_e_surface),
                            momentum_sample_with_duals,
                        )
                    } else {
                        let h_function = h(
                            &solution.solution,
                            None,
                            None,
                            &context.settings.lu_h_function,
                        );
                        let rescaled_momenta = momentum_sample
                            .rescaled_loop_momenta(&solution.solution, Subspace::None);

                        (
                            DualOrNot::NonDual(solution.solution.clone()),
                            DualOrNot::NonDual(h_function),
                            DualOrNot::NonDual(solution.derivative_at_solution.clone()),
                            rescaled_momenta,
                        )
                    };

                debug!("tstar: {}", tstar);
                debug!("h(tstar): {}", h_function);
                debug!("esurface derivative at tstar: {}", esurface_derivatives);

                let lu_params = LUParams { h_function, tstar };

                if !context.settings.subtraction.disable_threshold_subtraction {
                    kinematic_point
                        .dualized_momentum_sample_cache
                        .push(rescaled_momenta.clone());
                    kinematic_point
                        .lu_cut_parameter_cache
                        .push(lu_params.clone());
                    kinematic_point
                        .lu_cut_esurface_values
                        .push(esurface_derivatives.clone());
                }

                let prefactor = Complex::new_re(
                    if let Some((_channel_index, _alpha)) = &context.channel_id {
                        if matches!(lu_params.tstar, DualOrNot::Dual(_)) {
                            panic!("multi channeling with duals not supported yet");
                        }

                        self.multi_channeling_setup.compute_prefactor_impl(
                            *_channel_index,
                            &rescaled_momenta,
                            context.model,
                            _alpha,
                        )
                    } else {
                        F::from_f64(1.0)
                    },
                );

                let params = T::get_parameters(
                    &mut self.param_builder,
                    (
                        context.settings.general.enable_cache,
                        context.settings.general.debug_cache,
                    ),
                    &self.graph,
                    &rescaled_momenta,
                    hel,
                    &context.settings.additional_params(),
                    None,
                    None,
                    Some(&lu_params),
                );

                let cut_index = CutCFFIndex {
                    lu_cut_order: Some(num_esurfaces),
                    left_threshold_order: None,
                    right_threshold_order: None,
                };

                let result = self.integrand[raised_cut]
                    .get_mut(&cut_index)
                    .unwrap()
                    .evaluate(
                        params,
                        orientations,
                        context.settings,
                        context.evaluation_metadata,
                        context.record_primary_timing,
                    )
                    .expect("evaluation failed")
                    .pop()
                    .unwrap();

                debug!("pass 1 result {}", result);

                let mut params_for_pass_two = vec![];
                match result {
                    DualOrNot::Dual(dual_result) => {
                        params_for_pass_two
                            .extend_from_slice(&extract_t_derivatives_complex(dual_result));
                    }
                    DualOrNot::NonDual(non_dual_result) => {
                        params_for_pass_two.push(non_dual_result);
                    }
                }

                match esurface_derivatives {
                    DualOrNot::Dual(dual_e_surface) => {
                        extract_t_derivatives(dual_e_surface)[1..]
                            .iter()
                            .for_each(|v| {
                                params_for_pass_two.push(Complex::new_re(v.clone()));
                            });
                    }
                    DualOrNot::NonDual(non_dual_e_surface) => {
                        params_for_pass_two.push(Complex::new_re(non_dual_e_surface));
                    }
                }

                let pass_two_evaluator =
                    &mut self.raised_data.pass_two_evaluators[num_esurfaces - 1];

                let pass_two_result = evaluate_evaluator_single(
                    pass_two_evaluator,
                    &params_for_pass_two,
                    context.evaluation_metadata,
                    context.record_primary_timing,
                );

                debug!("pass_two_result: {:+16e}", pass_two_result);
                //debug!("param builder for cut {}: \n{}", cut, self.param_builder);

                let bare_contribution = pass_two_result * prefactor;
                bare_cut_total += bare_contribution.clone();
                cut_results[raised_cut].push(bare_contribution);
            }

            let ct_result = if context.settings.subtraction.disable_threshold_subtraction {
                Complex::new_re(momentum_sample.zero())
            } else {
                self.counterterm.evaluate(
                    &kinematic_point,
                    raised_cut,
                    &self.reversed_edges[raised_cut],
                    &self.lmbs,
                    &self.graph,
                    &self.graph.get_real_mass_vector(context.model),
                    context.rotation,
                    context.settings,
                    &mut self.param_builder,
                    orientations,
                    context.evaluation_metadata,
                    context.record_primary_timing,
                )?
            };

            threshold_counterterm_weights.push(ct_result.clone());
            cut_threshold_counterterms.push(ct_result.clone());

            if let Some(mut event) = accepted_event {
                let threshold_counterterm_total = threshold_counterterm_weights
                    .iter()
                    .fold(Complex::new_re(momentum_sample.zero()), |acc, value| {
                        acc + value.clone()
                    });
                event.weight = bare_cut_total.clone() + threshold_counterterm_total;

                if context.settings.general.store_additional_weights_in_event {
                    event
                        .additional_weights
                        .weights
                        .insert(AdditionalWeightKey::Original, bare_cut_total);
                    for (subset_index, threshold_counterterm) in
                        threshold_counterterm_weights.into_iter().enumerate()
                    {
                        event.additional_weights.weights.insert(
                            AdditionalWeightKey::ThresholdCounterterm { subset_index },
                            threshold_counterterm,
                        );
                    }
                }

                accepted_event_group.push(event);
            }
        }

        let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        for ((cut_id, result), ct_result) in cut_results
            .iter_enumerated()
            .zip(cut_threshold_counterterms.iter())
        {
            let mut total_bare_contribution = Complex::new_re(momentum_sample.zero());
            for (i, bare_contribution) in result.iter().enumerate() {
                debug!(
                    "cut {} contribution with {} esurfaces: {:+16e}",
                    cut_id.0,
                    i + 1,
                    bare_contribution
                );

                total_bare_contribution += bare_contribution;
            }
            debug!(
                "total bare contribution for cut {}: {:+16e}",
                cut_id.0, total_bare_contribution
            );

            all_cut_result += total_bare_contribution;

            debug!(
                "threshold counterterm for cut {}: {:+16e}",
                cut_id.0, ct_result
            );
            all_cut_result += ct_result;
        }

        let resolved_integral_unit = context
            .settings
            .general
            .integral_unit
            .resolve_for_cross_section(
                self.graph.initial_state_cut.iter_edges(&self.graph).count(),
            );
        let flux_factor = if context.settings.general.disable_flux_factor {
            F::from_f64(1.0)
        } else {
            match momentum_sample.external_moms().len() {
                1 => {
                    momentum_sample.one()
                        / (F::from_f64(2.0)
                            * &momentum_sample
                                .external_moms()
                                .first()
                                .as_ref()
                                .unwrap()
                                .temporal
                                .value)
                }
                2 => {
                    let mom_1 = &(momentum_sample.external_moms()[ExternalIndex::from(0)]);
                    let mom_2 = &(momentum_sample.external_moms()[ExternalIndex::from(1)]);
                    let mass_factor = self
                        .graph
                        .initial_state_cut
                        .iter_edges(&self.graph)
                        .map(|(_, e)| {
                            e.data
                                .mass
                                .value(context.model, &self.param_builder)
                                .unwrap()
                        })
                        .fold(Complex::new_re(momentum_sample.one()), |acc, mass| {
                            acc * &mass * &mass
                        })
                        .re;

                    let f = F::from_f64(4.0) * (mom_1.dot(mom_2).square() - mass_factor).sqrt();

                    momentum_sample.one() / f
                        * barn_conversion_factor(resolved_integral_unit, momentum_sample.one())
                }
                _ => unimplemented!(
                    "Flux factor for more than 3 or more incoming particles not implemented yet"
                ),
            }
        };

        let final_result = all_cut_result * flux_factor.clone();

        if context.settings.should_buffer_generated_events() {
            let flux_factor = Complex::new_re(flux_factor);
            for event in accepted_event_group.iter_mut() {
                event.weight *= flux_factor.clone();
                if !event.additional_weights.weights.is_empty() {
                    event.additional_weights.weights.insert(
                        AdditionalWeightKey::FullMultiplicativeFactor,
                        flux_factor.clone(),
                    );
                }
            }

            if !accepted_event_group.is_empty() {
                differential_result.event_groups.push(accepted_event_group);
            }
        }

        debug!(
            "{}",
            format!(
                "final result for graph: {}, {:+16e}",
                self.graph.name, final_result
            )
            .red()
        );

        differential_result.integrand_result = final_result;
        Ok(differential_result)
    }
    fn name(&self) -> String {
        self.graph.name.clone()
    }

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_num_orientations(&self) -> usize {
        filtered_orientation_count(&self.orientation_filter, &self.orientations)
    }

    fn get_tropical_sampler(&self) -> &momtrop::SampleGenerator<3> {
        unimplemented!(
            "Don't know how to generate subgraph table for forward scattering graphs yet"
        )
    }

    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>> {
        todo!()
    }
}

impl HasIntegrand for CrossSectionIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        create_grid(self)
    }

    fn name(&self) -> String {
        self.data.name.clone()
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        model: &Model,
        wgt: F<f64>,
        _iter: usize,
        use_f128: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        let result = evaluate_sample(self, model, sample, wgt, _iter, use_f128, max_eval);

        debug!(result = ?result,"Evaluating");

        result
    }

    fn get_n_dim(&self) -> usize {
        assert!(
            self.settings
                .sampling
                .get_parameterization_settings()
                .is_some(),
            "Tropical smapling not implemented for cross sections yet"
        );

        assert!(
            self.data
                .graph_terms
                .iter()
                .map(|term| term.graph.get_loop_number())
                .all_equal()
        );

        self.data.graph_terms[0].graph.get_loop_number() * 3
    }
}
