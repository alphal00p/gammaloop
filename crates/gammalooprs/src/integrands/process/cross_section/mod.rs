use crate::{
    DependentMomentaConstructor, GammaLoopContext, GammaLoopContextContainer,
    cff::{
        CutCFFIndex, esurface::Esurface, expression::OrientationID, orientations::GraphOrientation,
        surface::HybridSurfaceID,
    },
    graph::{
        ExternalConnection, FeynmanGraph, Graph, GraphGroup, GroupId, LmbIndex, LoopMomentumBasis,
        parse::complete_group_parsing,
    },
    integrands::{
        HasIntegrand,
        evaluation::{EvaluationResult, GraphEvaluationResult},
        process::{
            GraphTermEvaluationContext, ParamBuilder, SamplingChannelBridge,
            SamplingChannelCompileContext, SamplingChannelId, SamplingCutSide,
            evaluators::{ActiveF64Backend, EvaluatorStack, evaluate_evaluator_single},
            graph_to_group_id_for_group_structure,
            param_builder::LUParams,
            prepare_buffered_event,
            sampling_context::SamplingMapContext,
            sampling_maps::{
                ImplicitSurfaceRadialMap, SamplingEvaluationError, SamplingMapAffine,
                SamplingMapComponent, SamplingMapComposition, SamplingMapEmbedding,
            },
            sampling_selection::CompiledSamplingMap,
            threshold_multiplier::{
                ThresholdMultiplierEvaluatorCollection, ThresholdMultiplierExpression,
                ThresholdMultiplierLayout,
            },
        },
    },
    model::Model,
    momentum::{
        Energy, FourMomentum, Rotation, RotationMethod, ThreeMomentum,
        sample::{
            ExternalFourMomenta, ExternalIndex, LoopIndex, LoopMomenta, MomentumSample, Subspace,
            SubspaceData,
        },
    },
    observables::{
        AdditionalWeightKey, EventProcessingRuntime, GenericEvent, GenericEventGroup,
        GenericThresholdCountertermEventInfo,
    },
    processes::{
        self, CrossSectionCut, CrossSectionGraph, CutGroupData, CutGroupId, CutId,
        CutThresholdCountertermAssociations, GraphGenerationStats, GraphGroupSelectionPlan,
        IteratedCtCollection, LUCounterTermData, LUThresholdHelperOutputs, LeftThresholdId,
        RightThresholdId, ThresholdCountertermMetadataRegistry, ThresholdCountertermVariantId,
        ThresholdCountertermVariantStatus, TopologicalThresholdId,
    },
    settings::{
        GlobalSettings, RuntimeSettings,
        global::{CompilationOptimizationLevel, FrozenCompilationMode},
        runtime::{IntegralUnit, ParameterizationSettings},
    },
    subtraction::{
        generate_rstar_t_dependence_evaluator,
        lu_counterterm::{
            LUCTKinematicPoint, LUCounterTerm, LUCounterTermEvaluators, LUCountertermEvaluation,
            LUThresholdHelperEvaluators, LUVariantSubspaces,
        },
    },
    utils::{
        F, FloatLike, RuntimeCache, h, h_dual,
        hyperdual_utils::{
            DualOrNot, extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            shape_from_cut_cff_index, simple_n_deriv_shape,
        },
        newton_solver::{NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity},
        serde_utils::SmartSerde,
    },
};
use bincode::Encode;
use bincode_trait_derive::Decode;
use color_eyre::{Result, owo_colors::OwoColorize};
use eyre::Context;
use eyre::eyre;
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    slice,
    sync::Arc,
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
use tracing::{debug, warn};
use typed_index_collections::{TiVec, ti_vec};

use super::{
    GraphTerm, LmbMultiChannelingSetup, ProcessIntegrandImpl, create_grid, evaluate_sample,
    filtered_orientation_count, format_orientation_label, format_sampling_channel_label,
    histogram_process_info_for_integrand, resolve_visible_orientation_id,
    validate_group_orientation_catalogs, validate_process_runtime_settings,
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
    pub explicit_orientation_sum_only: bool,
    /// Frozen CP optimization choice; its validity after model updates is user supplied.
    pub symmetrize_left_right_states: bool,
    // pub builder_cache: ParamBuilder<f64>,
}

impl CrossSectionIntegrand {
    pub(crate) fn clone_with_graph_group_selection(
        &self,
        plan: &GraphGroupSelectionPlan,
    ) -> Result<Self> {
        let mut old_graph_to_new_group = vec![None; self.data.graph_terms.len()];
        let mut old_graph_is_master = vec![false; self.data.graph_terms.len()];
        for &old_group_id in plan.retained_group_ids() {
            let new_group_id = plan.new_group_id_for_old(old_group_id).ok_or_else(|| {
                eyre!(
                    "Graph-group selection is missing a compact id for cross-section group {}.",
                    old_group_id.0
                )
            })?;
            let group = self
                .data
                .graph_group_structure
                .get(old_group_id)
                .ok_or_else(|| {
                    eyre!(
                        "Graph-group selection refers to missing cross-section group {}.",
                        old_group_id.0
                    )
                })?;
            for old_graph_id in group {
                if old_graph_id >= old_graph_to_new_group.len() {
                    return Err(eyre!(
                        "Cross-section graph group {} refers to missing graph id {}.",
                        old_group_id.0,
                        old_graph_id
                    ));
                }
                old_graph_to_new_group[old_graph_id] = Some(new_group_id);
            }
            old_graph_is_master[group.master()] = true;
        }

        let mut graph_terms = self
            .data
            .graph_terms
            .iter()
            .enumerate()
            .filter_map(|(old_graph_id, source_graph_term)| {
                let new_group_id = old_graph_to_new_group[old_graph_id]?;
                let mut graph_term = source_graph_term.clone();
                graph_term.graph.group_id = Some(new_group_id);
                graph_term.graph.is_group_master = old_graph_is_master[old_graph_id];
                graph_term.multi_channeling_setup.graph.group_id = Some(new_group_id);
                graph_term.multi_channeling_setup.graph.is_group_master =
                    old_graph_is_master[old_graph_id];
                Some(graph_term)
            })
            .collect::<Vec<_>>();

        let mut parsed_graphs = graph_terms
            .iter()
            .map(|graph_term| graph_term.graph.clone())
            .collect::<Vec<_>>();
        let graph_group_structure = complete_group_parsing(&mut parsed_graphs)?;
        for (graph_term, parsed_graph) in graph_terms.iter_mut().zip(parsed_graphs) {
            graph_term.graph.group_id = parsed_graph.group_id;
            graph_term.graph.is_group_master = parsed_graph.is_group_master;
            graph_term.multi_channeling_setup.graph.group_id = parsed_graph.group_id;
            graph_term.multi_channeling_setup.graph.is_group_master = parsed_graph.is_group_master;
        }
        let graph_to_group_id = graph_to_group_id_for_group_structure(&graph_group_structure);

        Ok(Self {
            settings: self.settings.clone(),
            data: CrossSectionIntegrandData {
                name: self.data.name.clone(),
                compilation: self.data.compilation.clone(),
                loop_cache_id: self.data.loop_cache_id,
                external_cache_id: self.data.external_cache_id,
                base_external_cache_id: self.data.base_external_cache_id,
                rotations: self.data.rotations.clone(),
                graph_terms,
                n_incoming: self.data.n_incoming,
                external_connections: self.data.external_connections.clone(),
                graph_group_structure,
                graph_to_group_id,
                explicit_orientation_sum_only: self.data.explicit_orientation_sum_only,
                symmetrize_left_right_states: self.data.symmetrize_left_right_states,
            },
            event_processing_runtime: RuntimeCache::default(),
            active_f64_backend: self.active_f64_backend.clone(),
        })
    }

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
            FrozenCompilationMode::Symjit(optimization_level) => {
                let mut compile_times = Vec::with_capacity(self.data.graph_terms.len());
                for graph_term in &mut self.data.graph_terms {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let compile_started = Instant::now();
                    graph_term.for_each_generic_evaluator_mut(|evaluator| {
                        evaluator.activate_symjit(optimization_level)
                    })?;
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
            FrozenCompilationMode::Symjit(optimization_level) => {
                self.for_each_generic_evaluator_mut(|evaluator| {
                    evaluator.activate_symjit(optimization_level)
                })?;
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(None)
            }
            FrozenCompilationMode::Cpp(options) => self.activate_external_after_load(
                ActiveF64Backend::Cpp,
                options.optimization_level,
                allow_symjit_fallback,
            ),
            FrozenCompilationMode::Assembly(options) => self.activate_external_after_load(
                ActiveF64Backend::Assembly,
                options.optimization_level,
                allow_symjit_fallback,
            ),
        }
    }

    fn activate_external_after_load(
        &mut self,
        backend: ActiveF64Backend,
        optimization_level: CompilationOptimizationLevel,
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
                self.for_each_generic_evaluator_mut(|evaluator| {
                    evaluator.activate_symjit(optimization_level)
                })?;
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
            .with_context(|| "Error saving settings.toml file for cross-section integrand")?;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _): (CrossSectionIntegrandData, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let settings = SmartSerde::from_file(
            path.as_ref().join("settings.toml"),
            "runtime settings for cross-section integrand",
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

    pub(crate) fn invalidate_runtime_caches(&mut self) {
        self.event_processing_runtime.invalidate();
        for term in &mut self.data.graph_terms {
            term.multi_channeling_setup.invalidate_sampling();
        }
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
        self.invalidate_runtime_caches();
        if self.data.symmetrize_left_right_states {
            warn!(
                "This integrand was generated with symmetrize_left_right_states=true, which assumes CP symmetry. Complex couplings or model updates can invalidate that assumption; verifying it at the current parameter point is the user's responsibility"
            );
        }
        validate_process_runtime_settings(&self.settings, self.data.explicit_orientation_sum_only)?;

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
        validate_group_orientation_catalogs(
            &self.settings,
            &self.data.graph_terms,
            &self.data.graph_group_structure,
        )?;
        self.event_processing_runtime.set(
            EventProcessingRuntime::from_settings_with_model_and_process_info(
                &self.settings,
                model,
                &histogram_process_info_for_integrand(self)?,
            )?,
        );
        self.warm_up_sampling()
    }

    fn uses_explicit_orientation_sum_only(&self) -> bool {
        self.data.explicit_orientation_sum_only
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
    pub integrand: TiVec<CutGroupId, BTreeMap<CutCFFIndex, EvaluatorStack>>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub covariant_cut_representatives: BTreeMap<isize, isize>,
    pub topological_threshold_esurfaces: TiVec<TopologicalThresholdId, Esurface>,
    pub cut_threshold_associations: TiVec<CutId, CutThresholdCountertermAssociations>,
    pub reversed_edges: TiVec<CutGroupId, Vec<EdgeIndex>>,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub estimated_scale: Option<F<f64>>,
    pub real_mass_vec: Option<EdgeVec<Option<F<f64>>>>,
    pub param_builder: ParamBuilder<f64>,
    pub orientations: TiVec<OrientationID, EdgeVec<Orientation>>,
    production_orientation_keys: Vec<String>,
    pub orientation_filter: SubSet<OrientationID>,
    pub explicit_orientation_sum_only: bool,
    #[allow(private_interfaces)]
    pub counterterm: LUCounterTerm,
    pub cut_group_data: CutGroupData,
}

struct CutEventGenerationContext<'a> {
    model: &'a Model,
    channel_id: Option<SamplingChannelId>,
}

struct DeferredCutEvaluation<T: FloatLike> {
    cut_group_id: CutGroupId,
    kinematic_point: LUCTKinematicPoint<T>,
    bare_cut_total: Complex<F<T>>,
    threshold_counterterm_weights: Vec<Complex<F<T>>>,
    accepted_event: Option<GenericEvent<T>>,
}

impl CrossSectionGraphTerm {
    pub fn threshold_counterterm_metadata(&self) -> Option<&ThresholdCountertermMetadataRegistry> {
        self.counterterm.metadata_registry.as_ref()
    }

    fn build_threshold_multiplier_collection(
        graph: &CrossSectionGraph,
        cut_group_id: CutGroupId,
        counterterm_data: &LUCounterTermData,
        settings: &GlobalSettings,
    ) -> Result<Option<ThresholdMultiplierEvaluatorCollection>> {
        let resolved = graph
            .derived_data
            .resolved_threshold_counterterms
            .as_ref()
            .ok_or_else(|| {
                eyre!(
                    "graph '{}' has LU counterterms but no resolved threshold-counterterm variants",
                    graph.graph.name,
                )
            })?;
        let variant_ids = counterterm_data
            .left_variant_ids
            .iter()
            .chain(counterterm_data.right_variant_ids.iter())
            .copied()
            .collect::<Vec<_>>();
        if variant_ids
            .iter()
            .all(|variant_id| resolved.variants[*variant_id].multiplier.is_none())
        {
            return Ok(None);
        }

        let cut_group = &graph.derived_data.cut_group_data.cut_groups[cut_group_id];
        let mut esurface_ids = cut_group
            .cuts
            .iter()
            .map(|cut_id| graph.cut_esurface_id_map[*cut_id])
            .chain(
                cut_group
                    .related_esurface_group
                    .esurface_ids
                    .iter()
                    .copied(),
            )
            .collect::<BTreeSet<_>>();
        for variant_id in variant_ids {
            let variant = &resolved.variants[variant_id];
            esurface_ids.extend(variant.threshold_esurface_ids.iter().copied());
            esurface_ids.extend(variant.raised_esurface_group.esurface_ids.iter().copied());
            esurface_ids.extend(
                variant
                    .associations
                    .iter()
                    .map(|association| association.esurface_id),
            );
        }
        let layout = ThresholdMultiplierLayout::from_graph_esurfaces(&graph.graph, esurface_ids)
            .with_context(|| {
                format!(
                    "Failed to construct threshold-multiplier inputs for graph '{}' cut group {}",
                    graph.graph.name, cut_group_id.0,
                )
            })?;

        let parse_variant = |variant_id: ThresholdCountertermVariantId| -> Result<(
            ThresholdCountertermVariantId,
            Option<ThresholdMultiplierExpression>,
        )> {
            let variant = &resolved.variants[variant_id];
            let expression = variant
                .multiplier
                .as_ref()
                .map(|multiplier| {
                    if multiplier.symmetrize {
                        unimplemented!(
                            "symmetrized threshold-counterterm multipliers are not implemented"
                        );
                    }
                    layout
                        .parse_expression(&multiplier.expression)
                        .with_context(|| {
                            format!(
                                "Invalid threshold multiplier for graph '{}' cut group {} variant '{}' ({})",
                                graph.graph.name,
                                cut_group_id.0,
                                variant.name,
                                variant_id.0,
                            )
                        })
                })
                .transpose()?;
            Ok((variant_id, expression))
        };
        let left = counterterm_data
            .left_variant_ids
            .iter()
            .copied()
            .map(&parse_variant)
            .collect::<Result<Vec<_>>>()?;
        let right = counterterm_data
            .right_variant_ids
            .iter()
            .copied()
            .map(&parse_variant)
            .collect::<Result<Vec<_>>>()?;

        ThresholdMultiplierEvaluatorCollection::build(
            layout,
            left,
            right,
            &settings.generation.evaluator,
        )
        .with_context(|| {
            format!(
                "Failed to build threshold multipliers for graph '{}' cut group {}",
                graph.graph.name, cut_group_id.0,
            )
        })
    }

    pub fn from_cross_section_graph(
        graph: &CrossSectionGraph,
        settings: &GlobalSettings,
    ) -> Result<(Self, GraphGenerationStats)> {
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let mut stats = GraphGenerationStats::default();
        let production_orientation_ids = graph
            .derived_data
            .global_cff_expression
            .as_ref()
            .unwrap()
            .expression
            .orientations
            .iter_enumerated()
            .filter_map(|(orientation_id, orientation)| {
                (settings.generation.explicit_orientation_sum_only
                    || settings.generation.orientation_pattern.filter(orientation))
                .then_some(orientation_id)
            })
            .collect_vec();
        let selected_generation_orientations = production_orientation_ids
            .iter()
            .map(|orientation_id| {
                &graph
                    .derived_data
                    .global_cff_expression
                    .as_ref()
                    .unwrap()
                    .expression
                    .orientations[*orientation_id]
            })
            .collect_vec();
        // Every generalized residue map is a separate runtime channel. Its
        // physical directions are metadata and therefore must not deduplicate
        // maps that differ only by numerator/M sampling data.
        let orientations: TiVec<OrientationID, EdgeVec<Orientation>> =
            selected_generation_orientations
                .iter()
                .map(|orientation| orientation.orientation().clone())
                .collect();
        let production_orientation_keys = selected_generation_orientations
            .iter()
            .map(|orientation| orientation.residue_map_key())
            .collect_vec();
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
                orientation
                    .iter_denominator_nodes()
                    .filter_map(|tree_node| {
                        if let HybridSurfaceID::Esurface(esurface_id) = tree_node.data {
                            Some(esurface_id)
                        } else {
                            None
                        }
                    })
            })
            .collect::<HashSet<_>>();

        let active_cut_groups: TiVec<CutGroupId, bool> = graph
            .derived_data
            .cut_group_data
            .cut_groups
            .iter()
            .map(|cut_group| {
                cut_group
                    .related_esurface_group
                    .esurface_ids
                    .iter()
                    .any(|esurface_id| selected_generation_esurfaces.contains(esurface_id))
            })
            .collect();

        let masked_cut_parametric_integrand: TiVec<CutGroupId, _> = graph
            .derived_data
            .cut_paramatric_integrand
            .iter_enumerated()
            .map(|(cut_group_id, integrands)| {
                if active_cut_groups[cut_group_id] {
                    integrands.clone()
                } else {
                    integrands.zero_like()
                }
            })
            .collect();

        let mut active_left_thresholds: TiVec<_, TiVec<_, bool>> = TiVec::new();
        let mut active_right_thresholds: TiVec<_, TiVec<_, bool>> = TiVec::new();
        let mut active_iterated_thresholds: TiVec<_, _> = TiVec::new();
        let mut masked_threshold_counterterms: TiVec<CutGroupId, _> = TiVec::new();
        for (cut_group_id, counterterm_data) in
            graph.derived_data.threshold_counterterms.iter_enumerated()
        {
            let left_active: TiVec<_, bool> = counterterm_data
                .left_thresholds
                .iter()
                .map(|raised_group| {
                    active_cut_groups[cut_group_id]
                        && (settings.generation.explicit_orientation_sum_only
                            || raised_group.esurface_ids.iter().any(|esurface_id| {
                                selected_generation_esurfaces.contains(esurface_id)
                            }))
                })
                .collect();
            let right_active: TiVec<_, bool> = counterterm_data
                .right_thresholds
                .iter()
                .map(|raised_group| {
                    active_cut_groups[cut_group_id]
                        && (settings.generation.explicit_orientation_sum_only
                            || raised_group.esurface_ids.iter().any(|esurface_id| {
                                selected_generation_esurfaces.contains(esurface_id)
                            }))
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
        for (cut_group_id, integrand_for_cut_group) in
            masked_cut_parametric_integrand.iter_enumerated()
        {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            let mut cut_group_integrands = BTreeMap::new();
            for (cut_cff_index, integrand_for_subset) in integrand_for_cut_group.integrands.iter() {
                if crate::is_interrupted() {
                    return Err(eyre!("Generation interrupted by user"));
                }
                let dual_shape = shape_from_cut_cff_index(cut_cff_index);

                let (evaluator_stack, evaluator_timings) =
                    EvaluatorStack::from_integrand_with_timings(
                        integrand_for_subset,
                        &graph.graph.param_builder,
                        integrand_for_cut_group.integrands.numerators(),
                        (!settings.generation.explicit_orientation_sum_only)
                            .then_some((&orientations.raw, &production_orientation_ids)),
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
                cut_group_integrands.insert(*cut_cff_index, evaluator_stack);
            }
            integrand.push(cut_group_integrands);
            processes::cut_finished(
                "",
                &graph.graph.name,
                graph.derived_data.cut_group_data.cut_groups[cut_group_id]
                    .cuts
                    .len(),
            );
        }

        let mut ct_evaluators = TiVec::<CutGroupId, LUCounterTermEvaluators>::new();
        let include_threshold_metadata = graph
            .derived_data
            .resolved_threshold_counterterms
            .as_ref()
            .is_some_and(|resolved| {
                !resolved.legacy_equivalent
                    || (!graph.graph.threshold_counterterms.autogenerated
                        && !graph.graph.threshold_counterterms.cuts.is_empty())
            });
        let threshold_helper_outputs =
            match graph.derived_data.resolved_threshold_counterterms.as_ref() {
                Some(resolved) if !resolved.legacy_equivalent => LUThresholdHelperOutputs::Pieces,
                Some(_) if include_threshold_metadata => LUThresholdHelperOutputs::LegacyAndPieces,
                _ => LUThresholdHelperOutputs::Legacy,
            };
        for (cut_group_id, ct_data) in masked_threshold_counterterms.iter_enumerated() {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }

            let include_integrated = !settings
                .generation
                .threshold_subtraction
                .disable_integrated_ct;
            let optimization_settings = settings.generation.evaluator.optimization_settings();

            let build_single_helpers = |integrands: &crate::uv::forest::ParametricIntegrands,
                                        is_on_right: bool,
                                        loop_count: usize| {
                integrands
                    .integrands
                    .iter()
                    .map(|(cut_cff_index, _)| {
                        let lu_order = cut_cff_index.lu_cut_order.ok_or_else(|| {
                            eyre!("LU threshold counterterm helper index is missing lu_cut_order")
                        })?;
                        let threshold_order = if is_on_right {
                            cut_cff_index.right_threshold_order.ok_or_else(|| {
                                eyre!(
                                    "Right LU threshold helper index is missing right_threshold_order"
                                )
                            })?
                        } else {
                            cut_cff_index.left_threshold_order.ok_or_else(|| {
                                eyre!(
                                    "Left LU threshold helper index is missing left_threshold_order"
                                )
                            })?
                        };

                        let dual_shape = if lu_order > 1 {
                            Some(simple_n_deriv_shape(lu_order - 1))
                        } else {
                            None
                        };

                        let evaluator = graph.single_th_helper(
                            threshold_order as u8,
                            loop_count,
                            is_on_right,
                            include_integrated,
                            threshold_helper_outputs,
                            dual_shape,
                            optimization_settings.clone(),
                            &settings.generation.evaluator,
                        )?;
                        Ok((*cut_cff_index, evaluator))
                    })
                    .collect::<Result<BTreeMap<_, _>>>()
            };

            let build_iterated_helpers =
                |integrands: &crate::uv::forest::ParametricIntegrands,
                 left_loop_count: usize,
                 right_loop_count: usize| {
                    integrands
                    .integrands
                    .iter()
                    .map(|(cut_cff_index, _)| {
                        let lu_order = cut_cff_index.lu_cut_order.ok_or_else(|| {
                            eyre!(
                                "Iterated LU threshold counterterm helper index is missing lu_cut_order"
                            )
                        })?;
                        let left_threshold_order =
                            cut_cff_index.left_threshold_order.ok_or_else(|| {
                                eyre!(
                                    "Iterated LU threshold helper index is missing left_threshold_order"
                                )
                            })?;
                        let right_threshold_order =
                            cut_cff_index.right_threshold_order.ok_or_else(|| {
                                eyre!(
                                    "Iterated LU threshold helper index is missing right_threshold_order"
                                )
                            })?;


                        let dual_shape = if lu_order > 1 {
                            Some(simple_n_deriv_shape(lu_order - 1))
                        } else {
                            None
                        };
                        let evaluator = graph.iterated_th_helper(
                            left_threshold_order as u8,
                            right_threshold_order as u8,
                            left_loop_count,
                            right_loop_count,
                            include_integrated,
                            threshold_helper_outputs,
                            dual_shape,
                            optimization_settings.clone(),
                            &settings.generation.evaluator,
                        )
                        .with_context(|| {
                            format!(
                                "Failed to build iterated threshold helper for graph '{}' cut group {} with index {:?} (left loops {}, right loops {})",
                                graph.graph.name,
                                cut_group_id.0,
                                cut_cff_index,
                                left_loop_count,
                                right_loop_count,
                            )
                        })?;
                        Ok((*cut_cff_index, evaluator))
                    })
                    .collect::<Result<BTreeMap<_, _>>>()
                };

            let left_thresholds = ct_data
                .left_atoms
                .iter()
                .zip(&ct_data.left_subspaces)
                .map(|(integrands, subspace)| {
                    build_single_helpers(integrands, false, subspace.loopcount())
                })
                .collect::<Result<TiVec<_, _>>>()?;
            let right_thresholds = ct_data
                .right_atoms
                .iter()
                .zip(&ct_data.right_subspaces)
                .map(|(integrands, subspace)| {
                    build_single_helpers(integrands, true, subspace.loopcount())
                })
                .collect::<Result<TiVec<_, _>>>()?;
            let num_right_thresholds = ct_data.iterated.num_right_thresholds();
            let iterated = IteratedCtCollection::new(
                ct_data
                    .iterated
                    .iter()
                    .enumerate()
                    .map(|(flat_index, integrands)| {
                        let left_id = LeftThresholdId::from(flat_index / num_right_thresholds);
                        let right_id = RightThresholdId::from(flat_index % num_right_thresholds);
                        build_iterated_helpers(
                            integrands,
                            ct_data.left_subspaces[left_id].loopcount(),
                            ct_data.right_subspaces[right_id].loopcount(),
                        )
                    })
                    .collect::<Result<Vec<_>>>()?,
                left_thresholds.len(),
                right_thresholds.len(),
            );
            let threshold_helpers = LUThresholdHelperEvaluators {
                left_thresholds,
                right_thresholds,
                iterated,
            };
            let threshold_multipliers = Self::build_threshold_multiplier_collection(
                graph,
                cut_group_id,
                ct_data,
                settings,
            )?;

            let (evaluators, evaluator_timings) = LUCounterTermEvaluators::from_atoms(
                ct_data,
                graph.derived_data.cut_group_data.cut_groups[cut_group_id]
                    .related_esurface_group
                    .max_occurence,
                threshold_helpers,
                threshold_multipliers,
                &graph.graph.param_builder,
                settings,
                &orientations,
                &production_orientation_ids,
            );
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            stats.add_evaluator_build_timings(evaluator_timings);
            stats.evaluator_count += evaluators.generic_compileable_evaluator_count();
            ct_evaluators.push(evaluators);
        }

        let expression_esurfaces = &graph
            .derived_data
            .global_cff_expression
            .as_ref()
            .expect("global CFF expression should have been created")
            .expression
            .surfaces
            .esurface_cache;
        let mut thresholds = TiVec::new();
        for ct_data in &graph.derived_data.threshold_counterterms {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            thresholds.push((
                ct_data
                    .left_thresholds
                    .iter()
                    .map(|raised_group| expression_esurfaces[raised_group.esurface_ids[0]].clone())
                    .collect(),
                ct_data
                    .right_thresholds
                    .iter()
                    .map(|raised_group| expression_esurfaces[raised_group.esurface_ids[0]].clone())
                    .collect(),
            ));
        }

        let rstar_dependence_calculator = graph
            .derived_data
            .cut_group_data
            .cut_groups
            .iter()
            .map(|cut_group| {
                generate_rstar_t_dependence_evaluator(
                    cut_group
                        .related_esurface_group
                        .max_occurence
                        .saturating_sub(1),
                )
            })
            .collect::<Result<TiVec<CutGroupId, _>>>()?;

        let (variant_subspaces, metadata_registry) = if let Some(resolved) =
            graph.derived_data.resolved_threshold_counterterms.as_ref()
        {
            let variant_subspaces = if resolved.legacy_equivalent {
                None
            } else {
                Some(
                    graph
                        .derived_data
                        .threshold_counterterms
                        .iter()
                        .map(|counterterm_data| LUVariantSubspaces {
                            left_variant_ids: counterterm_data.left_variant_ids.clone(),
                            right_variant_ids: counterterm_data.right_variant_ids.clone(),
                            left: counterterm_data.left_subspaces.clone(),
                            right: counterterm_data.right_subspaces.clone(),
                        })
                        .collect::<TiVec<CutGroupId, _>>(),
                )
            };

            let metadata_registry = if include_threshold_metadata {
                let mut variant_statuses = resolved
                    .variants
                    .iter()
                    .map(|variant| ThresholdCountertermVariantStatus {
                        generated: variant
                            .associations
                            .iter()
                            .any(|association| association.eligible),
                        active: false,
                    })
                    .collect::<Vec<_>>();
                for (cut_group_id, counterterm_data) in
                    graph.derived_data.threshold_counterterms.iter_enumerated()
                {
                    for (&variant_id, &active) in counterterm_data
                        .left_variant_ids
                        .iter()
                        .zip(&active_left_thresholds[cut_group_id])
                        .chain(
                            counterterm_data
                                .right_variant_ids
                                .iter()
                                .zip(&active_right_thresholds[cut_group_id]),
                        )
                    {
                        variant_statuses[variant_id.0].active |= active;
                    }
                }
                let evaluator_registrations = ct_evaluators
                    .iter_enumerated()
                    .flat_map(|(cut_group_id, evaluators)| {
                        evaluators
                            .threshold_multipliers
                            .as_ref()
                            .into_iter()
                            .flat_map(move |collection| {
                                collection.metadata_registrations(Some(cut_group_id.0))
                            })
                    })
                    .collect();
                Some(ThresholdCountertermMetadataRegistry::build(
                    graph.graph.name.clone(),
                    resolved,
                    graph.derived_data.lmbs.as_ref().ok_or_else(|| {
                        eyre!(
                            "graph '{}' has threshold metadata but no loop-momentum bases",
                            graph.graph.name,
                        )
                    })?,
                    &variant_statuses,
                    evaluator_registrations,
                )?)
            } else {
                None
            };
            (variant_subspaces, metadata_registry)
        } else {
            if !graph.derived_data.threshold_counterterms.is_empty()
                || !ct_evaluators.is_empty()
                || !thresholds.is_empty()
            {
                return Err(eyre!(
                    "graph '{}' has LU counterterms but no resolved threshold-counterterm variants",
                    graph.graph.name,
                ));
            }
            // Threshold generation may be disabled entirely. This is the pre-existing empty LU
            // representation and must not allocate generalized subspace or metadata state.
            (None, None)
        };

        let counterterm = LUCounterTerm {
            evaluators: ct_evaluators,
            thresholds,
            subspaces: graph.derived_data.subspace_data.clone(),
            variant_subspaces,
            metadata_registry,
            rstar_dependence_calculator,
            active_cut_groups,
            active_left_thresholds,
            active_right_thresholds,
            active_iterated_thresholds,
        };

        let reversed_edges = graph
            .derived_data
            .cut_group_data
            .cut_groups
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
                covariant_cut_representatives: graph
                    .derived_data
                    .covariant_cut_representatives
                    .clone(),
                topological_threshold_esurfaces: graph
                    .derived_data
                    .topological_threshold_esurfaces
                    .clone(),
                cut_threshold_associations: graph.derived_data.cut_threshold_associations.clone(),
                multi_channeling_setup: LmbMultiChannelingSetup {
                    sampling_bridge: Default::default(),
                    sampling_bridge_quad: Default::default(),
                    sampling_bridge_arb: Default::default(),
                    sampling_catalogue: Default::default(),
                    sampling_programs: Default::default(),
                    lmb_basis_ids: TiVec::new(),
                    graph: graph.graph.clone(), // will be overwritten later,
                    all_bases: TiVec::new(),
                },
                lmbs: graph.derived_data.lmbs.as_ref().unwrap().clone(),
                estimated_scale: None,
                real_mass_vec: None,
                param_builder: graph.graph.param_builder.clone(),
                orientation_filter: SubSet::full(orientations.len()),
                orientations,
                production_orientation_keys,
                explicit_orientation_sum_only: settings.generation.explicit_orientation_sum_only,
                counterterm,
                reversed_edges,
                cut_group_data: graph.derived_data.cut_group_data.clone(),
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

        for (cut_group_id, integrands) in self.integrand.iter_mut().enumerate() {
            for (cut_cff_index, integrand) in integrands.iter_mut() {
                let n_derivatives = cut_cff_index.lu_cut_order.unwrap_or(0);
                integrand.compile(
                    format!(
                        "integrand_zen_cut_group_{}_deriv_{}",
                        cut_group_id, n_derivatives
                    ),
                    graph_path.clone(),
                    frozen_mode,
                )?;
            }
        }

        self.counterterm.compile(&graph_path, frozen_mode)?;

        for (index, evaluator) in self
            .cut_group_data
            .pass_two_evaluators
            .iter_mut()
            .enumerate()
        {
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
        for cut_group_integrands in self.integrand.iter_mut() {
            for evaluator_stack in cut_group_integrands.values_mut() {
                evaluator_stack.for_each_generic_evaluator_mut(&mut f)?;
            }
        }

        self.counterterm.for_each_generic_evaluator_mut(&mut f)?;

        for evaluator in self.cut_group_data.pass_two_evaluators.iter_mut() {
            f(evaluator)?;
        }

        Ok(())
    }

    fn generate_event_for_cut<T: FloatLike>(
        &self,
        event_context: CutEventGenerationContext<'_>,
        t_scaling_solution: &NewtonIterationResult<T>,
        momentum_sample: &MomentumSample<T>,
        cut_id: CutId,
        cut: &CrossSectionCut,
    ) -> Result<GenericEvent<T>> {
        let rescaled_momenta =
            momentum_sample.rescaled_loop_momenta(&t_scaling_solution.solution, Subspace::None);

        let mut new_event = GenericEvent::<T>::default();
        new_event.cut_info.cut_id = cut_id.0;
        new_event.cut_info.orientation_id = if self.explicit_orientation_sum_only {
            Some(0)
        } else {
            momentum_sample.sample.orientation
        };
        new_event.cut_info.sampling_channel_id = event_context.channel_id.map(usize::from);
        new_event.cut_info.sampling_channel_edge_ids = event_context
            .channel_id
            .map(|channel_id| {
                let catalogue = self
                    .multi_channeling_setup
                    .sampling_catalogue
                    .as_ref()
                    .ok_or_else(|| {
                        eyre!(
                            "sampling event metadata for graph '{}' requires warm_up",
                            self.graph.name
                        )
                    })?;
                match catalogue.entries.get(channel_id.index()) {
                    Some(super::SamplingCatalogueEntry::Lmb { edges, .. }) => {
                        Ok(Some(edges.iter().copied().collect()))
                    }
                    Some(_) => Ok(None),
                    None => Err(eyre!(
                        "sampling event channel {} is absent from graph '{}'",
                        channel_id.index(),
                        self.graph.name
                    )),
                }
            })
            .transpose()?
            .flatten();
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
                    event_context
                        .model
                        .get_particle_from_pdg(edge_pdg)
                        .get_anti_particle(event_context.model)
                        .pdg_code
                }
            };

            let mass_value =
                if let Some(mass) = d.data.mass.value(event_context.model, &self.param_builder) {
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
            new_event.cut_info.particle_pdgs.1.push(
                self.covariant_cut_representatives
                    .get(&cut_pdg)
                    .copied()
                    .unwrap_or(cut_pdg),
            );
        }

        Ok(new_event)
    }
}

impl GraphTerm for CrossSectionGraphTerm {
    fn sampling_setup(&self) -> &LmbMultiChannelingSetup {
        &self.multi_channeling_setup
    }

    fn sampling_setup_mut(&mut self) -> &mut LmbMultiChannelingSetup {
        &mut self.multi_channeling_setup
    }

    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64> {
        &mut self.param_builder
    }

    fn selected_lmb_basis_id(
        &self,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<LmbIndex> {
        self.multi_channeling_setup.selected_lmb_basis_id(
            &self.multi_channeling_setup.graph.name,
            parameterization_settings,
        )
    }

    fn bind_sampling_bridge<T: FloatLike>(
        &self,
        catalogue: &super::SamplingChannelCatalogue,
        programs: &[super::sampling_selection::SamplingChannelPrograms],
        parameterization_settings: &ParameterizationSettings,
        settings: &RuntimeSettings,
        external_momenta: &[[T; 4]],
        orientation: Option<usize>,
    ) -> Result<SamplingChannelBridge<T>> {
        let e_cm = settings.kinematics.e_cm;
        let parent_lmb = self
            .multi_channeling_setup
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect();
        let zero = F::<T>::default().zero();
        let mut context = SamplingChannelCompileContext::new(
            self.multi_channeling_setup.graph.name.clone(),
            parent_lmb,
            parameterization_settings.clone(),
            e_cm,
            self.graph.get_loop_number(),
        );
        context.orientation = orientation;
        if catalogue.named_entries().any(|channel| {
            channel
                .blocks
                .iter()
                .any(|block| !block.target.energy_edge_sets().is_empty())
        }) {
            if self.graph.loop_momentum_basis
                != self.multi_channeling_setup.graph.loop_momentum_basis
            {
                return Err(eyre!(
                    "cross-section surface sampling for graph '{}' requires the master's complete loop/external routing; graph parent {:?}, master '{}' parent {:?}",
                    self.graph.name,
                    self.graph.loop_momentum_basis.loop_edges,
                    self.multi_channeling_setup.graph.name,
                    self.multi_channeling_setup
                        .graph
                        .loop_momentum_basis
                        .loop_edges,
                ));
            }
            if external_momenta.len() != self.graph.loop_momentum_basis.ext_edges.len()
                || external_momenta
                    .iter()
                    .flatten()
                    .any(|value| !value.is_finite())
            {
                return Err(eyre!(
                    "cross-section surface sampling for graph '{}' needs {} finite external four-momenta, received {:?}",
                    self.graph.name,
                    self.graph.loop_momentum_basis.ext_edges.len(),
                    external_momenta
                ));
            }
            let cached_masses = self.real_mass_vec.as_ref().ok_or_else(|| {
                eyre!(
                    "physical-cut sampling for graph '{}' requires warmup mass data",
                    self.graph.name
                )
            })?;
            let masses = self.graph.new_edgevec(|_, edge, _| {
                cached_masses[edge]
                    .map(F::<T>::from_ff64)
                    .unwrap_or_else(|| zero.clone())
            });
            let externals =
                ExternalFourMomenta::from_iter(external_momenta.iter().map(|momentum| {
                    FourMomentum::from_args(
                        F(momentum[0].clone()),
                        F(momentum[1].clone()),
                        F(momentum[2].clone()),
                        F(momentum[3].clone()),
                    )
                }));
            for (cut_group_id, cut_group) in self.cut_group_data.cut_groups.iter_enumerated() {
                if !self.counterterm.cut_group_is_active(cut_group_id) {
                    continue;
                }
                for cut_id in &cut_group.cuts {
                    let surface = &self.cut_esurface[*cut_id];
                    let edges = surface
                        .energies
                        .iter()
                        .map(|edge| edge.0)
                        .sorted()
                        .collect_vec();
                    let selected = catalogue.named_entries().any(|channel| {
                        channel
                            .blocks
                            .iter()
                            .any(|block| block.target.host_cut() == Some(edges.as_slice()))
                    });
                    let cut_ids = context.physical_cut_ids.entry(edges.clone()).or_default();
                    if selected && let Some(previous_id) = cut_ids.first() {
                        let previous = &self.cut_esurface[CutId(*previous_id)];
                        if previous.external_shift.iter().sorted().collect_vec()
                            != surface.external_shift.iter().sorted().collect_vec()
                        {
                            return Err(eyre!(
                                "physical-cut sampling for graph '{}' is ambiguous: cut IDs {} and {} have the same energy edges {:?} but incompatible external shifts {:?} and {:?}",
                                self.graph.name,
                                previous_id,
                                cut_id.0,
                                edges,
                                previous.external_shift,
                                surface.external_shift,
                            ));
                        }
                    }
                    cut_ids.push(cut_id.0);
                    // Equivalent geometries can belong to several raised groups.
                    // Keep their largest derivative packet before deduplicating
                    // the root evaluator; a named profile inherits this bound.
                    context
                        .physical_cut_max_occurrences
                        .entry(edges)
                        .and_modify(|order| {
                            *order = (*order).max(cut_group.related_esurface_group.max_occurence)
                        })
                        .or_insert(cut_group.related_esurface_group.max_occurence);
                }
            }
            for channel in catalogue.named_entries() {
                let parent = &channel.definition.parent_lmb;
                let lmbs = TiVec::from(vec![
                    self.multi_channeling_setup.sampling_parent_lmb(parent)?,
                ]);
                let basis_id = LmbIndex::from(0);
                let lmb = &lmbs[basis_id];
                let frame = self
                    .multi_channeling_setup
                    .lmb_frame_map(lmb, external_momenta)?;
                // LU fixes generation K=0. In an affine native frame L=A K+BQ,
                // its fixed point is BQ; scaling L itself would change the cut.
                let origin = frame
                    .inverse(&vec![zero.0.clone(); 3 * parent.len()], &[])?
                    .coordinates;
                let origin_loops = LoopMomenta::from_iter(origin.chunks_exact(3).map(|v| {
                    ThreeMomentum::new(F(v[0].clone()), F(v[1].clone()), F(v[2].clone()))
                }));
                for block in &channel.blocks {
                    let energy_sets = block.target.energy_edge_sets();
                    if energy_sets.is_empty() {
                        continue;
                    }
                    let [edges] = energy_sets.as_slice() else {
                        return Err(eyre!(
                            "cross-section joint sampling channel '{}' is not enabled: routed-energy matching, physical host reconciliation and proposal-policy retention across native retries must be bound before production intersect targets can be used",
                            channel.name
                        ));
                    };
                    if context
                        .geometry_maps
                        .contains_key(&block.geometry_key(parent))
                    {
                        continue;
                    }
                    let host_edges = block.target.host_cut().ok_or_else(|| eyre!(
                        "cross-section surface channel '{}' needs at_cut(cut(...), ...) or a preceding phase_space(cut(...)) host for target {:?}",
                        channel.name, block.target))?;
                    let cut_id = context.physical_cut_ids.get(host_edges)
                        .and_then(|ids| ids.iter().find(|id| channel.definition.on_cut.is_empty() || channel.definition.on_cut.contains(id)))
                        .copied().map(CutId).ok_or_else(|| eyre!(
                            "sampling channel '{}' requests inactive or excluded host cut {:?}; active cuts {:?}, on_cut {:?}",
                            channel.name, host_edges, context.physical_cut_ids, channel.definition.on_cut))?;
                    let host = &self.cut_esurface[cut_id];
                    let active = block
                        .active_lmb
                        .iter()
                        .map(|edge| {
                            LoopIndex(
                                parent
                                    .iter()
                                    .position(|candidate| candidate == edge)
                                    .expect("resolved active parent edge"),
                            )
                        })
                        .collect_vec();
                    let preceding = block
                        .preceding_lmb
                        .iter()
                        .map(|edge| {
                            LoopIndex(
                                parent
                                    .iter()
                                    .position(|candidate| candidate == edge)
                                    .expect("resolved prior parent edge"),
                            )
                        })
                        .collect_vec();
                    let remaining = block
                        .remaining_lmb
                        .iter()
                        .map(|edge| {
                            LoopIndex(
                                parent
                                    .iter()
                                    .position(|candidate| candidate == edge)
                                    .expect("resolved remaining parent edge"),
                            )
                        })
                        .collect_vec();
                    let certify_independent = |surface: &Esurface,
                                               omitted: &[LoopIndex],
                                               role: &str|
                     -> Result<()> {
                        let dependencies = surface
                            .energies
                            .iter()
                            .flat_map(|edge| {
                                let row = lmb.edge_signatures[*edge].internal.to_momtrop_format();
                                omitted.iter().filter_map(move |index| {
                                    (row[index.0] != 0).then_some((
                                        edge.0,
                                        parent[index.0],
                                        row[index.0],
                                    ))
                                })
                            })
                            .collect_vec();
                        if !dependencies.is_empty() {
                            return Err(eyre!(
                                "sampling channel '{}' target {:?}, parent {:?}: {} depends on omitted/active coordinates (energy edge, LMB edge, signed coefficient) {:?}; host {:?}, active {:?}, preceding {:?}, remaining {:?}",
                                channel.name,
                                block.target,
                                parent,
                                role,
                                dependencies,
                                host_edges,
                                block.active_lmb,
                                block.preceding_lmb,
                                block.remaining_lmb
                            ));
                        }
                        Ok(())
                    };
                    let map = if block.target.is_phase_space() {
                        let outside = preceding.iter().chain(&remaining).copied().collect_vec();
                        certify_independent(host, &outside, "phase-space cut")?;
                        // This ray is centered at the actual LU fixed point, not
                        // a generic SOCP center. Only then does t*=R(direction)/r.
                        let zero_velocity =
                            LoopMomenta::from_iter((0..parent.len()).map(|_| {
                                ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())
                            }));
                        let value = host
                            .compute_self_and_r_derivative(
                                &zero,
                                &zero_velocity,
                                &origin_loops,
                                &externals,
                                &masses,
                                lmb,
                            )
                            .0;
                        if !value.0.is_finite()
                            || value >= -F::<T>::from_f64(e_cm) * zero.epsilon() * zero.from_i64(64)
                        {
                            return Err(SamplingEvaluationError::UncertainGeometry { detail: format!(
                                "phase-space channel '{}' host {:?} is not certified strictly interior at its LU fixed point BQ={:?}: eta={}; a shifted geometric center cannot retain the LU R/r profile",
                                channel.name, host_edges, origin, value) }.into());
                        }
                        let surface = host.clone();
                        let lmb = lmb.clone();
                        let masses = masses.clone();
                        let externals = externals.clone();
                        let origin_loops = origin_loops.clone();
                        let active = active.clone();
                        let dimension = 3 * active.len();
                        let evaluator = Arc::new(move |direction: &[T], radius: T| {
                            let radius = F(radius);
                            let mut velocity =
                                LoopMomenta::from_iter((0..lmb.loop_edges.len()).map(|_| {
                                    ThreeMomentum::new(radius.zero(), radius.zero(), radius.zero())
                                }));
                            for (&index, v) in active.iter().zip(direction.chunks_exact(3)) {
                                velocity[index] = ThreeMomentum::new(
                                    F(v[0].clone()),
                                    F(v[1].clone()),
                                    F(v[2].clone()),
                                );
                            }
                            let (value, derivative) = surface.compute_self_and_r_derivative(
                                &radius,
                                &velocity,
                                &origin_loops,
                                &externals,
                                &masses,
                                &lmb,
                            );
                            Ok((value.0, derivative.0))
                        });
                        let map =
                            CompiledSamplingMap::ImplicitSurface(ImplicitSurfaceRadialMap::new(
                                dimension,
                                vec![zero.0.clone(); dimension],
                                e_cm * parameterization_settings.b,
                                parameterization_settings.power,
                                evaluator,
                            )?);
                        let translation = block
                            .active_lmb
                            .iter()
                            .flat_map(|edge| {
                                let index = parent
                                    .iter()
                                    .position(|candidate| candidate == edge)
                                    .unwrap();
                                origin[3 * index..3 * index + 3].iter().cloned()
                            })
                            .collect_vec();
                        let matrix = (0..dimension)
                            .map(|i| {
                                (0..dimension)
                                    .map(|j| if i == j { zero.one().0 } else { zero.0.clone() })
                                    .collect()
                            })
                            .collect();
                        CompiledSamplingMap::Affine {
                            map: Box::new(map),
                            frame: SamplingMapAffine::new(matrix, translation)?,
                        }
                    } else {
                        certify_independent(
                            host,
                            &active.iter().chain(&remaining).copied().collect_vec(),
                            "host cut",
                        )?;
                        let associations = match block.target.cut_side() {
                            Some(SamplingCutSide::Left) => {
                                Some(&self.cut_threshold_associations[cut_id].left)
                            }
                            Some(SamplingCutSide::Right) => {
                                Some(&self.cut_threshold_associations[cut_id].right)
                            }
                            None => None,
                        };
                        // The full-graph equation includes fixed other-side cut
                        // energies. Resolving only runtime CTs would wrongly omit
                        // direct sampling targets whose subtraction is inactive.
                        let candidates = self
                            .topological_threshold_esurfaces
                            .iter_enumerated()
                            .filter(|(id, surface)| {
                                surface
                                    .energies
                                    .iter()
                                    .map(|edge| edge.0)
                                    .sorted()
                                    .eq(edges.iter().copied())
                                    && associations.is_none_or(|associations| {
                                        associations
                                            .iter()
                                            .any(|entry| entry.topological_threshold_id == *id)
                                    })
                            })
                            .collect_vec();
                        let (_, surface) = candidates.first().ok_or_else(|| eyre!(
                            "sampling channel '{}' has no {:?} target {:?} on host cut {} {:?}; topological candidates {:?}",
                            channel.name, block.target.cut_side(), edges, cut_id.0, host_edges,
                            self.topological_threshold_esurfaces.iter_enumerated().map(|(id, surface)| (id, &surface.energies, &surface.external_shift)).collect_vec()))?;
                        if candidates.iter().any(|(_, candidate)| {
                            candidate.external_shift.iter().sorted().collect_vec()
                                != surface.external_shift.iter().sorted().collect_vec()
                        }) {
                            return Err(eyre!(
                                "sampling channel '{}' target {:?} is ambiguous across topological equations {:?}",
                                channel.name,
                                block.target,
                                candidates
                            ));
                        }
                        certify_independent(surface, &remaining, "threshold target")?;
                        let subspace = SubspaceData::new_from_parent_basis_edges(
                            &block
                                .active_lmb
                                .iter()
                                .copied()
                                .map(EdgeIndex)
                                .collect_vec(),
                            &self.graph.full_filter(),
                            basis_id,
                            &self.graph,
                            &lmbs,
                        )?;
                        let complement = (0..parent.len())
                            .map(LoopIndex)
                            .filter(|index| !active.contains(index))
                            .collect_vec();
                        let fiber = surface.sampling_radial_map_in_subspace(
                            &subspace,
                            &lmbs,
                            &self.graph,
                            &masses,
                            &externals,
                            &complement,
                            settings,
                            e_cm * parameterization_settings.b,
                            parameterization_settings.power,
                        )?;
                        let host = host.clone();
                        let lmb = self.graph.loop_momentum_basis.clone();
                        let masses = masses.clone();
                        let externals = externals.clone();
                        let origin = origin.clone();
                        let frame = frame.clone();
                        let active_dimension = 3 * active.len();
                        let active = active.clone();
                        let preceding = preceding.clone();
                        let identity = RadialRootIdentity::new(format!(
                            "sampling host graph '{}' cut {}",
                            self.graph.name, cut_id.0
                        ));
                        let transform = Arc::new(move |raw_prior: &[T]| {
                            if raw_prior.len() != 3 * preceding.len() {
                                return Err(eyre!(
                                    "sampling host requires {} declared prior components, received {}",
                                    3 * preceding.len(),
                                    raw_prior.len()
                                ));
                            }
                            // Only signature-certified null directions remain
                            // unspecified here. Their BQ representative cannot
                            // affect either the host root or the target equation.
                            let mut native = origin.clone();
                            for (&index, values) in preceding.iter().zip(raw_prior.chunks_exact(3))
                            {
                                native[3 * index.0..3 * index.0 + 3].clone_from_slice(values);
                            }
                            let master = frame
                                .forward(&native, &mut SamplingMapContext::detached(&[]))?
                                .point;
                            let loops = LoopMomenta::from_iter(master.chunks_exact(3).map(|v| {
                                ThreeMomentum::new(
                                    F(v[0].clone()),
                                    F(v[1].clone()),
                                    F(v[2].clone()),
                                )
                            }));
                            let zero = F(native[0].clone()).zero();
                            let solution = host
                                .solve_lu_cut(
                                    &loops,
                                    &externals,
                                    &masses,
                                    &lmb,
                                    &F::from_f64(e_cm),
                                    &mut RadialRootDiagnostics::default(),
                                    &identity,
                                )
                                .map_err(|error| SamplingEvaluationError::UncertifiedRoot {
                                    detail: format!("{identity:?}: {error:?}"),
                                })?;
                            let tau = solution.solution;
                            let physical_master = master
                                .into_iter()
                                .map(|value| (F(value) * &tau).0)
                                .collect_vec();
                            if physical_master.iter().any(|value| !value.is_finite()) {
                                return Err(SamplingEvaluationError::Unrepresentable {
                                    operation: "conditional LU kinematics",
                                    detail: format!(
                                        "nonfinite prepared generation coordinates at t={tau}"
                                    ),
                                }
                                .into());
                            }
                            let physical_native = frame.inverse(&physical_master, &[])?.coordinates;
                            let physical_prior = complement
                                .iter()
                                .flat_map(|index| {
                                    physical_native[3 * index.0..3 * index.0 + 3]
                                        .iter()
                                        .cloned()
                                })
                                .collect_vec();
                            let inverse_tau = tau.one() / &tau;
                            let dimension = 3 * active.len();
                            let matrix = (0..dimension)
                                .map(|i| {
                                    (0..dimension)
                                        .map(|j| {
                                            if i == j {
                                                inverse_tau.0.clone()
                                            } else {
                                                zero.0.clone()
                                            }
                                        })
                                        .collect()
                                })
                                .collect();
                            let translation = active
                                .iter()
                                .flat_map(|index| {
                                    origin[3 * index.0..3 * index.0 + 3].iter().map(|value| {
                                        ((tau.one() - &inverse_tau) * F(value.clone())).0
                                    })
                                })
                                .collect_vec();
                            if physical_prior
                                .iter()
                                .chain(&translation)
                                .any(|value| !value.is_finite())
                            {
                                return Err(SamplingEvaluationError::Unrepresentable { operation: "conditional LU frame", detail: format!("nonfinite prepared context or affine translation at t={tau}") }.into());
                            }
                            // This matrix is exactly I/t for a validated positive
                            // root: a constructor failure, including determinant
                            // underflow/overflow, requires native precision rescue.
                            let affine = SamplingMapAffine::new(matrix, translation).map_err(
                                |error| SamplingEvaluationError::Unrepresentable {
                                    operation: "conditional LU Jacobian",
                                    detail: format!(
                                        "could not represent t^(-{dimension}) at t={tau}: {error}"
                                    ),
                                },
                            )?;
                            Ok((physical_prior, affine))
                        });
                        CompiledSamplingMap::Embedded(
                            SamplingMapEmbedding::from_composition(
                                SamplingMapComposition::then(vec![Box::new(fiber)])?,
                                (0..active_dimension).collect(),
                            )?
                            .with_context_transform(transform),
                        )
                    };
                    context.insert_geometry_map(
                        block.target.clone(),
                        parent.clone(),
                        block.active_lmb.clone(),
                        block.preceding_lmb.clone(),
                        map,
                    )?;
                }
            }
        }
        self.multi_channeling_setup
            .compile_sampling_channel_bridge_with_external(
                catalogue,
                programs,
                &context,
                external_momenta,
            )
    }

    fn sampling_channel_is_lmb(
        &self,
        channel_id: SamplingChannelId,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<bool> {
        self.multi_channeling_setup.sampling_channel_is_lmb(
            channel_id,
            &self.multi_channeling_setup.graph.name,
            parameterization_settings,
        )
    }

    fn sampling_channel_ids(
        &self,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Vec<SamplingChannelId>> {
        self.multi_channeling_setup.sampling_channel_ids(
            &self.multi_channeling_setup.graph.name,
            parameterization_settings,
        )
    }

    fn orientation_label(&self, orientation_id: usize) -> Option<String> {
        if self.explicit_orientation_sum_only {
            return (orientation_id == 0).then(|| {
                let n_edges = self
                    .orientations
                    .first()
                    .map(|orientation| orientation.iter().count())
                    .unwrap_or(0);
                "x".repeat(n_edges)
            });
        }

        self.orientations
            .get(resolve_visible_orientation_id(
                &self.orientation_filter,
                orientation_id,
            )?)
            .map(format_orientation_label)
    }

    fn sampling_channel_label(
        &self,
        channel_id: SamplingChannelId,
        parameterization_settings: &ParameterizationSettings,
    ) -> Result<Option<String>> {
        if self.multi_channeling_setup.sampling_channel_is_lmb(
            channel_id,
            &self.multi_channeling_setup.graph.name,
            parameterization_settings,
        )? {
            Ok(Some(format_sampling_channel_label(
                &self.multi_channeling_setup.sampling_channel_edge_ids(
                    channel_id,
                    &self.multi_channeling_setup.graph.name,
                    parameterization_settings,
                )?,
            )))
        } else {
            Ok(Some(self.multi_channeling_setup.sampling_channel_label(
                channel_id,
                &self.multi_channeling_setup.graph.name,
                parameterization_settings,
            )?))
        }
    }

    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()> {
        self.multi_channeling_setup.invalidate_sampling();
        self.graph.validate_real_masses(model)?;
        self.estimated_scale = Some(
            self.graph
                .expected_scale(F(settings.kinematics.e_cm), model),
        );

        if self.explicit_orientation_sum_only {
            self.orientation_filter = SubSet::full(self.orientations.len());
        } else {
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
            .renormalization_localization_scale_value(Complex::new_re(F(settings
                .general
                .renormalization_localization_scale)));
        self.graph
            .param_builder
            .mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq())));
        self.graph
            .param_builder
            .numerator_sampling_scale_value(Complex::new_re(F(settings
                .general
                .numerator_sampling_scale)));
        self.graph.param_builder.update_model_values(model);

        self.param_builder = self.graph.param_builder.clone();
        self.real_mass_vec = Some(self.graph.new_edgevec(|edge, _, _| {
            edge.mass_value(model, &self.param_builder)
                .map(|mass| mass.re)
        }));

        Ok(())
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        mut context: GraphTermEvaluationContext<'_, '_>,
    ) -> Result<GraphEvaluationResult<T>> {
        let orientations =
            momentum_sample.orientations(&self.orientation_filter, &self.orientations);

        // let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        let masses = self.graph.get_real_mass_vector(context.model);
        let hel = context.settings.kinematics.externals.get_helicities();
        let mut cut_results: TiVec<CutGroupId, Vec<Complex<F<T>>>> =
            ti_vec![Vec::new(); self.cut_group_data.cut_groups.len()];
        let mut cut_threshold_counterterms = ti_vec![
            Complex::new_re(momentum_sample.zero());
            self.cut_group_data.cut_groups.len()
        ];
        let mut deferred_cut_evaluations = Vec::new();
        let mut differential_result = GraphEvaluationResult::zero(momentum_sample.zero());
        let mut accepted_event_group = GenericEventGroup::default();

        crate::debug_tags!(#integration, #sample, #inspect;
            "loop moms: {}",
            momentum_sample.loop_moms()
        );

        // Record every active cut root at each precision, even if another cut fails. A
        // later cut then retains its own lower-precision baseline for roundoff rescue.
        let mut lu_solutions = BTreeMap::new();
        let mut lu_root_errors = Vec::new();
        for (cut_group_id, cut_group) in self.cut_group_data.cut_groups.iter_enumerated() {
            if !self.counterterm.cut_group_is_active(cut_group_id) {
                continue;
            }
            let representative_esurface = &self.cut_esurface[cut_group.cuts[0]];

            crate::debug_tags!(#integration, #cut, #inspect;
                "representative esurface: {:#?}",
                representative_esurface
            );

            // A physical LU cut needs an isolated positive root with a finite, positive
            // Jacobian. Validate the bracket and root before constructing any cut kinematics.
            let identity = RadialRootIdentity::new(format!(
                "LU cut graph '{}' cut group {} probe rotation {}",
                self.graph.name, cut_group_id.0, context.rotation.method,
            ));
            let solution = match representative_esurface.solve_lu_cut(
                momentum_sample.loop_moms(),
                momentum_sample.external_moms(),
                &masses,
                &self.graph.loop_momentum_basis,
                &F::from_f64(context.settings.kinematics.e_cm),
                &mut context.evaluation_metadata.radial_root_diagnostics,
                &identity,
            ) {
                Ok(solution) => solution,
                Err(error) => {
                    crate::debug_tags!(#integration, #cut, #solver;
                        graph = %self.graph.name,
                        cut_group_id = cut_group_id.0,
                        edges = ?representative_esurface.energies,
                        error = ?error,
                        "LU radial root requires precision escalation"
                    );
                    // Use the existing residue-failure path so a recoverable numerical
                    // failure is retried at higher precision and remains fatal at the final
                    // level, without ever exposing invalid cut kinematics to the evaluator.
                    lu_root_errors.push(format!(
                        "Could not solve LU cut group {} of graph '{}', edges {:?}: {:?}",
                        cut_group_id.0, self.graph.name, representative_esurface.energies, error,
                    ));
                    continue;
                }
            };

            crate::debug_tags!(#integration, #cut, #solver;
                "solution: {:?}",
                solution
            );

            lu_solutions.insert(cut_group_id, solution);
        }
        if !lu_root_errors.is_empty() {
            context
                .evaluation_metadata
                .record_threshold_counterterm_error(lu_root_errors.join("\n"));
            differential_result.integrand_result = Complex::new_re(F::from_f64(f64::NAN));
            return Ok(differential_result);
        }

        for (cut_group_id, cut_group) in self.cut_group_data.cut_groups.iter_enumerated() {
            let max_occurrence = cut_group.related_esurface_group.max_occurence;
            if !self.counterterm.cut_group_is_active(cut_group_id) {
                let zero = Complex::new_re(momentum_sample.zero());
                for _ in 1..=max_occurrence {
                    cut_results[cut_group_id].push(zero.clone());
                }
                cut_threshold_counterterms[cut_group_id] = zero;
                continue;
            }
            crate::debug_tags!(#integration, #cut;
                "\n =====START EVALUATION FOR CUT GROUP {}=====",
                cut_group_id.0
            );
            let solution = lu_solutions
                .remove(&cut_group_id)
                .expect("all active LU cut roots were validated before evaluation");
            let representative_esurface = &self.cut_esurface[cut_group.cuts[0]];

            let prepared_event = prepare_buffered_event(
                context.settings,
                context.rotation,
                context.event_processing_runtime.as_deref_mut(),
                || {
                    let mut generated = self.generate_event_for_cut::<T>(
                        CutEventGenerationContext {
                            model: context.model,
                            // Mapped sampling channels have already been mapped into the
                            // parent frame; preserve their canonical id in event
                            // metadata without routing the momenta through the
                            // default-sampling LMB reinterpretation path.
                            channel_id: context.sampling_channel,
                        },
                        &solution,
                        momentum_sample,
                        self.cut_group_data.cut_groups[cut_group_id].cuts[0],
                        &self.cuts[self.cut_group_data.cut_groups[cut_group_id].cuts[0]],
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
                for _ in 1..=max_occurrence {
                    cut_results[cut_group_id].push(zero.clone());
                }
                cut_threshold_counterterms[cut_group_id] = zero;
                continue;
            }

            let accepted_event = prepared_event.buffered_event;
            let mut bare_cut_total = Complex::new_re(momentum_sample.zero());
            let threshold_counterterm_weights = Vec::with_capacity(max_occurrence);
            let mut kinematic_point = LUCTKinematicPoint::new(momentum_sample.clone());
            // The sampling partition multiplies the fully subtracted graph
            // outside this LU calculation, after raised-residue derivatives.
            // It is distinct from overlap-group multichanneling internal to
            // threshold subtraction and is common to every cut at this point.
            for num_esurfaces in 1..=max_occurrence {
                let dual_shape = if num_esurfaces > 1 {
                    Some(HyperDual::<F<T>>::new(
                        self.cut_group_data.dual_shapes[num_esurfaces - 2].clone(),
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

                let result = self.integrand[cut_group_id]
                    .get_mut(&cut_index)
                    .unwrap()
                    .evaluate(
                        params,
                        orientations,
                        context.settings,
                        context.evaluation_metadata,
                        context.record_primary_timing,
                    )?
                    .pop()
                    .ok_or_else(|| eyre!("Evaluator returned no cut result"))?;

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
                    &mut self.cut_group_data.pass_two_evaluators[num_esurfaces - 1];

                let pass_two_result = evaluate_evaluator_single(
                    pass_two_evaluator,
                    &params_for_pass_two,
                    context.evaluation_metadata,
                    context.record_primary_timing,
                );

                debug!("pass_two_result: {:+16e}", pass_two_result);
                //debug!("param builder for cut {}: \n{}", cut, self.param_builder);

                let bare_contribution = pass_two_result;
                bare_cut_total += bare_contribution.clone();
                cut_results[cut_group_id].push(bare_contribution);
            }

            deferred_cut_evaluations.push(DeferredCutEvaluation {
                cut_group_id,
                kinematic_point,
                bare_cut_total,
                threshold_counterterm_weights,
                accepted_event,
            });
        }

        let deferred_points = deferred_cut_evaluations
            .iter()
            .map(|deferred| {
                (
                    deferred.cut_group_id,
                    deferred.kinematic_point.representative_sample(),
                )
            })
            .collect_vec();
        let shared_overlaps = if context.settings.subtraction.disable_threshold_subtraction {
            ti_vec![None; self.cut_group_data.cut_groups.len()]
        } else {
            self.counterterm.prepare_shared_overlaps(
                &deferred_points,
                &self.graph,
                &self.graph.get_real_mass_vector(context.model),
                &self.reversed_edges,
                &self.lmbs,
                context.settings,
                context.rotation,
            )?
        };

        for deferred in deferred_cut_evaluations {
            let cut_group_id = deferred.cut_group_id;
            let record_threshold_decomposition = deferred.accepted_event.is_some()
                && context.settings.general.store_additional_weights_in_event
                && self.counterterm.metadata_registry.is_some();
            let counterterm_evaluation =
                if context.settings.subtraction.disable_threshold_subtraction {
                    LUCountertermEvaluation {
                        total: Complex::new_re(momentum_sample.zero()),
                        components: record_threshold_decomposition.then(Vec::new),
                    }
                } else {
                    self.counterterm.evaluate(
                        &deferred.kinematic_point,
                        cut_group_id,
                        &self.reversed_edges[cut_group_id],
                        &self.lmbs,
                        &self.graph,
                        &self.graph.get_real_mass_vector(context.model),
                        context.rotation,
                        context.settings,
                        &mut self.param_builder,
                        orientations,
                        context.evaluation_metadata,
                        context.record_primary_timing,
                        record_threshold_decomposition,
                        shared_overlaps[cut_group_id].as_ref(),
                    )?
                };
            let threshold_decomposition = counterterm_evaluation.components.map(|components| {
                let mut decomposition = GenericThresholdCountertermEventInfo {
                    original: Complex::new_re(momentum_sample.zero()),
                    components,
                };
                decomposition.original = deferred.bare_cut_total.clone();
                decomposition
            });
            let ct_result = if let Some(decomposition) = &threshold_decomposition {
                decomposition.components.iter().fold(
                    Complex::new_re(momentum_sample.zero()),
                    |total, component| total + &component.weighted,
                )
            } else {
                counterterm_evaluation.total
            };

            let mut threshold_counterterm_weights = deferred.threshold_counterterm_weights;
            threshold_counterterm_weights.push(ct_result.clone());
            cut_threshold_counterterms[cut_group_id] = ct_result;

            if let Some(mut event) = deferred.accepted_event {
                let threshold_counterterm_total = threshold_counterterm_weights
                    .iter()
                    .fold(Complex::new_re(momentum_sample.zero()), |acc, value| {
                        acc + value.clone()
                    });
                event.weight = deferred.bare_cut_total.clone() + threshold_counterterm_total;

                if context.settings.general.store_additional_weights_in_event {
                    event
                        .additional_weights
                        .weights
                        .insert(AdditionalWeightKey::Original, deferred.bare_cut_total);
                    for (subset_index, threshold_counterterm) in
                        threshold_counterterm_weights.into_iter().enumerate()
                    {
                        event.additional_weights.weights.insert(
                            AdditionalWeightKey::ThresholdCounterterm { subset_index },
                            threshold_counterterm,
                        );
                    }
                    if let Some(decomposition) = threshold_decomposition {
                        event.weight = decomposition.total();
                        event.additional_weights.threshold_counterterms = Some(decomposition);
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
                event.apply_multiplicative_factor(&flux_factor);
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
        if self.explicit_orientation_sum_only {
            return 1;
        }

        filtered_orientation_count(&self.orientation_filter, &self.orientations)
    }

    fn production_orientation_keys(&self) -> &[String] {
        &self.production_orientation_keys
    }

    fn selected_production_orientation_keys(&self) -> Vec<&str> {
        if self.orientation_filter.is_full() {
            self.production_orientation_keys
                .iter()
                .map(String::as_str)
                .collect()
        } else {
            self.orientation_filter
                .included_iter()
                .map(|id| self.production_orientation_keys[id.0].as_str())
                .collect()
        }
    }

    fn get_tropical_sampler(&self) -> &momtrop::SampleGenerator<3> {
        unimplemented!(
            "Don't know how to generate subgraph table for forward scattering graphs yet"
        )
    }

    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>> {
        self.real_mass_vec
            .as_ref()
            .expect("real mass vector should be set during warmup")
            .clone()
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
