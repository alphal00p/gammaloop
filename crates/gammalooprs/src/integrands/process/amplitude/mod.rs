use std::{
    collections::HashSet,
    fs::{self},
    path::Path,
    time::Instant,
};

use bincode_trait_derive::{Decode, Encode};

use color_eyre::{Help, Result};

use colored::Colorize;
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{ModifySubSet, SubSetLike, subset::SubSet},
};
use momtrop::SampleGenerator;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use spenso::algebra::complex::Complex;
use symbolica::{
    atom::AtomCore,
    evaluate::OptimizationSettings,
    numerical_integration::{Grid, Sample},
};
use three_dimensional_reps::utils::rank_i64;
use tracing::{debug, info, instrument, warn};
use typed_index_collections::{TiVec, ti_vec};

use crate::{
    DependentMomentaConstructor, F, FloatLike, GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::{
            Esurface, EsurfaceCollection, EsurfaceRay, ExistingEsurfaces, GroupEsurfaceId,
            RaisedEsurfaceId, get_representative,
        },
        expression::OrientationID,
        surface::HybridSurfaceID,
    },
    graph::{
        FeynmanGraph, Graph, GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex,
        LoopMomentumBasis, parse::complete_group_parsing,
    },
    integrands::{
        HasIntegrand,
        evaluation::{EvaluationResult, GraphEvaluationResult},
        process::{
            CompiledSamplingMap, ParamBuilder, SamplingChannelBridge,
            SamplingChannelCompileContext, SamplingChannelId, SamplingMapDefinition,
            SurfaceRadialMap,
            evaluators::{ActiveF64Backend, EvaluatorStack},
            graph_to_group_id_for_group_structure,
            threshold_multiplier::ThresholdMultiplierEvaluatorCollection,
        },
    },
    model::Model,
    momentum::{
        FourMomentum, Helicity, Rotation, RotationMethod, SignOrZero, ThreeMomentum,
        sample::{ExternalFourMomenta, ExternalIndex, LoopMomenta, MomentumSample},
        signature::SignatureLike,
    },
    observables::{
        AdditionalWeightKey, EventProcessingRuntime, GenericEvent,
        GenericThresholdCountertermComponentWeight, GenericThresholdCountertermEventInfo,
        ThresholdCountertermComponentOccurrence,
    },
    processes::{
        AmplitudeGraph, GraphGenerationStats, GraphGroupSelectionPlan, GroupDerivedData,
        ThresholdCountertermMetadataRegistry, ThresholdCountertermVariantStatus,
    },
    settings::{
        GlobalSettings, RuntimeSettings,
        global::{CompilationOptimizationLevel, FrozenCompilationMode},
        runtime::{DiscreteGraphSamplingType, ParameterizationSettings, SamplingSettings},
    },
    subtraction::{
        amplitude_counterterm::{
            AmplitudeCountertermAtom, AmplitudeCountertermComponentEvaluation,
            AmplitudeCountertermData, AmplitudeCountertermEvaluation,
            AmplitudeLocalCountertermEvaluation, OverlapStructureWithKinematics,
        },
        overlap::{OverlapInput, SingleGraphOverlapData, find_maximal_overlap},
    },
    utils::{
        ArbPrec, DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, RuntimeCache, W_, compute_shift_part,
        serde_utils::SmartSerde, symbolica_ext::LOGPRINTOPTS,
    },
};

use super::{
    GraphTerm, GraphTermEvaluationContext, LmbMultiChannelingSetup, ProcessIntegrandImpl,
    create_grid, evaluate_sample, filtered_orientation_count, format_orientation_label,
    format_sampling_channel_label, histogram_process_info_for_integrand, prepare_buffered_event,
    resolve_visible_orientation_id, sampling_context::SamplingMapContext,
    sampling_selection::SamplingCatalogueEntry, validate_group_orientation_catalogs,
    validate_process_runtime_settings,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeGraphTerm {
    pub original_integrand: EvaluatorStack,
    pub orientations: TiVec<OrientationID, EdgeVec<Orientation>>,
    production_orientation_keys: Vec<String>,
    pub orientation_filter: SubSet<OrientationID>,
    pub explicit_orientation_sum_only: bool,
    pub esurfaces: EsurfaceCollection,
    pub threshold_counterterm: AmplitudeCountertermData,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub graph: Graph,
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder,
    pub real_mass_vec: Option<EdgeVec<Option<F<f64>>>>,
    pub master_external_signature: SignatureLike<ExternalIndex>,
    pub master_external_pdgs: Vec<isize>,
}

struct AmplitudeGraphTermEvaluation<T: FloatLike> {
    integrand_result: Complex<F<T>>,
    detailed_original: Option<Complex<F<T>>>,
    counterterms: AmplitudeCountertermEvaluation<T>,
}

pub(crate) fn amplitude_threshold_event_info<T: FloatLike>(
    registry: &ThresholdCountertermMetadataRegistry,
    original: Complex<F<T>>,
    components: Vec<AmplitudeCountertermComponentEvaluation<T>>,
) -> Result<GenericThresholdCountertermEventInfo<T>> {
    let components = components
        .into_iter()
        .map(|component| {
            let component_id =
                registry.component_id(None, component.kind, &[component.variant_id])?;
            Ok(GenericThresholdCountertermComponentWeight {
                component_id,
                occurrence: ThresholdCountertermComponentOccurrence::Amplitude {
                    raised_esurface_id: component.esurface_id.0,
                    overlap_group: component.overlap_group,
                },
                multiplier_values: std::iter::once(component.multiplier_value.clone()).collect(),
                effective_multiplier: component.multiplier_value,
                bare: component.bare,
                weighted: component.weighted,
                evaluation_skipped: component.evaluation_skipped,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(GenericThresholdCountertermEventInfo {
        original,
        components,
    })
}

/// Num(sigma_1,sigma_2,...)*(CFF_1 delta(edge(1),1) delta_(1,1,1,-1,1)+CFF_3 delta_(1,1,1,-1,1)+CFF_2 delta_(1,1,1,-1,1))
impl AmplitudeGraphTerm {
    pub fn threshold_counterterm_metadata(&self) -> Option<&ThresholdCountertermMetadataRegistry> {
        self.threshold_counterterm.metadata_registry.as_ref()
    }

    pub fn kinematics_for_threshold_approach(
        &mut self,
        settings: &RuntimeSettings,
        model: &Model,
        momentum_sample: &MomentumSample<ArbPrec>,
    ) -> Result<OverlapStructureWithKinematics<ArbPrec>> {
        self.threshold_counterterm.kinematics_for_approach(
            momentum_sample,
            &self.graph,
            model,
            &self.esurfaces,
            &Rotation::new(RotationMethod::Identity),
            settings,
        )
    }
    pub fn from_amplitude_graph(
        graph: &AmplitudeGraph,
        own_group_position: GraphGroupPosition,
        esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<(Self, GraphGenerationStats)> {
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_from_graph_start",
            graph = %graph.graph.name,
            "Generation timing milestone"
        );
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let mut stats = GraphGenerationStats::default();
        let production_orientation_ids = graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .expression
            .orientations
            .iter_enumerated()
            .filter_map(|(orientation_id, orientation)| {
                settings
                    .generation
                    .explicit_orientation_sum_only
                    .then_some(orientation_id)
                    .or_else(|| {
                        settings
                            .generation
                            .orientation_pattern
                            .filter(orientation)
                            .then_some(orientation_id)
                    })
            })
            .collect_vec();
        let selected_generation_orientations = production_orientation_ids
            .iter()
            .map(|orientation_id| {
                &graph
                    .derived_data
                    .cff_expression
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
                .map(|orientation| orientation.data.orientation.clone())
                .collect();
        let production_orientation_keys = selected_generation_orientations
            .iter()
            .map(|orientation| orientation.residue_map_key())
            .collect_vec();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #orientation, #summary;
            stage = "amplitude_graph_term_orientations_done",
            graph = %graph.graph.name,
            orientation_count = orientations.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
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
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_esurfaces_done",
            graph = %graph.graph.name,
            selected_esurface_count = selected_generation_esurfaces.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        crate::debug_tags!(#generation, #graph, #orientation, #compile, #dump;
            orientation_parametric_integrand = %graph.derived_data.all_mighty_integrand.printer(LOGPRINTOPTS.clone()),
            "Building evaluator for all orientations \n{}",
            graph.graph.param_builder.table()
        );

        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let original_started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_original_evaluator_start",
            graph = %graph.graph.name,
            orientation_count = orientations.len(),
            "Generation timing milestone"
        );
        let (original_integrand, evaluator_timings) =
            if settings.generation.explicit_orientation_sum_only {
                EvaluatorStack::new_explicit_sum_with_timings(
                    &[&graph.derived_data.all_mighty_integrand],
                    &graph.graph.param_builder,
                    None,
                    &settings.generation.evaluator,
                )?
            } else {
                EvaluatorStack::new_with_timings(
                    &[&graph.derived_data.all_mighty_integrand],
                    &graph.graph.param_builder,
                    orientations.as_slice().as_ref(),
                    &production_orientation_ids,
                    None,
                    &settings.generation.evaluator,
                )?
            };
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_original_evaluator_done",
            graph = %graph.graph.name,
            evaluator_count = original_integrand.generic_evaluator_count(),
            elapsed_ms = original_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            spenso_ms = evaluator_timings.spenso_time.as_secs_f64() * 1000.0,
            symbolica_ms = evaluator_timings.symbolica_time.as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        stats.add_evaluator_build_timings(evaluator_timings);
        stats.evaluator_count += original_integrand.generic_evaluator_count();

        let mut threshold_counterterm = AmplitudeCountertermData::new_empty(own_group_position);
        let resolved = graph.derived_data.resolved_threshold_counterterms.as_ref();
        let include_threshold_metadata = resolved.is_some_and(|resolved| {
            !resolved.legacy_equivalent
                || (!graph.graph.threshold_counterterms.autogenerated
                    && !graph.graph.threshold_counterterms.cuts.is_empty())
        });
        threshold_counterterm.legacy_equivalent =
            resolved.is_none_or(|resolved| resolved.legacy_equivalent);
        let selected_generation_raised_esurfaces: HashSet<RaisedEsurfaceId> =
            selected_generation_esurfaces
                .iter()
                .filter_map(|esurface_id| {
                    graph
                        .derived_data
                        .raised_esurface_ids
                        .get(*esurface_id)
                        .copied()
                })
                .collect();

        if threshold_counterterm.legacy_equivalent {
            let mut threshold_evaluators =
                Vec::with_capacity(graph.derived_data.threshold_counterterms.len());
            let active_mask: TiVec<RaisedEsurfaceId, bool> = graph
                .derived_data
                .threshold_counterterms
                .iter_enumerated()
                .map(|(raised_esurface_id, _)| {
                    selected_generation_raised_esurfaces.contains(&raised_esurface_id)
                })
                .collect();
            crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
                stage = "amplitude_graph_term_threshold_setup_done",
                graph = %graph.graph.name,
                threshold_count = graph.derived_data.threshold_counterterms.len(),
                active_threshold_count = active_mask.iter().filter(|active| **active).count(),
                "Generation timing milestone"
            );
            for (raised_esurface_id, ct) in
                graph.derived_data.threshold_counterterms.iter_enumerated()
            {
                if crate::is_interrupted() {
                    return Err(eyre!("Generation interrupted by user"));
                }
                let masked_counterterm = if active_mask[raised_esurface_id] {
                    ct.clone()
                } else {
                    ct.zero_like()
                };
                let (evaluator, evaluator_timings) = masked_counterterm.to_evaluator_with_timings(
                    &graph.graph.param_builder,
                    &orientations,
                    &production_orientation_ids,
                    settings,
                );
                stats.add_evaluator_build_timings(evaluator_timings);
                stats.evaluator_count += evaluator.generic_evaluator_count();
                threshold_evaluators.push(evaluator);
            }
            threshold_counterterm.evaluators = threshold_evaluators.into();
            threshold_counterterm.generated_mask = graph
                .derived_data
                .threshold_counterterms
                .iter()
                .map(AmplitudeCountertermAtom::is_generated)
                .collect();
            threshold_counterterm.active_mask = active_mask;
            threshold_counterterm.helper_evaluators = if include_threshold_metadata {
                let max_order = graph
                    .derived_data
                    .raised_data
                    .pass_two_evaluator
                    .as_ref()
                    .map_or(0, Vec::len);
                (1..=max_order)
                    .map(|order| {
                        crate::processes::threshold_counterterm_recording_helper(
                            order as u8,
                            graph.graph.get_loop_number(),
                            &settings.generation.evaluator,
                        )
                    })
                    .collect()
            } else {
                graph
                    .derived_data
                    .raised_data
                    .pass_two_evaluator
                    .clone()
                    .unwrap_or_default()
            };
            stats.evaluator_count += threshold_counterterm.helper_evaluators.len();
        } else {
            let resolved = resolved.expect("generalized amplitude thresholds must be resolved");
            threshold_counterterm.variant_metadata = resolved.variants.clone();
            if resolved.variants.len() != graph.derived_data.threshold_counterterm_variants.len() {
                return Err(eyre!(
                    "Graph '{}' has {} resolved threshold variants but {} symbolic variant counterterms",
                    graph.graph.name,
                    resolved.variants.len(),
                    graph.derived_data.threshold_counterterm_variants.len(),
                ));
            }

            if resolved
                .variants
                .iter()
                .any(|variant| variant.multiplier.is_some())
            {
                let esurface_ids = resolved.variants.iter().flat_map(|variant| {
                    variant
                        .threshold_esurface_ids
                        .iter()
                        .chain(&variant.raised_esurface_group.esurface_ids)
                        .copied()
                        .chain(
                            variant
                                .associations
                                .iter()
                                .map(|association| association.esurface_id),
                        )
                });
                let layout = super::threshold_multiplier::ThresholdMultiplierLayout::from_graph_esurfaces(
                    &graph.graph,
                    esurface_ids,
                )
                .with_context(|| {
                    format!(
                        "Failed to construct amplitude threshold-multiplier inputs for graph '{}'",
                        graph.graph.name,
                    )
                })?;
                let variants = resolved
                    .variants
                    .iter_enumerated()
                    .map(|(variant_id, variant)| {
                        let expression = variant
                            .multiplier
                            .as_ref()
                            .map(|multiplier| {
                                if multiplier.symmetrize {
                                    unimplemented!(
                                        "symmetrized threshold-counterterm multipliers are not implemented",
                                    );
                                }
                                layout
                                    .parse_expression(&multiplier.expression)
                                    .with_context(|| {
                                        format!(
                                            "Invalid threshold multiplier for amplitude graph '{}' variant '{}' ({})",
                                            graph.graph.name,
                                            variant.name,
                                            variant_id.0,
                                        )
                                    })
                            })
                            .transpose()?;
                        Ok((variant_id, expression))
                    })
                    .collect::<Result<Vec<_>>>()?;
                threshold_counterterm.threshold_multipliers =
                    ThresholdMultiplierEvaluatorCollection::build(
                        layout,
                        variants,
                        Vec::new(),
                        &settings.generation.evaluator,
                    )?;
            }

            for (variant_id, symbolic) in graph
                .derived_data
                .threshold_counterterm_variants
                .iter_enumerated()
            {
                let variant = &resolved.variants[variant_id];
                if symbolic.raised_esurface_id
                    != graph.derived_data.raised_esurface_ids
                        [variant.raised_esurface_group.esurface_ids[0]]
                {
                    return Err(eyre!(
                        "Graph '{}' threshold variant {} symbolic/resolved raised-surface mapping disagrees",
                        graph.graph.name,
                        variant_id.0,
                    ));
                }
                let active =
                    selected_generation_raised_esurfaces.contains(&symbolic.raised_esurface_id);
                let masked = if active {
                    symbolic.atom.clone()
                } else {
                    symbolic.atom.zero_like()
                };
                let (evaluator, evaluator_timings) = masked.to_evaluator_with_timings(
                    &graph.graph.param_builder,
                    &orientations,
                    &production_orientation_ids,
                    settings,
                );
                stats.add_evaluator_build_timings(evaluator_timings);
                stats.evaluator_count += evaluator.generic_evaluator_count();
                threshold_counterterm.variant_evaluators.push(evaluator);
                threshold_counterterm
                    .variant_generated_mask
                    .push(symbolic.atom.is_generated());
                threshold_counterterm.variant_active_mask.push(active);
                threshold_counterterm
                    .variant_raised_esurfaces
                    .push(symbolic.raised_esurface_id);
                threshold_counterterm
                    .variant_subspaces
                    .push(variant.subspace.clone());

                let max_order = variant.raised_esurface_group.max_occurence;
                let helpers = (1..=max_order)
                    .map(|order| {
                        crate::processes::threshold_counterterm_pieces_helper(
                            order as u8,
                            variant.subspace_loop_count,
                            &settings.generation.evaluator,
                        )
                    })
                    .collect_vec();
                stats.evaluator_count += helpers.len();
                threshold_counterterm
                    .variant_helper_evaluators
                    .push(helpers);
            }
        }
        threshold_counterterm.lmbs = graph
            .derived_data
            .threshold_topology
            .as_ref()
            .map(|(_, lmbs)| lmbs.clone())
            .or_else(|| graph.derived_data.lmbs.clone())
            .expect("amplitude graph term requires generated LMBs");
        threshold_counterterm.raised_data = graph.derived_data.raised_data.clone();
        if include_threshold_metadata {
            let resolved = resolved.expect("threshold metadata requires resolved variants");
            let variant_statuses = resolved
                .variants
                .iter_enumerated()
                .map(|(variant_id, variant)| {
                    if threshold_counterterm.legacy_equivalent {
                        let esurface_id = variant.raised_esurface_group.esurface_ids[0];
                        let raised_esurface_id =
                            graph.derived_data.raised_esurface_ids[esurface_id];
                        ThresholdCountertermVariantStatus {
                            generated: threshold_counterterm.generated_mask[raised_esurface_id],
                            active: threshold_counterterm.active_mask[raised_esurface_id],
                        }
                    } else {
                        ThresholdCountertermVariantStatus {
                            generated: threshold_counterterm.variant_generated_mask[variant_id],
                            active: threshold_counterterm.variant_active_mask[variant_id],
                        }
                    }
                })
                .collect::<Vec<_>>();
            let evaluator_registrations = threshold_counterterm
                .threshold_multipliers
                .as_ref()
                .map(|multipliers| multipliers.metadata_registrations(None))
                .unwrap_or_default();
            threshold_counterterm.metadata_registry =
                Some(ThresholdCountertermMetadataRegistry::build(
                    graph.graph.name.clone(),
                    resolved,
                    &threshold_counterterm.lmbs,
                    &variant_statuses,
                    evaluator_registrations,
                )?);
        }

        threshold_counterterm.local_esurface_exists = ti_vec![true; esurface_map.len()];
        threshold_counterterm.esurface_map = esurface_map;
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_from_graph_done",
            graph = %graph.graph.name,
            evaluator_count = stats.evaluator_count,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            spenso_ms = stats.evaluator_spenso_time.as_secs_f64() * 1000.0,
            symbolica_ms = stats.evaluator_symbolica_time.as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        Ok((
            AmplitudeGraphTerm {
                orientation_filter: SubSet::full(orientations.len()),
                orientations,
                production_orientation_keys,
                explicit_orientation_sum_only: settings.generation.explicit_orientation_sum_only,
                original_integrand,
                tropical_sampler: graph.derived_data.tropical_sampler.clone(),
                graph: graph.graph.clone(),
                multi_channeling_setup: LmbMultiChannelingSetup {
                    master_edge_masses: Default::default(),
                    sampling_bridge: Default::default(),
                    sampling_bridge_quad: Default::default(),
                    sampling_bridge_fixed256: Default::default(),
                    sampling_bridge_arb: Default::default(),
                    sampling_source: Default::default(),
                    sampling_catalogue: Default::default(),
                    sampling_programs: Default::default(),
                    lmb_basis_ids: TiVec::new(),
                    graph: graph.graph.clone(), // will be overwritten later,
                    all_bases: TiVec::new(),
                }, // to be taken from froup master
                lmbs: graph
                    .derived_data
                    .lmbs
                    .clone()
                    .expect("lmbs should have been created"),
                threshold_counterterm,
                estimated_scale: None,
                esurfaces: graph
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .expect("cff_expression should have been created")
                    .expression
                    .surfaces
                    .esurface_cache
                    .clone(),
                param_builder: graph.graph.param_builder.clone(),
                real_mass_vec: None,
                master_external_signature: graph.graph.get_external_signature(),
                master_external_pdgs: graph
                    .graph
                    .get_external_partcles()
                    .into_iter()
                    .map(|particle| particle.pdg_code)
                    .collect(),
            },
            stats,
        ))
    }

    #[instrument(
          name = "compile",
          level = "info",
          skip(self, path, override_existing, frozen_mode),
          fields(
              graph.name = %self.graph.name,
              path = %path.as_ref().display(),
          )
      )]
    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        frozen_mode: &FrozenCompilationMode,
    ) -> Result<std::time::Duration> {
        let compile_started = std::time::Instant::now();
        let graph_path = path.as_ref().join(&self.graph.name);

        fs::create_dir_all(&graph_path).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                graph_path.display()
            )
        })?;

        self.original_integrand.compile(
            "orientation_parametric_integrand",
            &graph_path,
            frozen_mode,
        )?;

        self.threshold_counterterm
            .compile(&graph_path, override_existing, frozen_mode)?;

        Ok(compile_started.elapsed())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        self.original_integrand
            .for_each_generic_evaluator_mut(&mut f)?;
        self.threshold_counterterm
            .for_each_generic_evaluator_mut(&mut f)?;
        Ok(())
    }

    pub(crate) fn generic_evaluator_count(&self) -> usize {
        self.original_integrand.generic_evaluator_count()
            + self.threshold_counterterm.generic_evaluator_count()
    }

    #[instrument(
          level = "debug",
          skip_all,
          fields(
              term.name = %self.name(),
          )
    )]
    fn generate_event<T: FloatLike>(
        &self,
        settings: &RuntimeSettings,
        orientation_id: Option<usize>,
        channel_id: Option<SamplingChannelId>,
    ) -> Result<GenericEvent<T>> {
        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::Amplitude(
                &self.master_external_signature,
            ))
            .with_context(|| {
                format!(
                    "when getting master externals to build amplitude event for graph: {}",
                    self.graph.name
                )
            })?;

        if externals.len() != self.master_external_pdgs.len() {
            return Err(eyre!(
                "Amplitude graph '{}' has inconsistent master external metadata: {} momenta vs {} PDGs.",
                self.graph.name,
                externals.len(),
                self.master_external_pdgs.len()
            ));
        }

        let mut event = GenericEvent::default();
        event.cut_info.cut_id = 0;
        event.cut_info.orientation_id = if self.explicit_orientation_sum_only {
            Some(0)
        } else {
            orientation_id
        };
        event.cut_info.sampling_channel_id = channel_id.map(usize::from);
        event.cut_info.sampling_channel_edge_ids = channel_id
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

        for ((sign, momentum), pdg) in self
            .master_external_signature
            .iter()
            .zip(externals.iter())
            .zip(self.master_external_pdgs.iter().copied())
        {
            match sign {
                SignOrZero::Plus => {
                    event.kinematic_configuration.0.push(momentum.clone());
                    event.cut_info.particle_pdgs.0.push(pdg);
                }
                SignOrZero::Minus => {
                    event.kinematic_configuration.1.push(momentum.clone());
                    event.cut_info.particle_pdgs.1.push(pdg);
                }
                SignOrZero::Zero => {
                    return Err(eyre!(
                        "Amplitude graph '{}' has an invalid zero-sign external momentum in its master signature.",
                        self.graph.name
                    ));
                }
            }
        }

        Ok(event)
    }

    fn sampling_target_surface<T: FloatLike>(
        &self,
        channel_name: &str,
        edges: &[usize],
        externals: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Result<(&Esurface, Option<&Esurface>)> {
        let zero = F::<T>::default().zero();
        let candidates = self
            .esurfaces
            .iter_enumerated()
            .filter(|(_, surface)| {
                surface
                    .energies
                    .iter()
                    .map(|edge| edge.0)
                    .sorted()
                    .eq(edges.iter().copied())
            })
            .collect_vec();
        let Some((_, representative)) = candidates.first() else {
            return Err(eyre!(
                "amplitude surface channel '{}' for graph '{}' requests unknown energy edges {:?}; catalogue candidates (ID, edges, external shift): {:?}",
                channel_name,
                self.graph.name,
                edges,
                self.esurfaces
                    .iter_enumerated()
                    .map(|(id, surface)| (id.0, &surface.energies, &surface.external_shift))
                    .collect_vec(),
            ));
        };
        // The catalogue also includes positive-shift counterparts, which
        // cannot bound a threshold at these kinematics. Keep the named
        // channel when none is eligible, using a normalized rootless map.
        let eligible = candidates
                    .iter()
                    .copied()
                    .map(|(id, surface)| {
                        let shift = surface.compute_shift_part_from_momenta(externals, lmb);
                        if !shift.0.is_finite() {
                            return Err(super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                                operation: "amplitude sampling external shift",
                                detail: format!("candidate {} with energy edges {:?} has nonfinite derived shift {shift}", id.0, surface.energies),
                            }.into());
                        }
                        Ok((shift < zero).then_some((id, surface)))
                    })
                    .collect::<Result<Vec<_>>>()?
                    .into_iter()
                    .flatten()
                    .unique_by(|(_, surface)| surface.external_shift.iter().sorted().collect_vec())
                    .collect_vec();
        if eligible.len() > 1 {
            return Err(eyre!(
                "amplitude surface channel '{}' for graph '{}' is ambiguous for energy edges {:?}: eligible catalogue candidates (ID, external shift, evaluated shift, equation) {:?}; distinct shifts cannot share an edge-only selector",
                channel_name,
                self.graph.name,
                edges,
                eligible
                    .iter()
                    .map(|(id, surface)| (
                        id.0,
                        &surface.external_shift,
                        surface.compute_shift_part_from_momenta(externals, lmb).0,
                        surface.to_atom(&[]).to_string(),
                    ))
                    .collect_vec(),
            ));
        }
        Ok((
            *representative,
            eligible.first().map(|(_, surface)| *surface),
        ))
    }

    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        context: &mut GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<AmplitudeGraphTermEvaluation<T>> {
        let prefactor = momentum_sample.one();

        let hel = context.settings.kinematics.externals.get_helicities();
        let orientations =
            momentum_sample.orientations(&self.orientation_filter, &self.orientations);

        debug!("loop_moms: {}", momentum_sample.loop_moms());
        debug!("jacobian: {:16e}", momentum_sample.jacobian());
        // debug!("Og paramBuilder: \n{}", self.param_builder.table());

        let input = T::get_parameters(
            &mut self.param_builder,
            (
                context.settings.general.enable_cache,
                context.settings.general.debug_cache,
            ),
            &self.graph,
            momentum_sample,
            hel,
            &context.settings.additional_params(),
            None,
            None,
            None,
        );
        let result = self
            .original_integrand
            .evaluate(
                input,
                orientations,
                context.settings,
                context.evaluation_metadata,
            )?
            .pop()
            .unwrap()
            .unwrap_real();
        // debug!("parambuilder 244: {}", self.param_builder);
        let counterterm_evaluation = self.threshold_counterterm.evaluate(
            momentum_sample,
            &self.graph,
            context.model,
            &self.esurfaces,
            context.rotation,
            context.settings,
            &mut self.param_builder,
            orientations,
            context.evaluation_metadata,
            context.settings.general.store_additional_weights_in_event
                && self.threshold_counterterm.metadata_registry.is_some(),
        )?;
        let sum_of_cts = counterterm_evaluation.total.clone();

        crate::debug_tags!(#integration, #subtraction;
            bare_cff = format!("{result:16e}"),
            "{}: {result:16e}",
            self.graph.name
        );
        crate::debug_tags!(#integration, #subtraction;
            cts = format!("{sum_of_cts:16e}"),
            "{}",
            self.graph.name
        );
        crate::debug_tags!(#integration, #subtraction; "result: {result:16e}");
        crate::debug_tags!(#integration, #subtraction; "sum_of_cts: {sum_of_cts:16e}");

        crate::debug_tags!(#integration, #subtraction;
            value = format!("{sum_of_cts:16e}"),
            "evaluated sum of threshold counterterms"
        );

        let diff = result.clone() - sum_of_cts.clone();
        let result_is_nonfinite = result.re.is_nan()
            || result.re.is_infinite()
            || result.im.is_nan()
            || result.im.is_infinite();
        let sum_of_cts_is_nonfinite = sum_of_cts.re.is_nan()
            || sum_of_cts.re.is_infinite()
            || sum_of_cts.im.is_nan()
            || sum_of_cts.im.is_infinite();
        let diff_is_nonfinite =
            diff.re.is_nan() || diff.re.is_infinite() || diff.im.is_nan() || diff.im.is_infinite();
        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
            stage = "amplitude_threshold_subtraction",
            graph = %self.graph.name,
            original = %format!("{:+16e}", result),
            cts = %format!("{:+16e}", sum_of_cts),
            diff = %format!("{:+16e}", diff),
            prefactor = %format!("{:+16e}", prefactor),
            original_nonfinite = result_is_nonfinite,
            cts_nonfinite = sum_of_cts_is_nonfinite,
            diff_nonfinite = diff_is_nonfinite,
            "amplitude threshold subtraction"
        );

        let component_factor = -prefactor.clone();
        let components = counterterm_evaluation.components.map(|components| {
            components
                .into_iter()
                .map(|mut component| {
                    if let Some(bare) = &mut component.bare {
                        *bare *= &component_factor;
                    }
                    component.weighted *= &component_factor;
                    component
                })
                .collect::<Vec<_>>()
        });
        let detailed_original = components
            .as_ref()
            .map(|_| result.clone() * prefactor.clone());
        let normalized_result = if let (Some(original), Some(components)) =
            (detailed_original.as_ref(), components.as_ref())
        {
            components
                .iter()
                .fold(original.clone(), |total, component| {
                    total + &component.weighted
                })
        } else {
            diff * prefactor.clone()
        };

        Ok(AmplitudeGraphTermEvaluation {
            integrand_result: normalized_result,
            detailed_original,
            counterterms: AmplitudeCountertermEvaluation {
                total: sum_of_cts * prefactor.clone(),
                local_counterterms: counterterm_evaluation
                    .local_counterterms
                    .into_iter()
                    .map(
                        |AmplitudeLocalCountertermEvaluation {
                             variant_id,
                             esurface_id,
                             overlap_group,
                             value,
                         }| AmplitudeLocalCountertermEvaluation {
                            variant_id,
                            esurface_id,
                            overlap_group,
                            value: value * prefactor.clone(),
                        },
                    )
                    .collect(),
                components,
            },
        })
    }
}

impl GraphTerm for AmplitudeGraphTerm {
    #[instrument(
          skip_all,
          fields(
              term.name = %self.name(),
          ),
          err
    )]
    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()> {
        self.multi_channeling_setup.master_edge_masses.invalidate();
        self.multi_channeling_setup.invalidate_sampling();
        if self.explicit_orientation_sum_only {
            self.orientation_filter = SubSet::full(self.orientations.len());
        } else {
            self.orientation_filter = SubSet::empty(self.orientations.len());
            for (id, o) in self.orientations.iter_enumerated() {
                if settings.general.orientation_pat.filter(o) {
                    self.orientation_filter.add(id);
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

        self.estimated_scale = Some(
            self.graph
                .expected_scale(F(settings.kinematics.e_cm), model),
        );

        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::Amplitude(
                &self.graph.get_external_signature(),
            ))
            .with_context(|| {
                format!("when getting externals to build amplitude graph term for integrand for graph: {}", self.graph.name)
            })?;

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
            let mut polarization_start = self
                .graph
                .param_builder
                .pairs
                .polarizations
                .value_range
                .start
                * multiplicative_offset;

            for pol_value in pols.iter() {
                values[polarization_start] = *pol_value;
                polarization_start += multiplicative_offset;
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
        self.multi_channeling_setup.warm_up_masses(settings, model);

        if matches!(&settings.sampling,
            SamplingSettings::DiscreteGraphs(discrete)
                if matches!(discrete.sampling_type, DiscreteGraphSamplingType::TropicalSampling(_)))
        {
            // Tropical compensation previously used get_energy_cache, which
            // rejects complex masses on paired edges. Preserve that contract
            // before caching only the real parts for canonical preparation.
            for (pair, edge_id, edge) in self.graph.iter_edges() {
                if pair.is_paired()
                    && let Some(mass) = edge.data.mass_value::<f64>(model, &self.param_builder)
                    && mass.im != mass.im.zero()
                {
                    return Err(eyre!(
                        "tropical sampling graph '{}' requires real masses; edge {} has mass {}",
                        self.graph.name,
                        edge_id,
                        mass
                    ));
                }
            }
        }
        let masses = self
            .graph
            .new_edgevec(|e, _, _| e.mass_value(model, &self.param_builder).map(|c| c.re));

        self.real_mass_vec = Some(masses);

        Ok(())
    }

    fn name(&self) -> String {
        self.graph.name.clone()
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

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn sampling_setup(&self) -> &LmbMultiChannelingSetup {
        &self.multi_channeling_setup
    }

    fn sampling_setup_mut(&mut self) -> &mut LmbMultiChannelingSetup {
        &mut self.multi_channeling_setup
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
        if let Some(channel) = catalogue
            .named_entries()
            .find(|channel| channel.definition.radial_profile.is_some())
        {
            return Err(eyre!(
                "sampling channel '{}' for amplitude '{}' requests a LU h profile, but amplitudes have no auxiliary Cutkosky-cut LU scale",
                channel.name,
                self.graph.name
            ));
        }
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
        let surface_channels = catalogue
            .entries
            .iter()
            .zip(programs)
            .filter_map(|(entry, programs)| match entry {
                super::sampling_selection::SamplingCatalogueEntry::Named(channel) => {
                    Some((channel, &programs.2))
                }
                _ => None,
            })
            .flat_map(|(channel, joint_program)| {
                channel
                    .blocks
                    .iter()
                    .filter(|block| !block.target.energy_edge_sets().is_empty())
                    .map(move |block| (channel, block, joint_program))
            })
            .collect_vec();
        if !surface_channels.is_empty() {
            let lmb = &self.graph.loop_momentum_basis;
            let master_lmb = &self.multi_channeling_setup.graph.loop_momentum_basis;
            if lmb != master_lmb {
                return Err(eyre!(
                    "amplitude surface sampling for graph '{}' requires the master's complete loop/external routing; graph parent {:?}, master '{}' parent {:?}; identical edge names alone do not identify the same cycles",
                    self.graph.name,
                    lmb.loop_edges,
                    self.multi_channeling_setup.graph.name,
                    master_lmb.loop_edges,
                ));
            }
            if external_momenta.len() != lmb.ext_edges.len()
                || external_momenta
                    .iter()
                    .flatten()
                    .any(|value| !value.is_finite())
            {
                return Err(eyre!(
                    "amplitude surface sampling for graph '{}' needs {} finite external four-momenta (including the dependent momentum), received {:?}",
                    self.graph.name,
                    lmb.ext_edges.len(),
                    external_momenta,
                ));
            }
            let cached_masses = self.real_mass_vec.as_ref().ok_or_else(|| {
                eyre!(
                    "amplitude surface sampling for graph '{}' requires warmup mass data",
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
            let origin = LoopMomenta::from_iter(
                (0..context.n_loop_momenta)
                    .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            for (channel, block, joint_program) in surface_channels {
                if !matches!(
                    block.target,
                    SamplingMapDefinition::Surface(_) | SamplingMapDefinition::Intersect(_)
                ) {
                    return Err(eyre!(
                        "amplitude sampling channel '{}' cannot use a Cutkosky host or side qualifier: {:?}",
                        channel.name,
                        block.target,
                    ));
                }
                let energy_sets = block.target.energy_edge_sets();
                let lmbs = TiVec::from(vec![
                    self.multi_channeling_setup
                        .sampling_parent_lmb(&channel.definition.parent_lmb)?,
                ]);
                let parent_id = LmbIndex::from(0);
                let lmb = &lmbs[parent_id];
                let subspace = crate::momentum::sample::SubspaceData::new_from_parent_basis_edges(
                    &block
                        .active_lmb
                        .iter()
                        .copied()
                        .map(EdgeIndex)
                        .collect_vec(),
                    &self.graph.full_filter(),
                    parent_id,
                    &self.graph,
                    &lmbs,
                )?;
                let active = subspace.iter_lmb_indices().collect_vec();
                let active_edges = subspace
                    .iter_basis_edges(&lmbs)
                    .map(|edge| edge.0)
                    .collect_vec();
                let complement = block
                    .preceding_lmb
                    .iter()
                    .map(|edge| {
                        lmb.loop_edges
                            .iter_enumerated()
                            .find(|(_, candidate)| candidate.0 == *edge)
                            .map(|(index, _)| index)
                            .ok_or_else(|| {
                                eyre!(
                                    "sampling prerequisite edge {edge} is absent from parent {:?}",
                                    channel.definition.parent_lmb
                                )
                            })
                    })
                    .collect::<Result<Vec<_>>>()?;
                let resolve_surface = |edges: &[usize]| -> Result<_> {
                    let (representative, eligible) =
                        self.sampling_target_surface(&channel.name, edges, &externals, lmb)?;
                    let unsampled_dependencies = block
                        .remaining_lmb
                        .iter()
                        .filter(|edge| {
                            let index = lmb
                                .loop_edges
                                .iter_enumerated()
                                .find(|(_, candidate)| candidate.0 == **edge)
                                .unwrap()
                                .0;
                            representative.energies.iter().any(|energy| {
                                lmb.edge_signatures[*energy].internal[index] != SignOrZero::Zero
                            })
                        })
                        .copied()
                        .collect_vec();
                    if !unsampled_dependencies.is_empty() {
                        return Err(eyre!(
                            "sampling product/ordered block '{}' target {:?} depends on unsampled parent edges {:?}; sample these prerequisites in a preceding block",
                            channel.name,
                            block.target,
                            unsampled_dependencies,
                        ));
                    }
                    let rows = representative
                        .energies
                        .iter()
                        .map(|edge| {
                            active
                                .iter()
                                .map(|index| match lmb.edge_signatures[*edge].internal[*index] {
                                    SignOrZero::Minus => -1_i64,
                                    SignOrZero::Zero => 0,
                                    SignOrZero::Plus => 1,
                                })
                                .collect_vec()
                        })
                        .collect_vec();
                    let rank = rank_i64(&rows);
                    if rank != active.len() {
                        return Err(eyre!(
                            "amplitude surface channel '{}' for graph '{}' has energy edges {:?} with routing rows {:?}, rank {}, but active subspace {:?} requires rank {}; spectator directions belong in a complement block",
                            channel.name,
                            self.graph.name,
                            edges,
                            rows,
                            rank,
                            active_edges,
                            active.len(),
                        ));
                    }
                    Ok((representative, eligible))
                };
                let beta = e_cm * parameterization_settings.b;
                if energy_sets.len() == 2 {
                    use super::{
                        SamplingMapAffine, SamplingMapComposition, SamplingMapEmbedding,
                        SharedEnergyJointMap,
                    };
                    use std::sync::Arc;
                    let surfaces = energy_sets
                        .iter()
                        .map(|edges| {
                            let (representative, eligible) = resolve_surface(edges)?;
                            Ok(eligible.unwrap_or(representative))
                        })
                        .collect::<Result<Vec<_>>>()?;
                    let (geometry, common) = surfaces[0].sampling_joint_geometry_in_subspace(
                        surfaces[1],
                        &subspace,
                        &lmbs,
                        &self.graph,
                        &masses,
                        &externals,
                        &complement,
                    )?;
                    // Both normals are energy residuals. The existing b setting
                    // steers the trial radius and ordinary fallback scale; power
                    // does not introduce anisotropy or change the uniform-R law.
                    let joint = SharedEnergyJointMap::new(
                        geometry,
                        3 * complement.len(),
                        F::<T>::from_f64(beta).0,
                        zero.one().0,
                        beta,
                        joint_program
                            .as_ref()
                            .ok_or_else(|| {
                                eyre!(
                                    "joint channel '{}' has no cached compiled program",
                                    channel.name,
                                )
                            })?
                            .clone(),
                    )?;
                    let prior_indices = complement.clone();
                    let loop_count = lmb.loop_edges.len();
                    let spatial = externals
                        .iter()
                        .map(|p| p.spatial.clone())
                        .collect::<crate::momentum::sample::ExternalThreeMomenta<F<T>>>();
                    let native_zero = zero.clone();
                    let transform = Arc::new(move |context: &mut SamplingMapContext<'_, T>| {
                        let previous = context.previous;
                        if previous.len() != 3 * prior_indices.len() {
                            return Err(eyre!(
                                "joint affine context has {} components, expected {}",
                                previous.len(),
                                3 * prior_indices.len()
                            ));
                        }
                        let mut loops = LoopMomenta::from_iter((0..loop_count).map(|_| {
                            ThreeMomentum::new(
                                native_zero.clone(),
                                native_zero.clone(),
                                native_zero.clone(),
                            )
                        }));
                        for (&index, point) in prior_indices.iter().zip(previous.chunks_exact(3)) {
                            loops[index] = ThreeMomentum::new(
                                F(point[0].clone()),
                                F(point[1].clone()),
                                F(point[2].clone()),
                            );
                        }
                        let offset: ThreeMomentum<F<T>> = common.compute_momentum(&loops, &spatial);
                        let translation =
                            [-offset.px, -offset.py, -offset.pz].map(|x| x.0).to_vec();
                        if translation.iter().any(|x| !x.is_finite()) {
                            return Err(super::sampling_maps::SamplingEvaluationError::Unrepresentable {
                                operation: "joint affine offset", detail: "native shared-energy routing produced a nonfinite translation".into(),
                            }.into());
                        }
                        let matrix = (0..3)
                            .map(|row| {
                                (0..3)
                                    .map(|column| {
                                        if row == column {
                                            native_zero.one().0
                                        } else {
                                            native_zero.0.clone()
                                        }
                                    })
                                    .collect()
                            })
                            .collect();
                        Ok((
                            previous.to_vec(),
                            SamplingMapAffine::new(matrix, translation)?,
                        ))
                    });
                    let map = if complement.is_empty() {
                        // With no prior coordinates the complete common-energy
                        // shift is immutable: bind the existing affine owner once.
                        let (_, frame) = transform(&mut SamplingMapContext::detached(&[]))?;
                        CompiledSamplingMap::Affine {
                            map: Box::new(CompiledSamplingMap::Joint(joint)),
                            frame,
                        }
                    } else {
                        CompiledSamplingMap::Embedded(
                            SamplingMapEmbedding::from_composition(
                                SamplingMapComposition::then(vec![Box::new(joint)])?,
                                (0..3).collect(),
                            )?
                            .with_context_transform(transform),
                        )
                    };
                    context.insert_geometry_map(
                        block.target.clone(),
                        channel.definition.parent_lmb.clone(),
                        active_edges,
                        block.preceding_lmb.clone(),
                        map,
                    )?;
                    continue;
                }
                let edges = energy_sets[0];
                let (_, eligible) = resolve_surface(edges)?;
                if active.len() < context.n_loop_momenta {
                    // Eligibility fixes the routed equation once; its body and
                    // interior center are classified for each sampled complement.
                    let Some(surface) = eligible else {
                        // Nonnegative external shifts certify no open interior,
                        // including exact massless pinches, without a numerical
                        // rediscovery of the minimum at every complement point.
                        context.insert_geometry_map(
                            block.target.clone(),
                            channel.definition.parent_lmb.clone(),
                            active_edges,
                            block.preceding_lmb.clone(),
                            CompiledSamplingMap::Surface(SurfaceRadialMap::new(
                                3 * active.len(),
                                vec![zero.0.clone(); 3 * active.len()],
                                None,
                                beta,
                                parameterization_settings.power,
                            )?),
                        )?;
                        continue;
                    };
                    let map = surface.sampling_radial_map_in_subspace(
                        &subspace,
                        &lmbs,
                        &self.graph,
                        &masses,
                        &externals,
                        &complement,
                        settings,
                        beta,
                        parameterization_settings.power,
                    )?;
                    context.insert_geometry_map(
                        block.target.clone(),
                        channel.definition.parent_lmb.clone(),
                        active_edges,
                        block.preceding_lmb.clone(),
                        CompiledSamplingMap::ImplicitSurface(map),
                    )?;
                    continue;
                }
                let existing = eligible.filter(|surface| {
                    surface
                        .classify_existence(
                            &externals,
                            lmb,
                            &masses,
                            &F::<T>::from_f64(e_cm),
                            &F::<T>::from_f64(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
                        )
                        .is_existing()
                });
                if let Some(surface) = existing {
                    let origin_value =
                        surface.compute_from_momenta(lmb, &masses, &origin, &externals);
                    let map = if origin_value.0.is_finite() && origin_value < zero {
                        surface.sampling_radial_map(
                            lmb,
                            &masses,
                            &externals,
                            beta,
                            parameterization_settings.power,
                        )?
                    } else {
                        // The complete-space chart has no sampled complement:
                        // its generic interior center is prepared once at binding.
                        surface.sampling_radial_map_in_subspace(
                            &subspace,
                            &lmbs,
                            &self.graph,
                            &masses,
                            &externals,
                            &[],
                            settings,
                            beta,
                            parameterization_settings.power,
                        )?
                    };
                    context.insert_geometry_map(
                        block.target.clone(),
                        channel.definition.parent_lmb.clone(),
                        active_edges,
                        block.preceding_lmb.clone(),
                        CompiledSamplingMap::ImplicitSurface(map),
                    )?;
                } else {
                    context.insert_geometry_map(
                        block.target.clone(),
                        channel.definition.parent_lmb.clone(),
                        active_edges,
                        block.preceding_lmb.clone(),
                        CompiledSamplingMap::Surface(SurfaceRadialMap::new(
                            3 * context.n_loop_momenta,
                            vec![zero.0.clone(); 3 * context.n_loop_momenta],
                            None,
                            beta,
                            parameterization_settings.power,
                        )?),
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

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        mut context: GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<GraphEvaluationResult<T>> {
        // The canonical row fixes the proposal. Check its original joint
        // equations before physical evaluation; this work and failures belong
        // to sampling time, including every native lane and probe rotation.
        let alignment_start = Instant::now();
        let alignment =
            (|| -> Result<()> {
                let Some(SamplingCatalogueEntry::Named(channel)) =
                    context.sampling_channel.and_then(|id| {
                        self.multi_channeling_setup
                            .sampling_catalogue
                            .as_ref()
                            .and_then(|catalogue| catalogue.entries.get(id.index()))
                    })
                else {
                    return Ok(());
                };
                let mut blocks = channel
                    .blocks
                    .iter()
                    .filter_map(|block| {
                        let sets = block.target.energy_edge_sets();
                        (sets.len() == 2).then_some((block, sets))
                    })
                    .peekable();
                if blocks.peek().is_none() {
                    return Ok(());
                }
                let canonical = context.canonical_sample.ok_or_else(|| {
                    eyre!("amplitude joint alignment requires its retained canonical row")
                })?;
                if canonical.graph_id != context.graph_id
                    || canonical.channel_id != context.sampling_channel
                {
                    return Err(eyre!(
                        "amplitude joint alignment has an incompatible canonical graph/channel row"
                    ));
                }
                let parent = self
                    .multi_channeling_setup
                    .sampling_parent_lmb(&channel.definition.parent_lmb)?;
                let canonical_point = canonical.sample.rotate(context.rotation, 0, 0);
                let cached_masses = self
                    .real_mass_vec
                    .as_ref()
                    .ok_or_else(|| eyre!("amplitude joint alignment requires warmup mass data"))?;
                let canonical_masses = self.graph.new_edgevec(|_, edge, _| {
                    cached_masses[edge]
                        .map(F::<ArbPrec>::from_ff64)
                        .unwrap_or_else(|| canonical_point.zero())
                });
                let native_masses = self.graph.get_real_mass_vector(context.model);
                for (block, sets) in blocks {
                    let targets = sets
                        .iter()
                        .map(|edges| {
                            // Resolve with the original canonical external data so a
                            // native retry cannot choose another shift counterpart.
                            let (representative, eligible) = self.sampling_target_surface(
                                &channel.name,
                                edges,
                                canonical.sample.external_moms(),
                                &parent,
                            )?;
                            Ok(eligible.unwrap_or(representative))
                        })
                        .collect::<Result<Vec<_>>>()?;
                    let original = |target: &Esurface| {
                        target.evaluate_routed_enclosed(
                            &canonical_point.one(),
                            canonical_point.loop_moms(),
                            canonical_point.external_moms(),
                            &canonical_masses,
                            &self.graph.loop_momentum_basis,
                        )
                    };
                    let materialized = |target: &Esurface| {
                        target.evaluate_routed_enclosed(
                            &momentum_sample.one(),
                            momentum_sample.loop_moms(),
                            momentum_sample.external_moms(),
                            &native_masses,
                            &self.graph.loop_momentum_basis,
                        )
                    };
                    // Amplitudes have no LU host constraint. The exact zero defect
                    // leaves the shared normal half-budget intact.
                    EsurfaceRay::<T>::verify_normal_alignment(
                    [original(targets[0])?, original(targets[1])?],
                    [materialized(targets[0])?, materialized(targets[1])?],
                    std::array::from_fn(|_| rug::Float::with_val(2048, 0)),
                    context.sampling_accuracy_budget,
                ).wrap_err_with(|| format!(
                    "amplitude graph '{}' channel '{}' target {:?}, parent {:?}, rotation {}",
                    self.graph.name, channel.name, block.target, channel.definition.parent_lmb,
                    context.rotation.method,
                ))?;
                }
                Ok(())
            })();
        context.evaluation_metadata.parameterization_time += alignment_start.elapsed();
        alignment?;

        let event_channel_id = context.sampling_channel;
        let prepared_event = prepare_buffered_event(
            context.settings,
            context.rotation,
            context.event_processing_runtime.take(),
            || {
                self.generate_event(
                    context.settings,
                    momentum_sample.sample.orientation,
                    event_channel_id,
                )
            },
        )?;
        if !prepared_event.selectors_pass {
            return Ok(GraphEvaluationResult {
                reference_moments: None,
                absolute_integrand_result: None,
                integrand_result: Complex::new_re(momentum_sample.zero()),
                event_groups: crate::observables::GenericEventGroupList::default(),
                event_processing_time: prepared_event.event_processing_time,
                generated_event_count: prepared_event.generated_event_count,
                accepted_event_count: prepared_event.accepted_event_count,
            });
        }

        let AmplitudeGraphTermEvaluation {
            integrand_result,
            detailed_original,
            counterterms: mut counterterm_evaluation,
        } = self.evaluate_impl(momentum_sample, &mut context)?;

        let mut event_groups = crate::observables::GenericEventGroupList::default();
        let generated_event_count = prepared_event.generated_event_count;
        let accepted_event_count = prepared_event.accepted_event_count;
        if let Some(mut event) = prepared_event.buffered_event {
            event.weight = integrand_result.clone();

            if context.settings.general.store_additional_weights_in_event {
                let original = detailed_original
                    .unwrap_or_else(|| integrand_result.clone() + &counterterm_evaluation.total);
                event
                    .additional_weights
                    .weights
                    .insert(AdditionalWeightKey::Original, original.clone());
                if let Some(components) = counterterm_evaluation.components.take() {
                    let registry = self
                        .threshold_counterterm
                        .metadata_registry
                        .as_ref()
                        .expect("recorded amplitude components require static metadata");
                    let decomposition =
                        amplitude_threshold_event_info(registry, original, components)?;
                    if decomposition.total() != integrand_result {
                        return Err(eyre!(
                            "Amplitude graph '{}' threshold decomposition does not reconcile with its event total",
                            self.graph.name,
                        ));
                    }
                    event.additional_weights.threshold_counterterms = Some(decomposition);
                }
                for threshold_counterterm in counterterm_evaluation.local_counterterms.into_iter() {
                    let key = if let Some(variant_id) = threshold_counterterm.variant_id {
                        AdditionalWeightKey::AmplitudeThresholdCountertermVariant {
                            variant_id: variant_id.0,
                            esurface_id: threshold_counterterm.esurface_id.0,
                            overlap_group: threshold_counterterm.overlap_group,
                        }
                    } else {
                        AdditionalWeightKey::AmplitudeThresholdCounterterm {
                            esurface_id: threshold_counterterm.esurface_id.0,
                            overlap_group: threshold_counterterm.overlap_group,
                        }
                    };
                    event
                        .additional_weights
                        .weights
                        .insert(key, -threshold_counterterm.value);
                }
            }

            event_groups.push_singleton(event);
        }

        Ok(GraphEvaluationResult {
            reference_moments: None,
            absolute_integrand_result: None,
            integrand_result,
            event_groups,
            event_processing_time: prepared_event.event_processing_time,
            generated_event_count,
            accepted_event_count,
        })
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

    fn get_tropical_sampler(&self) -> &SampleGenerator<3> {
        self.tropical_sampler
            .as_ref()
            .expect("Tropical sampler should be set.")
    }

    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64> {
        &mut self.param_builder
    }

    fn get_real_mass_vector(&self) -> Result<EdgeVec<Option<F<f64>>>> {
        self.real_mass_vec
            .as_ref()
            .cloned()
            .ok_or_else(|| eyre!("real mass vector is not initialized; call warm_up first"))
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeIntegrand {
    pub settings: RuntimeSettings,
    pub data: AmplitudeIntegrandData,
    pub(crate) event_processing_runtime: RuntimeCache<EventProcessingRuntime>,
    pub(crate) active_f64_backend: RuntimeCache<ActiveF64Backend>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeIntegrandData {
    pub rotations: Option<Vec<Rotation>>,
    pub name: String,
    pub compilation: FrozenCompilationMode,
    pub loop_cache_id: usize,
    pub external_cache_id: usize,
    /// Cache ID for the base (unrotated) external momentum configuration
    pub base_external_cache_id: usize,
    pub graph_terms: Vec<AmplitudeGraphTerm>,
    pub external_signature: SignatureLike<ExternalIndex>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub graph_to_group_id: Vec<usize>,
    pub group_derived_data: TiVec<GroupId, GroupDerivedData>,
    pub explicit_orientation_sum_only: bool,
}

pub mod export;
pub mod load;

impl AmplitudeIntegrand {
    pub(crate) fn clone_with_graph_group_selection(
        &self,
        plan: &GraphGroupSelectionPlan,
    ) -> Result<Self> {
        if !self.data.group_derived_data.is_empty()
            && self.data.group_derived_data.len() != self.data.graph_group_structure.len()
        {
            return Err(eyre!(
                "Amplitude integrand '{}' has {} graph groups but {} group-derived metadata entries.",
                self.data.name,
                self.data.graph_group_structure.len(),
                self.data.group_derived_data.len(),
            ));
        }

        let mut old_graph_to_new_group = vec![None; self.data.graph_terms.len()];
        let mut old_graph_is_master = vec![false; self.data.graph_terms.len()];
        for &old_group_id in plan.retained_group_ids() {
            let new_group_id = plan.new_group_id_for_old(old_group_id).ok_or_else(|| {
                eyre!(
                    "Graph-group selection is missing a compact id for amplitude group {}.",
                    old_group_id.0
                )
            })?;
            let group = self
                .data
                .graph_group_structure
                .get(old_group_id)
                .ok_or_else(|| {
                    eyre!(
                        "Graph-group selection refers to missing amplitude group {}.",
                        old_group_id.0
                    )
                })?;
            for old_graph_id in group {
                if old_graph_id >= old_graph_to_new_group.len() {
                    return Err(eyre!(
                        "Amplitude graph group {} refers to missing graph id {}.",
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

        let group_derived_data = if self.data.group_derived_data.is_empty() {
            TiVec::new()
        } else {
            plan.retained_group_ids()
                .iter()
                .map(|&old_group_id| self.data.group_derived_data[old_group_id].clone())
                .collect()
        };
        let graph_to_group_id = graph_to_group_id_for_group_structure(&graph_group_structure);

        Ok(Self {
            settings: self.settings.clone(),
            data: AmplitudeIntegrandData {
                rotations: self.data.rotations.clone(),
                name: self.data.name.clone(),
                compilation: self.data.compilation.clone(),
                loop_cache_id: self.data.loop_cache_id,
                external_cache_id: self.data.external_cache_id,
                base_external_cache_id: self.data.base_external_cache_id,
                graph_terms,
                external_signature: self.data.external_signature.clone(),
                graph_group_structure,
                graph_to_group_id,
                group_derived_data,
                explicit_orientation_sum_only: self.data.explicit_orientation_sum_only,
            },
            event_processing_runtime: RuntimeCache::default(),
            active_f64_backend: self.active_f64_backend.clone(),
        })
    }

    fn threshold_esurface_specifier(
        &self,
        group_id: GroupId,
        group_esurface_id: GroupEsurfaceId,
    ) -> String {
        let details = self.data.group_derived_data[group_id].esurface_map[group_esurface_id]
            .iter_enumerated()
            .filter_map(|(graph_group_pos, raised_esurface_id)| {
                raised_esurface_id.map(|raised_esurface_id| (graph_group_pos, raised_esurface_id))
            })
            .find(|(graph_group_pos, _)| {
                let graph_id = self.data.graph_group_structure[group_id][*graph_group_pos];
                self.data.graph_terms[graph_id]
                    .threshold_counterterm
                    .local_esurface_exists[group_esurface_id]
            })
            .or_else(|| {
                get_representative(
                    &self.data.group_derived_data[group_id].esurface_map[group_esurface_id],
                )
                .ok()
            })
            .map(|(graph_group_pos, raised_esurface_id)| {
                let graph_id = self.data.graph_group_structure[group_id][graph_group_pos];
                let graph_term = &self.data.graph_terms[graph_id];
                let esurface_id = graph_term.threshold_counterterm.raised_data.raised_groups
                    [raised_esurface_id]
                    .esurface_ids[0];
                let esurface = &graph_term.esurfaces[esurface_id];
                let loop_order = esurface.energies.len().saturating_sub(1);
                let edge_ids = esurface.energies.iter().map(|edge_id| edge_id.0).join(",");

                (loop_order, format!("{}({edge_ids})", group_esurface_id.0))
            });

        match details {
            Some((1, specifier)) => specifier.green().bold().to_string(),
            Some((2, specifier)) => specifier.blue().bold().to_string(),
            Some((_, specifier)) => specifier.cyan().bold().to_string(),
            None => format!("{}(?)", group_esurface_id.0).dimmed().to_string(),
        }
    }

    fn format_overlap_structure(
        &self,
        group_id: GroupId,
        overlap: &crate::subtraction::overlap::OverlapStructure,
    ) -> String {
        overlap
            .overlap_groups
            .iter()
            .map(|overlap_group| {
                let esurfaces = overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|existing_esurface_id| {
                        let group_esurface_id = overlap.existing_esurfaces[*existing_esurface_id];
                        self.threshold_esurface_specifier(group_id, group_esurface_id)
                    })
                    .join(", ");

                format!(
                    "{}: [{}]",
                    overlap_group.existing_esurfaces.len(),
                    esurfaces
                )
            })
            .join(", ")
    }

    pub(crate) fn kinematics_for_threshold_approach(
        &mut self,
        momentum_sample: &MomentumSample<ArbPrec>,
        model: &Model,
    ) -> Result<Vec<OverlapStructureWithKinematics<ArbPrec>>> {
        self.data
            .graph_terms
            .iter_mut()
            .map(|term| {
                term.kinematics_for_threshold_approach(&self.settings, model, momentum_sample)
            })
            .try_collect()
    }

    fn warn_on_off_shell_external_states(&self, model: &Model) -> Result<()> {
        let externals = self
            .settings
            .kinematics
            .externals
            .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                &self.data.external_signature,
            ))?;
        let masses = self.data.graph_terms[0].graph.get_external_masses(model);
        let particles = self.data.graph_terms[0].graph.get_external_partcles();
        let helicities = self.settings.kinematics.externals.get_helicities();

        if externals.len() != masses.len()
            || externals.len() != particles.len()
            || externals.len() != helicities.len()
        {
            return Ok(());
        }

        let one = F::<f64>::from_f64(1.0);
        let threshold = F::<f64>::from_f64(1.0e-8)
            * (F::<f64>::from_f64(self.settings.kinematics.e_cm).square() + one);

        for (external_index, (((momentum, mass), particle), helicity)) in externals
            .iter()
            .zip(masses.iter())
            .zip(particles.iter())
            .zip(helicities.iter())
            .enumerate()
        {
            if particle.is_scalar() || !matches!(helicity, Helicity::Signed(_)) {
                continue;
            }

            let expected_mass_sq = mass.square();
            let actual_mass_sq = momentum.square();
            let mass_sq_difference = actual_mass_sq - expected_mass_sq;
            let off_shellness = if mass_sq_difference < mass_sq_difference.zero() {
                -mass_sq_difference
            } else {
                mass_sq_difference
            };

            if off_shellness > threshold {
                warn!(
                    "External state {external_index} ('{}', PDG {}, helicity {}) is off shell: p^2 = {:+e}, expected m^2 = {:+e} (|Δ| = {:+e}). Fixed-helicity spinor/vector/tensor external states may be ill-defined at this kinematic point.",
                    particle.name,
                    particle.pdg_code,
                    helicity,
                    actual_mass_sq,
                    expected_mass_sq,
                    off_shellness
                );
            }
        }

        Ok(())
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
    ) -> Result<Vec<std::time::Duration>> {
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
                    let compile_started = std::time::Instant::now();
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
                Ok(vec![std::time::Duration::ZERO; self.data.graph_terms.len()])
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

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path> + Sync,
        override_existing: bool,
        thread_pool: &rayon::ThreadPool,
    ) -> Result<Vec<(String, std::time::Duration)>> {
        let frozen_mode = self.data.compilation.clone();
        let compile_times = thread_pool.install(|| {
            self.data
                .graph_terms
                .par_iter_mut()
                .map(|a| {
                    a.compile(path.as_ref(), override_existing, &frozen_mode)
                        .map(|duration| (a.graph.name.clone(), duration))
                })
                .collect::<Result<Vec<_>>>()
        })?;

        self.active_f64_backend
            .set(ActiveF64Backend::from_frozen_mode(&self.data.compilation));
        Ok(compile_times)
    }

    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self.data, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;

        // debug!("HE3");
        //
        self.settings
            .to_file(path.as_ref().join("settings.toml"), override_existing)
            .with_context(|| "Error saving settings.toml file for amplitude integrand")?;
        // debug!("HE");

        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let settings = SmartSerde::from_file(
            path.as_ref().join("settings.toml"),
            "runtime settings for amplitude integrand",
        )?;

        Ok(AmplitudeIntegrand {
            settings,
            data,
            event_processing_runtime: RuntimeCache::default(),
            active_f64_backend: RuntimeCache::default(),
        })
    }

    pub(crate) fn invalidate_runtime_caches(&mut self) {
        self.event_processing_runtime.invalidate();
        for term in &mut self.data.graph_terms {
            term.multi_channeling_setup.invalidate_sampling();
        }
    }

    pub(crate) fn get_existing_esurfaces(
        &mut self,
        model: &Model,
    ) -> TiVec<GroupId, ExistingEsurfaces> {
        let external_moms = self
            .settings
            .kinematics
            .externals
            .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                &self.data.external_signature,
            ))
            .expect("could not get externals");
        let e_cm = self.settings.kinematics.e_cm;

        let mut all_existing_esurfaces = TiVec::new();

        for (group_id, graph_group) in self.data.graph_group_structure.clone().iter_enumerated() {
            let group_esurface_map = self.data.group_derived_data[group_id].esurface_map.clone();

            for (_, graph_id) in graph_group.iter_enumerated() {
                self.data.graph_terms[graph_id]
                    .threshold_counterterm
                    .local_esurface_exists = ti_vec![false; group_esurface_map.len()];
            }

            let mut group_existing_esurfaces = ExistingEsurfaces::new();

            for (group_esurface_id, raised_esurface_map) in group_esurface_map.iter_enumerated() {
                let mapped_esurfaces = raised_esurface_map
                    .iter_enumerated()
                    .filter_map(|(graph_group_pos, option_raised_esurface_id)| {
                        option_raised_esurface_id
                            .map(|raised_esurface_id| (graph_group_pos, raised_esurface_id))
                    })
                    .collect_vec();

                let (representative_graph_group_pos, representative_raised_esurface_id) =
                    mapped_esurfaces
                        .first()
                        .copied()
                        .expect("no graph in group has this esurface, map corrupted");

                let mut representative_exists = false;
                let mut any_candidate_exists = false;

                for (graph_group_pos, raised_esurface_id) in mapped_esurfaces {
                    let graph_id = graph_group[graph_group_pos];
                    let candidate_exists = {
                        let graph_term = &self.data.graph_terms[graph_id];
                        let graph = &graph_term.graph;
                        let raised_group =
                            &graph_term.threshold_counterterm.raised_data.raised_groups
                                [raised_esurface_id];
                        let esurface_id = raised_group.esurface_ids[0];
                        let esurface = &graph_term.esurfaces[esurface_id];
                        let lmb = &graph.loop_momentum_basis;
                        let real_mass_vector = graph.get_real_mass_vector(model);
                        let candidate_existence = esurface.classify_existence(
                            &external_moms,
                            lmb,
                            &real_mass_vector,
                            &F(e_cm),
                            &F(self.settings.subtraction.esurface_existence_threshold),
                        );
                        let candidate_exists = candidate_existence.is_existing();
                        if tracing::event_enabled!(tracing::Level::DEBUG) {
                            let shift_part =
                                esurface.compute_shift_part_from_momenta(&external_moms, lmb);
                            let mass_sum: F<f64> = esurface
                                .energies
                                .iter()
                                .map(|index| &real_mass_vector[*index])
                                .fold(F::from_f64(0.0), |acc, x| acc + x);
                            let zero_vector = ThreeMomentum::new(F(0.0), F(0.0), F(0.0));
                            let shift_vector = esurface
                                .external_shift
                                .iter()
                                .map(|(index, sign)| {
                                    let external_signature = &lmb.edge_signatures[*index].external;
                                    compute_shift_part(external_signature, &external_moms).spatial
                                        * F::from_f64(*sign as f64)
                                })
                                .reduce(|acc, x| acc + x)
                                .unwrap_or(zero_vector);
                            let shift_vector_sq = shift_vector.norm_squared();
                            let existence_margin =
                                shift_part * shift_part - shift_vector_sq - mass_sum * mass_sum;
                            let lmb_reps = graph.integrand_replacement(
                                &graph.full_filter(),
                                &graph.loop_momentum_basis,
                                &[W_.x___],
                            );
                            let atom = esurface.lmb_atom_simplified(graph, &lmb_reps);
                            let raw_atom = esurface.to_atom(&[]);
                            let edge_ids = esurface
                                .energies
                                .iter()
                                .map(|edge_id| edge_id.0)
                                .collect_vec();
                            let local_esurface_ids = raised_group
                                .esurface_ids
                                .iter()
                                .map(|esurface_id| esurface_id.0)
                                .collect_vec();
                            let is_representative = graph_group_pos
                                == representative_graph_group_pos
                                && raised_esurface_id == representative_raised_esurface_id;
                            let counterterm = &graph_term.threshold_counterterm;
                            let (generated, active) = if counterterm.legacy_equivalent {
                                (
                                    counterterm.generated_mask.get(raised_esurface_id).copied(),
                                    counterterm.active_mask.get(raised_esurface_id).copied(),
                                )
                            } else {
                                let variant_ids = counterterm
                                    .variant_raised_esurfaces
                                    .iter_enumerated()
                                    .filter_map(|(variant_id, &candidate)| {
                                        (candidate == raised_esurface_id).then_some(variant_id)
                                    })
                                    .collect_vec();
                                (
                                    Some(variant_ids.iter().any(|&variant_id| {
                                        counterterm.variant_generated_mask[variant_id]
                                    })),
                                    Some(variant_ids.iter().any(|&variant_id| {
                                        counterterm.variant_active_mask[variant_id]
                                    })),
                                )
                            };

                            crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                                stage = "amplitude_threshold_esurface_candidate",
                                group_id = group_id.0,
                                group_esurface_id = group_esurface_id.0,
                                graph = %graph.name,
                                graph_group_pos = graph_group_pos.0,
                                raised_esurface_id = raised_esurface_id.0,
                                esurface_id = esurface_id.0,
                                local_esurface_ids = ?local_esurface_ids,
                                edges = ?edge_ids,
                                representative = is_representative,
                                representative_graph_group_pos = representative_graph_group_pos.0,
                                representative_raised_esurface_id = representative_raised_esurface_id.0,
                                candidate_exists,
                                candidate_status = candidate_existence.label(),
                                normalized_existence_margin = ?candidate_existence.normalized_margin(),
                                non_existing_reason = ?candidate_existence.non_existing_reason(),
                                classification = ?candidate_existence,
                                generated = ?generated,
                                active = ?active,
                                max_occurrence = raised_group.max_occurence,
                                shift_part = %format!("{:+16e}", shift_part),
                                shift_vector_sq = %format!("{:+16e}", shift_vector_sq),
                                mass_sum = %format!("{:+16e}", mass_sum),
                                existence_margin = %format!("{:+16e}", existence_margin),
                                file.atom = %atom,
                                file.raw_atom = %raw_atom,
                                "amplitude threshold esurface candidate"
                            );
                        }

                        candidate_exists
                    };

                    self.data.graph_terms[graph_id]
                        .threshold_counterterm
                        .local_esurface_exists[group_esurface_id] = candidate_exists;

                    if graph_group_pos == representative_graph_group_pos
                        && raised_esurface_id == representative_raised_esurface_id
                    {
                        representative_exists = candidate_exists;
                    }
                    any_candidate_exists |= candidate_exists;
                }

                crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                    stage = "amplitude_threshold_group_esurface_existence",
                    group_id = group_id.0,
                    group_esurface_id = group_esurface_id.0,
                    representative_exists,
                    any_candidate_exists,
                    "amplitude threshold group esurface existence"
                );

                if any_candidate_exists {
                    group_existing_esurfaces.push(group_esurface_id);
                }
            }

            all_existing_esurfaces.push(group_existing_esurfaces);
        }

        let groups_above_threshold = all_existing_esurfaces
            .iter()
            .filter(|existing_esurfaces| !existing_esurfaces.is_empty())
            .count();
        let existing_esurface_count: usize = all_existing_esurfaces
            .iter()
            .map(ExistingEsurfaces::len)
            .sum();
        let threshold_status = if existing_esurface_count == 0 {
            "below threshold"
        } else {
            "above threshold"
        };
        info!(
            integrand = %self.data.name,
            threshold_status,
            groups_above_threshold,
            group_count = all_existing_esurfaces.len(),
            existing_esurfaces = existing_esurface_count,
            "Input is {threshold_status}: {existing_esurface_count} existing threshold E-surfaces across {groups_above_threshold}/{} graph groups",
            all_existing_esurfaces.len()
        );

        all_existing_esurfaces
    }

    fn validate_runtime_threshold_counterterms(
        &self,
        existing_esurfaces: &TiVec<GroupId, ExistingEsurfaces>,
    ) -> Result<()> {
        // This is the runtime safety boundary for every generation-time threshold filter:
        // a surface required by the active masses and kinematics must have persisted symbolic
        // content, irrespective of why it was omitted while generating the integrand.
        for (group_id, group_existing_esurfaces) in existing_esurfaces.iter_enumerated() {
            for group_esurface_id in group_existing_esurfaces.iter().copied() {
                for (graph_group_pos, raised_esurface_id) in self.data.group_derived_data[group_id]
                    .esurface_map[group_esurface_id]
                    .iter_enumerated()
                    .filter_map(|(graph_group_pos, raised_esurface_id)| {
                        raised_esurface_id
                            .map(|raised_esurface_id| (graph_group_pos, raised_esurface_id))
                    })
                {
                    let graph_id = self.data.graph_group_structure[group_id][graph_group_pos];
                    let graph_term = &self.data.graph_terms[graph_id];
                    if !graph_term.threshold_counterterm.legacy_equivalent {
                        // Generalized directives are checked per projected variant at the actual
                        // evaluation sample. A geometry with no variants can be explicitly
                        // disabled and must not receive a legacy fallback here.
                        continue;
                    }
                    if !graph_term
                        .threshold_counterterm
                        .local_esurface_exists
                        .get(group_esurface_id)
                        .copied()
                        .unwrap_or(true)
                    {
                        continue;
                    }
                    let is_generated = graph_term
                        .threshold_counterterm
                        .generated_mask
                        .get(raised_esurface_id)
                        .copied()
                        .ok_or_else(|| {
                            eyre!(
                                "Threshold counterterm generation mask is inconsistent for graph '{}' and raised e-surface {}",
                                graph_term.graph.name,
                                raised_esurface_id.0
                            )
                        })?;

                    if !is_generated {
                        return Err(eyre!(
                            "Amplitude integrand '{}' was generated with specialized threshold-subtraction assumptions, but the current runtime model parameters or external kinematics require a trimmed threshold counterterm for graph '{}' and group e-surface {} ({})",
                            self.name(),
                            graph_term.graph.name,
                            group_esurface_id.0,
                            self.data.group_derived_data[group_id].esurface_atoms[group_esurface_id]
                        ))
                        .with_note(|| {
                            "Regenerate the integrand or restore compatible shared/per-integrand model parameters and external kinematics.".to_string()
                        });
                    }
                }
            }
        }

        Ok(())
    }
}

impl ProcessIntegrandImpl for AmplitudeIntegrand {
    type G = AmplitudeGraphTerm;

    fn external_cache_id(&self) -> usize {
        // info!("Getting cache id {}", self.data.external_cache_id);
        self.data.external_cache_id
    }

    fn increment_external_cache_id(&mut self, val: usize) {
        // info!(
        //     "Incrementing cache id {} by {val}",
        //     self.data.external_cache_id
        // );
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
        // info!(
        //     "Reverting external cache id from {} to base {}",
        //     self.data.external_cache_id,
        //     self.data.base_external_cache_id
        // );
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

    #[instrument(
          skip_all,
          fields(
              integrand.name = %self.name(),
          )
    )]
    fn warm_up(&mut self, model: &Model) -> Result<()> {
        self.invalidate_runtime_caches();
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
        for group in &self.data.graph_group_structure {
            if group.into_iter().all(|id| {
                self.data.graph_terms[id]
                    .threshold_counterterm
                    .legacy_equivalent
            }) {
                continue;
            }
            let mut catalogue =
                crate::subtraction::amplitude_counterterm::AmplitudeOverlapCatalogue {
                    graphs: group
                        .into_iter()
                        .map(|id| self.data.graph_terms[id].graph.clone())
                        .collect(),
                    lmbs: self.data.graph_terms[group.master()]
                        .threshold_counterterm
                        .lmbs
                        .clone(),
                    variants: Vec::new(),
                };
            for (position, graph_id) in group.iter_enumerated() {
                let term = &self.data.graph_terms[graph_id];
                let counterterm = &term.threshold_counterterm;
                for (variant_id, metadata) in counterterm.variant_metadata.iter_enumerated() {
                    if counterterm.variant_generated_mask[variant_id]
                        && counterterm.variant_active_mask[variant_id]
                    {
                        let raised = counterterm.variant_raised_esurfaces[variant_id];
                        let surface = term.esurfaces
                            [counterterm.raised_data.raised_groups[raised].esurface_ids[0]]
                            .clone();
                        catalogue
                            .variants
                            .push((position, variant_id, metadata.clone(), surface));
                    }
                }
            }
            let catalogue = std::sync::Arc::new(catalogue);
            for graph_id in group {
                self.data.graph_terms[graph_id]
                    .threshold_counterterm
                    .group_catalogue
                    .set(catalogue.clone());
            }
        }
        validate_group_orientation_catalogs(
            &self.settings,
            &self.data.graph_terms,
            &self.data.graph_group_structure,
        )?;
        let e_cm = F(self.settings.kinematics.e_cm);
        let constructor = DependentMomentaConstructor::Amplitude(&self.data.external_signature);
        let masses = self.data.graph_terms[0].graph.get_external_masses(model);

        self.settings
            .kinematics
            .externals
            .improve_and_cache(constructor, &masses, &e_cm)?;
        self.warn_on_off_shell_external_states(model)?;

        let thresholds_generated = self.data.graph_terms.iter().all(|term| {
            !term.threshold_counterterm.legacy_equivalent
                || !term.threshold_counterterm.evaluators.is_empty()
        });

        let is_tree_level = self.data.graph_terms[0].graph.get_loop_number() == 0;
        let existing_esurfaces = if thresholds_generated && !is_tree_level {
            debug!("esurface existence check");
            Some(self.get_existing_esurfaces(model))
        } else {
            None
        };

        if !thresholds_generated && !self.settings.subtraction.disable_threshold_subtraction {
            warn!(
                "Not all graphs have threshold counterterms generated, but threshold subtraction is not disabled. disable runtime threshold subtraction to remove this warning"
            );
            self.settings.subtraction.disable_threshold_subtraction = true;
        }

        if let Some(existing_esurfaces) = &existing_esurfaces {
            let existing_esurface_count: usize =
                existing_esurfaces.iter().map(ExistingEsurfaces::len).sum();
            if self.settings.subtraction.disable_threshold_subtraction
                && existing_esurface_count > 0
            {
                let groups_above_threshold = existing_esurfaces
                    .iter()
                    .filter(|existing_esurfaces| !existing_esurfaces.is_empty())
                    .count();
                warn!(
                    integrand = %self.data.name,
                    groups_above_threshold,
                    group_count = existing_esurfaces.len(),
                    existing_esurfaces = existing_esurface_count,
                    "Input is above threshold, but threshold subtraction is disabled. Turn threshold subtraction on for this input by setting subtraction.disable_threshold_subtraction=false."
                );
            }
        }

        if !self.settings.subtraction.disable_threshold_subtraction && !is_tree_level {
            let existing_esurfaces = existing_esurfaces
                .expect("threshold existence should be checked before runtime threshold setup");
            self.validate_runtime_threshold_counterterms(&existing_esurfaces)?;
            for (group_id, existing_esurfaces) in existing_esurfaces.iter_enumerated() {
                debug!(
                    "solving overlap for group {}, number of thresholds: {}",
                    group_id.0,
                    existing_esurfaces.len()
                );
                if tracing::event_enabled!(tracing::Level::DEBUG) {
                    for group_esurface_id in existing_esurfaces.iter() {
                        let Some((graph_group_pos, raised_esurface_id)) =
                            self.data.group_derived_data[group_id].esurface_map[*group_esurface_id]
                                .iter_enumerated()
                                .filter_map(|(graph_group_pos, raised_esurface_id)| {
                                    raised_esurface_id.map(|raised_esurface_id| {
                                        (graph_group_pos, raised_esurface_id)
                                    })
                                })
                                .find(|(graph_group_pos, _)| {
                                    let graph_id =
                                        self.data.graph_group_structure[group_id][*graph_group_pos];
                                    self.data.graph_terms[graph_id]
                                        .threshold_counterterm
                                        .local_esurface_exists[*group_esurface_id]
                                })
                        else {
                            continue;
                        };
                        let graph_id = self.data.graph_group_structure[group_id][graph_group_pos];
                        let graph_term = &self.data.graph_terms[graph_id];
                        let graph = &graph_term.graph;
                        let esurface_id =
                            graph_term.threshold_counterterm.raised_data.raised_groups
                                [raised_esurface_id]
                                .esurface_ids[0];
                        let lmb_reps = graph.integrand_replacement(
                            &graph.full_filter(),
                            &graph.loop_momentum_basis,
                            &[W_.x___],
                        );
                        let atom =
                            graph_term.esurfaces[esurface_id].lmb_atom_simplified(graph, &lmb_reps);
                        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                            stage = "amplitude_threshold_existing_esurface",
                            group_id = group_id.0,
                            group_esurface_id = group_esurface_id.0,
                            graph = %graph.name,
                            graph_group_pos = graph_group_pos.0,
                            raised_esurface_id = raised_esurface_id.0,
                            esurface_id = esurface_id.0,
                            file.atom = %atom,
                            "amplitude threshold existing esurface"
                        );
                    }
                }

                let legacy_existing_esurfaces: ExistingEsurfaces = existing_esurfaces
                    .iter()
                    .copied()
                    .filter(|&group_esurface_id| {
                        self.data.graph_group_structure[group_id]
                            .into_iter()
                            .any(|graph_id| {
                                let counterterm =
                                    &self.data.graph_terms[graph_id].threshold_counterterm;
                                counterterm.legacy_equivalent
                                    && counterterm.local_esurface_exists[group_esurface_id]
                            })
                    })
                    .collect();
                if legacy_existing_esurfaces.is_empty() {
                    for graph_id in self.data.graph_group_structure[group_id].into_iter() {
                        if !self.data.graph_terms[graph_id]
                            .threshold_counterterm
                            .legacy_equivalent
                        {
                            self.data.graph_terms[graph_id]
                                .threshold_counterterm
                                .overlap =
                                crate::subtraction::overlap::OverlapStructure::new_empty();
                        }
                    }
                    continue;
                }

                let graph_data = self.data.graph_group_structure[group_id]
                    .into_iter()
                    .map(|graph_id| {
                        let graph = &self.data.graph_terms[graph_id];
                        SingleGraphOverlapData {
                            lmb: &graph.graph.loop_momentum_basis,
                            esurfaces: &graph.esurfaces,
                            raised_data: &graph.threshold_counterterm.raised_data,
                            edge_masses: graph.graph.get_real_mass_vector::<f64>(model),
                        }
                    })
                    .collect();
                let local_esurface_exists = self.data.graph_group_structure[group_id]
                    .into_iter()
                    .map(|graph_id| {
                        let counterterm = &self.data.graph_terms[graph_id].threshold_counterterm;
                        if counterterm.legacy_equivalent {
                            counterterm.local_esurface_exists.clone()
                        } else {
                            ti_vec![false; counterterm.local_esurface_exists.len()]
                        }
                    })
                    .collect();

                let overlap_input = OverlapInput {
                    graph_data,
                    settings: &self.settings,
                    group_esurface_map: self.data.group_derived_data[group_id].esurface_map.clone(),
                    local_esurface_exists,
                };

                let external_moms = self
                    .settings
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                        &self.data.external_signature,
                    ))
                    .expect("could not get externals");

                let overlap = find_maximal_overlap(
                    &overlap_input,
                    &legacy_existing_esurfaces,
                    &external_moms,
                )
                .with_context(|| {
                    let readable_esurfaces = legacy_existing_esurfaces
                        .iter()
                        .map(|group_esurface_id| {
                            let (graph_group_pos, raised_esurface_id) = get_representative(
                                &self.data.group_derived_data[group_id].esurface_map
                                    [*group_esurface_id],
                            )
                            .unwrap();
                            let graph_id =
                                self.data.graph_group_structure[group_id][graph_group_pos];
                            let graph_term = &self.data.graph_terms[graph_id];
                            let graph = &graph_term.graph;
                            let esurface_id =
                                graph_term.threshold_counterterm.raised_data.raised_groups
                                    [raised_esurface_id]
                                    .esurface_ids[0];
                            let lmb_reps = graph.integrand_replacement(
                                &graph.full_filter(),
                                &graph.loop_momentum_basis,
                                &[W_.x___],
                            );

                            let esurface = &graph_term.esurfaces[esurface_id];
                            let atom = esurface.lmb_atom_simplified(graph, &lmb_reps);
                            (esurface_id, atom)
                        })
                        .collect_vec();

                    let mut msg = format!(
                        "finding overlap for group: {}, existing esurfaces:\n",
                        group_id.0
                    );

                    for readable_esurface in readable_esurfaces {
                        msg += &format!(
                            "esurface id: {}, atom: {}\n",
                            readable_esurface.0.0, readable_esurface.1
                        );
                    }

                    msg
                })?;

                info!(
                    "overlap structure of group {}: [{}]",
                    group_id.0,
                    self.format_overlap_structure(group_id, &overlap)
                );

                let loop_number = self
                    .data
                    .graph_terms
                    .first()
                    .unwrap()
                    .graph
                    .get_loop_number();
                let model_params = self
                    .data
                    .graph_terms
                    .first()
                    .unwrap()
                    .param_builder
                    .pairs
                    .model_parameters
                    .params
                    .clone();

                for graph_id in self.data.graph_group_structure[group_id].into_iter() {
                    let graph_term = &mut self.data.graph_terms[graph_id];
                    if !graph_term.threshold_counterterm.legacy_equivalent {
                        graph_term.threshold_counterterm.overlap =
                            crate::subtraction::overlap::OverlapStructure::new_empty();
                        continue;
                    }
                    let max_required_order = graph_term
                        .threshold_counterterm
                        .raised_data
                        .raised_groups
                        .iter()
                        .map(|raised_group| raised_group.max_occurence)
                        .max()
                        .unwrap_or(0);
                    let max_required_power = (2 * (max_required_order / 2 + 1)) as i32;

                    let mut localized_overlap = overlap.localized_to_existing_surfaces(
                        &graph_term.threshold_counterterm.local_esurface_exists,
                    );
                    localized_overlap
                        .build_evaluators(
                            &self.data.group_derived_data[group_id].esurface_atoms,
                            &OptimizationSettings::default(),
                            loop_number,
                            external_moms.len(),
                            model_params.clone(),
                            max_required_power,
                        )
                        .with_context(|| {
                            format!(
                                "Failed to build graph-local multi-channeling evaluators for group {} graph {}",
                                group_id.0, graph_term.graph.name
                            )
                        })?;

                    crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #overlap;
                        stage = "amplitude_threshold_localized_overlap",
                        group_id = group_id.0,
                        graph = %graph_term.graph.name,
                        existing_esurfaces = ?localized_overlap
                            .existing_esurfaces
                            .iter()
                            .map(|group_esurface_id| group_esurface_id.0)
                            .collect_vec(),
                        overlap_group_sizes = ?localized_overlap
                            .overlap_groups
                            .iter()
                            .map(|group| group.existing_esurfaces.len())
                            .collect_vec(),
                        "amplitude threshold localized overlap"
                    );

                    graph_term.threshold_counterterm.overlap = localized_overlap;
                }
            }
        }

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

    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.as_ref().expect("forgot warmup").iter()
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
    }

    fn graph_count(&self) -> usize {
        self.data.graph_terms.len()
    }

    fn get_master_graph(&self, group_id: GroupId) -> &Self::G {
        let group_master = self.data.graph_group_structure[group_id].master();
        &self.data.graph_terms[group_master]
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

    fn get_group(&self, group_id: GroupId) -> &GraphGroup {
        &self.data.graph_group_structure[group_id]
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor<'_> {
        DependentMomentaConstructor::Amplitude(&self.data.external_signature)
    }

    fn get_group_structure(&self) -> &TiVec<GroupId, GraphGroup> {
        &self.data.graph_group_structure
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

    fn groups_default_sample_events_by_graph_group(&self) -> bool {
        true
    }
    // fn get_builder_cache(&self) -> &ParamBuilder<f64> {
    //     &self.data.builder_cache
    // }
}

impl HasIntegrand for AmplitudeIntegrand {
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
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        evaluate_sample(self, model, sample, wgt, iter, use_arb_prec, max_eval)
    }

    fn get_n_dim(&self) -> usize {
        if self
            .settings
            .sampling
            .get_parameterization_settings()
            .is_some()
        {
            self.data.graph_terms[0].graph.get_loop_number() * 3
        } else {
            // let dimensions = self
            //     .data
            //     .graph_terms
            //     .iter()
            //     .map(|term| term.get_tropical_sampler().get_dimension())
            //     .sorted()
            //     .collect_vec();

            tracing::warn!(
                "get n dim called for tropical sampling, if groups are enabled this function panics, returning bs value to avoid this"
            );
            69

            //let median_dimension = dimensions[dimensions.len() / 2];
            //median_dimension
        }
    }
}

#[cfg(test)]
mod sampling_tests {
    use super::*;
    use crate::integrands::evaluation::EvaluationMetaData;
    use crate::{
        initialisation::test_initialise,
        integrands::process::{
            CompiledSamplingMap, ProcessIntegrand, SamplingMapAffine, SamplingSupport,
        },
        processes::Amplitude,
        utils::load_generic_model,
    };
    use linnet::half_edge::involution::EdgeIndex;

    #[test]
    fn generated_triangle_joint_binds_without_prerequisites() -> Result<()> {
        use crate::integrands::process::GaussianReferenceFunction;
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graphs = Graph::from_string(
            r#"digraph joint_triangle {
            node [num=1]; edge [num=1];
            incoming_a [style=invis]; incoming_b [style=invis]; outgoing [style=invis];
            incoming_a -> A:0 [id=0,particle=scalar_0];
            incoming_b -> B:1 [id=1,particle=scalar_0];
            C:2 -> outgoing [id=2,particle=scalar_0];
            A -> B [id=3,particle=scalar_1];
            B -> C [id=4,particle=scalar_1,lmb_id=0];
            C -> A [id=5,particle=scalar_1];
        }"#,
            &model,
        )?;
        let mut amplitude = Amplitude::from_graph_list("joint_triangle", graphs)?;
        // This generation checks the actual amplitude binder and reference
        // owner; it makes no physical threshold-subtraction claim.
        let global: GlobalSettings = toml::from_str(
            r#"
[generation]
override_lmb_heuristics = true
[generation.uv]
subtract_uv = false
generate_integrated = false
[generation.threshold_subtraction]
enable_thresholds = false
[generation.evaluator]
compile = false
summed = false
summed_function_map = true
iterative_orientation_optimization = false
"#,
        )?;
        let settings: RuntimeSettings = toml::from_str(
            r#"
[general]
evaluator_method = "SummedFunctionMap"
[kinematics]
e_cm = 21.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[11.0,-1.0,-4.0,0.0], [10.0,3.0,0.0,0.0], "dependent"]
helicities = [0,0,0]
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "summed"
default_channel_selection = ["joint", "ordinary"]
[sampling.channel_definitions.joint_triangle.joint]
around = "intersect(surface(3,4),surface(3,5))"
parent_lmb = [4]
subspace_lmb = [4]
[sampling.channel_definitions.joint_triangle.ordinary]
around = "lmb(4)"
parent_lmb = [4]
"#,
        )?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .stack_size(256 * 1024 * 1024)
            .build()?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(&model, "joint_triangle", &global, (&settings).into(), &pool)?;
        let runtime = amplitude.integrand.as_mut().unwrap();
        runtime.warm_up(&model)?;
        let ProcessIntegrand::Amplitude(generated) = &*runtime else {
            unreachable!()
        };
        let bridge = generated.data.graph_terms[0]
            .sampling_setup()
            .sampling_bridge::<f64>()?;
        let channel = &bridge.channels()[0];
        assert!(!channel.contract().requires_context);
        assert!(channel.contract().requires_proposal_policy);
        assert!(matches!(channel.map, CompiledSamplingMap::Affine { .. }));
        let cube = [0.19, 0.31, 0.67];
        let forward = bridge.forward(SamplingChannelId(0), &cube)?;
        assert!(
            forward
                .map
                .diagnostics
                .iter()
                .any(|d| d.contains("certified normal disk"))
        );
        let inverse = bridge
            .inverse(SamplingChannelId(0), &forward.raw_coordinates)?
            .expect("selected compact point");
        assert!((forward.map.jacobian * inverse.map.inverse_jacobian - 1.0).abs() < 1e-9);
        for (actual, expected) in inverse.map.coordinates.iter().zip(cube) {
            assert!((actual - expected).abs() < 1e-9);
        }
        let reference = GaussianReferenceFunction::new(5.0, vec![0.2, -0.1, 0.3])?;
        let source = Sample::Continuous(F(1.0), cube.map(F).to_vec());
        let value = runtime.evaluate_reference_sample_detailed(&source, &reference)?;
        assert!(!value.evaluation.evaluation_metadata.is_nan);
        assert!(value.evaluation.integrand_result.re.0 > 0.0);
        assert!(
            !value
                .evaluation
                .evaluation_metadata
                .sampling_proposal_policies
                .is_empty()
        );
        {
            use crate::{
                integrands::{
                    evaluation::PreciseEvaluationResult,
                    process::{
                        EvaluationTarget, evaluate_single,
                        gammaloop_sample::{GammaLoopSample, parameterize},
                        sampling_maps::SamplingEvaluationError,
                    },
                },
                settings::runtime::{Precision, StabilityLevelSetting},
                utils::QuadFloat,
            };
            // The original physical triangle (without threshold CTs) now also
            // exercises normal materialization through the production boundary.
            let mut physical = runtime.clone();
            physical.get_mut_settings().stability.levels = vec![
                StabilityLevelSetting::default_double(),
                StabilityLevelSetting::default_quad(),
            ];
            physical.get_mut_settings().stability.rotation_axis.clear();
            physical.warm_up(&model)?;
            let regular = physical
                .evaluate_sample_precise(&source, &model, F(1.0), false, Complex::new_zero())?
                .try_into_f64()?;
            assert!(!regular.evaluation_metadata.is_nan);
            let near = Sample::Continuous(F(1.0), vec![F(1.0e-12), F(cube[1]), F(cube[2])]);
            {
                let ProcessIntegrand::Amplitude(generated) = &mut physical else {
                    unreachable!()
                };
                let mut metadata = EvaluationMetaData::new_empty();
                metadata.sampling_proposal_policies.begin_collection();
                let canonical = parameterize::<ArbPrec, _>(&near, generated, &mut metadata)?;
                metadata.sampling_proposal_policies.seal();
                let rotation = Rotation::new(RotationMethod::EulerAngles(0.31, -0.17, 0.23));
                let double = canonical
                    .materialize::<f64>(GammaLoopSample::<f64>::relative_accuracy_budget(
                        generated.get_settings(),
                    ))?
                    .rotate(&rotation, 100, 101);
                let before = metadata.parameterization_time;
                let error = evaluate_single(
                    generated,
                    EvaluationTarget::Physical(&model),
                    &double,
                    &rotation,
                    &mut metadata,
                    Some(&canonical),
                )
                .unwrap_err();
                assert!(matches!(
                    error.downcast_ref::<SamplingEvaluationError>(),
                    Some(SamplingEvaluationError::UncertainGeometry { .. })
                ));
                assert!(
                    format!("{error:#}").contains("normal displacement"),
                    "{error:#}"
                );
                assert!(metadata.parameterization_time > before);
                let quad = canonical
                    .materialize::<QuadFloat>(
                        GammaLoopSample::<QuadFloat>::relative_accuracy_budget(
                            generated.get_settings(),
                        ),
                    )?
                    .rotate(&rotation, 102, 103);
                let value = evaluate_single(
                    generated,
                    EvaluationTarget::Physical(&model),
                    &quad,
                    &rotation,
                    &mut metadata,
                    Some(&canonical),
                )?;
                assert!(
                    !value.integrand_result.re.is_nan() && !value.integrand_result.re.is_infinite()
                );
                assert!(
                    !value.integrand_result.im.is_nan() && !value.integrand_result.im.is_infinite()
                );
            }
            let PreciseEvaluationResult::Quad(rescued) = physical.evaluate_sample_precise(
                &near,
                &model,
                F(1.0),
                false,
                Complex::new_zero(),
            )?
            else {
                panic!("the original joint normals must rescue this draw at Quad");
            };
            assert_eq!(
                rescued.evaluation_metadata.final_precision(),
                Some(Precision::Quad)
            );
            assert_eq!(rescued.evaluation_metadata.stability_results.len(), 2);
            assert!(!rescued.evaluation_metadata.is_nan);
            physical.get_mut_settings().stability.levels =
                vec![StabilityLevelSetting::default_quad()];
            physical.warm_up(&model)?;
            let PreciseEvaluationResult::Quad(forced) = physical.evaluate_sample_precise(
                &near,
                &model,
                F(1.0),
                false,
                Complex::new_zero(),
            )?
            else {
                unreachable!()
            };
            assert_eq!(rescued.integrand_result, forced.integrand_result);
            assert_eq!(
                rescued.parameterization_jacobian,
                forced.parameterization_jacobian
            );
            assert_eq!(
                rescued.evaluation_metadata.sampling_proposal_policies,
                forced.evaluation_metadata.sampling_proposal_policies
            );
        }
        {
            use crate::settings::runtime::DiscreteGraphSamplingSettings;
            let ProcessIntegrand::Amplitude(generated) = &*runtime else {
                unreachable!()
            };
            let mut term = generated.data.graph_terms[0].clone();
            let loops = LoopMomenta::from_iter(
                forward
                    .raw_coordinates
                    .chunks_exact(3)
                    .map(|p| ThreeMomentum::new(F(p[0]), F(p[1]), F(p[2]))),
            );
            let externals = settings
                .kinematics
                .externals
                .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                    &generated.data.external_signature,
                ))?;
            let old = term.graph.get_energy_cache(
                &model,
                &loops,
                &externals,
                &term.graph.loop_momentum_basis,
            );
            let masses = term.get_real_mass_vector()?;
            for (_, edge, _) in term.graph.iter_loop_edges() {
                let momentum = term.graph.loop_momentum_basis.edge_signatures[edge]
                    .compute_four_momentum_from_three(&loops, &externals);
                let current = momentum.spatial.on_shell_energy(masses[edge]).value;
                assert_eq!(
                    current, old[edge],
                    "tropical compensation energy must preserve the old cache equation"
                );
            }
            let mut tropical = settings.clone();
            tropical.sampling = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                sampling_type: DiscreteGraphSamplingType::TropicalSampling(Default::default()),
                ..Default::default()
            });
            term.warm_up(&tropical, &model)?;
            let mut complex_model = model.clone();
            complex_model.get_parameter_mut("mass_scalar_1")?.value =
                Some(Complex::new(F(1.0), F(0.1)));
            let error = term.warm_up(&tropical, &complex_model).unwrap_err();
            assert!(
                error.to_string().contains("tropical sampling")
                    && error.to_string().contains("requires real masses"),
                "{error:#}"
            );
        }

        Ok(())
    }

    #[test]
    fn generated_kite_reference_normalizes_unequal_groups_and_orientations() -> Result<()> {
        use crate::{
            integrands::process::{
                EvaluationTarget, GaussianReferenceFunction, evaluate_single,
                gammaloop_sample::{GammaLoopSample, parameterize},
            },
            settings::runtime::{
                DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, SamplingSettings,
            },
            utils::{ArbPrec, QuadFloat},
        };
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = Graph::from_string(
            include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/resources/graphs/massive_kite.dot"
            )),
            &model,
        )?
        .remove(0);
        // These identical graph copies form deliberately unequal groups. This
        // tests the normalized reference target, not a new physical graph sum.
        let graphs = (0..3)
            .map(|index| {
                let mut graph = graph.clone();
                graph.name = format!("reference_kite_{index}");
                graph.group_id = Some(GroupId(usize::from(index > 0)));
                graph.is_group_master = index != 2;
                graph
            })
            .collect();
        let mut amplitude = Amplitude::from_graph_list("reference_groups", graphs)?;
        let global: GlobalSettings = toml::from_str(
            r#"
[generation]
override_lmb_heuristics = true
[generation.uv]
subtract_uv = false
generate_integrated = false
[generation.threshold_subtraction]
enable_thresholds = false
[generation.evaluator]
compile = false
summed = false
summed_function_map = true
iterative_orientation_optimization = false
"#,
        )?;
        let settings: RuntimeSettings = toml::from_str(
            r#"
[general]
evaluator_method = "SummedFunctionMap"
[kinematics]
e_cm = 5.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0], "dependent"]
helicities = [0, 0]
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = false
"#,
        )?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .stack_size(256 * 1024 * 1024)
            .build()?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(
            &model,
            "reference_groups",
            &global,
            (&settings).into(),
            &pool,
        )?;
        let integrand = amplitude.integrand.as_mut().unwrap();
        integrand.warm_up(&model)?;
        assert_eq!(integrand.graph_count(), 3);
        assert_eq!(integrand.graph_group_count(), 2);
        let parameterization = settings.sampling.get_parameterization_settings().unwrap();
        let reference = GaussianReferenceFunction::new(2.0, vec![0.2, -0.1, 0.3, -0.2, 0.1, 0.4])?;
        let continuous = Sample::Continuous(
            F(1.0),
            vec![F(0.31), F(0.27), F(0.61), F(0.39), F(0.72), F(0.58)],
        );
        let summed = integrand.evaluate_reference_sample_detailed(&continuous, &reference)?;
        let expected = summed.evaluation.integrand_result.re.0
            * summed.evaluation.parameterization_jacobian.unwrap().0;
        let expected_moment =
            summed.moments.second_moment.0 * summed.evaluation.parameterization_jacobian.unwrap().0;
        for select_orientation in [false, true] {
            integrand.get_mut_settings().general.evaluator_method =
                crate::integrands::process::evaluators::EvaluatorMethod::SingleParametric;
            integrand.get_mut_settings().sampling =
                SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                    sample_orientations: select_orientation,
                    sampling_type: DiscreteGraphSamplingType::Default(parameterization.clone()),
                    ..Default::default()
                });
            integrand.warm_up(&model)?;
            let (mut total, mut total_moment) = (0.0, 0.0);
            for (group_id, group_probability) in [(0, 0.2), (1, 0.8)] {
                let count = if select_orientation {
                    integrand
                        .group_orientation_count(GroupId(group_id))
                        .unwrap()
                } else {
                    1
                };
                if select_orientation {
                    assert!(count > 1);
                }
                for orientation in 0..count {
                    let orientation_probability =
                        2.0 * (orientation + 1) as f64 / (count * (count + 1)) as f64;
                    let probability = group_probability * orientation_probability;
                    let inner = if select_orientation {
                        Sample::Discrete(F(1.0), orientation, Some(Box::new(continuous.clone())))
                    } else {
                        continuous.clone()
                    };
                    let sample =
                        Sample::Discrete(F(1.0 / probability), group_id, Some(Box::new(inner)));
                    let selected =
                        integrand.evaluate_reference_sample_detailed(&sample, &reference)?;
                    let factor = selected.evaluation.parameterization_jacobian.unwrap().0
                        * selected.evaluation.integrator_weight.0;
                    let value = selected.evaluation.integrand_result.re.0 * factor;
                    let group_fraction = if group_id == 0 { 1.0 / 3.0 } else { 2.0 / 3.0 };
                    assert!(
                        (value * probability / expected - group_fraction / count as f64).abs()
                            < 1.0e-12
                    );
                    total += probability * value;
                    total_moment += probability * selected.moments.second_moment.0 * factor;
                }
            }
            assert!((total / expected - 1.0).abs() < 1.0e-12);
            assert!((total_moment / expected_moment - 1.0).abs() < 1.0e-12);
        }
        // Exercise the precision-generic reference owner from the canonical
        // draw; the public acceptance API reports f64 after its physical stack.
        integrand.get_mut_settings().sampling = settings.sampling;
        integrand.warm_up(&model)?;
        let ProcessIntegrand::Amplitude(integrand) = integrand else {
            unreachable!()
        };
        let mut metadata = EvaluationMetaData::new_empty();
        metadata.sampling_proposal_policies.begin_collection();
        let canonical = parameterize::<ArbPrec, _>(&continuous, integrand, &mut metadata)?;
        metadata.sampling_proposal_policies.seal();
        let sample = canonical.materialize::<QuadFloat>(
            GammaLoopSample::<QuadFloat>::relative_accuracy_budget(integrand.get_settings()),
        )?;
        let rotation = Rotation::new(RotationMethod::Pi2X);
        let rotated = sample.rotate(&rotation, 100, 101);
        let result = evaluate_single(
            integrand,
            EvaluationTarget::Reference(&reference),
            &rotated,
            &rotation,
            &mut crate::integrands::evaluation::EvaluationMetaData::new_empty(),
            Some(&canonical),
        )?;
        assert!(
            (result.integrand_result.re.into_ff64().0 / summed.evaluation.integrand_result.re.0
                - 1.0)
                .abs()
                < 1.0e-12
        );
        assert!(
            (result
                .reference_moments
                .unwrap()
                .second_moment
                .into_ff64()
                .0
                / summed.moments.second_moment.0
                - 1.0)
                .abs()
                < 1.0e-12
        );
        Ok(())
    }

    #[test]
    fn generated_kite_sampling_uses_physical_six_dimensional_surfaces() -> Result<()> {
        test_initialise()?;
        let mut model = load_generic_model("scalars");
        let graphs = Graph::from_string(
            include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/resources/graphs/massive_kite.dot"
            )),
            &model,
        )?;
        let mut amplitude = Amplitude::from_graph_list("sampling_kite", graphs)?;
        let global: GlobalSettings = toml::from_str(
            r#"
[generation]
override_lmb_heuristics = true
[generation.uv]
subtract_uv = false
generate_integrated = false
[generation.evaluator]
compile = false
summed = false
summed_function_map = true
iterative_orientation_optimization = false
"#,
        )?;
        let settings: RuntimeSettings = toml::from_str(
            r#"
[general]
evaluator_method = "SummedFunctionMap"
integral_unit = "none"
[kinematics]
e_cm = 5.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0], "dependent"]
helicities = [0, 0]
[sampling]
graphs = "summed"
orientations = "summed"
sampling_multichanneling = true
sampling_channels = "summed"
power = 2.0
default_channel_selection = ["C", "D", "ordinary"]
[sampling.channel_definitions.massive_kite.C]
around = "surface(2,4,6)"
parent_lmb = [4,6]
subspace_lmb = [4,6]
[sampling.channel_definitions.massive_kite.D]
around = "surface(3,5,6)"
parent_lmb = [4,6]
subspace_lmb = [4,6]
[sampling.channel_definitions.massive_kite.ordinary]
around = "lmb(4,6)"
parent_lmb = [4,6]
"#,
        )?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .stack_size(256 * 1024 * 1024)
            .build()?;
        amplitude.preprocess(&model, &global.generation, &(&settings).into(), &pool)?;
        amplitude.build_integrand(&model, "sampling_kite", &global, (&settings).into(), &pool)?;
        let integrand = amplitude.integrand.as_mut().unwrap();
        integrand.warm_up(&model)?;
        // Exercise the warmed production owner on a clone, leaving the fresh
        // construction/geometry checks below independent of runtime caching.
        {
            use crate::{
                integrands::process::GaussianReferenceFunction,
                momentum::ExternalMomenta,
                settings::runtime::{SamplingSettings, kinematic::Externals},
                utils::QuadFloat,
            };
            let mut runtime = integrand.clone();
            let reference = GaussianReferenceFunction::new(2.0, vec![0.0; 6])?;
            let point = Sample::Continuous(
                F(1.0),
                vec![F(0.19), F(0.27), F(0.61), F(0.39), F(0.72), F(0.58)],
            );
            let ProcessIntegrand::Amplitude(amplitude) = &runtime else {
                unreachable!()
            };
            let setup = &amplitude.data.graph_terms[0].multi_channeling_setup;
            let address = setup.sampling_bridge::<f64>()? as *const _ as usize;
            let original = setup
                .sampling_bridge::<f64>()?
                .forward(SamplingChannelId(0), &[0.19, 0.27, 0.61, 0.39, 0.72, 0.58])?;
            let first = runtime.evaluate_reference_sample_detailed(&point, &reference)?;
            for _ in 0..3 {
                let repeated = runtime.evaluate_reference_sample_detailed(&point, &reference)?;
                assert_eq!(
                    repeated.evaluation.integrand_result,
                    first.evaluation.integrand_result
                );
                let ProcessIntegrand::Amplitude(amplitude) = &runtime else {
                    unreachable!()
                };
                assert_eq!(
                    amplitude.data.graph_terms[0]
                        .multi_channeling_setup
                        .sampling_bridge::<f64>()? as *const _ as usize,
                    address
                );
            }
            // In-place external edits keep e_cm fixed and must replace both
            // improved numeric caches before the next bridge captures them.
            let Externals::Constant { momenta, .. } =
                &mut runtime.get_mut_settings().kinematics.externals;
            momenta[0] = ExternalMomenta::Independent([26.0_f64.sqrt(), 0.0, 0.0, 1.0].map(F));
            let error = runtime
                .evaluate_reference_sample_detailed(&point, &reference)
                .unwrap_err();
            assert!(format!("{error:#}").contains("call warm_up"), "{error:#}");
            runtime.warm_up(&model)?;
            let ProcessIntegrand::Amplitude(amplitude) = &runtime else {
                unreachable!()
            };
            let constructor =
                DependentMomentaConstructor::Amplitude(&amplitude.data.external_signature);
            let external = amplitude
                .settings
                .kinematics
                .externals
                .get_dependent_externals::<f64>(constructor)?;
            let external_quad = amplitude
                .settings
                .kinematics
                .externals
                .get_dependent_externals::<QuadFloat>(constructor)?;
            assert!((external[ExternalIndex(0)].spatial.pz.0 - 1.0).abs() < 1.0e-12);
            assert!(
                (external_quad[ExternalIndex(0)]
                    .spatial
                    .pz
                    .clone()
                    .into_ff64()
                    .0
                    - 1.0)
                    .abs()
                    < 1.0e-12
            );
            let moved = amplitude.data.graph_terms[0]
                .multi_channeling_setup
                .sampling_bridge::<f64>()?
                .forward(SamplingChannelId(0), &[0.19, 0.27, 0.61, 0.39, 0.72, 0.58])?;
            assert_ne!(original.raw_coordinates, moved.raw_coordinates);
            // The public shared warmup helper also clears a previous bridge
            // when called directly and a replacement cannot be compiled.
            let ProcessIntegrand::Amplitude(amplitude) = &mut runtime else {
                unreachable!()
            };
            let masses = amplitude.data.graph_terms[0].real_mass_vec.take();
            assert!(amplitude.warm_up_sampling().is_err());
            assert!(
                amplitude.data.graph_terms[0]
                    .multi_channeling_setup
                    .sampling_bridge::<f64>()
                    .is_err()
            );
            amplitude.data.graph_terms[0].real_mass_vec = masses;
            amplitude.warm_up_sampling()?;
            // A failed warmup cannot keep the previous geometry usable.
            let SamplingSettings::MultiChanneling(channels) =
                &mut runtime.get_mut_settings().sampling
            else {
                unreachable!()
            };
            channels
                .parameterization_settings
                .sampling_channels
                .default_channel_selection = vec!["missing_channel".to_string()];
            assert!(runtime.warm_up(&model).is_err());
            let error = runtime
                .evaluate_reference_sample_detailed(&point, &reference)
                .unwrap_err();
            assert!(format!("{error:#}").contains("call warm_up"), "{error:#}");
            *runtime.get_mut_settings() = settings.clone();
            runtime.warm_up(&model)?;
            let restored = runtime.evaluate_reference_sample_detailed(&point, &reference)?;
            assert_eq!(
                restored.evaluation.integrand_result,
                first.evaluation.integrand_result
            );
        }
        // A depends only on q4, so its independent product block freezes the
        // analytic two-energy center during binding. Its boosted invariant is
        // above threshold, but the interior depth is below Double's native
        // cancellation budget. Quad binds the same physical equation.
        {
            use crate::{
                integrands::process::{
                    GaussianReferenceFunction, sampling_maps::SamplingEvaluationError,
                },
                momentum::ExternalMomenta,
                settings::runtime::{
                    Precision, SamplingSettings, StabilityLevelSetting, kinematic::Externals,
                },
                utils::QuadFloat,
            };
            let mut runtime = integrand.clone();
            let runtime_settings = runtime.get_mut_settings();
            runtime_settings.stability.rotation_axis.clear();
            runtime_settings.stability.levels = vec![
                StabilityLevelSetting::default_double(),
                StabilityLevelSetting::default_quad(),
            ];
            let SamplingSettings::MultiChanneling(channels) = &mut runtime_settings.sampling else {
                unreachable!()
            };
            channels.parameterization_settings.power = 1.0;
            channels
                .parameterization_settings
                .sampling_channels
                .default_channel_selection = vec!["frozen_A".into()];
            let definitions = channels
                .parameterization_settings
                .sampling_channels
                .channel_definitions
                .get_mut("massive_kite")
                .unwrap();
            let mut definition = definitions["C"].clone();
            definition.around = "product(block(lmb(4),surface(3,4)),complement(6))".into();
            definition.subspace_lmb = vec![4];
            definitions.insert("frozen_A".into(), definition);
            let Externals::Constant { momenta, .. } = &mut runtime_settings.kinematics.externals;
            momenta[0] = ExternalMomenta::Independent([1.0e6 + 2.1e-6, 0.0, 0.0, 1.0e6].map(F));
            runtime.warm_up(&model)?;
            let ProcessIntegrand::Amplitude(amplitude) = &mut runtime else {
                unreachable!()
            };
            let setup = &amplitude.data.graph_terms[0].multi_channeling_setup;
            let failure = setup.sampling_bridge::<f64>().unwrap_err();
            assert!(
                matches!(
                    failure.downcast_ref::<SamplingEvaluationError>(),
                    Some(SamplingEvaluationError::UncertainGeometry { .. })
                ),
                "{failure:#}"
            );
            assert!(setup.sampling_bridge::<QuadFloat>().is_ok());
            assert!(setup.sampling_catalogue.as_ref().is_some());
            // Removing an input which a fresh bind requires proves that repeated
            // requests reuse the cached result, without another center solve.
            let masses = amplitude.data.graph_terms[0].real_mass_vec.take();
            for _ in 0..3 {
                assert_eq!(
                    format!(
                        "{:#}",
                        amplitude.prepare_sampling_precision::<f64>().unwrap_err()
                    ),
                    format!("{failure:#}")
                );
                amplitude.prepare_sampling_precision::<QuadFloat>()?;
            }
            amplitude.data.graph_terms[0].real_mass_vec = masses;
            let source =
                Sample::Continuous(F(1.0), [0.19, 0.27, 0.61, 0.39, 0.72, 0.58].map(F).to_vec());
            // q3 = q4 + Q, so the equal-mass A center is q4 = -Q/2.
            // Keep this probe near that center instead of testing a negligible tail.
            let reference =
                GaussianReferenceFunction::new(5.0, vec![0.0, 0.0, -5.0e5, 0.0, 0.0, 0.0])?;
            let rescued = runtime.evaluate_reference_sample_detailed(&source, &reference)?;
            assert_eq!(
                rescued.evaluation.evaluation_metadata.final_precision(),
                Some(Precision::Double)
            );
            assert_eq!(
                rescued
                    .evaluation
                    .evaluation_metadata
                    .stability_results
                    .len(),
                1
            );
            assert!(rescued.evaluation.integrand_result.re.0 > 0.0);
            assert!(rescued.moments.second_moment.0 > 0.0);
            let mut double_only = runtime.clone();
            double_only.get_mut_settings().stability.levels =
                vec![StabilityLevelSetting::default_double()];
            double_only.warm_up(&model)?;
            let direct = double_only.evaluate_reference_sample_detailed(&source, &reference)?;
            assert_eq!(
                rescued.evaluation.integrand_result,
                direct.evaluation.integrand_result
            );
            assert_eq!(rescued.moments.second_moment, direct.moments.second_moment);
            // A valid canonical map remains usable with a Double-only physical
            // stack, even when its optional Double component bridge cannot bind.
            let ProcessIntegrand::Amplitude(amplitude) = &double_only else {
                unreachable!()
            };
            assert!(
                amplitude.data.graph_terms[0]
                    .multi_channeling_setup
                    .sampling_catalogue
                    .as_ref()
                    .is_some()
            );
            assert!(
                amplitude.data.graph_terms[0]
                    .multi_channeling_setup
                    .sampling_bridge::<f64>()
                    .is_err()
            );
            assert!(
                amplitude.data.graph_terms[0]
                    .multi_channeling_setup
                    .sampling_bridge::<ArbPrec>()
                    .is_ok()
            );
            // A structural bind failure still invalidates all precisions even
            // when the previous epoch contained a usable Quad binding.
            let ProcessIntegrand::Amplitude(amplitude) = &mut runtime else {
                unreachable!()
            };
            let masses = amplitude.data.graph_terms[0].real_mass_vec.take();
            let error = amplitude.warm_up_sampling().unwrap_err();
            assert!(error.downcast_ref::<SamplingEvaluationError>().is_none());
            assert!(
                format!("{error:#}").contains("warmup mass data"),
                "{error:#}"
            );
            assert!(
                amplitude.data.graph_terms[0]
                    .multi_channeling_setup
                    .sampling_catalogue
                    .as_ref()
                    .is_none()
            );
            amplitude.data.graph_terms[0].real_mass_vec = masses;
            let Externals::Constant { momenta, .. } =
                &mut runtime.get_mut_settings().kinematics.externals;
            momenta[0] = ExternalMomenta::Independent([5.0, 0.0, 0.0, 0.0].map(F));
            runtime.warm_up(&model)?;
            let ProcessIntegrand::Amplitude(amplitude) = &runtime else {
                unreachable!()
            };
            assert!(
                amplitude.data.graph_terms[0]
                    .multi_channeling_setup
                    .sampling_bridge::<f64>()
                    .is_ok(),
                "changed externals must invalidate the cached Double failure"
            );
        }
        let mut runtime_for_rescue = integrand.clone();
        let mut seam_coordinates = None;
        let ProcessIntegrand::Amplitude(integrand) = integrand else {
            unreachable!("generated an amplitude")
        };
        let term = integrand.data.graph_terms.first_mut().unwrap();
        assert_eq!(term.production_orientation_keys.len(), 18);
        // Requested parent order changes coordinates, never the generated
        // catalogue or its exact external routing. Compare the same native
        // point after swapping its two loop blocks, then exercise both compiled
        // physical channels in the common master frame.
        {
            use crate::integrands::process::SamplingMapComponent;
            let setup = &term.multi_channeling_setup;
            let generated = setup.all_bases.clone();
            let first = setup.sampling_parent_lmb(&[5, 6])?;
            let reversed = setup.sampling_parent_lmb(&[6, 5])?;
            assert_eq!(
                first.loop_edges.iter().copied().collect_vec(),
                [EdgeIndex(5), EdgeIndex(6)]
            );
            assert_eq!(first.tree, reversed.tree);
            assert_eq!(first.ext_edges, reversed.ext_edges);
            let external = [[5.0, 0.0, 0.0, 1.0]; 2];
            let first_point = setup.lmb_frame_map(&first, &external)?.forward(
                &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                &mut SamplingMapContext::detached(&[]),
            )?;
            let reversed_point = setup.lmb_frame_map(&reversed, &external)?.forward(
                &[4.0, 5.0, 6.0, 1.0, 2.0, 3.0],
                &mut SamplingMapContext::detached(&[]),
            )?;
            assert_eq!(first_point.point, reversed_point.point);
            assert_eq!(first_point.jacobian, reversed_point.jacobian);
            assert_eq!(setup.all_bases, generated);
            for invalid in [&[5][..], &[5, 5], &[5, 999]] {
                assert!(setup.sampling_parent_lmb(invalid).is_err());
            }
            for parent in [vec![5, 6], vec![6, 5]] {
                let mut parameterization =
                    settings.sampling.get_parameterization_settings().unwrap();
                parameterization.sampling_channels.default_channel_selection = vec!["C".into()];
                let definition = parameterization
                    .sampling_channels
                    .channel_definitions
                    .get_mut("massive_kite")
                    .unwrap()
                    .get_mut("C")
                    .unwrap();
                definition.parent_lmb = parent.clone();
                definition.subspace_lmb = parent;
                let bridge =
                    term.compile_sampling_bridge(&parameterization, &settings, &external, None)?;
                assert!(matches!(
                    bridge.channels()[0].map,
                    CompiledSamplingMap::Affine { .. }
                ));
                let coordinates = [0.19, 0.27, 0.61, 0.39, 0.72, 0.58];
                let mapped = bridge.forward(SamplingChannelId(0), &coordinates)?;
                let inverse = bridge
                    .inverse(SamplingChannelId(0), &mapped.raw_coordinates)?
                    .expect("full-support map must contain the supplied point");
                assert!((mapped.map.jacobian * inverse.map.inverse_jacobian - 1.0).abs() < 1.0e-9);
                for (actual, expected) in inverse.map.coordinates.iter().zip(coordinates) {
                    assert!((actual - expected).abs() < 1.0e-10);
                }
            }
        }
        // This routed factory/component gate supplies the regular point for
        // the production joint-card and original-source policy checks below.
        fn check_joint_matcher<T: FloatLike>(
            term: &AmplitudeGraphTerm,
            program: crate::integrands::process::SamplingExpressionEvaluator,
        ) -> Result<Vec<T>> {
            use crate::integrands::process::{SamplingMapComponent, SharedEnergyJointMap};
            use crate::momentum::sample::SubspaceData;
            let one = F::<T>::default().one();
            let zero = one.zero();
            let half = &one / one.from_usize(2);
            // A norm-preserving rotation about the boosted z axis keeps the
            // energy geometry while avoiding the ordinary inverse's azimuth seam.
            let transverse = [3, 4].map(|value| one.from_usize(value) / one.from_usize(5));
            let energy = one.from_usize(26).sqrt();
            let externals: ExternalFourMomenta<F<T>> = (0..2)
                .map(|_| {
                    FourMomentum::from_args(energy.clone(), zero.clone(), zero.clone(), one.clone())
                })
                .collect();
            let spatial = externals
                .iter()
                .map(|p| p.spatial.clone())
                .collect::<crate::momentum::sample::ExternalThreeMomenta<F<T>>>();
            let masses: EdgeVec<F<T>> = term
                .real_mass_vec
                .as_ref()
                .unwrap()
                .iter()
                .map(|(_, m)| m.map(F::<T>::from_ff64).unwrap_or_else(|| zero.clone()))
                .collect();
            let mut regular_raw = None;
            for (parent, active_edge) in [(vec![4, 6], 6), (vec![5, 4], 5)] {
                let lmbs = TiVec::from(vec![
                    term.multi_channeling_setup.sampling_parent_lmb(&parent)?,
                ]);
                let id = LmbIndex::from(0);
                let lmb = &lmbs[id];
                let subspace = SubspaceData::new_from_parent_basis_edges(
                    &[EdgeIndex(active_edge)],
                    &term.graph.full_filter(),
                    id,
                    &term.graph,
                    &lmbs,
                )?;
                let active = subspace.iter_lmb_indices().next().unwrap();
                let prior = lmb
                    .loop_edges
                    .iter_enumerated()
                    .find_map(|(i, e)| (*e == EdgeIndex(4)).then_some(i))
                    .unwrap();
                let select = |edges: &[usize]| {
                    term.esurfaces
                        .iter()
                        .find(|s| {
                            s.energies
                                .iter()
                                .map(|e| e.0)
                                .sorted()
                                .eq(edges.iter().copied())
                                && s.compute_shift_part_from_momenta(&externals, lmb) < zero
                        })
                        .unwrap()
                        .clone()
                };
                let c = select(&[2, 4, 6]);
                let d = select(&[3, 5, 6]);
                let (prepare, common) = c.sampling_joint_geometry_in_subspace(
                    &d,
                    &subspace,
                    &lmbs,
                    &term.graph,
                    &masses,
                    &externals,
                    &[prior],
                )?;
                let context = [transverse[0].0.clone(), transverse[1].0.clone(), (-&half).0];
                let geometry = prepare(&context)?;
                assert_eq!(
                    geometry.shifts[0],
                    [
                        transverse[0].0.clone(),
                        transverse[1].0.clone(),
                        half.0.clone()
                    ]
                );
                assert_eq!(geometry.shifts[1], context);
                for sum in &geometry.energy_sums {
                    assert!(
                        (F(sum.clone()) - (&energy - one.from_usize(3) * &half)).abs()
                            < one.epsilon().sqrt()
                    );
                }
                let mut loops = LoopMomenta::from_iter(
                    (0..2).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                );
                loops[prior] = ThreeMomentum::new(
                    F(context[0].clone()),
                    F(context[1].clone()),
                    F(context[2].clone()),
                );
                let offset: ThreeMomentum<F<T>> = common.compute_momentum(&loops, &spatial);
                assert_eq!(
                    offset.norm_squared() > zero,
                    active_edge == 5,
                    "alternate parent must exercise a nonzero shared-energy affine offset"
                );
                let kernel = SharedEnergyJointMap::new(
                    prepare.clone(),
                    3,
                    (&one / one.from_usize(8)).0,
                    one.0.clone(),
                    5.0,
                    program.clone(),
                )?;
                let cube = [19, 31, 67].map(|v| (one.from_usize(v) / one.from_usize(100)).0);
                let forward = kernel.forward(&cube, &mut SamplingMapContext::detached(&context))?;
                if active_edge == 6 {
                    regular_raw = Some(
                        context
                            .iter()
                            .cloned()
                            .chain(forward.point.iter().cloned())
                            .collect(),
                    );
                }
                assert!(
                    forward
                        .diagnostics
                        .iter()
                        .any(|d| d.contains("certified normal disk"))
                );
                let inverse = kernel
                    .inverse(&forward.point, &mut SamplingMapContext::detached(&context))?
                    .expect("selected compact point");
                assert!(
                    (F(forward.jacobian) * F(inverse.inverse_jacobian) - &one).abs()
                        < one.epsilon().sqrt() * one.from_usize(100)
                );
                for point in [
                    forward.point,
                    vec![half.0.clone(), (-&one).0, one.0.clone()],
                ] {
                    let x = ThreeMomentum::new(
                        F(point[0].clone()),
                        F(point[1].clone()),
                        F(point[2].clone()),
                    );
                    loops[active] = &x - &offset;
                    let e0 = (x.norm_squared() + F(geometry.masses[0].clone()).square()).sqrt();
                    for (i, surface) in [&c, &d].into_iter().enumerate() {
                        let shift = ThreeMomentum::new(
                            F(geometry.shifts[i][0].clone()),
                            F(geometry.shifts[i][1].clone()),
                            F(geometry.shifts[i][2].clone()),
                        );
                        let partner = &x + shift;
                        let expected = &e0
                            + (partner.norm_squared() + F(geometry.masses[i + 1].clone()).square())
                                .sqrt()
                            - F(geometry.energy_sums[i].clone());
                        let original =
                            surface.compute_from_momenta(lmb, &masses, &loops, &externals);
                        assert!(
                            (original - expected).abs()
                                < one.epsilon().sqrt() * one.from_usize(100)
                        );
                    }
                }
                assert!(
                    c.sampling_joint_geometry_in_subspace(
                        &d,
                        &subspace,
                        &lmbs,
                        &term.graph,
                        &masses,
                        &externals,
                        &[]
                    )
                    .err()
                    .unwrap()
                    .to_string()
                    .contains("unsampled")
                );
                assert!(
                    c.sampling_joint_geometry_in_subspace(
                        &c,
                        &subspace,
                        &lmbs,
                        &term.graph,
                        &masses,
                        &externals,
                        &[prior]
                    )
                    .err()
                    .unwrap()
                    .to_string()
                    .contains("common full signed routing")
                );
                let mut repeated = c.clone();
                repeated.energies.push(EdgeIndex(2));
                assert!(
                    repeated
                        .sampling_joint_geometry_in_subspace(
                            &d,
                            &subspace,
                            &lmbs,
                            &term.graph,
                            &masses,
                            &externals,
                            &[prior]
                        )
                        .err()
                        .unwrap()
                        .to_string()
                        .contains("two varying energy occurrences")
                );
                // Diagnostic mass input isolates a fixed-energy range failure;
                // it does not claim this is a regenerated physical kite model.
                let mut tiny_masses = masses.clone();
                tiny_masses[EdgeIndex(4)] = &one / one.from_usize(10).powi(200);
                let mut tiny_externals = externals.clone();
                for momentum in &mut tiny_externals {
                    momentum.temporal.value = &tiny_masses[EdgeIndex(4)] * one.from_usize(2);
                }
                let (tiny, _) = c.sampling_joint_geometry_in_subspace(
                    &d,
                    &subspace,
                    &lmbs,
                    &term.graph,
                    &tiny_masses,
                    &tiny_externals,
                    &[prior],
                )?;
                let tiny_result = tiny(&[zero.0.clone(), zero.0.clone(), zero.0.clone()]);
                if matches!(
                    T::sampling_precision(),
                    crate::utils::SamplingPrecision::Fixed256
                        | crate::utils::SamplingPrecision::Arb
                ) {
                    let sum = F(tiny_result?.energy_sums[0].clone());
                    assert_eq!(sum, tiny_masses[EdgeIndex(4)]);
                    assert!(sum > zero);
                } else {
                    assert!(tiny_result.err().unwrap().downcast_ref::<crate::integrands::process::sampling_maps::SamplingEvaluationError>().is_some());
                }
                let short = ExternalFourMomenta::from_iter(externals.iter().take(1).cloned());
                assert!(
                    c.sampling_joint_geometry_in_subspace(
                        &d,
                        &subspace,
                        &lmbs,
                        &term.graph,
                        &masses,
                        &short,
                        &[prior]
                    )
                    .err()
                    .unwrap()
                    .to_string()
                    .contains("external ports")
                );
            }
            Ok(regular_raw.unwrap())
        }
        let joint_program =
            crate::integrands::process::SharedEnergyJointMap::<f64>::compile_program()?;
        let regular_raw = check_joint_matcher::<f64>(term, joint_program.clone())?;
        check_joint_matcher::<crate::utils::QuadFloat>(term, joint_program.clone())?;
        check_joint_matcher::<crate::utils::ArbPrec>(term, joint_program)?;
        {
            use crate::{
                integrands::process::GaussianReferenceFunction,
                momentum::ExternalMomenta,
                settings::runtime::{
                    SamplingSettingsParser, StabilityLevelSetting, kinematic::Externals,
                },
            };
            let mut joint_runtime = runtime_for_rescue.clone();
            let mut parser: SamplingSettingsParser =
                toml::from_str(&toml::to_string(&settings.sampling)?)?;
            parser.default_channel_selection = vec!["joint".into(), "ordinary".into()];
            let definitions = parser.channel_definitions.get_mut("massive_kite").unwrap();
            let mut definition = definitions["C"].clone();
            definition.around =
                "then(complement(4),block(lmb(6),intersect(surface(2,4,6),surface(3,5,6))))".into();
            definition.subspace_lmb = vec![6];
            definitions.insert("joint".into(), definition);
            let native_settings = joint_runtime.get_mut_settings();
            native_settings.sampling = toml::from_str(&toml::to_string(&parser)?)?;
            let Externals::Constant { momenta, .. } = &mut native_settings.kinematics.externals;
            momenta[0] = ExternalMomenta::Independent([26_f64.sqrt(), 0.0, 0.0, 1.0].map(F));
            native_settings.stability.rotation_axis.clear();
            joint_runtime.warm_up(&model)?;
            assert_eq!(joint_runtime.group_orientation_count(GroupId(0)), Some(18));
            let ProcessIntegrand::Amplitude(amplitude) = &joint_runtime else {
                unreachable!()
            };
            assert_eq!(
                amplitude.data.graph_terms[0]
                    .production_orientation_keys
                    .len(),
                18
            );
            let bridge = amplitude.data.graph_terms[0]
                .sampling_setup()
                .sampling_bridge::<f64>()?;
            assert!(bridge.channels()[0].map.contract().requires_proposal_policy);
            assert_eq!(
                bridge.channels()[0].map.contract().support,
                SamplingSupport::Restricted
            );
            // The factory point comes from a certified regular complement. Use
            // the complete compiled inverse to recover the original six cube coordinates.
            let inverse = bridge
                .inverse(SamplingChannelId(0), &regular_raw)?
                .expect("regular joint point");
            let cube = inverse.map.coordinates;
            let mapped = bridge.forward(SamplingChannelId(0), &cube)?;
            assert!(
                mapped
                    .map
                    .diagnostics
                    .iter()
                    .any(|d| d.contains("certified normal disk"))
            );
            assert!((mapped.map.jacobian * inverse.map.inverse_jacobian - 1.0).abs() < 1e-9);
            assert!(
                bridge
                    .inverse(SamplingChannelId(1), &mapped.raw_coordinates)?
                    .is_some()
            );
            let reference =
                GaussianReferenceFunction::new(1.2, vec![0.4, -0.3, 0.2, -0.2, 0.1, 0.35])?;
            let point = Sample::Continuous(F(1.0), cube.iter().copied().map(F).collect());
            let standard = joint_runtime.evaluate_reference_sample_detailed(&point, &reference)?;
            assert!(!standard.evaluation.evaluation_metadata.is_nan);
            assert!(standard.evaluation.integrand_result.re.0 > 0.0);
            for level in [
                StabilityLevelSetting::default_quad(),
                StabilityLevelSetting::default_arb(),
            ] {
                let mut native = joint_runtime.clone();
                native.get_mut_settings().stability.levels = vec![level];
                native.warm_up(&model)?;
                let value = native.evaluate_reference_sample_detailed(&point, &reference)?;
                assert!(!value.evaluation.evaluation_metadata.is_nan);
                assert!(
                    (value.evaluation.integrand_result.re.0
                        / standard.evaluation.integrand_result.re.0
                        - 1.0)
                        .abs()
                        < 1e-8
                );
                assert!(
                    (value.moments.second_moment.0 / standard.moments.second_moment.0 - 1.0).abs()
                        < 1e-8
                );
            }
            // The actual physical body is distinct from the reference and
            // controlled replay targets. Retain all 18 orientations and CT terms.
            {
                use crate::{
                    integrands::process::MomentumSpaceEvaluationInput, settings::runtime::SumMode,
                };
                let mut physical = joint_runtime.clone();
                let mut selected_settings = parser.clone();
                selected_settings.graphs = SumMode::MonteCarlo;
                selected_settings.sampling_channels = SumMode::MonteCarlo;
                physical.get_mut_settings().sampling =
                    toml::from_str(&toml::to_string(&selected_settings)?)?;
                physical.get_mut_settings().stability.levels =
                    vec![StabilityLevelSetting::default_arb()];
                physical.warm_up(&model)?;
                assert_eq!(physical.group_orientation_count(GroupId(0)), Some(18));
                let selected_source = Sample::Discrete(
                    F(1.0),
                    0,
                    Some(Box::new(Sample::Discrete(
                        F(1.0),
                        0,
                        Some(Box::new(point.clone())),
                    ))),
                );
                let selected_x = physical
                    .evaluate_sample_precise(
                        &selected_source,
                        &model,
                        F(1.0),
                        false,
                        Complex::new_zero(),
                    )?
                    .try_into_f64()?;
                let mut raw_input = MomentumSpaceEvaluationInput {
                    loop_momenta: mapped
                        .raw_coordinates
                        .chunks_exact(3)
                        .map(|p| ThreeMomentum::new(F(p[0]), F(p[1]), F(p[2])))
                        .collect(),
                    integrator_weight: F(1.0),
                    graph_id: Some(0),
                    group_id: None,
                    orientation: None,
                    channel_id: None,
                };
                let raw = physical
                    .evaluate_momentum_configuration_precise(&model, &raw_input, false)?
                    .try_into_f64()?;
                raw_input.graph_id = None;
                raw_input.group_id = Some(GroupId(0));
                raw_input.channel_id = Some(SamplingChannelId(0));
                let selected_raw = physical
                    .evaluate_momentum_configuration_precise(&model, &raw_input, false)?
                    .try_into_f64()?;
                let norm = raw.integrand_result.re.0.hypot(raw.integrand_result.im.0);
                assert!(norm.is_finite() && norm > 0.0);
                assert!(!raw.evaluation_metadata.is_nan);
                // Raw momentum evaluation has no separate parameterization;
                // X-space reports unity after folding its map into the value.
                for (result, factor, parameterization) in [
                    (&selected_x, mapped.selected_factor()?, Some(F(1.0))),
                    (&selected_raw, mapped.partition.weight(0).unwrap(), None),
                ] {
                    assert!(!result.evaluation_metadata.is_nan);
                    assert_eq!(result.parameterization_jacobian, parameterization);
                    assert!(
                        !result
                            .evaluation_metadata
                            .sampling_proposal_policies
                            .is_empty()
                    );
                    for (actual, original) in [
                        (result.integrand_result.re.0, raw.integrand_result.re.0),
                        (result.integrand_result.im.0, raw.integrand_result.im.0),
                    ] {
                        assert!(
                            (actual - factor * original).abs() < 1e-7 * factor.abs() * norm,
                            "physical joint accounting: actual {actual}, factor {factor}, raw {original}"
                        );
                    }
                }
            }
            // Halton dispersion is diagnostic rather than a randomized error
            // bar. This summed-channel acceptance tests the complete estimator.
            let coordinates = (1..=8192)
                .map(|sample| {
                    [2, 3, 5, 7, 11, 13]
                        .map(|base| {
                            let (mut index, mut fraction, mut value) = (sample, 1.0, 0.0);
                            while index > 0 {
                                fraction /= base as f64;
                                value += (index % base) as f64 * fraction;
                                index /= base;
                            }
                            value
                        })
                        .to_vec()
                })
                .collect::<Vec<_>>();
            super::super::tests::check_joint_policy_source_replay(
                &joint_runtime,
                &model,
                SamplingChannelId(0),
                &cube,
            )?;
            let started = std::time::Instant::now();
            let probe =
                joint_runtime.evaluate_reference_coordinates(&coordinates[..8], &reference)?;
            assert_eq!(probe.finite_sample_count, 8);
            info!(
                stage = "joint_reference_cost_probe",
                sample_count = 8,
                elapsed_seconds = started.elapsed().as_secs_f64(),
                "all-orientation joint reference preparation and native evaluation cost"
            );
            let report = joint_runtime.evaluate_reference_coordinates(&coordinates, &reference)?;
            assert_eq!(report.finite_sample_count, coordinates.len());
            assert!((report.normalization - 1.0).abs() < 0.06, "{report:?}");
            assert!(
                (report.second_moment / report.expected_second_moment - 1.0).abs() < 0.08,
                "{report:?}"
            );
        }

        // The same generated C equation also defines a genuine three-dimensional
        // fiber: k=q4 varies while l=q6 was sampled by the preceding block.
        // This factory/composition gate complements the production command-card
        // binding below; automatic channel discovery remains separate work.
        fn check_fiber<T: FloatLike>(
            term: &AmplitudeGraphTerm,
            settings: &RuntimeSettings,
        ) -> Result<()> {
            use crate::graph::lmb::LMBwithEdges;
            use crate::integrands::process::sampling_maps::SamplingEvaluationError;
            use crate::integrands::process::{
                PreparedSurfaceStatus, SamplingMapComponent, SamplingMapComposition,
                SamplingMapEmbedding, SurfaceRadialMap,
            };
            use crate::momentum::sample::SubspaceData;
            let one = F::<T>::from_f64(1.0);
            let zero = one.zero();
            let lmbs = TiVec::from(vec![term.graph.loop_momentum_basis.clone()]);
            let subspace = SubspaceData::new_from_parent_basis_edges(
                &[EdgeIndex(4)],
                &term.graph.full_filter(),
                LmbIndex::from(0),
                &term.graph,
                &lmbs,
            )?;
            let k = subspace.iter_lmb_indices().next().unwrap();
            let l = lmbs[LmbIndex::from(0)]
                .loop_edges
                .iter_enumerated()
                .find_map(|(index, edge)| (*edge == EdgeIndex(6)).then_some(index))
                .unwrap();
            assert_eq!(lmbs[LmbIndex::from(0)].loop_edges[k], EdgeIndex(4));
            assert_eq!(lmbs[LmbIndex::from(0)].ext_edges.len(), 2);
            let masses: EdgeVec<F<T>> = term
                .real_mass_vec
                .as_ref()
                .unwrap()
                .iter()
                .map(|(_, mass)| mass.map(F::<T>::from_ff64).unwrap_or_else(|| zero.clone()))
                .collect();
            let externals: ExternalFourMomenta<F<T>> = (0..2)
                .map(|_| {
                    FourMomentum::from_args(
                        one.from_i64(5),
                        zero.clone(),
                        zero.clone(),
                        zero.clone(),
                    )
                })
                .collect();
            let surface = term
                .esurfaces
                .iter()
                .find(|surface| {
                    surface.energies == [EdgeIndex(2), EdgeIndex(4), EdgeIndex(6)]
                        && surface
                            .compute_shift_part_from_momenta(&externals, &lmbs[LmbIndex::from(0)])
                            < zero
                })
                .unwrap();
            let fiber = surface.sampling_radial_map_in_subspace(
                &subspace,
                &lmbs,
                &term.graph,
                &masses,
                &externals,
                &[l],
                settings,
                2.0,
                1.0,
            )?;
            let cube = [19, 31, 67].map(|value| (one.from_i64(value) / one.from_i64(100)).0);
            // A=[3,4] is independent of q6 in this parent, whereas C=[2,4,6]
            // depends on it. Only the former may omit that later coordinate.
            let independent = term
                .esurfaces
                .iter()
                .find(|candidate| {
                    candidate.energies == [EdgeIndex(3), EdgeIndex(4)]
                        && candidate
                            .compute_shift_part_from_momenta(&externals, &lmbs[LmbIndex::from(0)])
                            < zero
                })
                .unwrap()
                .sampling_radial_map_in_subspace(
                    &subspace,
                    &lmbs,
                    &term.graph,
                    &masses,
                    &externals,
                    &[],
                    settings,
                    2.0,
                    1.0,
                )?;
            assert_eq!(independent.contract().support, SamplingSupport::Full);
            let independent_point = independent.forward(&cube)?;
            assert!(
                independent
                    .inverse(&independent_point.point)?
                    .inverse_jacobian
                    .is_finite()
            );
            assert!(
                surface
                    .sampling_radial_map_in_subspace(
                        &subspace,
                        &lmbs,
                        &term.graph,
                        &masses,
                        &externals,
                        &[],
                        settings,
                        2.0,
                        1.0,
                    )
                    .unwrap_err()
                    .to_string()
                    .contains("depends on unsampled parent edge")
            );
            for (numerator, denominator, existing) in
                [(0, 1, true), (1, 1, true), (21, 10, false), (3, 1, false)]
            {
                let l = one.from_i64(numerator) / one.from_i64(denominator);
                let context = [l.0.clone(), zero.0.clone(), zero.0.clone()];
                let (center, status) = fiber.prepare_context(&context)?;
                let minimum = (l.square() + one.from_i64(4)).sqrt() + (l.square() + &one).sqrt()
                    - one.from_i64(5);
                assert_eq!(status.as_ref().unwrap().is_existing(), existing);
                assert_eq!(F(center[0].clone()), -&l / one.from_i64(2));
                if let Some(PreparedSurfaceStatus::Existing { normalized_margin }) = status {
                    assert!(
                        (F(normalized_margin.unwrap()) + minimum / one.from_i64(5)).abs()
                            < one.epsilon() * one.from_i64(100)
                    );
                }
                let mapped = fiber.forward_with_context(&cube, &context)?;
                let inverse = fiber.inverse_with_context(&mapped.point, &context)?;
                let tolerance = one.epsilon().sqrt() * one.from_i64(100);
                assert!(
                    (F(mapped.jacobian.clone()) * F(inverse.inverse_jacobian) - &one).abs()
                        < tolerance
                );
                for (actual, expected) in inverse.coordinates.iter().zip(&cube) {
                    assert!((F(actual.clone()) - F(expected.clone())).abs() < tolerance);
                }
                assert_eq!(
                    mapped
                        .diagnostics
                        .iter()
                        .any(|label| label.contains("regular_root")),
                    existing
                );
            }
            // The algebraic boundary l²=96/25 is irrational in the stored
            // coordinates. Rounded equality must request precision, while
            // both separated neighboring fibers are classified correctly.
            let boundary = (one.from_i64(96) / one.from_i64(25)).sqrt();
            let context = [boundary.0.clone(), zero.0.clone(), zero.0.clone()];
            assert!(
                fiber
                    .prepare_context(&context)
                    .unwrap_err()
                    .downcast_ref::<SamplingEvaluationError>()
                    .is_some()
            );
            let offset = one.epsilon().sqrt();
            for (l, existing) in [
                (&boundary * (&one - &offset), true),
                (&boundary * (&one + &offset), false),
            ] {
                assert_eq!(
                    fiber
                        .prepare_context(&[l.0, zero.0.clone(), zero.0.clone()])?
                        .1
                        .unwrap()
                        .is_existing(),
                    existing
                );
            }
            assert!(fiber.prepare_context(&[]).is_err());
            // Diagnostic mass overrides exercise the same routed geometry, not
            // a differently generated physical amplitude. An asymmetric zero
            // mass puts one energy exactly at its minimizing endpoint.
            for (first_mass, second_mass) in [(0, 2), (0, 0)] {
                let mut diagnostic_masses = masses.clone();
                diagnostic_masses[EdgeIndex(2)] = one.from_i64(first_mass);
                diagnostic_masses[EdgeIndex(4)] = one.from_i64(second_mass);
                let massless = surface.sampling_radial_map_in_subspace(
                    &subspace,
                    &lmbs,
                    &term.graph,
                    &diagnostic_masses,
                    &externals,
                    &[l],
                    settings,
                    2.0,
                    1.0,
                )?;
                let context = [one.0.clone(), zero.0.clone(), zero.0.clone()];
                let (center, status) = massless.prepare_context(&context)?;
                assert!(status.unwrap().is_existing());
                let mut loops = LoopMomenta::from_iter(
                    (0..2).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                );
                loops[k] = ThreeMomentum::new(
                    F(center[0].clone()),
                    F(center[1].clone()),
                    F(center[2].clone()),
                );
                loops[l] = ThreeMomentum::new(one.clone(), zero.clone(), zero.clone());
                let velocity = LoopMomenta::from_iter(
                    (0..2).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                );
                let value = surface
                    .compute_self_and_r_derivative(
                        &zero,
                        &velocity,
                        &loops,
                        &externals,
                        &diagnostic_masses,
                        &lmbs[LmbIndex::from(0)],
                    )
                    .0;
                assert!(value < zero);
                for sign in [-1, 1] {
                    let mut moved = loops.clone();
                    moved[k].py += one.from_i64(sign) / one.from_i64(10);
                    let shifted = surface
                        .compute_self_and_r_derivative(
                            &zero,
                            &velocity,
                            &moved,
                            &externals,
                            &diagnostic_masses,
                            &lmbs[LmbIndex::from(0)],
                        )
                        .0;
                    assert!(
                        shifted > value,
                        "transverse perturbation must increase the minimum"
                    );
                }
                let forward = massless.forward_with_context(&cube, &context)?;
                let inverse = massless.inverse_with_context(&forward.point, &context)?;
                assert!(
                    (F(forward.jacobian) * F(inverse.inverse_jacobian) - &one).abs()
                        < one.epsilon().sqrt() * one.from_i64(100)
                );
            }
            // The alternate parent p=q5=k+l makes both A energies carry a
            // large common complement shift. Double loses their unit spatial
            // separation; higher precision must recover the absent fiber.
            let alternate = TiVec::from(vec![
                term.lmb_with_loop_edges([EdgeIndex(5), EdgeIndex(6)].as_slice())?,
            ]);
            assert_eq!(alternate[LmbIndex::from(0)].ext_edges.len(), 2);
            let active_p = SubspaceData::new_from_parent_basis_edges(
                &[EdgeIndex(5)],
                &term.graph.full_filter(),
                LmbIndex::from(0),
                &term.graph,
                &alternate,
            )?;
            let complement = alternate[LmbIndex::from(0)]
                .loop_edges
                .iter_enumerated()
                .filter_map(|(index, edge)| (*edge == EdgeIndex(6)).then_some(index))
                .collect_vec();
            let external: ExternalFourMomenta<F<T>> = (0..2)
                .map(|_| {
                    FourMomentum::from_args(
                        one.from_i64(21) / one.from_i64(10),
                        one.clone(),
                        zero.clone(),
                        zero.clone(),
                    )
                })
                .collect();
            let a = term
                .esurfaces
                .iter()
                .find(|surface| {
                    surface.energies == [EdgeIndex(3), EdgeIndex(4)]
                        && surface.compute_shift_part_from_momenta(
                            &external,
                            &alternate[LmbIndex::from(0)],
                        ) < zero
                })
                .unwrap();
            let cancellation = a.sampling_radial_map_in_subspace(
                &active_p,
                &alternate,
                &term.graph,
                &masses,
                &external,
                &complement,
                settings,
                2.0,
                1.0,
            )?;
            let large = one.from_i64(2).powi(55);
            let prepared =
                cancellation.prepare_context(&[large.0.clone(), zero.0.clone(), zero.0.clone()]);
            if one.epsilon() * &large > one.clone() / one.from_i64(100) {
                assert!(
                    prepared
                        .unwrap_err()
                        .downcast_ref::<SamplingEvaluationError>()
                        .is_some()
                );
            } else {
                assert!(matches!(
                    prepared?.1,
                    Some(PreparedSurfaceStatus::Absent { .. })
                ));
            }
            // Three varying energies exercise the generic routed/SOCP branch.
            // A large boost makes the origin exterior although the invariant
            // energy still exceeds the three-particle threshold.
            let full = SubspaceData::new_from_parent_basis_edges(
                &[EdgeIndex(6), EdgeIndex(4)],
                &term.graph.full_filter(),
                LmbIndex::from(0),
                &term.graph,
                &lmbs,
            )?;
            assert_eq!(
                full.iter_basis_edges(&lmbs).collect_vec(),
                [EdgeIndex(4), EdgeIndex(6)]
            );
            let origin = LoopMomenta::from_iter(
                (0..2).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            for (energy, boost, existing) in [
                (one.from_i64(89).sqrt(), one.from_i64(8), true),
                (one.from_i64(2), zero.clone(), false),
            ] {
                let external: ExternalFourMomenta<F<T>> = (0..2)
                    .map(|_| {
                        FourMomentum::from_args(
                            energy.clone(),
                            zero.clone(),
                            zero.clone(),
                            boost.clone(),
                        )
                    })
                    .collect();
                let general = surface.sampling_radial_map_in_subspace(
                    &full,
                    &lmbs,
                    &term.graph,
                    &masses,
                    &external,
                    &[],
                    settings,
                    2.0,
                    1.0,
                )?;
                let (center, status) = general.prepare_context(&[])?;
                assert_eq!(status.unwrap().is_existing(), existing);
                let origin_value = surface
                    .compute_self_and_r_derivative(
                        &zero,
                        &origin,
                        &origin,
                        &external,
                        &masses,
                        &lmbs[LmbIndex::from(0)],
                    )
                    .0;
                assert!(origin_value > zero);
                if existing {
                    let center = LoopMomenta::from_iter(center.chunks_exact(3).map(|p| {
                        ThreeMomentum::new(F(p[0].clone()), F(p[1].clone()), F(p[2].clone()))
                    }));
                    let value = surface
                        .compute_self_and_r_derivative(
                            &zero,
                            &origin,
                            &center,
                            &external,
                            &masses,
                            &lmbs[LmbIndex::from(0)],
                        )
                        .0;
                    assert!(value < zero, "native SOCP center value={value}");
                }
            }
            // The ordered determinant includes complement-dependent centers
            // and roots. Its off-diagonal entries need no extra Jacobian.
            let ordinary = SurfaceRadialMap::new(3, vec![zero.0.clone(); 3], None, 2.0, 1.0)?;
            let map = SamplingMapEmbedding::from_composition(
                SamplingMapComposition::then(vec![Box::new(ordinary), Box::new(fiber)])?,
                [l, k]
                    .into_iter()
                    .flat_map(|index| [3 * index.0, 3 * index.0 + 1, 3 * index.0 + 2])
                    .collect(),
            )?;
            let coordinates =
                [13, 23, 47, 19, 31, 67].map(|value| (one.from_i64(value) / one.from_i64(100)).0);
            let forward = map.forward(&coordinates, &mut SamplingMapContext::detached(&[]))?;
            let inverse = map
                .inverse(&forward.point, &mut SamplingMapContext::detached(&[]))?
                .expect("full-support map must contain the supplied point");
            assert!(
                (F(forward.jacobian.clone()) * F(inverse.inverse_jacobian) - &one).abs()
                    < one.epsilon().sqrt() * one.from_i64(100)
            );
            let step = one / zero.from_i64(1000000);
            let mut matrix = vec![vec![0.0; 6]; 6];
            for column in 0..6 {
                let mut plus = coordinates.clone();
                let mut minus = coordinates.clone();
                plus[column] = (F(plus[column].clone()) + &step).0;
                minus[column] = (F(minus[column].clone()) - &step).0;
                let plus = map.forward(&plus, &mut SamplingMapContext::detached(&[]))?;
                let minus = map.forward(&minus, &mut SamplingMapContext::detached(&[]))?;
                for (row, derivatives) in matrix.iter_mut().enumerate() {
                    derivatives[column] = ((F(plus.point[row].clone())
                        - F(minus.point[row].clone()))
                        / (&step * step.from_i64(2)))
                    .into_f64();
                }
            }
            assert!(
                (SamplingMapAffine::new(matrix, vec![0.0; 6])?.determinant()
                    / F(forward.jacobian).into_f64()
                    - 1.0)
                    .abs()
                    < 2.0e-6
            );
            Ok(())
        }
        check_fiber::<f64>(term, &settings)?;
        check_fiber::<crate::utils::QuadFloat>(term, &settings)?;
        check_fiber::<crate::utils::ArbPrec>(term, &settings)?;
        // Compile the user's shortcut and explicit ordering through the real
        // warmed graph, retaining all 18 orientations and the ordinary channel.
        {
            use crate::integrands::process::GaussianReferenceFunction;
            use crate::settings::runtime::{SamplingSettingsParser, SumMode};
            let reference =
                GaussianReferenceFunction::new(1.2, vec![0.4, -0.3, 0.2, -0.2, 0.1, 0.35])?;
            let cube = vec![0.19, 0.27, 0.61, 0.39, 0.72, 0.58];
            let point = Sample::Continuous(F(1.0), cube.iter().copied().map(F).collect());
            let mut previous = None;
            for around in ["surface(2,4,6)", "then(complement(6),surface(2,4,6))"] {
                let mut runtime = runtime_for_rescue.clone();
                let mut parser: SamplingSettingsParser =
                    toml::from_str(&toml::to_string(&settings.sampling)?)?;
                parser.power = 1.0;
                parser.default_channel_selection = vec!["C".to_owned(), "ordinary".to_owned()];
                let definition = parser
                    .channel_definitions
                    .get_mut("massive_kite")
                    .unwrap()
                    .get_mut("C")
                    .unwrap();
                definition.around = around.to_owned();
                definition.subspace_lmb = vec![4];
                runtime.get_mut_settings().sampling = toml::from_str(&toml::to_string(&parser)?)?;
                runtime.warm_up(&model)?;
                assert_eq!(runtime.group_orientation_count(GroupId(0)), Some(18));
                let summed = runtime.evaluate_reference_sample_detailed(&point, &reference)?;
                if let Some(previous) = previous {
                    assert_eq!(summed.evaluation.integrand_result, previous);
                }
                previous = Some(summed.evaluation.integrand_result);
                if around.starts_with("surface") {
                    // Reuse the saved-state acceptance quadrature and physical
                    // raw-momentum moment, including both existing/absent fibers.
                    let coordinates = (1..=8192)
                        .map(|sample| {
                            [2, 3, 5, 7, 11, 13]
                                .map(|base| {
                                    let (mut index, mut fraction, mut value) = (sample, 1.0, 0.0);
                                    while index > 0 {
                                        fraction /= base as f64;
                                        value += (index % base) as f64 * fraction;
                                        index /= base;
                                    }
                                    value
                                })
                                .to_vec()
                        })
                        .collect::<Vec<_>>();
                    let report =
                        runtime.evaluate_reference_coordinates(&coordinates, &reference)?;
                    assert_eq!(report.finite_sample_count, coordinates.len());
                    assert!((report.normalization - 1.0).abs() < 0.02, "{report:?}");
                    assert!(
                        (report.second_moment - report.expected_second_moment).abs() < 0.2,
                        "{report:?}"
                    );
                }
                parser.graphs = SumMode::MonteCarlo;
                parser.sampling_channels = SumMode::MonteCarlo;
                runtime.get_mut_settings().sampling = toml::from_str(&toml::to_string(&parser)?)?;
                runtime.warm_up(&model)?;
                let mut value = 0.0;
                let mut moment = 0.0;
                for id in 0..2 {
                    let report = runtime.evaluate_reference_discrete_coordinates(
                        &[0, id],
                        std::slice::from_ref(&cube),
                        &reference,
                    )?;
                    value += report.normalization;
                    moment += report.second_moment;
                }
                assert!((value / summed.evaluation.integrand_result.re.0 - 1.0).abs() < 1.0e-12);
                assert!((moment / summed.moments.second_moment.0 - 1.0).abs() < 1.0e-12);
                let definition = parser
                    .channel_definitions
                    .get_mut("massive_kite")
                    .unwrap()
                    .get_mut("C")
                    .unwrap();
                definition.around = "product(surface(2,4,6),complement(6))".to_owned();
                runtime.get_mut_settings().sampling = toml::from_str(&toml::to_string(&parser)?)?;
                let error = runtime.warm_up(&model).unwrap_err();
                assert!(
                    format!("{error:#}").contains("depends on unsampled parent edges [6]"),
                    "{error:#}"
                );
            }
        }
        let parameterization = settings.sampling.get_parameterization_settings().unwrap();
        // Amplitudes have no auxiliary LU variable even if supplied an h
        // configuration; the same profile spelling must reject here explicitly.
        for (around, diagnostic) in [
            (
                "surface(2,4,6)",
                "radial_profile=lu_h requires exactly one phase_space(cut(...)) block",
            ),
            (
                "phase_space(cut(2,4,6))",
                "amplitudes have no auxiliary Cutkosky-cut LU scale",
            ),
        ] {
            let mut invalid_lu = parameterization.clone();
            let definition = invalid_lu
                .sampling_channels
                .channel_definitions
                .get_mut("massive_kite")
                .unwrap()
                .values_mut()
                .next()
                .unwrap();
            definition.around = around.to_owned();
            definition.radial_profile = Some(Default::default());
            // This rejection compiles the eager profile before binding its
            // physical host, so use the same development-build stack allowance
            // as the isolated LU-profile compilation tests.
            let error = std::thread::scope(|scope| {
                let term = &mut *term;
                let runtime_settings = &settings;
                std::thread::Builder::new()
                    .stack_size(64 * 1024 * 1024)
                    .spawn_scoped(scope, move || {
                        term.compile_sampling_bridge(
                            &invalid_lu,
                            runtime_settings,
                            &[[5.0, 0.0, 0.0, 0.0]; 2],
                            None,
                        )
                        .unwrap_err()
                        .to_string()
                    })
                    .unwrap()
                    .join()
                    .unwrap()
            });
            assert!(error.contains(diagnostic), "{error}");
        }
        let coordinates = [0.19, 0.27, 0.61, 0.39, 0.72, 0.58];
        let mut rest_point = Vec::new();

        for external in [[5.0, 0.0, 0.0, 0.0], [26.0_f64.sqrt(), 0.0, 0.0, 1.0]] {
            let externals = [external; 2];
            let bridge =
                term.compile_sampling_bridge(&parameterization, &settings, &externals, None)?;
            let cached = term
                .multi_channeling_setup
                .sampling_bridge::<f64>()?
                .forward(SamplingChannelId(0), &coordinates)?;
            let fresh = bridge.forward(SamplingChannelId(0), &coordinates)?;
            if external[3] == 0.0 {
                assert_eq!(cached.raw_coordinates, fresh.raw_coordinates);
                assert_eq!(cached.map.jacobian, fresh.map.jacobian);
                // Physical active cube coordinates follow canonical parent
                // order. Reordering metadata identifies the same subspace;
                // it does not request a permutation of input cube axes.
                let mut reordered = parameterization.clone();
                reordered
                    .sampling_channels
                    .channel_definitions
                    .get_mut("massive_kite")
                    .unwrap()
                    .get_mut("C")
                    .unwrap()
                    .subspace_lmb = vec![6, 4];
                let reordered =
                    term.compile_sampling_bridge(&reordered, &settings, &externals, None)?;
                let point = reordered.forward(SamplingChannelId(0), &coordinates)?;
                assert_eq!(point.raw_coordinates, fresh.raw_coordinates);
                assert_eq!(point.map.jacobian, fresh.map.jacobian);
                let inverse = reordered
                    .inverse(SamplingChannelId(0), &point.raw_coordinates)?
                    .expect("full-support map must contain the supplied point");
                for (actual, expected) in inverse.map.coordinates.iter().zip(coordinates) {
                    assert!((actual - expected).abs() < 1.0e-10);
                }
            } else {
                assert_ne!(cached.raw_coordinates, fresh.raw_coordinates);
                assert_eq!(cached.raw_coordinates, rest_point);
            }

            assert_eq!(
                bridge
                    .channels()
                    .iter()
                    .map(|channel| channel.name.as_str())
                    .collect_vec(),
                ["C", "D", "ordinary"]
            );
            assert!(matches!(
                bridge.channels()[0].map,
                CompiledSamplingMap::ImplicitSurface(_)
            ));
            assert!(matches!(
                bridge.channels()[1].map,
                CompiledSamplingMap::ImplicitSurface(_)
            ));
            let forward = bridge.forward(SamplingChannelId(0), &coordinates)?;
            if external[3] == 0.0 {
                rest_point = forward.raw_coordinates.clone();
            } else {
                assert_ne!(
                    rest_point, forward.raw_coordinates,
                    "external spatial updates must reach the captured graph equation"
                );
            }
            assert!(forward.map.jacobian.is_finite() && forward.map.jacobian > 0.0);
            assert!((forward.partition.weights.iter().sum::<f64>() - 1.0).abs() < 1.0e-13);
            for channel_id in 0..3 {
                let inverse = bridge
                    .inverse(SamplingChannelId(channel_id), &forward.raw_coordinates)?
                    .expect("full-support map must contain the supplied point");
                assert_eq!(inverse.partition, forward.partition);
                let recovered =
                    bridge.forward(SamplingChannelId(channel_id), &inverse.map.coordinates)?;
                for (original, recovered) in forward
                    .raw_coordinates
                    .iter()
                    .zip(&recovered.raw_coordinates)
                {
                    assert!((original - recovered).abs() < 1.0e-9);
                }
            }

            // Independently evaluate C=E(k+l+Q)+E(k)+E(l)-Q0 in the
            // actual parent [q4=k,q6=l]. Its radius varies with direction.
            let norm = forward
                .raw_coordinates
                .iter()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt();
            let direction = forward
                .raw_coordinates
                .iter()
                .map(|x| x / norm)
                .collect_vec();
            let physical_surface = |radius: f64| {
                let k = [
                    radius * direction[0],
                    radius * direction[1],
                    radius * direction[2],
                ];
                let l = [
                    radius * direction[3],
                    radius * direction[4],
                    radius * direction[5],
                ];
                let energy_k = (1.0 + k.iter().map(|x| x * x).sum::<f64>()).sqrt();
                let energy_l = (1.0 + l.iter().map(|x| x * x).sum::<f64>()).sqrt();
                let energy_p = (1.0
                    + (0..3)
                        .map(|axis| (k[axis] + l[axis] + external[axis + 1]).powi(2))
                        .sum::<f64>())
                .sqrt();
                energy_k + energy_l + energy_p - external[0]
            };
            let (mut lower, mut upper) = (0.0, 16.0);
            for _ in 0..64 {
                let midpoint = 0.5 * (lower + upper);
                if physical_surface(midpoint) < 0.0 {
                    lower = midpoint;
                } else {
                    upper = midpoint;
                }
            }
            let root = 0.5 * (lower + upper);
            assert!(physical_surface(root).abs() < 1.0e-12);
            let split = root / (root + 5.0 * parameterization.b);
            if external[3] == 0.0 {
                let mut seam = coordinates;
                seam[0] = split + 1.0e-10;
                let error = bridge.forward(SamplingChannelId(0), &seam).unwrap_err();
                assert!(error.downcast_ref::<crate::integrands::process::sampling_maps::SamplingEvaluationError>().is_some(), "{error}");
                seam_coordinates = Some(seam);
            }
            let mut near_shell = coordinates;
            near_shell[0] = split + 1.0e-5;
            // This binary64 draw has a representable radius but insufficient
            // normal-density accuracy in Double. Keep the strict bridge guard
            // and establish the same physical localization at Quad precision.
            if external[3] == 0.0 {
                let error = bridge
                    .forward(SamplingChannelId(0), &near_shell)
                    .unwrap_err();
                assert!(matches!(error.downcast_ref::<crate::integrands::process::sampling_maps::SamplingEvaluationError>(),
                    Some(crate::integrands::process::sampling_maps::SamplingEvaluationError::Unrepresentable {
                        operation: "sampling inverse density consistency", ..
                    })));
            }
            use crate::utils::QuadFloat;
            let quad_bridge = term.compile_sampling_bridge(
                &parameterization,
                &settings,
                &[external.map(|value| F::<QuadFloat>::from_f64(value).0); 2],
                None,
            )?;
            let mapped = quad_bridge.forward(
                SamplingChannelId(0),
                &near_shell.map(QuadFloat::from_f64_exact_binary),
            )?;
            let mapped_radius = mapped
                .raw_coordinates
                .iter()
                .map(|x| F(*x).square())
                .fold(F::<QuadFloat>::default(), |sum, value| sum + value)
                .sqrt()
                .into_ff64()
                .0;
            assert!(
                physical_surface(mapped_radius).abs() < 1.0e-7,
                "production surface profile must focus on the real C equation"
            );

            // Both sides must approach the independently routed physical
            // threshold with the intended square-root radial determinant.
            let angular = coordinates[2..].iter().enumerate().fold(
                std::f64::consts::TAU,
                |jacobian, (i, u)| {
                    jacobian * 2.0 * (1.0 - (2.0 * u - 1.0).powi(2)).sqrt().powi(3 - i as i32)
                },
            );
            let beta = 5.0 * parameterization.b;
            for side in [-1.0, 1.0] {
                for offset in [1.0e-2, 2.0e-3, 4.0e-4] {
                    let mut cube = coordinates;
                    cube[0] = split + side * offset;
                    let point = bridge.forward(SamplingChannelId(0), &cube)?;
                    let r = point
                        .raw_coordinates
                        .iter()
                        .map(|v| v * v)
                        .sum::<f64>()
                        .sqrt();
                    let delta = r - root;
                    assert!(side * delta > 0.0 && side * physical_surface(r) > 0.0);
                    let radial = point.map.jacobian / (angular * r.powi(5));
                    let expected = if delta < 0.0 {
                        2.0 * (root + beta) * (-delta / root).sqrt()
                    } else {
                        2.0 * (root + beta)
                            * (delta / beta).sqrt()
                            * (1.0 + (delta / beta).sqrt()).powi(2)
                    };
                    assert!((radial / expected - 1.0).abs() < 1.0e-7);
                    let inverse = bridge
                        .inverse(SamplingChannelId(0), &point.raw_coordinates)?
                        .expect("full-support map must contain the supplied point");
                    assert!(
                        (point.map.jacobian * inverse.map.inverse_jacobian - 1.0).abs() < 1.0e-7
                    );
                }
            }

            let step = 1.0e-5;
            let mut matrix = vec![vec![0.0; 6]; 6];
            for axis in 0..6 {
                let mut plus = coordinates;
                let mut minus = coordinates;
                plus[axis] += step;
                minus[axis] -= step;
                let plus = bridge.forward(SamplingChannelId(0), &plus)?;
                let minus = bridge.forward(SamplingChannelId(0), &minus)?;
                for (component, row) in matrix.iter_mut().enumerate() {
                    row[axis] = (plus.raw_coordinates[component]
                        - minus.raw_coordinates[component])
                        / (2.0 * step);
                }
            }
            let determinant = SamplingMapAffine::new(matrix, vec![0.0; 6])?.determinant();
            assert!((determinant / forward.map.jacobian - 1.0).abs() < 1.0e-4);
        }

        // Below threshold, at the pinch, and with no negative-shift member,
        // the same canonical channel retains its normalized full-space map.
        for external in [[2.0, 0.0, 0.0, 0.0], [3.0, 0.0, 0.0, 0.0], [0.0; 4]] {
            let bridge =
                term.compile_sampling_bridge(&parameterization, &settings, &[external; 2], None)?;
            assert!(matches!(
                bridge.channels()[0].map,
                CompiledSamplingMap::Surface(_)
            ));
            let forward = bridge.forward(SamplingChannelId(0), &coordinates)?;
            let inverse = bridge
                .inverse(SamplingChannelId(0), &forward.raw_coordinates)?
                .expect("full-support map must contain the supplied point");
            for (left, right) in coordinates.iter().zip(&inverse.map.coordinates) {
                assert!((left - right).abs() < 1.0e-10);
            }
        }
        let externals = [[5.0, 0.0, 0.0, 0.0]; 2];
        model.get_parameter_mut("mass_scalar_1")?.value = Some(Complex::new_re(F(2.0)));
        term.warm_up(&settings, &model)?;
        let bridge =
            term.compile_sampling_bridge(&parameterization, &settings, &externals, None)?;
        assert!(
            matches!(bridge.channels()[0].map, CompiledSamplingMap::Surface(_)),
            "rewarming after a mass update must invalidate the old physical surface"
        );
        model.get_parameter_mut("mass_scalar_1")?.value = Some(Complex::new_re(F(1.0)));
        term.warm_up(&settings, &model)?;

        for (momenta, diagnostic) in [
            (vec![externals[0]], "including the dependent"),
            (vec![[f64::NAN; 4]; 2], "finite external"),
        ] {
            assert!(
                term.compile_sampling_bridge(&parameterization, &settings, &momenta, None)
                    .unwrap_err()
                    .to_string()
                    .contains(diagnostic)
            );
        }
        let cached_masses = term.real_mass_vec.take();
        assert!(
            term.compile_sampling_bridge(&parameterization, &settings, &externals, None)
                .unwrap_err()
                .to_string()
                .contains("warmup mass data")
        );
        term.real_mass_vec = cached_masses;
        let boosted = term.compile_sampling_bridge(
            &parameterization,
            &settings,
            &[[125.0_f64.sqrt(), 0.0, 0.0, 10.0]; 2],
            None,
        )?;
        let CompiledSamplingMap::ImplicitSurface(map) = &boosted.channels()[0].map else {
            panic!("boosted existing surface needs a physical interior-center chart")
        };
        assert!(
            map.center()
                .iter()
                .any(|component| component.abs() > 1.0e-6)
        );
        for _ in 0..2 {
            let point = boosted.forward(SamplingChannelId(0), &coordinates)?;
            let inverse = boosted
                .inverse(SamplingChannelId(0), &point.raw_coordinates)?
                .expect("full-support map must contain the supplied point");
            assert!((point.map.jacobian * inverse.map.inverse_jacobian - 1.0).abs() < 1.0e-8);
        }
        let master_lmb = term
            .multi_channeling_setup
            .graph
            .loop_momentum_basis
            .clone();
        term.multi_channeling_setup
            .graph
            .loop_momentum_basis
            .edge_signatures[EdgeIndex(6)]
        .internal = [1_isize, 1].into_iter().collect();
        assert!(
            term.compile_sampling_bridge(&parameterization, &settings, &externals, None)
                .unwrap_err()
                .to_string()
                .contains("same cycles")
        );
        term.multi_channeling_setup.graph.loop_momentum_basis = master_lmb;
        for (around, subspace, diagnostic) in [
            ("surface(2,4,999)", vec![4, 6], "unknown energy edges"),
            ("surface(3,4)", vec![4, 6], "rank 1"),
        ] {
            let mut invalid = parameterization.clone();
            let definition = invalid
                .sampling_channels
                .channel_definitions
                .get_mut("massive_kite")
                .unwrap()
                .get_mut("C")
                .unwrap();
            definition.around = around.into();
            definition.subspace_lmb = subspace;
            assert!(
                term.compile_sampling_bridge(&invalid, &settings, &externals, None)
                    .unwrap_err()
                    .to_string()
                    .contains(diagnostic)
            );
        }
        // The CFF cache includes positive and negative shifts; only conflicting
        // eligible equations cause an ambiguity. Include both in diagnostics.
        let mut conflicting = term
            .esurfaces
            .iter()
            .find(|surface| {
                surface.energies.iter().map(|edge| edge.0).eq([2, 4, 6])
                    && surface.compute_shift_part_from_momenta(
                        &ExternalFourMomenta::from_iter(
                            [FourMomentum::from_args(F(5.0), F(0.0), F(0.0), F(0.0)); 2],
                        ),
                        &term.graph.loop_momentum_basis,
                    ) < F(0.0)
            })
            .unwrap()
            .clone();
        for (_, coefficient) in &mut conflicting.external_shift {
            *coefficient *= 2;
        }
        term.esurfaces.push(conflicting);
        let error = term
            .compile_sampling_bridge(&parameterization, &settings, &externals, None)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("ambiguous")
                && error.contains("evaluated shift")
                && error.contains("equation")
        );
        // Finite external inputs can overflow a derived temporal shift.
        // This is a native precision failure, never an absence classification.
        let error = term
            .compile_sampling_bridge(
                &parameterization,
                &settings,
                &[[1.0e308, 0.0, 0.0, 0.0]; 2],
                None,
            )
            .unwrap_err();
        assert!(matches!(
            error.downcast_ref::<crate::integrands::process::sampling_maps::SamplingEvaluationError>(),
            Some(crate::integrands::process::sampling_maps::SamplingEvaluationError::Unrepresentable {
                operation: "amplitude sampling external shift", ..
            })
        ), "{error:#}");
        term.esurfaces.pop();
        // Prepare the original binary64 cube once, including every foreign
        // C/D/LMB density and threshold CT. Double physical arithmetic can
        // collapse the focused distance; Quad evaluates the same retained draw.
        {
            use crate::{
                integrands::evaluation::{PreciseEvaluationResult, StabilityStatus},
                settings::runtime::{Precision, StabilityLevelSetting},
                utils::QuadFloat,
            };
            use symbolica::prelude::SingleFloat;
            let stability = &mut runtime_for_rescue.get_mut_settings().stability;
            stability.rotation_axis.clear();
            stability.levels = vec![
                StabilityLevelSetting::default_double(),
                StabilityLevelSetting::default_quad(),
            ];
            runtime_for_rescue.warm_up(&model)?;
            let source = Sample::Continuous(F(1.0), seam_coordinates.unwrap().map(F).to_vec());
            let rescued = runtime_for_rescue.evaluate_sample_precise(
                &source,
                &model,
                F(1.0),
                false,
                Complex::new_zero(),
            )?;
            let PreciseEvaluationResult::Quad(rescued) = rescued else {
                panic!("the actual kite surface must reconstruct at Quad precision");
            };
            assert_eq!(
                rescued.evaluation_metadata.final_precision(),
                Some(Precision::Quad)
            );
            assert_eq!(rescued.evaluation_metadata.stability_results.len(), 2);
            assert_eq!(
                rescued.evaluation_metadata.stability_results[0].precision,
                Precision::Double
            );
            // The Double body is evaluated once without rotation probes, so its
            // nonfinite result has Unknown status; a pre-body map failure used
            // to record Unstable(0). The second level rescues this same draw.
            assert_eq!(
                rescued.evaluation_metadata.stability_results[0].status,
                StabilityStatus::Unknown,
                "{:?}",
                rescued.evaluation_metadata.stability_results
            );
            assert!(
                rescued.integrand_result.re.0.is_finite()
                    && rescued.integrand_result.im.0.is_finite()
            );
            assert_ne!(
                rescued.integrand_result,
                Complex::new_re(rescued.integrand_result.re.zero())
            );
            assert_eq!(
                rescued.parameterization_jacobian,
                Some(F::<QuadFloat>::default().one())
            );
            let reference = crate::integrands::process::GaussianReferenceFunction::new(
                1.2,
                vec![0.4, -0.3, 0.2, -0.2, 0.1, 0.35],
            )?;
            let reference_rescued =
                runtime_for_rescue.evaluate_reference_sample_detailed(&source, &reference)?;
            assert_eq!(
                reference_rescued
                    .evaluation
                    .evaluation_metadata
                    .final_precision(),
                Some(Precision::Double)
            );
            assert_eq!(
                reference_rescued
                    .evaluation
                    .evaluation_metadata
                    .stability_results
                    .len(),
                1
            );
            runtime_for_rescue.get_mut_settings().stability.levels =
                vec![StabilityLevelSetting::default_quad()];
            runtime_for_rescue.warm_up(&model)?;
            let direct = runtime_for_rescue.evaluate_sample_precise(
                &source,
                &model,
                F(1.0),
                false,
                Complex::new_zero(),
            )?;
            let PreciseEvaluationResult::Quad(direct) = direct else {
                unreachable!()
            };
            assert_eq!(rescued.integrand_result, direct.integrand_result);
            // The smooth reference needs no native map reconstruction or
            // physical threshold rescue. Compare identical Double evaluations.
            let mut reference_control = runtime_for_rescue.clone();
            reference_control.get_mut_settings().stability.levels =
                vec![StabilityLevelSetting::default_double()];
            reference_control.warm_up(&model)?;
            let reference_direct =
                reference_control.evaluate_reference_sample_detailed(&source, &reference)?;
            assert_eq!(
                reference_rescued.evaluation.integrand_result,
                reference_direct.evaluation.integrand_result
            );
            assert_eq!(
                reference_rescued.moments.second_moment,
                reference_direct.moments.second_moment
            );

            let ProcessIntegrand::Amplitude(runtime) = &mut runtime_for_rescue else {
                unreachable!()
            };
            let setup = &runtime.data.graph_terms[0].multi_channeling_setup;
            assert!(
                setup.sampling_bridge::<f64>().is_err(),
                "Quad-only warmup must not require Double geometry"
            );
            assert!(setup.sampling_bridge::<QuadFloat>().is_ok());
            assert!(setup.sampling_bridge::<crate::utils::ArbPrec>().is_ok());
            let catalogue = setup.sampling_catalogue.as_ref().unwrap() as *const _ as usize;
            let programs = setup.sampling_programs.as_ref().unwrap() as *const _ as usize;
            runtime.prepare_sampling_precision::<crate::utils::ArbPrec>()?;
            let setup = &runtime.data.graph_terms[0].multi_channeling_setup;
            assert_eq!(
                setup.sampling_catalogue.as_ref().unwrap() as *const _ as usize,
                catalogue
            );
            assert_eq!(
                setup.sampling_programs.as_ref().unwrap() as *const _ as usize,
                programs
            );
            assert!(setup.sampling_bridge::<crate::utils::ArbPrec>().is_ok());
            // Compare the very same binary draw in Quad and Arb, evaluating
            // C independently at both resulting momenta before any narrowing.
            // This checks localization/density numerically, not a root enclosure.
            use crate::utils::ArbPrec;
            let one = F::<ArbPrec>::default().one();
            let eta = |point: &[ArbPrec]| {
                let energy = |vectors: &[&[ArbPrec]]| {
                    (one.clone()
                        + (0..3)
                            .map(|i| {
                                vectors
                                    .iter()
                                    .fold(one.zero(), |sum, vector| sum + F(vector[i].clone()))
                                    .square()
                            })
                            .fold(one.zero(), |sum, square| sum + square))
                    .sqrt()
                };
                energy(&[&point[..3]])
                    + energy(&[&point[3..]])
                    + energy(&[&point[..3], &point[3..]])
                    - one.from_i64(5)
            };
            for sign in [-1.0, 1.0] {
                let mut cube = seam_coordinates.unwrap();
                cube[0] -= if sign < 0.0 { 2.0e-10 } else { 0.0 };
                let quad = setup.sampling_bridge::<QuadFloat>()?.forward(
                    SamplingChannelId(0),
                    &cube.map(QuadFloat::from_f64_exact_binary),
                )?;
                let arb = setup.sampling_bridge::<ArbPrec>()?.forward(
                    SamplingChannelId(0),
                    &cube.map(ArbPrec::from_f64_exact_binary),
                )?;
                let quad_point = quad
                    .raw_coordinates
                    .iter()
                    .map(|v| F(*v).higher().0)
                    .collect_vec();
                let eta_quad = eta(&quad_point);
                let eta_arb = eta(&arb.raw_coordinates);
                assert!(eta_arb.clone().into_ff64().0 * sign > 0.0);
                assert!(((&eta_quad / &eta_arb - &one).abs()).into_ff64().0 < 1.0e-6);
                // Foreign scores must be evaluated at this common supplied
                // Quad point, rather than at the nearby Arb forward point.
                let arb_inverse = setup
                    .sampling_bridge::<ArbPrec>()?
                    .inverse(SamplingChannelId(0), &quad_point)?
                    .expect("full-support map must contain the supplied point");
                let ratio = F(quad.map.jacobian).higher() * F(arb_inverse.map.inverse_jacobian);
                assert!((ratio - &one).abs().into_ff64().0 < 1.0e-6);
            }
        }
        // Reuse the generated amplitude at fixed raw momenta, so this range
        // regression is independent of a focused map or root certificate.
        // The mass input and its square fit binary64; the full amplitude does not.
        {
            use crate::{
                integrands::evaluation::PreciseEvaluationResult,
                integrands::process::MomentumSpaceEvaluationInput,
                settings::runtime::StabilityLevelSetting,
            };
            use symbolica::prelude::SingleFloat;
            let mut heavy_model = model.clone();
            heavy_model.get_parameter_mut("mass_scalar_1")?.value =
                Some(Complex::new_re(F(1.0e60)));
            runtime_for_rescue.get_mut_settings().stability.levels =
                vec![StabilityLevelSetting::default_arb()];
            runtime_for_rescue.warm_up(&heavy_model)?;
            let input = MomentumSpaceEvaluationInput {
                loop_momenta: vec![
                    ThreeMomentum::new(F(0.7), F(0.2), F(-0.4)),
                    ThreeMomentum::new(F(-0.3), F(0.8), F(0.1)),
                ],
                integrator_weight: F(1.0),
                graph_id: Some(0),
                group_id: None,
                orientation: None,
                channel_id: None,
            };
            let result = runtime_for_rescue.evaluate_momentum_configuration_precise(
                &heavy_model,
                &input,
                true,
            )?;
            let PreciseEvaluationResult::Arb(result) = result else {
                panic!("explicit Arb output was not retained")
            };
            assert!(
                result.integrand_result.re.0.is_finite()
                    && result.integrand_result.im.0.is_finite()
            );
            assert_ne!(
                result.integrand_result,
                Complex::new_re(result.integrand_result.re.zero())
            );
            assert!(
                [&result.integrand_result.re, &result.integrand_result.im]
                    .into_iter()
                    .any(|value| value != &value.zero() && value.clone().into_ff64().0 == 0.0)
            );
            assert!(!result.evaluation_metadata.is_nan);
            assert!(
                runtime_for_rescue
                    .evaluate_momentum_configuration(&heavy_model, &input, true)
                    .unwrap_err()
                    .to_string()
                    .contains("f64 integration/reporting boundary")
            );
        }
        Ok(())
    }
}
