use std::{
    collections::HashSet,
    fs::{self},
    path::Path,
};

use bincode_trait_derive::{Decode, Encode};

use color_eyre::{Help, Result};

use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
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
use tracing::{debug, info, instrument, warn};
use typed_index_collections::TiVec;

use crate::{
    DependentMomentaConstructor, F, FloatLike, GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::{
            EsurfaceCollection, ExistingEsurfaces, GroupEsurfaceId, RaisedEsurfaceId,
            get_representative,
        },
        expression::OrientationID,
        surface::HybridSurfaceID,
    },
    graph::{
        FeynmanGraph, Graph, GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex,
        LoopMomentumBasis,
    },
    integrands::{
        HasIntegrand,
        evaluation::{EvaluationResult, GraphEvaluationResult},
        process::{
            ChannelIndex, ParamBuilder,
            evaluators::{ActiveF64Backend, EvaluatorStack},
        },
    },
    model::Model,
    momentum::{
        Helicity, Rotation, RotationMethod, SignOrZero,
        sample::{ExternalIndex, MomentumSample},
        signature::SignatureLike,
    },
    observables::{AdditionalWeightKey, EventProcessingRuntime, GenericEvent},
    processes::{AmplitudeGraph, GraphGenerationStats, GroupDerivedData},
    settings::{GlobalSettings, RuntimeSettings, global::FrozenCompilationMode},
    subtraction::{
        amplitude_counterterm::{
            AmplitudeCountertermAtom, AmplitudeCountertermData, AmplitudeCountertermEvaluation,
            OverlapStructureWithKinematics,
        },
        overlap::{OverlapInput, SingleGraphOverlapData, find_maximal_overlap},
    },
    utils::{ArbPrec, W_, serde_utils::SmartSerde, symbolica_ext::LOGPRINTOPTS},
};

use super::{
    GraphTerm, GraphTermEvaluationContext, LmbMultiChannelingSetup, ProcessIntegrandImpl,
    RuntimeCache, create_grid, evaluate_sample, filtered_orientation_count,
    format_lmb_channel_label, format_orientation_label, histogram_process_info_for_integrand,
    prepare_buffered_event, resolve_visible_orientation_id,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeGraphTerm {
    pub original_integrand: EvaluatorStack,
    pub orientations: TiVec<OrientationID, EdgeVec<Orientation>>,
    pub orientation_filter: SubSet<OrientationID>,
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

/// Num(sigma_1,sigma_2,...)*(CFF_1 delta(edge(1),1) delta_(1,1,1,-1,1)+CFF_3 delta_(1,1,1,-1,1)+CFF_2 delta_(1,1,1,-1,1))
impl AmplitudeGraphTerm {
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
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let mut stats = GraphGenerationStats::default();
        let selected_generation_orientations = graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .filter(|orientation| settings.generation.orientation_pattern.filter(*orientation))
            .collect_vec();
        let orientations: TiVec<OrientationID, EdgeVec<Orientation>> =
            selected_generation_orientations
                .iter()
                .map(|a| a.data.orientation.clone())
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

        crate::debug_tags!(#generation, #graph, #orientation, #compile, #dump;
            orientation_parametric_integrand = %graph.derived_data.all_mighty_integrand.printer(LOGPRINTOPTS),
            "Building evaluator for all orientations \n{}",
            graph.graph.param_builder.table()
        );

        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let (original_integrand, evaluator_timings) = EvaluatorStack::new_with_timings(
            &[&graph.derived_data.all_mighty_integrand],
            &graph.graph.param_builder,
            orientations.as_slice().as_ref(),
            None,
            &settings.generation.evaluator,
        )?;
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        stats.add_evaluator_build_timings(evaluator_timings);
        stats.evaluator_count += original_integrand.generic_evaluator_count();

        let mut threshold_counterterm = AmplitudeCountertermData::new_empty(own_group_position);
        let mut threshold_evaluators =
            Vec::with_capacity(graph.derived_data.threshold_counterterms.len());
        let selected_generation_raised_esurfaces: HashSet<RaisedEsurfaceId> =
            if graph.derived_data.threshold_counterterms.is_empty() {
                HashSet::new()
            } else {
                selected_generation_esurfaces
                    .iter()
                    .map(|esurface_id| graph.derived_data.raised_esurface_ids[*esurface_id])
                    .collect()
            };
        let active_mask: TiVec<RaisedEsurfaceId, bool> = graph
            .derived_data
            .threshold_counterterms
            .iter_enumerated()
            .map(|(raised_esurface_id, _)| {
                selected_generation_raised_esurfaces.contains(&raised_esurface_id)
            })
            .collect();
        for (raised_esurface_id, ct) in graph.derived_data.threshold_counterterms.iter_enumerated()
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
                settings,
            );
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
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
        // `generated_mask` answers "does this threshold slot exist symbolically at all?",
        // while `active_mask` answers "is that generated slot compatible with the
        // selected generation-time orientation subset for this evaluator build?"
        threshold_counterterm.active_mask = active_mask;
        threshold_counterterm.raised_data = graph.derived_data.raised_data.clone();
        threshold_counterterm.helper_evaluators = graph
            .derived_data
            .raised_data
            .pass_two_evaluator
            .clone()
            .unwrap_or_default();
        stats.evaluator_count += threshold_counterterm.helper_evaluators.len();

        threshold_counterterm.esurface_map = esurface_map;

        Ok((
            AmplitudeGraphTerm {
                orientation_filter: SubSet::full(orientations.len()),
                orientations,
                original_integrand,
                tropical_sampler: graph.derived_data.tropical_sampler.clone(),
                graph: graph.graph.clone(),
                multi_channeling_setup: LmbMultiChannelingSetup {
                    channels: TiVec::new(),
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
        channel_id: Option<ChannelIndex>,
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
        event.cut_info.orientation_id = orientation_id;
        event.cut_info.lmb_channel_id = channel_id.map(usize::from);
        event.cut_info.lmb_channel_edge_ids =
            channel_id.map(|channel_id| self.multi_channeling_setup.channel_edge_ids(channel_id));

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

    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        context: &mut GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<(Complex<F<T>>, AmplitudeCountertermEvaluation<T>)> {
        let (momentum_sample, prefactor) = if let Some((channel_id, alpha)) = &context.channel_id {
            self.multi_channeling_setup
                .reinterpret_loop_momenta_and_compute_prefactor(
                    *channel_id,
                    momentum_sample,
                    0,
                    context.model,
                    alpha,
                )
        } else {
            (momentum_sample.clone(), momentum_sample.one())
        };

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
            &momentum_sample,
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
                context.record_primary_timing,
            )?
            .pop()
            .unwrap()
            .unwrap_real();
        // debug!("parambuilder 244: {}", self.param_builder);
        let counterterm_evaluation = self.threshold_counterterm.evaluate(
            &momentum_sample,
            &self.graph,
            context.model,
            &self.esurfaces,
            context.rotation,
            context.settings,
            &mut self.param_builder,
            orientations,
            context.evaluation_metadata,
            context.record_primary_timing,
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

        let diff = result - sum_of_cts.clone();

        Ok((
            diff * prefactor.clone(),
            AmplitudeCountertermEvaluation {
                total: sum_of_cts * prefactor.clone(),
                local_counterterms: counterterm_evaluation
                    .local_counterterms
                    .into_iter()
                    .map(|counterterm| counterterm * prefactor.clone())
                    .collect(),
            },
        ))
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
            .mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq())));
        self.graph.param_builder.update_model_values(model);

        self.param_builder = self.graph.param_builder.clone();

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

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_num_channels(&self) -> usize {
        self.multi_channeling_setup.channels.len()
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        mut context: GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<GraphEvaluationResult<T>> {
        let event_channel_id = context
            .channel_id
            .as_ref()
            .map(|(channel_id, _)| *channel_id);
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
                integrand_result: Complex::new_re(momentum_sample.zero()),
                event_groups: crate::observables::GenericEventGroupList::default(),
                event_processing_time: prepared_event.event_processing_time,
                generated_event_count: prepared_event.generated_event_count,
                accepted_event_count: prepared_event.accepted_event_count,
            });
        }

        let (integrand_result, counterterm_evaluation) =
            self.evaluate_impl(momentum_sample, &mut context)?;

        let mut event_groups = crate::observables::GenericEventGroupList::default();
        let generated_event_count = prepared_event.generated_event_count;
        let accepted_event_count = prepared_event.accepted_event_count;
        if let Some(mut event) = prepared_event.buffered_event {
            event.weight = integrand_result.clone();

            if context.settings.general.store_additional_weights_in_event {
                let original = integrand_result.clone() + counterterm_evaluation.total;
                event
                    .additional_weights
                    .weights
                    .insert(AdditionalWeightKey::Original, original);
                for (subset_index, threshold_counterterm) in counterterm_evaluation
                    .local_counterterms
                    .into_iter()
                    .enumerate()
                {
                    event.additional_weights.weights.insert(
                        AdditionalWeightKey::ThresholdCounterterm { subset_index },
                        -threshold_counterterm,
                    );
                }
            }

            event_groups.push_singleton(event);
        }

        Ok(GraphEvaluationResult {
            integrand_result,
            event_groups,
            event_processing_time: prepared_event.event_processing_time,
            generated_event_count,
            accepted_event_count,
        })
    }

    fn get_num_orientations(&self) -> usize {
        filtered_orientation_count(&self.orientation_filter, &self.orientations)
    }

    fn get_tropical_sampler(&self) -> &SampleGenerator<3> {
        self.tropical_sampler
            .as_ref()
            .expect("Tropical sampler should be set.")
    }

    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64> {
        &mut self.param_builder
    }

    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>> {
        self.real_mass_vec
            .as_ref()
            .expect("real mass vector should be set")
            .clone()
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
}

pub mod export;
pub mod load;

impl AmplitudeIntegrand {
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
            FrozenCompilationMode::Symjit => {
                let mut compile_times = Vec::with_capacity(self.data.graph_terms.len());
                for graph_term in &mut self.data.graph_terms {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let compile_started = std::time::Instant::now();
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

    pub(crate) fn invalidate_event_processing_runtime(&mut self) {
        self.event_processing_runtime.invalidate();
    }

    pub(crate) fn get_existing_esurfaces(
        &self,
        model: &Model,
    ) -> TiVec<GroupId, ExistingEsurfaces> {
        self.data
            .group_derived_data
            .iter_enumerated()
            .map(|(group_id, group_derived_data)| {
                let group_esurface_map = &group_derived_data.esurface_map;
                let external_moms = self
                    .settings
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                        &self.data.external_signature,
                    ))
                    .expect("could not get externals");

                let e_cm = self.settings.kinematics.e_cm;

                group_esurface_map
                    .iter_enumerated()
                    .filter_map(|(group_esurface_id, raised_esurface_map)| {
                        let esurface_exists = raised_esurface_map
                            .iter_enumerated()
                            .find_map(|(graph_group_pos, option_raised_esurface_id)| {
                                option_raised_esurface_id.map(|raised_esurface_id| {
                                    // extreme indexing
                                    let graph = &self.data.graph_terms[self
                                        .data
                                        .graph_group_structure[group_id][graph_group_pos]];
                                    let esurface_id =
                                        graph.threshold_counterterm.raised_data.raised_groups
                                            [raised_esurface_id]
                                            .esurface_ids[0];

                                    let esurface = &graph.esurfaces[esurface_id];
                                    esurface.exists(
                                        &external_moms,
                                        &graph.graph.loop_momentum_basis,
                                        &graph.graph.get_real_mass_vector(model),
                                        &F(e_cm),
                                    )
                                })
                            })
                            .expect("no graph in group has this esurface, map corrupted");

                        if esurface_exists {
                            Some(group_esurface_id)
                        } else {
                            None
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn validate_runtime_threshold_counterterms(
        &self,
        existing_esurfaces: &TiVec<GroupId, ExistingEsurfaces>,
    ) -> Result<()> {
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
                            "Amplitude integrand '{}' was generated with specialized threshold-subtraction assumptions, but the current runtime model parameters require a trimmed threshold counterterm for graph '{}' and group e-surface {} ({})",
                            self.name(),
                            graph_term.graph.name,
                            group_esurface_id.0,
                            self.data.group_derived_data[group_id].esurface_atoms[group_esurface_id]
                        ))
                        .with_note(|| {
                            "Regenerate the integrand or restore compatible shared/per-integrand model parameters.".to_string()
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
        let e_cm = F(self.settings.kinematics.e_cm);
        let constructor = DependentMomentaConstructor::Amplitude(&self.data.external_signature);
        let masses = self.data.graph_terms[0].graph.get_external_masses(model);

        self.settings
            .kinematics
            .externals
            .improve_and_cache(constructor, &masses, &e_cm)?;
        self.warn_on_off_shell_external_states(model)?;

        let thresholds_generated = self
            .data
            .graph_terms
            .iter()
            .all(|term| !term.threshold_counterterm.evaluators.is_empty());

        if !thresholds_generated && !self.settings.subtraction.disable_threshold_subtraction {
            warn!(
                "Not all graphs have threshold counterterms generated, but threshold subtraction is not disabled. disable runtime threshold subtraction to remove this warning"
            );
            self.settings.subtraction.disable_threshold_subtraction = true;
        }

        let is_tree_level = self.data.graph_terms[0].graph.get_loop_number() == 0;

        if !self.settings.subtraction.disable_threshold_subtraction && !is_tree_level {
            debug!("esurface existence check");
            let existing_esurfaces = self.get_existing_esurfaces(model);
            self.validate_runtime_threshold_counterterms(&existing_esurfaces)?;
            for (group_id, existing_esurfaces) in existing_esurfaces.iter_enumerated() {
                debug!(
                    "solving overlap for group {}, number of thresholds: {}",
                    group_id.0,
                    existing_esurfaces.len()
                );

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

                let overlap_input = OverlapInput {
                    graph_data,
                    settings: &self.settings,
                    group_esurface_map: self.data.group_derived_data[group_id].esurface_map.clone(),
                };

                let external_moms = self
                    .settings
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                        &self.data.external_signature,
                    ))
                    .expect("could not get externals");

                let mut overlap =
                    find_maximal_overlap(&overlap_input, existing_esurfaces, &external_moms)
                        .with_context(|| {
                            let readable_esurfaces = existing_esurfaces
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

                overlap
                    .build_evaluators(
                        &self.data.group_derived_data[group_id].esurface_atoms,
                        &OptimizationSettings::default(),
                        self.data
                            .graph_terms
                            .first()
                            .unwrap()
                            .graph
                            .get_loop_number(),
                        external_moms.len(),
                        self.data
                            .graph_terms
                            .first()
                            .unwrap()
                            .param_builder
                            .pairs
                            .model_parameters
                            .params
                            .clone(),
                    )
                    .with_context(|| {
                        format!(
                            "Failed to build multi-channeling evaluators for group {}",
                            group_id.0
                        )
                    })?;

                info!(
                    "overlap structure of group {}: {:?}",
                    group_id.0,
                    overlap
                        .overlap_groups
                        .iter()
                        .map(|group| group.existing_esurfaces.len())
                        .collect_vec()
                );

                for graph_id in self.data.graph_group_structure[group_id].into_iter() {
                    self.data.graph_terms[graph_id]
                        .threshold_counterterm
                        .overlap = overlap.clone();
                }
            }
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
