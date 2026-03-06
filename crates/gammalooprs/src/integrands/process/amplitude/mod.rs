use std::{
    fs::{self},
    path::Path,
};

use bincode_trait_derive::{Decode, Encode};

use color_eyre::Result;

use eyre::Context;
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
            EsurfaceCollection, EsurfaceID, ExistingEsurfaces, GroupEsurfaceId, get_representative,
        },
        expression::AmplitudeOrientationID,
    },
    graph::{
        FeynmanGraph, Graph, GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex,
        LoopMomentumBasis,
    },
    integrands::HasIntegrand,
    integrands::evaluation::EvaluationResult,
    integrands::process::{ChannelIndex, ParamBuilder, evaluators::EvaluatorStack},
    model::Model,
    momentum::sample::{ExternalIndex, MomentumSample},
    momentum::signature::SignatureLike,
    momentum::{Rotation, RotationMethod},
    processes::{AmplitudeGraph, GroupDerivedData},
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        amplitude_counterterm::AmplitudeCountertermData,
        overlap::{OverlapInput, SingleGraphOverlapData, find_maximal_overlap},
    },
    utils::{W_, serde_utils::SmartSerde, symbolica_ext::LOGPRINTOPTS},
};

use super::{
    GraphTerm, LmbMultiChannelingSetup, ProcessIntegrandImpl, create_grid, evaluate_sample,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeGraphTerm {
    pub original_integrand: EvaluatorStack,
    pub orientations: TiVec<AmplitudeOrientationID, EdgeVec<Orientation>>,
    pub orientation_filter: SubSet<AmplitudeOrientationID>,
    pub esurfaces: EsurfaceCollection,
    pub threshold_counterterm: AmplitudeCountertermData,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub graph: Graph,
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder,
    pub real_mass_vec: Option<EdgeVec<Option<F<f64>>>>,
}

/// Num(sigma_1,sigma_2,...)*(CFF_1 delta(edge(1),1) delta_(1,1,1,-1,1)+CFF_3 delta_(1,1,1,-1,1)+CFF_2 delta_(1,1,1,-1,1))
impl AmplitudeGraphTerm {
    pub fn from_amplitude_graph(
        graph: &AmplitudeGraph,
        own_group_position: GraphGroupPosition,
        esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<Self> {
        let orientations: TiVec<AmplitudeOrientationID, EdgeVec<Orientation>> = graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|a| a.data.orientation.clone())
            .collect();

        debug!(orientation_parametric_integrand = %graph.derived_data.all_mighty_integrand.printer(LOGPRINTOPTS), "Building evaluator for all orientations \n{}",graph.graph.param_builder.table());

        let original_integrand = EvaluatorStack::new(
            &[&graph.derived_data.all_mighty_integrand],
            &graph.graph.param_builder,
            orientations.as_slice().as_ref(),
            &settings.generation.evaluator,
        )?;

        let mut threshold_counterterm = AmplitudeCountertermData::new_empty(own_group_position);

        threshold_counterterm.evaluators = graph
            .derived_data
            .threshold_counterterms
            .iter()
            .map(|ct| ct.to_evaluator(&graph.graph.param_builder, &orientations, settings))
            .collect();

        threshold_counterterm.esurface_map = esurface_map;

        Ok(AmplitudeGraphTerm {
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
        })
    }

    #[instrument(
          name = "compile",
          level = "info",
          skip(self, path, override_existing, settings),
          fields(
              graph.name = %self.graph.name,
              path = %path.as_ref().display(),
          )
      )]
    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
    ) -> Result<()> {
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
            settings,
        )?;

        self.threshold_counterterm
            .compile(&graph_path, override_existing, settings);

        Ok(())
    }

    #[instrument(
          skip_all,
          fields(
              term.name = %self.name(),
          )
    )]
    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        settings: &RuntimeSettings,
        rotation: &Rotation,
        channel_id: Option<(ChannelIndex, F<T>)>,
    ) -> Result<Complex<F<T>>> {
        let (momentum_sample, prefactor) = if let Some((channel_id, alpha)) = channel_id {
            self.multi_channeling_setup
                .reinterpret_loop_momenta_and_compute_prefactor(
                    channel_id,
                    momentum_sample,
                    0,
                    model,
                    &alpha,
                )
        } else {
            (momentum_sample.clone(), momentum_sample.one())
        };

        let hel = settings.kinematics.externals.get_helicities();
        let orientations =
            momentum_sample.orientations(&self.orientation_filter, &self.orientations);

        debug!("loop_moms: {}", momentum_sample.loop_moms());
        debug!("jacobian: {:16e}", momentum_sample.jacobian());
        // debug!("Og paramBuilder: \n{}", self.param_builder.table());

        let input = T::get_parameters(
            &mut self.param_builder,
            (settings.general.enable_cache, settings.general.debug_cache),
            &self.graph,
            &momentum_sample,
            hel,
            &settings.additional_params(),
            None,
            None,
            None,
        );

        let result = self
            .original_integrand
            .evaluate(input, orientations, settings)?
            .pop()
            .unwrap();
        // debug!("parambuilder 244: {}", self.param_builder);
        let sum_of_cts = self.threshold_counterterm.evaluate(
            &momentum_sample,
            &self.graph,
            model,
            &self.esurfaces,
            rotation,
            settings,
            &mut self.param_builder,
            orientations,
        );

        debug!(
            bare_cff = format!("{result:16e}"),
            "{}: {result:16e}", self.graph.name
        );
        debug!(cts = format!("{sum_of_cts:16e}"), "{}", self.graph.name);
        debug!("result: {result:16e}");
        debug!("sum_of_cts: {sum_of_cts:16e}");

        debug!(
            value = format!("{sum_of_cts:16e}"),
            "evaluated sum of threshold counterterms"
        );

        let diff = result - sum_of_cts;

        Ok(diff * prefactor)
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
                values[polarization_start] = pol_value.clone();
                polarization_start += multiplicative_offset;
            }
        }

        self.graph
            .param_builder
            .m_uv_value(Complex::new_re(F(settings.general.m_uv)));
        self.graph
            .param_builder
            .mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq)));
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

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_num_channels(&self) -> usize {
        self.multi_channeling_setup.channels.len()
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        settings: &RuntimeSettings,
        rotation: &Rotation,
        channel_id: Option<(ChannelIndex, F<T>)>,
    ) -> Result<Complex<F<T>>> {
        self.evaluate_impl(momentum_sample, model, settings, rotation, channel_id)
    }

    fn get_num_orientations(&self) -> usize {
        self.orientations.len()
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
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeIntegrandData {
    pub rotations: Option<Vec<Rotation>>,
    pub name: String,
    pub loop_cache_id: usize,
    pub external_cache_id: usize,
    /// Cache ID for the base (unrotated) external momentum configuration
    pub base_external_cache_id: usize,
    pub graph_terms: Vec<AmplitudeGraphTerm>,
    pub external_signature: SignatureLike<ExternalIndex>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub group_derived_data: TiVec<GroupId, GroupDerivedData>,
}

pub mod export;
pub mod load;

impl AmplitudeIntegrand {
    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path> + Sync,
        override_existing: bool,
        settings: &GlobalSettings,
        thread_pool: &rayon::ThreadPool,
    ) -> Result<()> {
        thread_pool.install(|| {
            self.data
                .graph_terms
                .par_iter_mut()
                .try_for_each(|a| a.compile(path.as_ref(), override_existing, settings))
        })?;

        Ok(())
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

        Ok(AmplitudeIntegrand { settings, data })
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
                    .filter_map(|(group_esurface_id, esurface_map)| {
                        let esurface_exists = esurface_map
                            .iter_enumerated()
                            .find_map(|(graph_group_pos, option_esurface_id)| {
                                option_esurface_id.map(|esurface_id| {
                                    // extreme indexing
                                    let graph = &self.data.graph_terms[self
                                        .data
                                        .graph_group_structure[group_id][graph_group_pos]];

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
                                    let (graph_group_pos, esurface_id) = get_representative(
                                        &self.data.group_derived_data[group_id].esurface_map
                                            [*group_esurface_id],
                                    )
                                    .unwrap();
                                    let graph_id =
                                        self.data.graph_group_structure[group_id][graph_group_pos];
                                    let graph = &self.data.graph_terms[graph_id].graph;
                                    let lmb_reps = graph.integrand_replacement(
                                        &graph.full_filter(),
                                        &graph.loop_momentum_basis,
                                        &[W_.x___],
                                    );

                                    let esurface =
                                        &self.data.graph_terms[graph_id].esurfaces[esurface_id];
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

        Ok(())
    }

    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.as_ref().expect("forgot warmup").iter()
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
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

    fn get_group(&self, group_id: GroupId) -> &GraphGroup {
        &self.data.graph_group_structure[group_id]
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor<'_> {
        DependentMomentaConstructor::Amplitude(&self.data.external_signature)
    }

    fn get_group_structure(&self) -> &TiVec<GroupId, GraphGroup> {
        &self.data.graph_group_structure
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
