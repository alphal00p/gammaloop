use std::{
    fs::{self, File},
    path::Path,
};

use bincode_trait_derive::{Decode, Encode};

use bitvec::vec::BitVec;

use color_eyre::Result;

use eyre::{eyre, Context};
use itertools::Itertools;
use libc::group;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use momtrop::SampleGenerator;
use spenso::algebra::complex::Complex;
use symbolica::{
    evaluate::OptimizationSettings,
    numerical_integration::{Grid, Sample},
};
use tracing::{debug, instrument};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        esurface::{self, get_existing_esurfaces, EsurfaceCollection, ExistingEsurfaces},
        expression::{AmplitudeOrientationID, GraphOrientation},
    },
    evaluation_result::EvaluationResult,
    gammaloop_integrand::ParamBuilder,
    graph::{
        FeynmanGraph, Graph, GraphGroup, GraphGroupPosition, GroupId, LmbIndex, LoopMomentumBasis,
    },
    integrands::HasIntegrand,
    model::Model,
    momentum::{Rotation, RotationMethod},
    momentum_sample::{ExternalIndex, MomentumSample},
    processes::{AmplitudeDerivedData, AmplitudeGraph, GroupDerivedData},
    settings::{GlobalSettings, RuntimeSettings},
    signature::SignatureLike,
    status_debug, status_info, status_warn,
    subtraction::{
        amplitude_counterterm::AmplitudeCountertermData,
        overlap::{find_maximal_overlap, OverlapInput, SingleGraphOverlapData},
    },
    utils::bitvec_ext::BinVec,
    DependentMomentaConstructor, FloatLike, GammaLoopContext, GammaLoopContextContainer, F,
};

use super::{
    create_grid, evaluate_sample, GammaloopIntegrand, GenericEvaluator, GenericEvaluatorFloat,
    GraphTerm, LmbMultiChannelingSetup,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeGraphTerm {
    pub orientation_parametric_integrand: GenericEvaluator,
    pub iterative_integrand_evaluator: Option<GenericEvaluator>,
    pub orientations: TiVec<AmplitudeOrientationID, EdgeVec<Orientation>>,
    pub orientation_filter: BinVec,
    pub esurfaces: EsurfaceCollection,
    pub threshold_counterterm: AmplitudeCountertermData,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub graph: Graph,
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder,
}

/// Num(sigma_1,sigma_2,...)*(CFF_1 delta(edge(1),1) delta_(1,1,1,-1,1)+CFF_3 delta_(1,1,1,-1,1)+CFF_2 delta_(1,1,1,-1,1))

impl AmplitudeGraphTerm {
    pub fn from_amplitude_graph(
        graph: &AmplitudeGraph,
        own_group_position: GraphGroupPosition,
        model: &Model,
        settings: &GlobalSettings,
    ) -> Result<Self> {
        let param_builder = ParamBuilder::new(&graph.graph, model);

        let orientations: TiVec<AmplitudeOrientationID, EdgeVec<Orientation>> = graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|a| a.data.orientation.clone())
            .collect();

        let orientation_parametric_integrand = GenericEvaluator::new_from_builder(
            [graph.derived_data.all_mighty_integrand.clone()],
            &param_builder,
            OptimizationSettings::default(),
        )
        .unwrap();

        let iterative_integrand_evaluator = if settings
            .generation
            .evaluator_settings
            .iterative_orientation_optimization
        {
            GenericEvaluator::new_from_builder(
                orientations
                    .iter()
                    .map(|a| a.select(&graph.derived_data.all_mighty_integrand)),
                &param_builder,
                OptimizationSettings::default(),
            )
        } else {
            None
        };

        let mut threshold_counterterm = AmplitudeCountertermData::new_empty(own_group_position);

        threshold_counterterm.evaluators = graph
            .derived_data
            .threshold_counterterms
            .iter()
            .map(|ct| ct.to_evaluator(&param_builder, OptimizationSettings::default()))
            .collect();

        Ok(AmplitudeGraphTerm {
            orientation_filter: BinVec(BitVec::repeat(true, orientations.len())),
            orientations,
            iterative_integrand_evaluator,
            orientation_parametric_integrand,
            tropical_sampler: graph.derived_data.tropical_sampler.clone(),
            graph: graph.graph.clone(),
            multi_channeling_setup: graph
                .derived_data
                .multi_channeling_setup
                .clone()
                .expect("multi_channeling_setup should have been created"),
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

            param_builder,
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
    ) {
        let graph_path = path.as_ref().join(&self.graph.name);

        let r = fs::create_dir_all(&graph_path)
            .with_context(|| {
                format!(
                    "Trying to create directory to save amplitude {}",
                    graph_path.display()
                )
            })
            .unwrap();
        self.orientation_parametric_integrand.compile(
            graph_path.join("orientation_parametric_integrand"),
            format!("{}_orientation_parametric_integrand", &self.graph.name,),
            graph_path.join("orientation_parametric_integrand"),
            settings.generation.gammaloop_compile_options.inline_asm(),
        );

        self.threshold_counterterm
            .compile(&graph_path, override_existing, settings);

        self.iterative_integrand_evaluator.as_mut().map(|e| {
            e.compile(
                graph_path.join("iterative"),
                format!("{}_iterative", &self.graph.name,),
                graph_path.join("iterative"),
                settings.generation.gammaloop_compile_options.inline_asm(),
            )
        });
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
        settings: &RuntimeSettings,
        rotation: &Rotation,
    ) -> Complex<F<T>> {
        let hel = settings.kinematics.externals.get_helicities();
        let orientations =
            momentum_sample.orientations(&self.orientation_filter.0, &self.orientations);

        let evaluator = &self.orientation_parametric_integrand;
        let mut result = Complex::new_re(F(T::from_f64(0.)));

        {
            let a = T::get_parameters(
                &mut self.param_builder,
                &self.graph,
                momentum_sample,
                hel,
                None,
            );

            let iterative = self
                .iterative_integrand_evaluator
                .as_ref()
                .map(|ev| <T as GenericEvaluatorFloat>::get_evaluator(ev)(&a));

            for (i, e) in orientations.iter() {
                if let Some(iterative) = &iterative {
                    result += &iterative[i.0]
                } else {
                    self.param_builder.orientation_value(e);
                    let a = T::get_parameters(
                        &mut self.param_builder,
                        &self.graph,
                        momentum_sample,
                        hel,
                        None,
                    );
                    result += <T as GenericEvaluatorFloat>::get_evaluator_single(evaluator)(&a)
                }
            }
        }
        // status_debug!("last_params"; data = self.param_builder);
        debug!("evaluated integrand: {:16e}", result);

        let sum_of_cts = self.threshold_counterterm.evaluate(
            momentum_sample,
            &self.graph,
            &self.esurfaces,
            rotation,
            settings,
            &mut self.param_builder,
            orientations,
        );

        debug!("evaluated threshold counterterm: {:16e}", sum_of_cts);
        result - sum_of_cts
    }
}

impl GraphTerm for AmplitudeGraphTerm {
    type DerivedData = AmplitudeDerivedData;

    #[instrument(
          skip_all,
          fields(
              term.name = %self.name(),
          ),
          err
    )]
    fn warm_up(
        &mut self,
        derived_data: &AmplitudeDerivedData,
        settings: &RuntimeSettings,
    ) -> Result<()> {
        let a: BitVec = self
            .orientations
            .iter()
            .map(|a| settings.general.orientation_pat.filter(a))
            .collect();

        self.orientation_filter = BinVec(a);
        self.estimated_scale = Some(
            self.graph
                .underlying
                .expected_scale(settings.kinematics.e_cm),
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

        self.param_builder.add_external_four_mom(&externals);
        self.param_builder.polarizations_values(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );

        self.param_builder
            .m_uv_value(Complex::new_re(settings.general.m_uv));
        self.param_builder
            .mu_r_sq_value(Complex::new_re(settings.general.mu_r_sq));

        Ok(())
    }

    fn name(&self) -> String {
        self.graph.name.clone()
    }

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_multi_channeling_setup(&self) -> &LmbMultiChannelingSetup {
        &self.multi_channeling_setup
    }

    fn get_lmbs(&self) -> &TiVec<LmbIndex, LoopMomentumBasis> {
        &self.lmbs
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        settings: &RuntimeSettings,
        rotation: &Rotation,
    ) -> Complex<F<T>> {
        self.evaluate_impl(momentum_sample, settings, rotation)
    }

    fn get_num_orientations(&self) -> usize {
        self.orientations.len()
    }

    fn get_tropical_sampler(&self) -> &SampleGenerator<3> {
        &self
            .tropical_sampler
            .as_ref()
            .expect("Tropical sampler should be set.")
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

    pub graph_terms: Vec<AmplitudeGraphTerm>,
    pub external_signature: SignatureLike<ExternalIndex>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub group_derived_data: TiVec<GroupId, GroupDerivedData>,
}

impl AmplitudeIntegrand {
    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
    ) -> Result<()> {
        for a in &mut self.data.graph_terms {
            a.compile(path.as_ref(), override_existing, settings);
        }

        Ok(())
    }

    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self.data, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;

        // debug!("HE3");
        serde_yaml::to_writer(
            File::create(path.as_ref().join("settings.yaml"))?,
            &self.settings,
        )?;
        // debug!("HE");

        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let settings = serde_yaml::from_reader(File::open(path.as_ref().join("settings.yaml"))?)?;

        Ok(AmplitudeIntegrand { settings, data })
    }

    pub(crate) fn get_existing_esurfaces(&self) -> TiVec<GroupId, ExistingEsurfaces> {
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
                                        &graph.graph.loop_momentum_basis,
                                        &graph.graph.get_real_mass_vector(),
                                        &external_moms,
                                        &e_cm,
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

impl GammaloopIntegrand for AmplitudeIntegrand {
    type G = AmplitudeGraphTerm;

    #[instrument(
          skip_all,
          fields(
              integrand.name = %self.name(),
          )
    )]
    fn warm_up(&mut self, derived_data: &[&AmplitudeDerivedData]) -> Result<()> {
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

        for (a, derived_data) in self.data.graph_terms.iter_mut().zip(derived_data) {
            a.warm_up(derived_data, &self.settings)?;
        }

        let thresholds_generated = self
            .data
            .graph_terms
            .iter()
            .all(|term| !term.threshold_counterterm.evaluators.is_empty());

        if !thresholds_generated && !self.settings.subtraction.disable_threshold_subtraction {
            status_warn!("Not all graphs have threshold counterterms generated, but threshold subtraction is not disabled. disable runtime threshold subtraction to remove this warning");
            self.settings.subtraction.disable_threshold_subtraction = true;
        }

        if !self.settings.subtraction.disable_threshold_subtraction {
            let existing_esurfaces = self.get_existing_esurfaces();
            for (group_id, existing_esurfaces) in existing_esurfaces.iter_enumerated() {
                let graph_data = self.data.graph_group_structure[group_id]
                    .into_iter()
                    .map(|graph_id| {
                        let graph = &self.data.graph_terms[graph_id];
                        SingleGraphOverlapData {
                            lmb: &graph.graph.loop_momentum_basis,
                            esurfaces: &graph.esurfaces,
                            edge_masses: graph.graph.get_real_mass_vector(),
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
                    find_maximal_overlap(&overlap_input, existing_esurfaces, &external_moms)?;

                overlap.build_evaluators(
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
                        .model_parameters
                        .params
                        .clone(),
                )?;

                status_info!("overlap: {}", overlap);
            }

            todo!("store overlap,  and fix evaluation")
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

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor {
        DependentMomentaConstructor::Amplitude(&self.data.external_signature)
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
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: Complex<F<f64>>,
    ) -> EvaluationResult {
        evaluate_sample(self, sample, wgt, iter, use_f128, max_eval)
    }

    fn get_n_dim(&self) -> usize {
        if self
            .settings
            .sampling
            .get_parameterization_settings()
            .is_some()
        {
            self.data.graph_terms[0].graph.underlying.get_loop_number() * 3
        } else {
            let dimensions = self
                .data
                .graph_terms
                .iter()
                .map(|term| term.get_tropical_sampler().get_dimension())
                .sorted()
                .collect_vec();

            let median_dimension = dimensions[dimensions.len() / 2];
            median_dimension
        }
    }
}
