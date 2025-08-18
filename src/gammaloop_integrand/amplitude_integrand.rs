use std::{
    fs::{self, File},
    path::Path,
};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use log::debug;
use momtrop::SampleGenerator;
use spenso::algebra::complex::Complex;
use symbolica::{
    evaluate::OptimizationSettings,
    numerical_integration::{Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        esurface::{get_existing_esurfaces, EsurfaceCollection},
        expression::AmplitudeOrientationID,
    },
    evaluation_result::EvaluationResult,
    gammaloop_integrand::ParamBuilder,
    graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    integrands::HasIntegrand,
    model::Model,
    momentum::Rotation,
    momentum_sample::{ExternalIndex, MomentumSample},
    processes::AmplitudeGraph,
    signature::SignatureLike,
    subtraction::{amplitude_counterterm::AmplitudeCountertermData, overlap::find_maximal_overlap},
    DependentMomentaConstructor, FloatLike, GammaLoopContext, GammaLoopContextContainer,
    RuntimeSettings, F,
};

use super::{
    create_grid, evaluate_sample, GammaloopIntegrand, GenericEvaluator, GenericEvaluatorFloat,
    GraphTerm, LmbMultiChannelingSetup,
};

const HARD_CODED_M_UV: F<f64> = F(1000.0);
const HARD_CODED_M_R_SQ: F<f64> = F(1000.0);

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeGraphTerm {
    pub orientation_parametric_integrand: GenericEvaluator,
    pub orientations: TiVec<AmplitudeOrientationID, EdgeVec<Orientation>>,
    pub esurfaces: EsurfaceCollection,
    pub threshold_counterterm: AmplitudeCountertermData,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub graph: Graph,
    pub estimated_scale: F<f64>,
    pub param_builder: ParamBuilder,
}

/// Num(sigma_1,sigma_2,...)*(CFF_1 delta(edge(1),1) delta_(1,1,1,-1,1)+CFF_3 delta_(1,1,1,-1,1)+CFF_2 delta_(1,1,1,-1,1))

impl AmplitudeGraphTerm {
    pub fn from_amplitude_graph(
        graph: &AmplitudeGraph,
        settings: &RuntimeSettings,
        model: &Model,
    ) -> Self {
        let estimated_scale = graph
            .graph
            .underlying
            .expected_scale(settings.kinematics.e_cm);

        let esurfaces = &graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache;

        let esurface_data = graph.derived_data.esurface_data.as_ref().unwrap();
        let externals = settings.kinematics.externals.get_dependent_externals(
            DependentMomentaConstructor::Amplitude(&graph.graph.get_external_signature()),
        );

        let existing_esurfaces = get_existing_esurfaces(
            esurfaces,
            esurface_data,
            &externals,
            &graph.graph.loop_momentum_basis,
            settings.general.debug,
            settings.kinematics.e_cm,
        );

        let mut param_builder = graph.param_builder_core(model);
        param_builder.add_external_four_mom(&externals);
        param_builder.polarizations_values(
            &graph.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );
        param_builder.model_parameters_value(model);
        param_builder.mu_r_sq_value(Complex::new_zero());
        param_builder.m_uv_value(Complex::new_zero());

        // (momentum_sample);

        let edge_masses = graph.graph.new_edgevec(|edge, _, _| edge.mass_value());

        let overlap = find_maximal_overlap(
            &graph.graph.loop_momentum_basis,
            &existing_esurfaces,
            esurfaces,
            &edge_masses,
            &externals,
            &settings,
        );

        let mut threshold_counterterm = graph.derived_data.threshold_counterterm.clone();

        threshold_counterterm.overlap = overlap;

        threshold_counterterm
            .param_builder
            .add_external_four_mom(&externals);
        threshold_counterterm.param_builder.polarizations_values(
            &graph.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );
        threshold_counterterm
            .param_builder
            .model_parameters_value(model);
        threshold_counterterm
            .param_builder
            .mu_r_sq_value(Complex::new_zero());
        threshold_counterterm
            .param_builder
            .m_uv_value(Complex::new_zero());

        AmplitudeGraphTerm {
            orientations: graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations
                .iter()
                .map(|a| a.data.orientation.clone())
                .collect(),
            orientation_parametric_integrand: GenericEvaluator::new_from_builder(
                [graph.derived_data.all_mighty_integrand.clone()],
                &param_builder,
                OptimizationSettings::default(),
            )
            .unwrap(),
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
            estimated_scale,
            esurfaces: graph
                .derived_data
                .cff_expression
                .as_ref()
                .expect("cff_expression should have been created")
                .surfaces
                .esurface_cache
                .clone(),

            param_builder,
        }
    }

    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        settings: &RuntimeSettings,
        rotation: &Rotation,
    ) -> Complex<F<T>> {
        if let Some(forced_orientations) = &settings.general.force_orientations {
            if momentum_sample.sample.orientation.is_none() {
                return forced_orientations
                    .iter()
                    .map(|orientation_usize| {
                        let mut new_sample = momentum_sample.clone();
                        new_sample.sample.orientation = Some(*orientation_usize);
                        self.evaluate(&new_sample, settings, rotation)
                    })
                    .fold(
                        Complex::new_re(momentum_sample.zero()),
                        |sum, orientation_result| sum + orientation_result,
                    );
            }
        }

        if settings.general.debug > 4 {
            println!("Evaluating graph: {}", self.graph.name);
        }
        let hel = settings.kinematics.externals.get_helicities();

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let orientation = &self.orientations[AmplitudeOrientationID::from(orientation_id)];
                self.param_builder.orientation_value(orientation);
                let a = T::get_parameters(
                    &mut self.param_builder,
                    &self.graph,
                    momentum_sample,
                    hel,
                    None,
                );
                let evaluator = &self.orientation_parametric_integrand;
                <T as GenericEvaluatorFloat>::get_evaluator_single(evaluator)(&a)
            }
            None => {
                let evaluator = &self.orientation_parametric_integrand;
                let mut res = Complex::new_re(F(T::from_f64(0.)));

                for e in &self.orientations {
                    self.param_builder.orientation_value(e);
                    let a = T::get_parameters(
                        &mut self.param_builder,
                        &self.graph,
                        momentum_sample,
                        hel,
                        None,
                    );
                    res += <T as GenericEvaluatorFloat>::get_evaluator_single(evaluator)(&a)
                }
                res
            }
        };
        debug!("evaluated integrand: {:16e}", result);

        if !settings.general.disable_threshold_subtraction {
            let sum_of_cts = self.threshold_counterterm.evaluate(
                momentum_sample,
                &self.graph,
                &self.esurfaces,
                rotation,
                settings,
                // param_builder,
            );
            debug!("evaluated threshold counterterm: {:16e}", sum_of_cts);
            result - sum_of_cts
        } else {
            result
        }
    }
}

impl GraphTerm for AmplitudeGraphTerm {
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
    pub rotations: Vec<Rotation>,
    pub name: String,

    pub graph_terms: Vec<AmplitudeGraphTerm>,
    pub external_signature: SignatureLike<ExternalIndex>,
    // pub builder_cache: ParamBuilder<f64>,
    // pub model_parameter_cache: Vec<Complex<F<f64>>>,
}

impl AmplitudeIntegrand {
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
}

impl GammaloopIntegrand for AmplitudeIntegrand {
    type G = AmplitudeGraphTerm;

    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.iter()
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
    }

    fn get_graph(&self, graph_id: usize) -> &Self::G {
        &self.data.graph_terms[graph_id]
    }

    fn get_terms(&self) -> impl Iterator<Item = &Self::G> {
        self.data.graph_terms.iter()
    }

    fn get_settings(&self) -> &RuntimeSettings {
        &self.settings
    }

    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G {
        &mut self.data.graph_terms[graph_id]
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
