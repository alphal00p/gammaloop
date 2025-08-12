use std::{
    fs::{self, File},
    path::Path,
};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::Context;
use itertools::Itertools;
use log::debug;
use momtrop::SampleGenerator;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::AtomCore,
    numerical_integration::{Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::{esurface::EsurfaceID, expression::AmplitudeOrientationID},
    evaluation_result::EvaluationResult,
    integrands::HasIntegrand,
    momentum::Rotation,
    momentum_sample::{ExternalIndex, MomentumSample},
    new_gammaloop_integrand::ParamBuilder,
    new_graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    signature::SignatureLike,
    subtraction::overlap::OverlapStructure,
    DependentMomentaConstructor, FloatLike, GammaLoopContext, GammaLoopContextContainer,
    Polarizations, Settings, F,
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
    pub integrand_evaluator_all_orientations: GenericEvaluator,
    pub integrand_evaluators: TiVec<AmplitudeOrientationID, GenericEvaluator>,
    pub counterterm_evaluators: TiVec<EsurfaceID, GenericEvaluator>,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: SampleGenerator<3>,
    pub graph: Graph,
    pub estimated_scale: F<f64>,
    pub overlap: OverlapStructure,
    pub param_builder: ParamBuilder<f64>,
}

impl AmplitudeGraphTerm {
    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
        // mut param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        if let Some(forced_orientations) = &settings.general.force_orientations {
            if momentum_sample.sample.orientation.is_none() {
                return forced_orientations
                    .iter()
                    .map(|orientation_usize| {
                        let mut new_sample = momentum_sample.clone();
                        new_sample.sample.orientation = Some(*orientation_usize);
                        self.evaluate(&new_sample, settings)
                    })
                    .fold(
                        Complex::new_re(momentum_sample.zero()),
                        |sum, orientation_result| sum + orientation_result,
                    );
            }
        }

        debug!("Sample: \n{:>40}", momentum_sample);
        let hel = settings.kinematics.externals.get_helicities();

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let a =
                    T::get_parameters(&mut self.param_builder, &self.graph, momentum_sample, hel);
                let orientation_id = AmplitudeOrientationID::from(orientation_id);
                let orientation_evaluator = &self.integrand_evaluators[orientation_id];
                <T as GenericEvaluatorFloat>::get_evaluator(orientation_evaluator)(&a)
            }
            None => {
                let evaluator = &self.integrand_evaluator_all_orientations;

                // let replaced = self.param_builder.replace_non_emr(evaluator.expr.clone());
                // println!("replaced: {:+>}", replaced.collect_num());
                // evaluator.validate(&param_builder);
                let a =
                    T::get_parameters(&mut self.param_builder, &self.graph, momentum_sample, hel);
                // self.param_builder.validate();

                <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&a)
            }
        };

        // if settings.general.debug > 0 {
        debug!("result of graph {}: {:16e}", self.graph.name, result);
        // }

        result
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
        settings: &Settings,
        // param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        self.evaluate_impl(momentum_sample, settings)
    }

    fn get_num_orientations(&self) -> usize {
        self.integrand_evaluators.len()
    }

    fn get_tropical_sampler(&self) -> &SampleGenerator<3> {
        &self.tropical_sampler
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeIntegrand {
    pub settings: Settings,
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

    fn get_settings(&self) -> &Settings {
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
