use itertools::Itertools;
use momtrop::SampleGenerator;
use spenso::complex::Complex;
use symbolica::numerical_integration::{Grid, Sample};
use typed_index_collections::TiVec;

use crate::{
    cff::cut_expression::OrientationID,
    evaluation_result::EvaluationResult,
    integrands::HasIntegrand,
    momentum::Rotation,
    momentum_sample::{ExternalIndex, MomentumSample},
    new_graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    signature::SignatureLike,
    utils::{FloatLike, F},
    DependentMomentaConstructor, Polarizations, Settings,
};

use super::{
    create_grid, evaluate_sample, GammaloopIntegrand, GenericEvaluator, GenericEvaluatorFloat,
    GraphTerm, LmbMultiChannelingSetup,
};

#[derive(Clone)]
pub struct AmplitudeGraphTerm {
    pub bare_cff_evaluator: GenericEvaluator,
    pub bare_cff_orientation_evaluatos: TiVec<OrientationID, GenericEvaluator>,
    pub graph: Graph,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: SampleGenerator<3>,
}

impl AmplitudeGraphTerm {
    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        if let Some(forced_orientations) = &settings.general.force_orientations {
            if momentum_sample.sample.orientation.is_none() {
                return forced_orientations
                    .iter()
                    .map(|orientation_usize| {
                        let mut new_sample = momentum_sample.clone();
                        new_sample.sample.orientation =
                            Some(OrientationID::from(*orientation_usize));
                        self.evaluate(&new_sample, settings)
                    })
                    .fold(
                        Complex::new_re(momentum_sample.zero()),
                        |sum, orientation_result| sum + orientation_result,
                    );
            }
        }

        let params = self
            .graph
            .underlying
            .get_energy_cache(
                momentum_sample.loop_moms(),
                momentum_sample.external_moms(),
                &self.graph.loop_momentum_basis,
            )
            .into_iter()
            .map(|x| x.1)
            .collect_vec();

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let orientation_evaluator = &self.bare_cff_orientation_evaluatos[orientation_id];
                <T as GenericEvaluatorFloat>::get_evaluator(orientation_evaluator)(&params)
            }
            None => {
                let evaluator = &self.bare_cff_evaluator;
                <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params)
            }
        };

        Complex::new_re(result)
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
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.evaluate(momentum_sample, settings)
    }

    fn get_num_orientations(&self) -> usize {
        self.bare_cff_orientation_evaluatos.len()
    }

    fn get_tropical_sampler(&self) -> &SampleGenerator<3> {
        &self.tropical_sampler
    }
}

#[derive(Clone)]
pub struct AmplitudeIntegrand {
    pub settings: Settings,
    pub rotations: Vec<Rotation>,
    pub polarizations: Vec<Polarizations>,
    pub graph_terms: Vec<AmplitudeGraphTerm>,
    pub external_signature: SignatureLike<ExternalIndex>,
}

impl GammaloopIntegrand for AmplitudeIntegrand {
    type G = AmplitudeGraphTerm;

    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.rotations.iter()
    }

    fn get_terms(&self) -> impl Iterator<Item = &Self::G> {
        self.graph_terms.iter()
    }

    fn get_settings(&self) -> &Settings {
        &self.settings
    }

    fn get_graph(&self, graph_id: usize) -> &Self::G {
        &self.graph_terms[graph_id]
    }

    fn get_polarizations(&self) -> &[Polarizations] {
        &self.polarizations
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor {
        DependentMomentaConstructor::Amplitude(&self.external_signature)
    }
}

impl HasIntegrand for AmplitudeIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        create_grid(self)
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
        todo!()
    }
}
