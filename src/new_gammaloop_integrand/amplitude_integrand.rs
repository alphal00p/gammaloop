use color_eyre::Result;
use itertools::Itertools;
use momtrop::SampleGenerator;
use serde::de::value;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::Atom,
    numerical_integration::{Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        esurface::{Esurface, EsurfaceID},
        expression::AmplitudeOrientationID,
    },
    evaluation_result::EvaluationResult,
    integrands::HasIntegrand,
    momentum::Rotation,
    momentum_sample::{ExternalIndex, MomentumSample},
    new_gammaloop_integrand::ParamBuilder,
    new_graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    signature::SignatureLike,
    subtraction::overlap::OverlapStructure,
    DependentMomentaConstructor, FloatLike, Polarizations, Settings, F,
};

use super::{
    create_grid, evaluate_sample, GammaloopIntegrand, GenericEvaluator, GenericEvaluatorFloat,
    GraphTerm, LmbMultiChannelingSetup,
};

const HARD_CODED_M_UV: F<f64> = F(1000.0);
const HARD_CODED_M_R_SQ: F<f64> = F(1000.0);

#[derive(Clone)]
pub struct AmplitudeGraphTerm {
    pub bare_cff_evaluator: GenericEvaluator,
    pub bare_cff_orientation_evaluators: TiVec<AmplitudeOrientationID, GenericEvaluator>,
    pub counterterm_evaluators: TiVec<EsurfaceID, GenericEvaluator>,
    pub graph: Graph,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: SampleGenerator<3>,
    pub estimated_scale: F<f64>,
    pub overlap: OverlapStructure,
}

impl AmplitudeGraphTerm {
    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
        mut param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        if let Some(forced_orientations) = &settings.general.force_orientations {
            if momentum_sample.sample.orientation.is_none() {
                return forced_orientations
                    .iter()
                    .map(|orientation_usize| {
                        let mut new_sample = momentum_sample.clone();
                        new_sample.sample.orientation = Some(*orientation_usize);
                        self.evaluate(&new_sample, settings, param_builder.clone())
                    })
                    .fold(
                        Complex::new_re(momentum_sample.zero()),
                        |sum, orientation_result| sum + orientation_result,
                    );
            }
        }

        if settings.general.debug > 0 {
            println!("loop_momenta: {:?}", momentum_sample.loop_moms());
            println!("external_momenta: {:?}", momentum_sample.external_moms());
        }

        param_builder.emr_spatial_value(
            self.graph
                .underlying
                .get_emr_vec_cache(
                    momentum_sample.loop_moms(),
                    momentum_sample.external_moms(),
                    &self.graph.loop_momentum_basis,
                )
                .into_iter()
                .map(|q| Complex::new_re(q))
                .collect(),
        );

        let params = param_builder.build_values_amplitude().unwrap();

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let orientation_id = AmplitudeOrientationID::from(orientation_id);
                let orientation_evaluator = &self.bare_cff_orientation_evaluators[orientation_id];
                <T as GenericEvaluatorFloat>::get_evaluator(orientation_evaluator)(&params)
            }
            None => {
                let evaluator = &self.bare_cff_evaluator;
                <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params)
            }
        };

        if settings.general.debug > 0 {
            println!("result of graph {}: {:16e}", self.graph.name, result);
        }

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
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
        param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        self.evaluate(momentum_sample, settings, param_builder)
    }

    fn get_num_orientations(&self) -> usize {
        self.bare_cff_orientation_evaluators.len()
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
    pub model_parameter_cache: Vec<Complex<F<f64>>>,
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

    fn get_model_parameter_cache<T: FloatLike>(&self) -> Vec<Complex<F<T>>> {
        self.model_parameter_cache
            .iter()
            .map(|x| Complex::new(F::from_ff64(x.re), F::from_ff64(x.im)))
            .collect()
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
        if self
            .settings
            .sampling
            .get_parameterization_settings()
            .is_some()
        {
            self.graph_terms[0].graph.underlying.get_loop_number() * 3
        } else {
            let dimensions = self
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
