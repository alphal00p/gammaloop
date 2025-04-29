use spenso::complex::Complex;
use symbolica::numerical_integration::{Grid, Sample};
use typed_index_collections::TiVec;

use crate::{
    cff::cut_expression::OrientationID,
    evaluation_result::EvaluationResult,
    integrands::HasIntegrand,
    momentum::{Polarization, Rotation},
    momentum_sample::MomentumSample,
    new_graph::{Graph, LmbIndex, LoopMomentumBasis},
    utils::{FloatLike, F},
    Polarizations, Settings,
};

use super::{
    create_grid, GammaloopIntegrand, GenericEvaluator, GraphTerm, LmbMultiChannelingSetup,
};

#[derive(Clone)]
pub struct AmplitudeGraphTerm {
    bare_cff_evaluator: GenericEvaluator,
    bare_cff_orientation_evaluatos: TiVec<OrientationID, GenericEvaluator>,
    graph: Graph,
    multi_channeling_setup: LmbMultiChannelingSetup,
    lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
}

impl AmplitudeGraphTerm {
    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        todo!()
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
        todo!()
    }
}

#[derive(Clone)]
pub struct AmplitudeIntegrand {
    pub settings: Settings,
    pub rotations: Vec<Rotation>,
    pub polarizations: Vec<Polarizations>,
    pub graph_terms: Vec<AmplitudeGraphTerm>,
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
        todo!()
    }

    fn get_n_dim(&self) -> usize {
        todo!()
    }
}
