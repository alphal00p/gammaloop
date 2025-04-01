use symbolica::numerical_integration::{Grid, Sample};

use crate::{evaluation_result::EvaluationResult, integrands::HasIntegrand, utils::F, Settings};

#[derive(Clone)]
pub struct AmplitudeIntegrand {
    pub settings: Settings,
}

impl HasIntegrand for AmplitudeIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        todo!()
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
    ) -> EvaluationResult {
        todo!()
    }

    fn get_n_dim(&self) -> usize {
        todo!()
    }
}
