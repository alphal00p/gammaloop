use std::time::Duration;

use log::info;
use serde::Serialize;
use spenso::complex::Complex;
use symbolica::{
    evaluate::ExpressionEvaluator,
    numerical_integration::{ContinuousGrid, Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::esurface::Esurface,
    evaluation_result::{EvaluationMetaData, EvaluationResult},
    graph,
    integrands::HasIntegrand,
    momentum::ThreeMomentum,
    momentum_sample::{ExternalFourMomenta, LoopMomenta},
    new_cs::{CrossSectionGraph, CutId},
    new_graph::{FeynmanGraph, Graph},
    signature::ExternalSignature,
    utils::{self, FloatLike, F},
    IntegratedCounterTermRange, Settings,
};

#[derive(Clone)]
pub struct CrossSectionIntegrand {
    pub settings: Settings,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
}

#[derive(Clone)]
pub struct CrossSectionGraphTerm {
    pub bare_cff_evaluators: TiVec<CutId, ExpressionEvaluator<F<f64>>>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
}

impl HasIntegrand for CrossSectionIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        Grid::Continuous(ContinuousGrid::new(
            3,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
            self.settings.integrator.bin_number_evolution.clone(),
            self.settings.integrator.train_on_avg,
        ))
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
    ) -> EvaluationResult {
        assert_eq!(self.graph_terms.len(), 1);

        let graph_term = &self.graph_terms[0];
        assert_eq!(graph_term.bare_cff_evaluators.len(), 1);
        assert_eq!(graph_term.cut_esurface.len(), 1);
        let cut_id = CutId::from(0);

        let esurface = &graph_term.cut_esurface[cut_id];

        // println!("esurface: {:?}", esurface);

        let signature = ExternalSignature::from_iter(vec![1i8, 0]);
        let externals = self
            .settings
            .kinematics
            .externals
            .get_dependent_externals(&signature);

        // println!("externals: {:?}", &self.settings.kinematics.externals);

        let externals = ExternalFourMomenta::from_iter(externals);

        let raw_sample = match sample {
            Sample::Continuous(_, xs) => xs,
            _ => panic!("Expected continuous sample"),
        };

        let e_cm = self.settings.kinematics.e_cm;
        let e_cm_squared = e_cm * e_cm;

        let (loop_momenta, jacobian) =
            utils::global_parameterize(raw_sample, e_cm_squared, &self.settings, false);

        let real_loop_momenta = ThreeMomentum::from(loop_momenta[0]);
        let real_real_loop_momnta = LoopMomenta::from_iter([real_loop_momenta]);
        let loop_momenta = real_real_loop_momnta;

        let factor: F<f64> = F::from_f64((2. * std::f64::consts::PI).powi(-3));

        // println!("loop_momenta: {:?}", loop_momenta);
        let (tstar_initial, _) = esurface.get_radius_guess(
            &loop_momenta,
            &externals,
            &graph_term.graph.loop_momentum_basis,
        );

        // println!("tstar_initial: {}", tstar_initial);
        // println!("lmb: {:?}", graph_term.graph.loop_momentum_basis);

        let center = LoopMomenta::from_iter([ThreeMomentum::from([F(0.0); 3])]);

        // println!("am here");
        let function = |t: &_| {
            esurface.compute_self_and_r_derivative(
                t,
                &loop_momenta,
                &center,
                &externals,
                &graph_term.graph.loop_momentum_basis,
                &graph_term.graph.underlying.get_real_mass_vector(),
            )
        };

        let newton_result =
            newton_iteration_and_derivative(&tstar_initial, function, &F(1.0), 20, &e_cm);

        let t_star = newton_result.solution;
        let grad_eta = newton_result.derivative_at_solution;

        // println!("jacobian: {:?}", jacobian);
        // println!("externals: {:?}", externals);
        // println!("t_star: {}", t_star);

        //println!("grad_eta: {}", grad_eta);

        let h_function_settings = match &self.settings.subtraction.integrated_ct_settings.range {
            IntegratedCounterTermRange::Compact => panic!(),
            IntegratedCounterTermRange::Infinite {
                h_function_settings,
            } => h_function_settings,
        };

        let h = utils::h(t_star, None, None, h_function_settings);

        let rescaled_loop_momenta = loop_momenta.rescale(&t_star, None);

        let mut params = graph_term
            .graph
            .underlying
            .get_energy_cache(
                &rescaled_loop_momenta,
                &externals,
                &graph_term.graph.loop_momentum_basis,
            )
            .get_raw();

        params.push(t_star);
        params.push(h);
        params.push(grad_eta);

        let mut out = [F(0.0)];
        self.graph_terms[0].bare_cff_evaluators[cut_id].evaluate(&params, &mut out);

        let mut res = out[0];
        let is_nan = if res.is_nan() {
            res = F(0.0);
            res.is_nan()
        } else {
            res.is_nan()
        };

        EvaluationResult {
            integrand_result: Complex::new(res * jacobian * factor, F(0.0)),
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata: EvaluationMetaData {
                total_timing: Duration::ZERO,
                rep3d_evaluation_time: Duration::ZERO,
                parameterization_time: Duration::ZERO,
                relative_instability_error: Complex::new(F(0.0), F(0.0)),
                is_nan,
                highest_precision: crate::Precision::Double,
            },
        }
    }

    fn get_n_dim(&self) -> usize {
        3
    }
}

/// root finding, returns the derivative at the root, so that we don't have to recompute it.
/// Also returns the value of the function whose root is being found and the number of iterations used for debug information
fn newton_iteration_and_derivative<T: FloatLike>(
    guess: &F<T>,
    f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
    tolerance: &F<T>,
    max_iterations: usize,
    e_cm: &F<T>,
) -> NewtonIterationResult<T> {
    let mut x = guess.clone();
    let (mut val_f_x, mut val_df_x) = f_x_and_df_x(&x);

    let mut iteration = 0;

    while iteration < max_iterations && val_f_x.abs() > guess.epsilon() * tolerance * e_cm {
        x -= val_f_x / val_df_x;
        (val_f_x, val_df_x) = f_x_and_df_x(&x);
        iteration += 1;
    }

    NewtonIterationResult {
        solution: x,
        derivative_at_solution: val_df_x,
        error_of_function: val_f_x,
        num_iterations_used: iteration,
    }
}

#[derive(Serialize, Clone)]
struct NewtonIterationResult<T: FloatLike> {
    solution: F<T>,
    derivative_at_solution: F<T>,
    error_of_function: F<T>,
    num_iterations_used: usize,
}
