use bincode::Encode;
use bincode_trait_derive::Decode;
use color_eyre::Result;
use colored::Colorize;
use itertools::Itertools;
use serde::Serialize;
use spenso::algebra::complex::Complex;
use std::{
    fs::{self, File},
    path::Path,
};
use symbolica::numerical_integration::{Grid, Sample};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CutOrientationData, SuperGraphOrientationID},
        esurface::Esurface,
    },
    evaluation_result::EvaluationResult,
    integrands::HasIntegrand,
    momentum::{Rotation, ThreeMomentum},
    momentum_sample::{LoopMomenta, MomentumSample},
    new_cs::{CrossSectionCut, CutId},
    new_gammaloop_integrand::ParamBuilder,
    new_graph::{ExternalConnection, FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    utils::{self, FloatLike, F},
    DependentMomentaConstructor, GammaLoopContext, GammaLoopContextContainer,
    IntegratedCounterTermRange, Polarizations, Settings,
};

use super::{
    create_grid, evaluate_sample, GammaloopIntegrand, GenericEvaluator, GenericEvaluatorFloat,
    GraphTerm, LmbMultiChannelingSetup,
};

const TOLERANCE: F<f64> = F(2.0);
const HARD_CODED_M_UV: F<f64> = F(1000.0);
const HARD_CODED_M_R_SQ: F<f64> = F(1000.0);

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionIntegrand {
    pub settings: Settings,
    pub data: CrossSectionIntegrandData,
}
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionIntegrandData {
    pub name: String,
    pub polarizations: Vec<Polarizations>,
    pub rotations: Vec<Rotation>,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
    pub n_incoming: usize,
    pub external_connections: Vec<ExternalConnection>,
    pub model_parameter_cache: Vec<Complex<F<f64>>>,
}

impl CrossSectionIntegrand {
    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;
        serde_yaml::to_writer(
            File::create(path.as_ref().join("settings.yaml"))?,
            &self.settings,
        )?;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;
        let settings = serde_yaml::from_reader(File::open(path.as_ref().join("settings.yaml"))?)?;

        Ok(CrossSectionIntegrand { settings, data })
    }
}

impl GammaloopIntegrand for CrossSectionIntegrand {
    type G = CrossSectionGraphTerm;
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.iter()
    }

    fn get_terms(&self) -> impl Iterator<Item = &Self::G> {
        self.data.graph_terms.iter()
    }

    fn get_settings(&self) -> &Settings {
        &self.settings
    }

    fn get_graph(&self, graph_id: usize) -> &Self::G {
        &self.data.graph_terms[graph_id]
    }

    fn get_polarizations(&self) -> &[Polarizations] {
        &self.data.polarizations
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor {
        DependentMomentaConstructor::CrossSection {
            external_connections: &self.data.external_connections,
        }
    }

    fn get_model_parameter_cache<T: FloatLike>(&self) -> Vec<Complex<F<T>>> {
        self.data
            .model_parameter_cache
            .iter()
            .map(|x| Complex::new(F::from_ff64(x.re), F::from_ff64(x.im)))
            .collect()
    }
}

#[derive(Clone, Encode, Decode)]
pub struct OrientationEvaluator {
    pub orientation_data: CutOrientationData,
    pub evaluators: Vec<GenericEvaluator>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionGraphTerm {
    pub bare_cff_evaluators: TiVec<CutId, GenericEvaluator>,
    pub bare_cff_orientation_evaluators: TiVec<SuperGraphOrientationID, OrientationEvaluator>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub estimated_scale: F<f64>,
}

impl GraphTerm for CrossSectionGraphTerm {
    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
        param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        self.evaluate(momentum_sample, settings, param_builder)
    }

    fn get_multi_channeling_setup(&self) -> &LmbMultiChannelingSetup {
        &self.multi_channeling_setup
    }

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_lmbs(&self) -> &TiVec<LmbIndex, LoopMomentumBasis> {
        &self.lmbs
    }

    fn get_num_orientations(&self) -> usize {
        self.bare_cff_orientation_evaluators.len()
    }

    fn get_tropical_sampler(&self) -> &momtrop::SampleGenerator<3> {
        unimplemented!(
            "Don't know how to generate subgraph table for forward scattering graphs yet"
        )
    }
}

impl CrossSectionGraphTerm {
    fn self_get_cuts_to_evaluate(
        &self,
        orientation: Option<SuperGraphOrientationID>,
    ) -> Vec<(CutId, &Esurface)> {
        if let Some(orientation_id) = orientation {
            self.bare_cff_orientation_evaluators[orientation_id]
                .orientation_data
                .cuts
                .iter()
                .map(|cut_id| (*cut_id, &self.cut_esurface[*cut_id]))
                .collect()
        } else {
            self.cut_esurface.iter_enumerated().collect()
        }
    }

    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
        param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        // implementation of forced orientations, only works with sample orientation disabled
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

        let center = LoopMomenta::from_iter(vec![
            ThreeMomentum::from([
                momentum_sample.zero(),
                momentum_sample.zero(),
                momentum_sample.zero()
            ]);
            self.graph.underlying.get_loop_number()
        ]);

        let params = self
            .self_get_cuts_to_evaluate(
                momentum_sample
                    .sample
                    .orientation
                    .map(SuperGraphOrientationID::from),
            )
            .into_iter()
            .map(|(_cut_id, esurface)| {
                let (tstar_initial, _tstar_initial_negative) = esurface.get_radius_guess(
                    momentum_sample.loop_moms(),
                    momentum_sample.external_moms(),
                    &self.graph.loop_momentum_basis,
                );

                let root_function = |t: &_| {
                    esurface.compute_self_and_r_derivative(
                        t,
                        momentum_sample.loop_moms(),
                        &center,
                        momentum_sample.external_moms(),
                        &self.graph.loop_momentum_basis,
                        &self.graph.underlying.get_real_mass_vector(),
                    )
                };

                let e_cm = F::from_ff64(settings.kinematics.e_cm);
                let tolerance = F::from_ff64(TOLERANCE);

                let newton_result = newton_iteration_and_derivative(
                    &tstar_initial,
                    root_function,
                    &tolerance,
                    20,
                    &e_cm,
                );

                let rescaled_sample =
                    momentum_sample.rescaled_loop_momenta(&newton_result.solution, None);

                let h_function_settings = match &settings.subtraction.integrated_ct_settings.range {
                    IntegratedCounterTermRange::Compact => panic!(),
                    IntegratedCounterTermRange::Infinite {
                        h_function_settings,
                    } => h_function_settings,
                };

                let h_function = utils::h(&newton_result.solution, None, None, h_function_settings);

                if settings.general.debug > 1 && _cut_id == CutId::from(2) {
                    println!(
                        "generated parameters for cut: {}",
                        format!("{}", _cut_id).green()
                    );

                    println!("solver: {:#?}", newton_result);
                    println!("h_function: {}", format!("{:16e}", h_function).green());
                    println!("rescaled loop momenta: {:?}", rescaled_sample.loop_moms());
                }

                let mut cut_param_builder = param_builder.clone();

                cut_param_builder.emr_spatial_value(
                    self.graph
                        .underlying
                        .get_emr_vec_cache(
                            rescaled_sample.loop_moms(),
                            rescaled_sample.external_moms(),
                            &self.graph.loop_momentum_basis,
                        )
                        .into_iter()
                        .map(|q| Complex::new_re(q))
                        .collect(),
                );

                cut_param_builder.tstar_value(Complex::new_re(newton_result.solution));
                cut_param_builder.h_function_value(Complex::new_re(h_function));
                cut_param_builder.derivative_at_tstar_value(Complex::new_re(
                    newton_result.derivative_at_solution,
                ));

                let params = cut_param_builder.build_values_cs().unwrap();
                params
            });

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let orientation_id = SuperGraphOrientationID::from(orientation_id);
                let orientation_evaluator = &self.bare_cff_orientation_evaluators[orientation_id];
                orientation_evaluator
                    .evaluators
                    .iter()
                    .zip(params)
                    .map(|(evaluator, params)| {
                        let cut_results =
                            <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params);
                        cut_results
                    })
                    .fold(
                        Complex::new_re(momentum_sample.zero()),
                        |sum, cut_result| sum + cut_result,
                    )
            }
            None => self
                .bare_cff_evaluators
                .iter()
                .zip_eq(params)
                .enumerate()
                .map(|(id, (evaluator, params))| {
                    let cut_results =
                        <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params);
                    if settings.general.debug > 0 && id == 2 {
                        println!(
                            "cut: {}, result: {}",
                            format!("{}", id).green(),
                            format!("{:16e}", cut_results).blue()
                        );
                    }
                    cut_results
                })
                .fold(
                    Complex::new_re(momentum_sample.zero()),
                    |sum, cut_result| sum + cut_result,
                ),
        };

        let final_result = result;
        if settings.general.debug > 0 {
            println!(
                "sum of all cuts: {}",
                format!("{:16e}", final_result).blue()
            );
        }
        final_result
    }
}

impl HasIntegrand for CrossSectionIntegrand {
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
        _iter: usize,
        use_f128: bool,
        max_eval: Complex<F<f64>>,
    ) -> EvaluationResult {
        let result = evaluate_sample(self, sample, wgt, _iter, use_f128, max_eval);
        if self.settings.general.debug > 0 {
            println!("result: {:?}", result.integrand_result);
        }

        result
    }

    fn get_n_dim(&self) -> usize {
        assert!(
            self.settings
                .sampling
                .get_parameterization_settings()
                .is_some(),
            "Tropical smapling not implemented for cross sections yet"
        );

        assert!(self
            .data
            .graph_terms
            .iter()
            .map(|term| term.graph.underlying.get_loop_number())
            .all_equal());

        self.data.graph_terms[0].graph.underlying.get_loop_number() * 3
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

#[derive(Serialize, Clone, Debug)]
struct NewtonIterationResult<T: FloatLike> {
    solution: F<T>,
    derivative_at_solution: F<T>,
    error_of_function: F<T>,
    num_iterations_used: usize,
}
