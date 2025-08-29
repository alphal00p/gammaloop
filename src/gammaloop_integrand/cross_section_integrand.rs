use bincode::Encode;
use bincode_trait_derive::Decode;
use color_eyre::Result;
use colored::Colorize;
use itertools::Itertools;
use log::{debug, info};
use spenso::algebra::complex::Complex;
use std::{
    fs::{self, File},
    io::{Read, Write},
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
    gammaloop_integrand::ParamBuilder,
    graph::{ExternalConnection, FeynmanGraph, Graph, GroupId, LmbIndex, LoopMomentumBasis},
    integrands::HasIntegrand,
    momentum::{Rotation, RotationMethod, ThreeMomentum},
    momentum_sample::{LoopMomenta, MomentumSample},
    processes::{CrossSectionCut, CrossSectionDerivedData, CutId},
    settings::{runtime::IntegratedCounterTermRange, RuntimeSettings},
    utils::{self, newton_solver::newton_iteration_and_derivative, FloatLike, F},
    DependentMomentaConstructor, GammaLoopContext, GammaLoopContextContainer,
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
    pub settings: RuntimeSettings,
    pub data: CrossSectionIntegrandData,
}
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionIntegrandData {
    pub name: String,
    // pub polarizations: Vec<Polarizations>,
    pub rotations: Option<Vec<Rotation>>,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
    pub n_incoming: usize,
    pub external_connections: Vec<ExternalConnection>,
    // pub builder_cache: ParamBuilder<f64>,
}

impl CrossSectionIntegrand {
    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;
        File::create(path.as_ref().join("settings.toml"))?
            .write(toml::to_string_pretty(&self.settings)?.as_bytes())?;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;
        let mut buf = vec![];
        File::open(path.as_ref().join("settings.toml"))?.read(&mut buf)?;
        let settings = toml::from_slice(&buf)?;

        Ok(CrossSectionIntegrand { settings, data })
    }
}

impl GammaloopIntegrand for CrossSectionIntegrand {
    type G = CrossSectionGraphTerm;
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.as_ref().expect("forgot warmup").iter()
    }

    fn warm_up(&mut self, derived_data: &[&CrossSectionDerivedData]) -> Result<()> {
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
        Ok(())
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
    }

    fn get_group_masters(&self) -> impl Iterator<Item = &Self::G> {
        todo!();
        std::iter::empty()
    }

    fn get_settings(&self) -> &RuntimeSettings {
        &self.settings
    }

    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G {
        &mut self.data.graph_terms[graph_id]
    }

    fn get_master_graph(&self, group_id: GroupId) -> &Self::G {
        todo!()
    }

    fn get_group(&self, group_id: GroupId) -> &crate::graph::GraphGroup {
        todo!()
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor {
        DependentMomentaConstructor::CrossSection {
            external_connections: &self.data.external_connections,
        }
    }

    // fn get_builder_cache(&self) -> &ParamBuilder<f64> {
    //     &self.data.builder_cache
    // }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
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
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder<f64>,
}

impl GraphTerm for CrossSectionGraphTerm {
    type DerivedData = CrossSectionDerivedData;

    fn warm_up(
        &mut self,
        derived_data: &CrossSectionDerivedData,
        settings: &RuntimeSettings,
    ) -> Result<()> {
        self.estimated_scale = Some(
            self.graph
                .underlying
                .expected_scale(settings.kinematics.e_cm),
        );

        Ok(())
    }
    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        settings: &RuntimeSettings,
        rotation: &Rotation,
    ) -> Complex<F<T>> {
        todo!()
        // self.evaluate_impl(momentum_sample, settings, rotation)
    }
    fn name(&self) -> String {
        self.graph.name.clone()
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

    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        settings: &RuntimeSettings,
        param_builder: ParamBuilder<T>,
    ) -> Complex<F<T>> {
        // implementation of forced orientations, only works with sample orientation disabled

        debug!("loop_momenta: {:?}", momentum_sample.loop_moms());
        debug!("external_momenta: {:?}", momentum_sample.external_moms());

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
                    IntegratedCounterTermRange::Compact {} => panic!(),
                    IntegratedCounterTermRange::Infinite {
                        h_function_settings,
                    } => h_function_settings,
                };

                let h_function = utils::h(&newton_result.solution, None, None, h_function_settings);

                debug!(
                    "generated parameters for cut: {}",
                    format!("{}", _cut_id).green()
                );

                debug!("solver: {:#?}", newton_result);
                debug!("h_function: {}", format!("{:16e}", h_function).green());
                debug!("rescaled loop momenta: {:?}", rescaled_sample.loop_moms());

                let mut cut_param_builder = param_builder.clone();

                // cut_param_builder.emr_spatial_value(
                //     self.graph
                //         .underlying
                //         .get_emr_vec_cache(
                //             rescaled_sample.loop_moms(),
                //             rescaled_sample.external_moms(),
                //             &self.graph.loop_momentum_basis,
                //         )
                //         .into_iter()
                //         .map(|q| Complex::new_re(q))
                //         .collect(),
                // );

                // cut_param_builder.tstar_value(Complex::new_re(newton_result.solution));
                // cut_param_builder.h_function_value(Complex::new_re(h_function));
                // cut_param_builder.derivative_at_tstar_value(Complex::new_re(
                //     newton_result.derivative_at_solution,
                // ));

                let params: Vec<_> = cut_param_builder.values;
                params
            });

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let orientation_id = SuperGraphOrientationID::from(orientation_id);
                todo!();
                // let orientation_evaluator =
                //     &mut self.bare_cff_orientation_evaluators[orientation_id];
                // orientation_evaluator
                //     .evaluators
                //     .iter_mut()
                //     .zip(params)
                //     .map(|(evaluator, params)| {
                //         let cut_results =
                //             <T as GenericEvaluatorFloat>::get_evaluator_single(evaluator)(&params);
                //         cut_results
                //     })
                //     .fold(
                //         Complex::new_re(momentum_sample.zero()),
                //         |sum, cut_result| sum + cut_result,
                //     )
                Complex::new_re(momentum_sample.zero())
            }
            None => {
                // self
                // .bare_cff_evaluators
                // .iter_mut()
                // .zip_eq(params)
                // .enumerate()
                // .map(|(id, (evaluator, params))| {
                //     let cut_results =
                //         <T as GenericEvaluatorFloat>::get_evaluator_single(evaluator)(&params);
                //     info!(
                //         "cut: {}, result: {}",
                //         format!("{}", id).green(),
                //         format!("{:16e}", cut_results).blue()
                //     );

                //     cut_results
                // })
                // .fold(
                //     Complex::new_re(momentum_sample.zero()),
                //     |sum, cut_result| sum + cut_result,
                // ),
                Complex::new_re(momentum_sample.zero())
            }
        };

        let final_result = result;
        info!(
            "sum of all cuts: {}",
            format!("{:16e}", final_result).blue()
        );

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
        info!("result: {:?}", result.integrand_result);

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
