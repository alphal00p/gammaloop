#![allow(dead_code)]

use bincode::Encode;
use bincode_trait_derive::Decode;
use bitvec::vec::BitVec;
use color_eyre::Result;
use colored::Colorize;
use eyre::Context;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use log::{debug, info};
use spenso::algebra::complex::Complex;
use std::{
    fs::{self, File},
    io::{Read, Write},
    path::Path,
};
use symbolica::{
    evaluate::OptimizationSettings,
    numerical_integration::{Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CutOrientationData, SuperGraphOrientationID},
        esurface::Esurface,
        expression::GraphOrientation,
    },
    evaluation_result::EvaluationResult,
    gammaloop_integrand::{
        param_builder::LUParams, GenericEvaluatorFloat, ParamBuilder, UpdateAndGetParams,
    },
    graph::{
        ExternalConnection, FeynmanGraph, Graph, GraphGroup, GroupId, LmbIndex, LoopMomentumBasis,
    },
    integrands::HasIntegrand,
    model::Model,
    momentum::{Rotation, RotationMethod, ThreeMomentum},
    momentum_sample::{LoopMomenta, MomentumSample},
    processes::{Amplitude, CrossSectionCut, CrossSectionGraph, CutId, GroupDerivedData},
    settings::{
        runtime::{HFunctionSettings, IntegratedCounterTermRange},
        GlobalSettings, RuntimeSettings,
    },
    utils::{
        self, bitvec_ext::BinVec, h, newton_solver::newton_iteration_and_derivative,
        serde_utils::SmartSerde, FloatLike, Length, F,
    },
    DependentMomentaConstructor, GammaLoopContext, GammaLoopContextContainer,
};

use super::{
    create_grid, evaluate_sample, GammaloopIntegrand, GenericEvaluator, GraphTerm,
    LmbMultiChannelingSetup,
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
    pub loop_cache_id: usize,
    pub external_cache_id: usize,
    /// Cache ID for the base (unrotated) external momentum configuration
    pub base_external_cache_id: usize,
    // pub polarizations: Vec<Polarizations>,
    pub rotations: Option<Vec<Rotation>>,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
    pub n_incoming: usize,
    pub external_connections: Vec<ExternalConnection>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    // pub builder_cache: ParamBuilder<f64>,
}

impl CrossSectionIntegrand {
    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self.data, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;

        self.settings
            .to_file(path.as_ref().join("settings.toml"), override_existing)
            .with_context(|| "Error saving settings.toml file for amplitude integrand")?;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _): (CrossSectionIntegrandData, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let settings = SmartSerde::from_file(
            path.as_ref().join("settings.toml"),
            "runtime settings for amplitude integrand",
        )?;

        Ok(CrossSectionIntegrand { settings, data })
    }
}

impl GammaloopIntegrand for CrossSectionIntegrand {
    type G = CrossSectionGraphTerm;

    fn external_cache_id(&self) -> usize {
        self.data.external_cache_id
    }

    fn increment_external_cache_id(&mut self, val: usize) {
        self.data.external_cache_id += val
    }

    fn signal_external_momenta_changed(&mut self) {
        self.increment_external_cache_id(1);
        // Update base cache ID when the fundamental configuration changes
        self.data.base_external_cache_id = self.data.external_cache_id;
    }

    fn get_current_external_cache_id(&self) -> usize {
        self.external_cache_id()
    }

    /// Revert to the base external cache ID for the current configuration
    fn revert_to_base_external_cache_id(&mut self) {
        self.data.external_cache_id = self.data.base_external_cache_id;
    }

    fn get_base_external_cache_id(&self) -> usize {
        self.data.base_external_cache_id
    }

    fn increment_loop_cache_id(&mut self, val: usize) {
        self.data.loop_cache_id += val
    }

    fn loop_cache_id(&self) -> usize {
        self.data.loop_cache_id
    }
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.as_ref().expect("forgot warmup").iter()
    }

    fn warm_up(&mut self, model: &Model) -> Result<()> {
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

        for a in self.data.graph_terms.iter_mut() {
            a.warm_up(&self.settings, model)?;
        }
        Ok(())
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
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

    fn get_master_graph(&self, group_id: GroupId) -> &Self::G {
        println!("group structure: {:?}", self.data.graph_group_structure);
        let group_master = self.data.graph_group_structure[group_id].master();

        &self.data.graph_terms[group_master]
    }

    fn get_group(&self, group_id: GroupId) -> &crate::graph::GraphGroup {
        &self.data.graph_group_structure[group_id]
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor {
        DependentMomentaConstructor::CrossSection
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
    pub iterative_integrand: Option<TiVec<CutId, GenericEvaluator>>,
    pub parametric_integrand: TiVec<CutId, GenericEvaluator>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder<f64>,
    pub orientations: TiVec<SuperGraphOrientationID, EdgeVec<Orientation>>,
    pub orientation_filter: BinVec,
}

impl CrossSectionGraphTerm {
    pub fn from_cross_section_graph(
        graph: &CrossSectionGraph,
        settings: &GlobalSettings,
    ) -> Result<Self> {
        let orientations: TiVec<SuperGraphOrientationID, EdgeVec<Orientation>> = graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientation_data
            .iter()
            .map(|data| data.orientation.clone())
            .collect();

        let parametric_integrand = graph
            .derived_data
            .cut_paramatric_integrand
            .iter()
            .map(|integrand_for_cut| {
                GenericEvaluator::new_from_builder(
                    [integrand_for_cut.clone()],
                    &graph.graph.param_builder,
                    OptimizationSettings::default(),
                )
                .unwrap()
            })
            .collect::<TiVec<CutId, GenericEvaluator>>();

        let iterative_integrand = if settings
            .generation
            .evaluator
            .iterative_orientation_optimization
        {
            Some(
                graph
                    .derived_data
                    .cut_paramatric_integrand
                    .iter()
                    .map(|integrand_for_cut| {
                        GenericEvaluator::new_from_builder(
                            orientations.iter().map(|or| or.select(integrand_for_cut)),
                            &graph.graph.param_builder,
                            OptimizationSettings::default(),
                        )
                        .unwrap()
                    })
                    .collect::<TiVec<CutId, GenericEvaluator>>(),
            )
        } else {
            None
        };

        Ok(Self {
            iterative_integrand,
            parametric_integrand,
            graph: graph.graph.clone(),
            cut_esurface: graph.cut_esurface.clone(),
            cuts: graph.cuts.clone(),
            multi_channeling_setup: graph
                .derived_data
                .multi_channeling_setup
                .as_ref()
                .unwrap()
                .clone(),
            lmbs: graph.derived_data.lmbs.as_ref().unwrap().clone(),
            estimated_scale: None,
            param_builder: graph.graph.param_builder.clone(),
            orientation_filter: BinVec(BitVec::repeat(true, orientations.len())),
            orientations,
        })
    }
}

impl GraphTerm for CrossSectionGraphTerm {
    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64> {
        &mut self.param_builder
    }

    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()> {
        self.estimated_scale = Some(
            self.graph
                .expected_scale(F(settings.kinematics.e_cm), model),
        );

        self.orientation_filter = BinVec(
            self.orientations
                .iter()
                .map(|or| settings.general.orientation_pat.filter(or))
                .collect(),
        );

        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::CrossSection)
            .with_context(|| {
                format!(
                    "Failed to get dependent external momenta for graph {}",
                    self.graph.name
                )
            })?;
        self.graph.param_builder.add_external_four_mom(&externals);

        self.graph.param_builder.add_external_four_mom(&externals);
        let pols = self.graph.param_builder.pairs.polarizations_values(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );

        self.graph.param_builder.values[self
            .graph
            .param_builder
            .pairs
            .polarizations
            .value_range
            .clone()]
        .clone_from_slice(&pols);

        self.graph
            .param_builder
            .m_uv_value(Complex::new_re(F(settings.general.m_uv)));
        self.graph
            .param_builder
            .mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq)));
        self.graph.param_builder.update_model_values(model);

        self.param_builder = self.graph.param_builder.clone();

        Ok(())
    }
    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        settings: &RuntimeSettings,
        rotation: &Rotation,
    ) -> Complex<F<T>> {
        let orientations =
            momentum_sample.orientations(&self.orientation_filter.0, &self.orientations);

        // let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        let center = LoopMomenta::from_iter(vec![
            ThreeMomentum {
                px: momentum_sample.zero(),
                py: momentum_sample.zero(),
                pz: momentum_sample.zero(),
            };
            momentum_sample.loop_moms().len()
        ]);
        let masses = self.graph.get_real_mass_vector(&model);
        let hel = settings.kinematics.externals.get_helicities();
        let mut cut_results = TiVec::<CutId, Complex<F<T>>>::new();

        debug!("loop moms: {}", momentum_sample.loop_moms());

        for (cut, esurface) in self.cut_esurface.iter_enumerated() {
            let function = |t: &F<T>| {
                esurface.compute_self_and_r_derivative(
                    t,
                    momentum_sample.loop_moms(),
                    &center,
                    momentum_sample.external_moms(),
                    &self.graph.loop_momentum_basis,
                    &masses,
                )
            };

            let (guess, _) = esurface.get_radius_guess(
                momentum_sample.loop_moms(),
                momentum_sample.external_moms(),
                &self.graph.loop_momentum_basis,
            );

            let solution = newton_iteration_and_derivative(
                &guess,
                function,
                &F::from_f64(1.0),
                20,
                &F::from_f64(settings.kinematics.e_cm),
            );

            let h_function = h(
                &solution.solution,
                None,
                None,
                &HFunctionSettings::default(),
            );

            let lu_params = LUParams {
                h_function,
                tstar: solution.solution.clone(),
                esurface_derivative: solution.derivative_at_solution.clone(),
            };

            let rescaled_momenta = MomentumSample {
                sample: momentum_sample
                    .sample
                    .rescaled_loop_momenta(&solution.solution, None),
            };

            debug!("rescaled loop moms: {}", rescaled_momenta.loop_moms());

            let mut result = Complex::new_re(momentum_sample.zero());
            let params = T::get_parameters(
                &mut self.param_builder,
                settings.general.enable_cache,
                &self.graph,
                &rescaled_momenta,
                hel,
                None,
                Some(&lu_params),
            );

            if cut == CutId::from(2) {
                for expr in &self.iterative_integrand.as_ref().unwrap()[cut].exprs {
                    debug!("expr: {}", expr);
                }
            }

            let iterative = self
                .iterative_integrand
                .as_mut()
                .map(|ev| <T as GenericEvaluatorFloat>::get_evaluator(&mut ev[cut])(&params));

            for (i, e) in orientations.iter() {
                if let Some(iterative) = &iterative {
                    result += &iterative[i.0]
                } else {
                    self.param_builder.orientation_value(e);
                    let a = T::get_parameters(
                        &mut self.param_builder,
                        settings.general.enable_cache,
                        &self.graph,
                        &rescaled_momenta,
                        hel,
                        None,
                        Some(&lu_params),
                    );
                    result += <T as GenericEvaluatorFloat>::get_evaluator_single(
                        &mut self.parametric_integrand[cut],
                    )(&a)
                }
            }

            debug!("param builder for cut {}: \n{}", cut, self.param_builder);

            cut_results.push(result);
        }

        let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        for (cut_id, result) in cut_results.iter_enumerated() {
            debug!("Result for cut {}: {}", cut_id, result);
            all_cut_result += result;
        }

        all_cut_result
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
        self.orientations.len()
    }

    fn get_tropical_sampler(&self) -> &momtrop::SampleGenerator<3> {
        unimplemented!(
            "Don't know how to generate subgraph table for forward scattering graphs yet"
        )
    }

    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>> {
        todo!()
    }
}

impl CrossSectionGraphTerm {
    fn evaluate_impl<T: FloatLike>(
        &mut self,
        model: &Model,
        momentum_sample: &MomentumSample<T>,
        settings: &RuntimeSettings,
        param_builder: ParamBuilder<T>,
    ) -> Result<Complex<F<T>>> {
        todo!()
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
        model: &Model,
        wgt: F<f64>,
        _iter: usize,
        use_f128: bool,
        max_eval: Complex<F<f64>>,
    ) -> EvaluationResult {
        let result = evaluate_sample(self, model, sample, wgt, _iter, use_f128, max_eval);
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
            .map(|term| term.graph.get_loop_number())
            .all_equal());

        self.data.graph_terms[0].graph.get_loop_number() * 3
    }
}
