use core::f64;

use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::algebra::complex::Complex;
use symbolica::{domains::float::RealNumberLike, evaluate::OptimizationSettings};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::SuperGraphOrientationID,
        esurface::{Esurface, EsurfaceID},
        expression::GraphOrientation,
    },
    gammaloop_integrand::{GenericEvaluator, ParamBuilder},
    graph::{Graph, LmbIndex, LoopMomentumBasis},
    momentum_sample::{MomentumSample, SubspaceData},
    processes::{
        CutId, IteratedCtCollection, LUCounterTermData, LeftThresholdId, RightThresholdId,
    },
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        overlap::find_maximal_overlap,
        overlap_subspace::{self, OverlapInput},
    },
    utils::{FloatLike, F},
};

pub struct LUCounterTermEvaluators {
    pub parametric_left_thresholds_evaluator: TiVec<LeftThresholdId, GenericEvaluator>,
    pub parametric_right_threshold_evaluator: TiVec<RightThresholdId, GenericEvaluator>,
    pub parametric_iterated_evaluator: IteratedCtCollection<GenericEvaluator>,
    pub iterative_left_thresholds_evaluator: Option<TiVec<LeftThresholdId, GenericEvaluator>>,
    pub iterative_right_threshold_evaluator: Option<TiVec<RightThresholdId, GenericEvaluator>>,
    pub iterative_iterated_evaluator: Option<IteratedCtCollection<GenericEvaluator>>,
}

impl LUCounterTermEvaluators {
    pub fn from_atoms(
        counterterm_data: &LUCounterTermData,
        param_builder: &ParamBuilder,
        settings: &GlobalSettings,
        orientations: &TiVec<SuperGraphOrientationID, EdgeVec<Orientation>>,
    ) -> Self {
        let parametric_left_thresholds_evaluator = counterterm_data
            .left_atoms
            .iter()
            .map(|atom| {
                GenericEvaluator::new_from_builder(
                    [atom.clone()],
                    &param_builder,
                    OptimizationSettings::default(),
                )
                .unwrap()
            })
            .collect();

        let parametric_right_threshold_evaluator = counterterm_data
            .right_atoms
            .iter()
            .map(|atom| {
                GenericEvaluator::new_from_builder(
                    [atom.clone()],
                    &param_builder,
                    OptimizationSettings::default(),
                )
                .unwrap()
            })
            .collect();

        let parametric_iterated_evaluator = counterterm_data.iterated.map_ref(|atom| {
            GenericEvaluator::new_from_builder(
                [atom.clone()],
                param_builder,
                OptimizationSettings::default(),
            )
            .unwrap()
        });

        let iterative_left_thresholds_evaluator = if settings
            .generation
            .evaluator
            .iterative_orientation_optimization
        {
            Some(
                counterterm_data
                    .left_atoms
                    .iter()
                    .map(|atom| {
                        GenericEvaluator::new_from_builder(
                            orientations.iter().map(|or| or.select(atom)),
                            &param_builder,
                            OptimizationSettings::default(),
                        )
                        .unwrap()
                    })
                    .collect(),
            )
        } else {
            None
        };

        let iterative_right_threshold_evaluator = if settings
            .generation
            .evaluator
            .iterative_orientation_optimization
        {
            Some(
                counterterm_data
                    .right_atoms
                    .iter()
                    .map(|atom| {
                        GenericEvaluator::new_from_builder(
                            orientations.iter().map(|or| or.select(atom)),
                            &param_builder,
                            OptimizationSettings::default(),
                        )
                        .unwrap()
                    })
                    .collect(),
            )
        } else {
            None
        };

        let iterative_iterated_evaluator = if settings
            .generation
            .evaluator
            .iterative_orientation_optimization
        {
            Some(counterterm_data.iterated.map_ref(|atom| {
                GenericEvaluator::new_from_builder(
                    orientations.iter().map(|or| or.select(atom)),
                    param_builder,
                    OptimizationSettings::default(),
                )
                .unwrap()
            }))
        } else {
            None
        };

        LUCounterTermEvaluators {
            parametric_left_thresholds_evaluator,
            parametric_right_threshold_evaluator,
            parametric_iterated_evaluator,
            iterative_left_thresholds_evaluator,
            iterative_right_threshold_evaluator,
            iterative_iterated_evaluator,
        }
    }
}

pub struct LUCounterTerm {
    pub evaluators: LUCounterTermEvaluators,
    pub thresholds: TiVec<
        CutId,
        (
            TiVec<LeftThresholdId, Esurface>,
            TiVec<RightThresholdId, Esurface>,
        ),
    >,
    pub subspaces: TiVec<CutId, (SubspaceData, SubspaceData)>,
}

impl LUCounterTerm {
    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        cut_id: CutId,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        settings: &RuntimeSettings,
    ) -> Complex<F<T>> {
        let (left_subspace, right_subspace) = &self.subspaces[cut_id];
        let (left_thresholds_typed, right_thresholds_typed) = &self.thresholds[cut_id];

        let left_thresholds = TiVec::from_ref(&left_thresholds_typed.raw);
        let right_thresholds = TiVec::from_ref(&right_thresholds_typed.raw);

        let masses_f64: EdgeVec<F<f64>> = masses.iter().map(|(_, m)| F(m.to_f64())).collect();
        let loop_moms_f64 = momentum_sample
            .loop_moms()
            .iter()
            .map(|lm| lm.to_f64())
            .collect();
        let external_moms_f64 = momentum_sample
            .external_moms()
            .iter()
            .map(|em| em.to_f64())
            .collect();

        let e_cm = F::from_f64(settings.kinematics.e_cm);

        let left_existing_esurfaces = self.thresholds[cut_id]
            .0
            .iter_enumerated()
            .filter_map(|(left_id, esurface)| {
                if esurface.exists_subspace(
                    momentum_sample.loop_moms(),
                    momentum_sample.external_moms(),
                    left_subspace,
                    all_lmbs,
                    graph,
                    masses,
                    &e_cm,
                ) {
                    Some(EsurfaceID::from(left_id.0))
                } else {
                    return None;
                }
            })
            .collect();

        let right_existing_esurfaces = self.thresholds[cut_id]
            .1
            .iter_enumerated()
            .filter_map(|(right_id, esurface)| {
                if esurface.exists_subspace(
                    momentum_sample.loop_moms(),
                    momentum_sample.external_moms(),
                    right_subspace,
                    all_lmbs,
                    graph,
                    masses,
                    &e_cm,
                ) {
                    Some(EsurfaceID::from(right_id.0))
                } else {
                    None
                }
            })
            .collect();

        let left_overlap_input = OverlapInput {
            graph,
            subspace: left_subspace,
            settings,
            edge_masses: masses_f64.clone(),
            lmbs: all_lmbs,
            thresholds: left_thresholds,
        };

        let right_overlap_input = OverlapInput {
            graph,
            subspace: right_subspace,
            settings,
            edge_masses: masses_f64,
            lmbs: all_lmbs,
            thresholds: right_thresholds,
        };

        let left_overlap = if let Ok(left_overlap) = overlap_subspace::find_maximal_overlap(
            &left_overlap_input,
            &left_existing_esurfaces,
            &loop_moms_f64,
            &external_moms_f64,
        ) {
            left_overlap
        } else {
            return Complex::new_re(F::from_f64(f64::NAN));
        };

        let right_overlap = if let Ok(right_overlap) = overlap_subspace::find_maximal_overlap(
            &right_overlap_input,
            &right_existing_esurfaces,
            &loop_moms_f64,
            &external_moms_f64,
        ) {
            right_overlap
        } else {
            return Complex::new_re(F::from_f64(f64::NAN));
        };

        todo!()
    }
}
