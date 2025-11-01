use linnet::half_edge::involution::{EdgeVec, Orientation};
use symbolica::evaluate::OptimizationSettings;
use typed_index_collections::TiVec;

use crate::{
    cff::{cut_expression::SuperGraphOrientationID, expression::GraphOrientation},
    gammaloop_integrand::{GenericEvaluator, ParamBuilder},
    processes::{IteratedCtCollection, LUCounterTermData, LeftThresholdId, RightThresholdId},
    settings::GlobalSettings,
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
