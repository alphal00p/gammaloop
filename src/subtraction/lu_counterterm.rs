use core::f64;

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::algebra::complex::{sub, Complex};
use symbolica::{
    domains::float::{NumericalFloatLike, Real, RealNumberLike},
    evaluate::OptimizationSettings,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::SuperGraphOrientationID,
        esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId},
        expression::GraphOrientation,
    },
    gammaloop_integrand::{GenericEvaluator, ParamBuilder},
    graph::{Graph, LmbIndex, LoopMomentumBasis},
    model::Model,
    momentum::Rotation,
    momentum_sample::{LoopMomenta, MomentumSample, SubspaceData},
    processes::{
        CutId, IteratedCtCollection, LUCounterTermData, LeftThresholdId, RightThresholdId,
    },
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        overlap::find_maximal_overlap,
        overlap_subspace::{self, OverlapGroup, OverlapInput, OverlapStructure},
    },
    utils::{
        newton_solver::{newton_iteration_and_derivative, NewtonIterationResult},
        FloatLike, F,
    },
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

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    real_mass_vector: &'a EdgeVec<F<T>>,
    e_cm: F<T>,
    graph: &'a Graph,
    subspace: &'a SubspaceData,
    all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
    rotation_for_overlap: &'a Rotation,
    settings: &'a RuntimeSettings,
    esurface_collection: &'a EsurfaceCollection,
    transformed_sample: MomentumSample<T>,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    fn new(
        graph: &'a Graph,
        rotation_for_overlap: &'a Rotation,
        settings: &'a RuntimeSettings,
        esurface_collection: &'a EsurfaceCollection,
        sample: &'a MomentumSample<T>,
        overlap_structure: &'a OverlapStructure,
        masses: &'a EdgeVec<F<T>>,
        all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
        subspace: &'a SubspaceData,
    ) -> Self {
        let e_cm = F::from_f64(settings.kinematics.e_cm);
        let transformed_sample =
            sample.lmb_transform(&graph.loop_momentum_basis, subspace.get_lmb(all_lmbs));

        Self {
            real_mass_vector: masses,
            e_cm,
            graph,
            rotation_for_overlap,
            settings,
            esurface_collection,
            overlap_structure,
            transformed_sample,
            all_lmbs,
            subspace,
        }
    }

    fn new_overlap_builder(&'a self, overlap_group: &'a OverlapGroup) -> OverlapBuilder<'a, T> {
        let center = &overlap_group.center;
        let subspace = self.subspace;

        let (unrotated_center, rotated_center) = (
            center.cast(),
            center.rotate(self.rotation_for_overlap).cast(),
        );

        let shifted_loop_momenta = self.transformed_sample.loop_moms() - &rotated_center;
        let radius = shifted_loop_momenta
            .hyper_radius_squared(Some(&subspace.iter_lmb_indices().collect_vec()))
            .sqrt();
        let unit_shifted_momenta = shifted_loop_momenta.rescale(
            &radius.inv(),
            Some(&subspace.iter_lmb_indices().collect_vec()),
        );

        OverlapBuilder {
            counterterm_builder: self,
            overlap_group,
            rotated_center,
            _unrotated_center: unrotated_center,
            unit_shifted_momenta,
            radius,
        }
    }
}

struct OverlapBuilder<'a, T: FloatLike> {
    counterterm_builder: &'a CounterTermBuilder<'a, T>,
    overlap_group: &'a OverlapGroup,
    rotated_center: LoopMomenta<F<T>>,
    _unrotated_center: LoopMomenta<F<T>>,
    unit_shifted_momenta: LoopMomenta<F<T>>,
    radius: F<T>,
}

impl<'a, T: FloatLike> OverlapBuilder<'a, T> {
    fn new_esurface_builder(
        &'a self,
        existing_esurface_id: ExistingEsurfaceId,
    ) -> EsurfaceCTBuilder<'a, T> {
        let esurface_id = self
            .counterterm_builder
            .overlap_structure
            .existing_esurfaces[existing_esurface_id];

        EsurfaceCTBuilder {
            overlap_builder: self,
            _existing_esurface_id: existing_esurface_id,
            esurface: &self.counterterm_builder.esurface_collection[esurface_id],
            esurface_id,
        }
    }
}

const MAX_ITERATIONS: usize = 40;
const TOLERANCE: f64 = 1.0;

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    _existing_esurface_id: ExistingEsurfaceId,
    esurface: &'a Esurface,
    esurface_id: EsurfaceID,
}

impl<'a, T: FloatLike> EsurfaceCTBuilder<'a, T> {
    fn solve_rstar(self) -> RstarSolution<'a, T> {
        let subspace = self.overlap_builder.counterterm_builder.subspace;
        let graph = self.overlap_builder.counterterm_builder.graph;
        let lmbs = self.overlap_builder.counterterm_builder.all_lmbs;
        let masses = self.overlap_builder.counterterm_builder.real_mass_vector;

        let (radius_guess, _) = self.esurface.get_radius_guess_subspace(
            &self.overlap_builder.unit_shifted_momenta,
            &self
                .overlap_builder
                .counterterm_builder
                .transformed_sample
                .external_moms(),
            subspace,
            lmbs,
            graph,
            masses,
        );

        let function = |r: &_| {
            self.esurface.compute_self_and_r_derivative_subspace(
                r,
                &self.overlap_builder.unit_shifted_momenta,
                &self.overlap_builder.rotated_center,
                self.overlap_builder
                    .counterterm_builder
                    .transformed_sample
                    .external_moms(),
                &self.overlap_builder.counterterm_builder.real_mass_vector,
                subspace,
                lmbs,
                graph,
            )
        };

        let solution = newton_iteration_and_derivative(
            &radius_guess,
            function,
            &F::from_f64(TOLERANCE),
            MAX_ITERATIONS,
            &self.overlap_builder.counterterm_builder.e_cm,
        );

        RstarSolution {
            esurface_ct_builder: self,
            solution,
        }
    }
}

struct RstarSolution<'a, T: FloatLike> {
    esurface_ct_builder: EsurfaceCTBuilder<'a, T>,
    solution: NewtonIterationResult<T>,
}

impl<'a, T: FloatLike> RstarSolution<'a, T> {
    fn rstar_samples(self) -> RstarSample<'a, T> {
        let subspace = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .subspace;

        let rstar_loop_momenta = &self
            .esurface_ct_builder
            .overlap_builder
            .unit_shifted_momenta
            .rescale(
                &self.solution.solution,
                Some(&subspace.iter_lmb_indices().collect_vec()),
            )
            + &self.esurface_ct_builder.overlap_builder.rotated_center;

        let mut rstar_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .transformed_sample
            .clone();

        rstar_sample.sample.loop_moms = rstar_loop_momenta;

        RstarSample {
            rstar_solution: self,
            rstar_sample,
        }
    }
}

struct RstarSample<'a, T: FloatLike> {
    rstar_solution: RstarSolution<'a, T>,
    rstar_sample: MomentumSample<T>,
}
