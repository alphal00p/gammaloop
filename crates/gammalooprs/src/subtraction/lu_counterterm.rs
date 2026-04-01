use core::f64;

use bincode_trait_derive::{Decode, Encode};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use spenso::algebra::complex::Complex;
use symbolica::domains::float::RealLike;
use tracing::debug;
use typed_index_collections::TiVec;

use crate::{
    GammaLoopContext,
    cff::{
        esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId},
        expression::OrientationID,
    },
    graph::{Graph, LmbIndex, LoopMomentumBasis},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            GenericEvaluator, ParamBuilder, ThresholdParams,
            evaluators::{EvaluatorStack, SingleOrAllOrientations, evaluate_evaluator_single},
            param_builder::LUParams,
        },
    },
    momentum::{
        Rotation,
        sample::{LoopMomenta, MomentumSample, SubspaceData},
    },
    processes::{
        IteratedCtCollection, LUCounterTermData, LeftThresholdId, RaisedCutId, RightThresholdId,
        build_derivative_structure,
    },
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        evaluate_integrated_ct_normalisation, evaluate_uv_damper,
        overlap_subspace::{self, OverlapGroup, OverlapInput, OverlapStructure},
    },
    utils::{
        F, FloatLike,
        hyperdual_utils::shape_for_t_derivatives,
        newton_solver::{NewtonIterationResult, newton_iteration_and_derivative},
    },
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LUCounterTermEvaluators {
    pub left_thresholds_evaluator: TiVec<LeftThresholdId, Vec<EvaluatorStack>>,
    pub right_thresholds_evaluator: TiVec<RightThresholdId, Vec<EvaluatorStack>>,
    pub iterated_evaluator: IteratedCtCollection<Vec<EvaluatorStack>>,
    pub pass_two_evaluator: GenericEvaluator,
}

impl LUCounterTermEvaluators {
    pub fn from_atoms(
        counterterm_data: &LUCounterTermData,
        param_builder: &ParamBuilder,
        settings: &GlobalSettings,
        orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    ) -> Self {
        let left_thresholds_evaluator = counterterm_data
            .left_atoms
            .iter()
            .map(|parametric_integrands| {
                parametric_integrands
                    .integrands
                    .iter()
                    .enumerate()
                    .map(|(num_esurface, atom)| {
                        let dual_shape = if num_esurface > 1 {
                            Some(shape_for_t_derivatives(num_esurface - 1))
                        } else {
                            None
                        };

                        EvaluatorStack::new(
                            &[atom.clone()],
                            param_builder,
                            &orientations.raw,
                            dual_shape,
                            &settings.generation.evaluator,
                        )
                        .unwrap()
                    })
                    .collect()
            })
            .collect();

        let right_thresholds_evaluator = counterterm_data
            .right_atoms
            .iter()
            .map(|parametric_integrands| {
                parametric_integrands
                    .integrands
                    .iter()
                    .enumerate()
                    .map(|(num_esurface, atom)| {
                        let dual_shape = if num_esurface > 1 {
                            Some(shape_for_t_derivatives(num_esurface - 1))
                        } else {
                            None
                        };

                        EvaluatorStack::new(
                            &[atom.clone()],
                            param_builder,
                            &orientations.raw,
                            dual_shape,
                            &settings.generation.evaluator,
                        )
                        .unwrap()
                    })
                    .collect()
            })
            .collect();

        let iterated_evaluator = counterterm_data.iterated.map_ref(|parametric_integrands| {
            parametric_integrands
                .integrands
                .iter()
                .enumerate()
                .map(|(num_esurface, atom)| {
                    let dual_shape = if num_esurface > 1 {
                        Some(shape_for_t_derivatives(num_esurface - 1))
                    } else {
                        None
                    };

                    EvaluatorStack::new(
                        &[atom.clone()],
                        param_builder,
                        &orientations.raw,
                        dual_shape,
                        &settings.generation.evaluator,
                    )
                    .unwrap()
                })
                .collect()
        });

        let pass_two_evaluator = build_derivative_structure(1);

        LUCounterTermEvaluators {
            left_thresholds_evaluator,
            right_thresholds_evaluator,
            iterated_evaluator,
            pass_two_evaluator,
        }
    }
}

type CutThresholds = (
    TiVec<LeftThresholdId, Esurface>,
    TiVec<RightThresholdId, Esurface>,
);

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct LUCounterTerm {
    pub evaluators: TiVec<RaisedCutId, LUCounterTermEvaluators>,
    pub thresholds: TiVec<RaisedCutId, CutThresholds>,
    pub subspaces: TiVec<RaisedCutId, (SubspaceData, SubspaceData)>,
}

pub struct LUCTKinematicPoint<T: FloatLike> {
    pub dualized_momentum_sample_cache: Vec<MomentumSample<T>>,
    pub lu_cut_parameter_cache: Vec<LUParams<T>>,
    pub lu_cut_esurface_values: Vec<F<T>>,
}

impl<T: FloatLike> LUCTKinematicPoint<T> {
    pub fn representative_sample(&self) -> &MomentumSample<T> {
        &self.dualized_momentum_sample_cache[0]
    }

    pub fn lmb_transform(&self, from: &LoopMomentumBasis, to: &LoopMomentumBasis) -> Self {
        let transformed_samples = self
            .dualized_momentum_sample_cache
            .iter()
            .map(|sample| sample.lmb_transform(from, to))
            .collect();

        Self {
            dualized_momentum_sample_cache: transformed_samples,
            lu_cut_parameter_cache: self.lu_cut_parameter_cache.clone(),
            lu_cut_esurface_values: self.lu_cut_esurface_values.clone(),
        }
    }

    pub fn non_dual_cut_params(&self) -> LUParams<T> {
        self.lu_cut_parameter_cache[0].clone()
    }

    pub fn new() -> Self {
        Self {
            dualized_momentum_sample_cache: Vec::new(),
            lu_cut_parameter_cache: Vec::new(),
            lu_cut_esurface_values: Vec::new(),
        }
    }
}

impl LUCounterTerm {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate<T: FloatLike>(
        &mut self,
        kinematic_point: &LUCTKinematicPoint<T>,
        cut_id: RaisedCutId,
        reversed_edges: &[EdgeIndex],
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        rotation: &Rotation,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder<f64>,
        orientations: SingleOrAllOrientations<'_, OrientationID>,
        evaluation_meta_data: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Complex<F<T>> {
        let (left_subspace, right_subspace) = &self.subspaces[cut_id];
        let (sample_left_transformed, sample_right_transformed) = (
            kinematic_point
                .lmb_transform(&graph.loop_momentum_basis, left_subspace.get_lmb(all_lmbs)),
            kinematic_point
                .lmb_transform(&graph.loop_momentum_basis, right_subspace.get_lmb(all_lmbs)),
        );
        let (left_thresholds_typed, right_thresholds_typed) = &self.thresholds[cut_id];

        let left_thresholds = TiVec::from_ref(&left_thresholds_typed.raw);
        let right_thresholds = TiVec::from_ref(&right_thresholds_typed.raw);

        debug!("possible left thresholds: {}", left_thresholds.len());
        debug!("possible right thresholds: {}", right_thresholds.len());

        let masses_f64: EdgeVec<F<f64>> = masses.iter().map(|(_, m)| F(m.to_f64())).collect();
        let sample_left_transformed_f64 = sample_left_transformed
            .representative_sample()
            .loop_moms()
            .iter()
            .map(|lm| lm.to_f64())
            .collect();
        let sample_right_transformed_f64 = sample_right_transformed
            .representative_sample()
            .loop_moms()
            .iter()
            .map(|lm| lm.to_f64())
            .collect();
        let external_moms_f64 = kinematic_point
            .representative_sample()
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
                    sample_left_transformed.representative_sample().loop_moms(),
                    kinematic_point.representative_sample().external_moms(),
                    left_subspace,
                    all_lmbs,
                    graph,
                    masses,
                    reversed_edges,
                    &e_cm,
                ) {
                    Some(EsurfaceID::from(left_id.0))
                } else {
                    None
                }
            })
            .collect::<TiVec<ExistingEsurfaceId, _>>();

        let right_existing_esurfaces = self.thresholds[cut_id]
            .1
            .iter_enumerated()
            .filter_map(|(right_id, esurface)| {
                if esurface.exists_subspace(
                    sample_right_transformed.representative_sample().loop_moms(),
                    kinematic_point.representative_sample().external_moms(),
                    right_subspace,
                    all_lmbs,
                    graph,
                    masses,
                    reversed_edges,
                    &e_cm,
                ) {
                    Some(EsurfaceID::from(right_id.0))
                } else {
                    None
                }
            })
            .collect::<TiVec<ExistingEsurfaceId, _>>();

        debug!(
            "number of thresholds on the left: {}",
            left_existing_esurfaces.len()
        );
        debug!(
            "number of thresholds on the right: {}",
            right_existing_esurfaces.len()
        );

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
            &sample_left_transformed_f64,
            &external_moms_f64,
        ) {
            left_overlap
        } else {
            return Complex::new_re(F::from_f64(f64::NAN));
        };

        let right_overlap = if let Ok(right_overlap) = overlap_subspace::find_maximal_overlap(
            &right_overlap_input,
            &right_existing_esurfaces,
            &sample_right_transformed_f64,
            &external_moms_f64,
        ) {
            right_overlap
        } else {
            return Complex::new_re(F::from_f64(f64::NAN));
        };

        debug!(
            "left overlap structure: {:?}",
            left_overlap
                .overlap_groups
                .iter()
                .map(|group| group.existing_esurfaces.len())
                .collect_vec()
        );

        debug!(
            "right overlap structure: {:?}",
            right_overlap
                .overlap_groups
                .iter()
                .map(|group| group.existing_esurfaces.len())
                .collect_vec()
        );

        let left_counterterm_builder = CounterTermBuilder::new(
            graph,
            rotation,
            settings,
            left_overlap_input.thresholds,
            sample_left_transformed,
            &left_overlap,
            masses,
            all_lmbs,
            left_subspace,
        );

        let left_overlap_builders = left_overlap
            .overlap_groups
            .iter()
            .map(|overlap_group| left_counterterm_builder.new_overlap_builder(overlap_group))
            .collect_vec();

        let left_overlap_samples = left_overlap_builders
            .iter()
            .map(|overlap_builder| {
                overlap_builder
                    .overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|esurface| {
                        overlap_builder
                            .new_esurface_builder(*esurface)
                            .solve_rstar()
                            .rstar_sample()
                    })
                    .collect_vec()
            })
            .collect_vec();

        let right_counterterm_builder = CounterTermBuilder::new(
            graph,
            rotation,
            settings,
            right_overlap_input.thresholds,
            sample_right_transformed,
            &right_overlap,
            masses,
            all_lmbs,
            right_subspace,
        );

        let right_overlap_builders = right_overlap
            .overlap_groups
            .iter()
            .map(|overlap_group| right_counterterm_builder.new_overlap_builder(overlap_group))
            .collect_vec();

        let right_overlap_samples = right_overlap_builders
            .iter()
            .map(|overlap_builder| {
                overlap_builder
                    .overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|esurface| {
                        overlap_builder
                            .new_esurface_builder(*esurface)
                            .solve_rstar()
                            .rstar_sample()
                    })
                    .collect_vec()
            })
            .collect_vec();

        let mut left_evaluations = Complex::new_re(kinematic_point.representative_sample().zero());

        let non_dual_lu_cut_params = kinematic_point.non_dual_cut_params();

        for samples_group in left_overlap_samples.iter() {
            for sample in samples_group {
                debug!("left threshold parameters");

                let left_threshold_params: ThresholdParams<T> =
                    sample.extract_threshold_parameters(true);
                let inverse_transformed_sample = sample.get_inverse_transformed_sample();
                let left_threshold_id = LeftThresholdId::from(sample.get_esurface_id().0);

                let params = T::get_parameters(
                    param_builder,
                    (false, false),
                    graph,
                    &inverse_transformed_sample,
                    settings.kinematics.externals.get_helicities(),
                    &settings.additional_params(),
                    Some(&left_threshold_params),
                    None,
                    Some(&non_dual_lu_cut_params),
                );

                let mut result_of_this_ct = self.evaluators[cut_id].left_thresholds_evaluator
                    [left_threshold_id][0]
                    .evaluate(
                        params,
                        orientations,
                        settings,
                        evaluation_meta_data,
                        record_primary_timing,
                    )
                    .unwrap()
                    .pop()
                    .unwrap()
                    .unwrap_real();

                result_of_this_ct *= &sample.value_of_multi_channeling_factor;
                left_evaluations += result_of_this_ct;
            }
        }

        let mut right_evaluations = Complex::new_re(kinematic_point.representative_sample().zero());

        for samples_group in right_overlap_samples.iter() {
            for sample in samples_group {
                debug!("right threshold parameters");
                let right_threshold_params: ThresholdParams<T> =
                    sample.extract_threshold_parameters(true);
                let inverse_transformed_sample = sample.get_inverse_transformed_sample();
                let right_threshold_id = RightThresholdId::from(sample.get_esurface_id().0);

                let params = T::get_parameters(
                    param_builder,
                    (false, false),
                    graph,
                    &inverse_transformed_sample,
                    settings.kinematics.externals.get_helicities(),
                    &settings.additional_params(),
                    None,
                    Some(&right_threshold_params),
                    Some(&non_dual_lu_cut_params),
                );

                let mut result_of_this_ct = self.evaluators[cut_id].right_thresholds_evaluator
                    [right_threshold_id][0]
                    .evaluate(
                        params,
                        orientations,
                        settings,
                        evaluation_meta_data,
                        record_primary_timing,
                    )
                    .unwrap()
                    .pop()
                    .unwrap()
                    .unwrap_real();
                result_of_this_ct *= &sample.value_of_multi_channeling_factor;

                debug!(
                    "evaluation of ct for esurface id {}: {:+16e}",
                    sample.get_esurface_id().0,
                    result_of_this_ct,
                );

                right_evaluations += result_of_this_ct;
            }
        }

        let flattened_left_iter = left_overlap_samples.iter().flatten();
        let flattened_right_iter = right_overlap_samples.iter().flatten();
        let cartesian_product_iter = flattened_left_iter.cartesian_product(flattened_right_iter);

        let mut cartesian_product_result =
            Complex::new_re(kinematic_point.representative_sample().zero());

        for (sample_left, sample_right) in cartesian_product_iter {
            let left_threshold_params: ThresholdParams<T> =
                sample_left.extract_threshold_parameters(false);
            let right_threshold_params: ThresholdParams<T> =
                sample_right.extract_threshold_parameters(false);
            let multi_channeling_factor = &sample_left.value_of_multi_channeling_factor
                * &sample_right.value_of_multi_channeling_factor;
            let iterated_index = (
                LeftThresholdId::from(sample_left.get_esurface_id().0),
                RightThresholdId::from(sample_right.get_esurface_id().0),
            );
            let inverse_transformed_momentum_sample =
                merge_and_inverse_transform(sample_left, sample_right);

            let params = T::get_parameters(
                param_builder,
                (false, false),
                graph,
                &inverse_transformed_momentum_sample,
                settings.kinematics.externals.get_helicities(),
                &settings.additional_params(),
                Some(&left_threshold_params),
                Some(&right_threshold_params),
                Some(&non_dual_lu_cut_params),
            );

            let mut result_of_this_ct = self.evaluators[cut_id].iterated_evaluator[iterated_index]
                [0]
            .evaluate(
                params,
                orientations,
                settings,
                evaluation_meta_data,
                record_primary_timing,
            )
            .unwrap()
            .pop()
            .unwrap()
            .unwrap_real();

            result_of_this_ct *= multi_channeling_factor;
            debug!(
                "evaluation of ct for esurfaces {}, {}: {:+16e}",
                sample_left.get_esurface_id().0,
                sample_right.get_esurface_id().0,
                result_of_this_ct,
            );
            cartesian_product_result += result_of_this_ct;
        }

        debug!("left ct evaluation: {:+16e}", left_evaluations);
        debug!("right ct evaluation: {:+16e}", right_evaluations);
        debug!(
            "cartesian product ct evaluation: {:+16e}",
            cartesian_product_result
        );

        let pass_one_result = -(left_evaluations + right_evaluations - cartesian_product_result);

        let params_for_pass_two = [
            pass_one_result,
            Complex::new_re(kinematic_point.lu_cut_esurface_values[0].clone()),
        ];
        let pass_two_result = evaluate_evaluator_single(
            &mut self.evaluators[cut_id].pass_two_evaluator,
            &params_for_pass_two,
            evaluation_meta_data,
            record_primary_timing,
        );

        pass_two_result
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
    transformed_kinematic_point: LUCTKinematicPoint<T>,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        graph: &'a Graph,
        rotation_for_overlap: &'a Rotation,
        settings: &'a RuntimeSettings,
        esurface_collection: &'a EsurfaceCollection,
        transformed_kinematic_point: LUCTKinematicPoint<T>,
        overlap_structure: &'a OverlapStructure,
        masses: &'a EdgeVec<F<T>>,
        all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
        subspace: &'a SubspaceData,
    ) -> Self {
        let e_cm = F::from_f64(settings.kinematics.e_cm);

        Self {
            real_mass_vector: masses,
            e_cm,
            graph,
            rotation_for_overlap,
            settings,
            esurface_collection,
            overlap_structure,
            transformed_kinematic_point,
            all_lmbs,
            subspace,
        }
    }

    fn new_overlap_builder(&'a self, overlap_group: &'a OverlapGroup) -> OverlapBuilder<'a, T> {
        let center = &overlap_group.center;
        let subspace = self.subspace;

        let rotated_center = center.rotate(self.rotation_for_overlap).cast();

        let shifted_loop_momenta = self
            .transformed_kinematic_point
            .representative_sample()
            .loop_moms()
            - &rotated_center;

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
            unit_shifted_momenta,
            radius,
        }
    }
}

struct OverlapBuilder<'a, T: FloatLike> {
    counterterm_builder: &'a CounterTermBuilder<'a, T>,
    overlap_group: &'a OverlapGroup,
    rotated_center: LoopMomenta<F<T>>,
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
            self.overlap_builder
                .counterterm_builder
                .transformed_kinematic_point
                .representative_sample()
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
                    .transformed_kinematic_point
                    .representative_sample()
                    .external_moms(),
                self.overlap_builder.counterterm_builder.real_mass_vector,
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
    fn rstar_sample(self) -> RstarSample<'a, T> {
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
            .transformed_kinematic_point
            .representative_sample()
            .clone();

        rstar_sample.sample.loop_moms = rstar_loop_momenta;

        let multi_channeling_denominator = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .overlap_structure
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .complement
                    .iter()
                    .map(|existing_esurface_id| {
                        let esurface_id = self
                            .esurface_ct_builder
                            .overlap_builder
                            .counterterm_builder
                            .overlap_structure
                            .existing_esurfaces[*existing_esurface_id];

                        let esurface = &self
                            .esurface_ct_builder
                            .overlap_builder
                            .counterterm_builder
                            .esurface_collection[esurface_id];

                        let esurface_value = esurface.compute_from_momenta(
                            subspace.get_lmb(
                                self.esurface_ct_builder
                                    .overlap_builder
                                    .counterterm_builder
                                    .all_lmbs,
                            ),
                            self.esurface_ct_builder
                                .overlap_builder
                                .counterterm_builder
                                .real_mass_vector,
                            rstar_sample.loop_moms(),
                            rstar_sample.external_moms(),
                        );

                        &esurface_value * &esurface_value
                    })
                    .fold(rstar_sample.one(), |acc, val| acc * val)
            })
            .fold(rstar_sample.zero(), |acc, val| acc + val);

        let multichanneling_numerator = self
            .esurface_ct_builder
            .overlap_builder
            .overlap_group
            .complement
            .iter()
            .map(|existing_esurface_id| {
                let esurface_id = self
                    .esurface_ct_builder
                    .overlap_builder
                    .counterterm_builder
                    .overlap_structure
                    .existing_esurfaces[*existing_esurface_id];

                let esurface = &self
                    .esurface_ct_builder
                    .overlap_builder
                    .counterterm_builder
                    .esurface_collection[esurface_id];

                let esurface_value = esurface.compute_from_momenta(
                    subspace.get_lmb(
                        self.esurface_ct_builder
                            .overlap_builder
                            .counterterm_builder
                            .all_lmbs,
                    ),
                    self.esurface_ct_builder
                        .overlap_builder
                        .counterterm_builder
                        .real_mass_vector,
                    rstar_sample.loop_moms(),
                    rstar_sample.external_moms(),
                );

                &esurface_value * &esurface_value
            })
            .fold(rstar_sample.one(), |acc, val| acc * val);

        let value_of_multi_channeling_factor =
            Complex::new_re(multichanneling_numerator / multi_channeling_denominator);

        debug!(
            "value of multi-channeling factor: {}",
            value_of_multi_channeling_factor
        );

        RstarSample {
            rstar_solution: self,
            rstar_sample,
            value_of_multi_channeling_factor,
        }
    }
}

struct RstarSample<'a, T: FloatLike> {
    rstar_solution: RstarSolution<'a, T>,
    rstar_sample: MomentumSample<T>,
    // this can already be computed here because of non-crossing cuts
    value_of_multi_channeling_factor: Complex<F<T>>,
}

impl<'a, T: FloatLike> RstarSample<'a, T> {
    fn extract_threshold_parameters(&self, is_first_call: bool) -> ThresholdParams<T> {
        let radius = &self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .radius;

        let e_cm = &self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .e_cm;

        let uv_localisation_settings = &self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .settings
            .subtraction
            .local_ct_settings
            .uv_localisation;

        let radius_star = &self.rstar_solution.solution.solution;
        let esurface_derivative = &self.rstar_solution.solution.derivative_at_solution;
        let uv_damp_plus = evaluate_uv_damper(radius, radius_star, e_cm, uv_localisation_settings);
        let uv_damp_minus =
            evaluate_uv_damper(&-radius, radius_star, e_cm, uv_localisation_settings);

        let h_function = evaluate_integrated_ct_normalisation(
            radius,
            radius_star,
            e_cm,
            &self
                .rstar_solution
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .settings
                .subtraction
                .integrated_ct_settings,
        );

        if is_first_call {
            let edges_in_esurface = self
                .rstar_solution
                .esurface_ct_builder
                .esurface
                .energies
                .iter()
                .map(|i| {
                    self.rstar_solution
                        .esurface_ct_builder
                        .overlap_builder
                        .counterterm_builder
                        .graph[i]
                        .0
                        .name
                        .clone()
                })
                .collect_vec();

            debug!("esurface_id: {}", self.get_esurface_id().0);
            debug!("edges in esurface: {:?}", edges_in_esurface);
            debug!("radius: {}", radius);
            debug!("radius_star: {}", radius_star);
            debug!("esurface_derivative: {}", esurface_derivative);
            debug!("uv_damp_plus: {}", uv_damp_plus);
            debug!("uv_damp_minus: {}", uv_damp_minus);
            debug!("h_function: {}", h_function);
        }

        ThresholdParams {
            radius: radius.clone(),
            radius_star: radius_star.clone(),
            esurface_derivative: esurface_derivative.clone(),
            uv_damp_plus,
            uv_damp_minus,
            h_function,
        }
    }

    fn get_inverse_transformed_sample(&self) -> MomentumSample<T> {
        let subspace = self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .subspace;
        let current_lmb = subspace.get_lmb(
            self.rstar_solution
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .all_lmbs,
        );
        let target_lmb = &self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .graph
            .loop_momentum_basis;

        self.rstar_sample.lmb_transform(current_lmb, target_lmb)
    }

    fn get_esurface_id(&self) -> EsurfaceID {
        self.rstar_solution.esurface_ct_builder.esurface_id
    }
}

fn merge_and_inverse_transform<T: FloatLike>(
    left_sample: &RstarSample<'_, T>,
    right_sample: &RstarSample<'_, T>,
) -> MomentumSample<T> {
    let left_subspace = left_sample
        .rstar_solution
        .esurface_ct_builder
        .overlap_builder
        .counterterm_builder
        .subspace;

    let right_subspace = right_sample
        .rstar_solution
        .esurface_ct_builder
        .overlap_builder
        .counterterm_builder
        .subspace;

    assert!(
        left_subspace.is_mergable_with(right_subspace),
        "incompatible subspaces for merging samples"
    );

    let mut merged_sample = left_sample.rstar_sample.clone();
    for lmb_index in right_subspace.iter_lmb_indices() {
        let right_momentum = right_sample.rstar_sample.loop_moms()[lmb_index].clone();
        merged_sample.sample.loop_moms[lmb_index] = right_momentum;
    }

    let current_lmb = left_subspace.get_lmb(
        left_sample
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .all_lmbs,
    );

    let target_lmb = &left_sample
        .rstar_solution
        .esurface_ct_builder
        .overlap_builder
        .counterterm_builder
        .graph
        .loop_momentum_basis;

    merged_sample.lmb_transform(current_lmb, target_lmb)
}
