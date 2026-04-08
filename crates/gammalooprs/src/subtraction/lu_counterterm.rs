use core::f64;
use std::path::Path;

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use spenso::algebra::complex::Complex;
use symbolica::domains::{
    dual::HyperDual,
    float::{Real, RealLike},
};
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
        Rotation, SignOrZero, ThreeMomentum,
        sample::{
            ExternalFourMomenta, ExternalIndex, LoopIndex, LoopMomenta, MomentumSample,
            SubspaceData,
        },
        signature::LoopSignature,
    },
    processes::{
        EvaluatorBuildTimings, IteratedCtCollection, LUCounterTermData, LeftThresholdId,
        RaisedCutId, RightThresholdId, build_derivative_structure,
    },
    settings::{GlobalSettings, RuntimeSettings, global::FrozenCompilationMode},
    subtraction::{
        RstarTDependenceEvaluator, evaluate_integrated_ct_normalisation,
        evaluate_integrated_ct_normalisation_dual, evaluate_uv_damper, evaluate_uv_damper_dual,
        overlap_subspace::{self, OverlapGroup, OverlapInput, OverlapStructure},
    },
    utils::{
        F, FloatLike,
        hyperdual_utils::{
            DualOrNot, extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            shape_for_t_derivatives,
        },
        newton_solver::{NewtonIterationResult, newton_iteration_and_derivative},
    },
};

fn zero_dual_or_not_complex<T: FloatLike>(order: usize, zero: &F<T>) -> DualOrNot<Complex<F<T>>> {
    if order == 0 {
        DualOrNot::NonDual(Complex::new_re(zero.clone()))
    } else {
        DualOrNot::Dual(HyperDual::new(shape_for_t_derivatives(order)))
    }
}

fn negate_dual_or_not_complex<T: FloatLike>(
    value: DualOrNot<Complex<F<T>>>,
) -> DualOrNot<Complex<F<T>>> {
    match value {
        DualOrNot::Dual(dual) => DualOrNot::Dual(-dual),
        DualOrNot::NonDual(non_dual) => DualOrNot::NonDual(-non_dual),
    }
}

fn multiply_dual_or_not_complex<T: FloatLike>(
    lhs: DualOrNot<Complex<F<T>>>,
    rhs: &DualOrNot<Complex<F<T>>>,
) -> DualOrNot<Complex<F<T>>> {
    match (lhs, rhs) {
        (DualOrNot::NonDual(lhs), DualOrNot::NonDual(rhs)) => DualOrNot::NonDual(lhs * rhs),
        (DualOrNot::Dual(lhs), DualOrNot::Dual(rhs)) => DualOrNot::Dual(lhs * rhs.clone()),
        (DualOrNot::Dual(lhs), DualOrNot::NonDual(rhs)) => {
            let rhs_dual = new_constant(&lhs, rhs);
            DualOrNot::Dual(lhs * rhs_dual)
        }
        (DualOrNot::NonDual(lhs), DualOrNot::Dual(rhs)) => {
            let lhs_dual = new_constant(rhs, &lhs);
            DualOrNot::Dual(lhs_dual * rhs.clone())
        }
    }
}

fn real_dual_to_complex<T: FloatLike>(dual: HyperDual<F<T>>) -> HyperDual<Complex<F<T>>> {
    let values = dual.values.into_iter().map(Complex::new_re).collect_vec();
    HyperDual::from_values(shape_for_t_derivatives(values.len() - 1), values)
}

fn dualize_external_momenta<T: FloatLike>(
    shape: &HyperDual<F<T>>,
    external_moms: &ExternalFourMomenta<F<T>>,
) -> ExternalFourMomenta<HyperDual<F<T>>> {
    external_moms
        .iter()
        .map(|momentum| momentum.map_ref(&|component| new_constant(shape, component)))
        .collect()
}

fn dualize_loop_momenta<T: FloatLike>(
    shape: &HyperDual<F<T>>,
    loop_moms: &LoopMomenta<F<T>>,
) -> LoopMomenta<HyperDual<F<T>>> {
    loop_moms
        .iter()
        .map(|momentum| momentum.map_ref(&|component| new_constant(shape, component)))
        .collect()
}

fn dual_shifted_radius<T: FloatLike>(
    shifted_momenta: &LoopMomenta<HyperDual<F<T>>>,
    subspace: &SubspaceData,
) -> HyperDual<F<T>> {
    let zero = new_constant(
        &shifted_momenta[LoopIndex(0)].px,
        &shifted_momenta[LoopIndex(0)].px.values[0].zero(),
    );

    subspace
        .iter_lmb_indices()
        .fold(zero, |acc, loop_index| {
            acc + shifted_momenta[loop_index].norm_squared()
        })
        .sqrt()
}

fn compute_shift_part_from_dual_momenta_in_subspace<T: FloatLike>(
    esurface: &Esurface,
    loop_moms: &LoopMomenta<HyperDual<F<T>>>,
    external_moms: &ExternalFourMomenta<HyperDual<F<T>>>,
    subspace: &SubspaceData,
    all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    graph: &Graph,
    masses: &EdgeVec<F<T>>,
) -> HyperDual<F<T>> {
    let lmb = subspace.get_lmb(all_lmbs);
    let zero = new_constant(
        &external_moms[ExternalIndex(0)].temporal.value,
        &external_moms[ExternalIndex(0)].temporal.value.values[0].zero(),
    );

    let full_external_shift = esurface
        .external_shift
        .iter()
        .map(|(index, sign)| {
            let external_signature = &lmb.edge_signatures[*index].external;
            let sign = new_constant(&zero, &F::from_f64(*sign as f64));
            let external_energy = external_signature
                .try_apply(&external_moms.raw)
                .map(|momentum| momentum.temporal.value)
                .unwrap_or_else(|| zero.clone());

            sign * external_energy
        })
        .reduce(|acc, value| acc + value)
        .unwrap_or_else(|| zero.clone());

    let spatial_externals = external_moms
        .iter()
        .map(|momentum| momentum.spatial.clone())
        .collect::<TiVec<ExternalIndex, _>>();

    let remaining_shift = subspace
        .does_not_contain(&esurface.energies, graph)
        .map(|index| {
            let signature = &lmb.edge_signatures[index];
            let momentum = signature
                .try_compute_momentum(&loop_moms.0, &spatial_externals.raw)
                .unwrap_or_else(|| unreachable!());
            let mass = &masses[index];
            let lifted_mass = new_constant(&momentum.px, mass);

            (momentum.norm_squared() + lifted_mass.clone() * lifted_mass).sqrt()
        })
        .reduce(|acc, value| acc + value)
        .unwrap_or_else(|| zero.clone());

    full_external_shift + remaining_shift
}

#[allow(clippy::too_many_arguments)]
fn compute_self_and_r_derivative_subspace_dual<T: FloatLike>(
    esurface: &Esurface,
    radius: &HyperDual<F<T>>,
    shifted_unit_loops_in_subspace: &LoopMomenta<HyperDual<F<T>>>,
    center_in_subspace: &LoopMomenta<HyperDual<F<T>>>,
    external_moms: &ExternalFourMomenta<HyperDual<F<T>>>,
    real_mass_vector: &EdgeVec<F<T>>,
    subspace: &SubspaceData,
    all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    graph: &Graph,
) -> (HyperDual<F<T>>, HyperDual<F<T>>) {
    let spatial_part_of_externals = external_moms
        .iter()
        .map(|momentum| momentum.spatial.clone())
        .collect::<TiVec<ExternalIndex, _>>();

    let loops: LoopMomenta<HyperDual<F<T>>> = shifted_unit_loops_in_subspace
        .iter_enumerated()
        .map(|(loop_index, shifted_unit_momenta)| {
            if subspace.contains_loop_index(loop_index) {
                shifted_unit_momenta * radius + &center_in_subspace[loop_index]
            } else {
                shifted_unit_momenta.clone()
            }
        })
        .collect();

    let shift = compute_shift_part_from_dual_momenta_in_subspace(
        esurface,
        shifted_unit_loops_in_subspace,
        external_moms,
        subspace,
        all_lmbs,
        graph,
        real_mass_vector,
    );

    let lmb = subspace.get_lmb(all_lmbs);
    let zero = new_constant(radius, &radius.values[0].zero());
    let (derivative, energy_sum) = subspace
        .contains(&esurface.energies, graph)
        .map(|index| {
            let signature = &lmb.edge_signatures[index];
            let momentum = signature
                .try_compute_momentum(&loops.0, &spatial_part_of_externals.raw)
                .unwrap_or_else(|| unreachable!());
            let unit_loop_part = compute_loop_part_subspace_dual(
                &signature.internal,
                shifted_unit_loops_in_subspace,
                subspace,
            );
            let mass = &real_mass_vector[index];
            let lifted_mass = new_constant(&momentum.px, mass);
            let energy = (momentum.norm_squared() + lifted_mass.clone() * lifted_mass).sqrt();
            let numerator = momentum * &unit_loop_part;

            (numerator / &energy, energy)
        })
        .fold((zero.clone(), zero), |(der_sum, en_sum), (der, en)| {
            (der_sum + der, en_sum + en)
        });

    (energy_sum + shift, derivative)
}

fn compute_loop_part_subspace_dual<T: FloatLike>(
    loop_signature: &LoopSignature,
    loop_moms: &LoopMomenta<HyperDual<F<T>>>,
    subspace: &SubspaceData,
) -> ThreeMomentum<HyperDual<F<T>>> {
    let projected: LoopSignature = subspace.project_loop_signature(loop_signature).collect();
    let zero = new_constant(
        &loop_moms[LoopIndex(0)].px,
        &loop_moms[LoopIndex(0)].px.values[0].zero(),
    );
    let mut result = ThreeMomentum::new(zero.clone(), zero.clone(), zero);

    for (loop_index, sign) in projected.iter_enumerated() {
        match sign {
            SignOrZero::Zero => {}
            SignOrZero::Plus => result += &loop_moms[loop_index],
            SignOrZero::Minus => result -= &loop_moms[loop_index],
        }
    }

    result
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LUCounterTermEvaluators {
    pub left_thresholds_evaluator: TiVec<LeftThresholdId, Vec<EvaluatorStack>>,
    pub right_thresholds_evaluator: TiVec<RightThresholdId, Vec<EvaluatorStack>>,
    pub iterated_evaluator: IteratedCtCollection<Vec<EvaluatorStack>>,
    pub pass_two_evaluator: Vec<GenericEvaluator>,
}

impl LUCounterTermEvaluators {
    pub(crate) fn generic_evaluator_count(&self) -> usize {
        let left = self
            .left_thresholds_evaluator
            .iter()
            .flat_map(|evaluators| evaluators.iter())
            .map(EvaluatorStack::generic_evaluator_count)
            .sum::<usize>();
        let right = self
            .right_thresholds_evaluator
            .iter()
            .flat_map(|evaluators| evaluators.iter())
            .map(EvaluatorStack::generic_evaluator_count)
            .sum::<usize>();
        let iterated = self
            .iterated_evaluator
            .iter()
            .flat_map(|evaluators| evaluators.iter())
            .map(EvaluatorStack::generic_evaluator_count)
            .sum::<usize>();

        left + right + iterated + 1
    }

    pub fn from_atoms(
        counterterm_data: &LUCounterTermData,
        param_builder: &ParamBuilder,
        settings: &GlobalSettings,
        orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    ) -> (Self, EvaluatorBuildTimings) {
        let mut timings = EvaluatorBuildTimings::default();
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

                        let (evaluator, evaluator_timings) = EvaluatorStack::new_with_timings(
                            std::slice::from_ref(atom),
                            param_builder,
                            &orientations.raw,
                            dual_shape,
                            &settings.generation.evaluator,
                        )
                        .unwrap();
                        timings += evaluator_timings;
                        evaluator
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

                        let (evaluator, evaluator_timings) = EvaluatorStack::new_with_timings(
                            std::slice::from_ref(atom),
                            param_builder,
                            &orientations.raw,
                            dual_shape,
                            &settings.generation.evaluator,
                        )
                        .unwrap();
                        timings += evaluator_timings;
                        evaluator
                    })
                    .collect()
            })
            .collect();

        let iterated_timings = std::cell::Cell::new(EvaluatorBuildTimings::default());
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

                    let (evaluator, evaluator_timings) = EvaluatorStack::new_with_timings(
                        std::slice::from_ref(atom),
                        param_builder,
                        &orientations.raw,
                        dual_shape,
                        &settings.generation.evaluator,
                    )
                    .unwrap();
                    let mut timings = iterated_timings.get();
                    timings += evaluator_timings;
                    iterated_timings.set(timings);
                    evaluator
                })
                .collect()
        });
        timings += iterated_timings.get();

        let symbolica_started = std::time::Instant::now();
        let max_num_esurfaces = counterterm_data
            .left_atoms
            .iter()
            .map(|integrands| integrands.integrands.len())
            .chain(
                counterterm_data
                    .right_atoms
                    .iter()
                    .map(|integrands| integrands.integrands.len()),
            )
            .chain(
                counterterm_data
                    .iterated
                    .iter()
                    .map(|integrands| integrands.integrands.len()),
            )
            .max()
            .unwrap_or(1);

        let pass_two_evaluator = (1..=max_num_esurfaces)
            .map(|order| build_derivative_structure(order as u8))
            .collect(, &settings.generation.evaluator);
        timings.symbolica_time += symbolica_started.elapsed();

        (
            LUCounterTermEvaluators {
                left_thresholds_evaluator,
                right_thresholds_evaluator,
                iterated_evaluator,
                pass_two_evaluator,
            },
            timings,
        )
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        cut_id: RaisedCutId,
        frozen_mode: &FrozenCompilationMode,
    ) -> color_eyre::Result<()> {
        for (threshold_id, evaluators) in self.left_thresholds_evaluator.iter_mut_enumerated() {
            for (num_esurface, evaluator) in evaluators.iter_mut().enumerate() {
                let name = format!(
                    "cut_{}_left_threshold_{}_{}",
                    cut_id.0, threshold_id.0, num_esurface
                );
                evaluator.compile(&name, path.as_ref(), frozen_mode)?;
            }
        }

        for (threshold_id, evaluators) in self.right_thresholds_evaluator.iter_mut_enumerated() {
            for (num_esurface, evaluator) in evaluators.iter_mut().enumerate() {
                let name = format!(
                    "cut_{}_right_threshold_{}_{}",
                    cut_id.0, threshold_id.0, num_esurface
                );
                evaluator.compile(&name, path.as_ref(), frozen_mode)?;
            }
        }

        for (iterated_index, evaluators) in self.iterated_evaluator.iter_mut().enumerate() {
            for (num_esurface, evaluator) in evaluators.iter_mut().enumerate() {
                let name = format!(
                    "cut_{}_iterated_{}_{}",
                    cut_id.0, iterated_index, num_esurface
                );
                evaluator.compile(&name, path.as_ref(), frozen_mode)?;
            }
        }

        self.pass_two_evaluator.compile_external(
            path.as_ref()
                .join(format!("cut_{}_pass_two", cut_id.0))
                .with_extension("cpp"),
            format!("cut_{}_pass_two", cut_id.0),
            path.as_ref()
                .join(format!("cut_{}_pass_two", cut_id.0))
                .with_extension("so"),
            frozen_mode,
        )?;

        Ok(())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> color_eyre::Result<()>,
    ) -> color_eyre::Result<()> {
        for evaluators in self.left_thresholds_evaluator.iter_mut() {
            for evaluator in evaluators.iter_mut() {
                evaluator.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        for evaluators in self.right_thresholds_evaluator.iter_mut() {
            for evaluator in evaluators.iter_mut() {
                evaluator.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        for evaluators in self.iterated_evaluator.iter_mut() {
            for evaluator in evaluators.iter_mut() {
                evaluator.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        f(&mut self.pass_two_evaluator)?;

        Ok(())
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
    pub active_cuts: TiVec<RaisedCutId, bool>,
    pub active_left_thresholds: TiVec<RaisedCutId, TiVec<LeftThresholdId, bool>>,
    pub active_right_thresholds: TiVec<RaisedCutId, TiVec<RightThresholdId, bool>>,
    pub active_iterated_thresholds: TiVec<RaisedCutId, IteratedCtCollection<bool>>,
    pub rstar_dependence_calculator: TiVec<RaisedCutId, RstarTDependenceEvaluator>,
}

pub struct LUCTKinematicPoint<T: FloatLike> {
    pub unrescaled_sample: MomentumSample<T>,
    pub dualized_momentum_sample_cache: Vec<MomentumSample<T>>,
    pub lu_cut_parameter_cache: Vec<LUParams<T>>,
    pub lu_cut_esurface_values: Vec<DualOrNot<F<T>>>,
}

impl<T: FloatLike> LUCTKinematicPoint<T> {
    pub fn unrescaled_sample(&self) -> &MomentumSample<T> {
        &self.unrescaled_sample
    }

    pub fn representative_sample(&self) -> &MomentumSample<T> {
        &self.dualized_momentum_sample_cache[0]
    }

    pub fn sample_for_order(&self, order: usize) -> &MomentumSample<T> {
        &self.dualized_momentum_sample_cache[order]
    }

    pub fn lmb_transform(&self, from: &LoopMomentumBasis, to: &LoopMomentumBasis) -> Self {
        let transformed_samples = self
            .dualized_momentum_sample_cache
            .iter()
            .map(|sample| sample.lmb_transform(from, to))
            .collect();

        Self {
            unrescaled_sample: self.unrescaled_sample.clone(),
            dualized_momentum_sample_cache: transformed_samples,
            lu_cut_parameter_cache: self.lu_cut_parameter_cache.clone(),
            lu_cut_esurface_values: self.lu_cut_esurface_values.clone(),
        }
    }

    pub fn non_dual_cut_params(&self) -> LUParams<T> {
        self.lu_cut_parameter_cache[0].clone()
    }

    pub fn cut_params_for_order(&self, order: usize) -> &LUParams<T> {
        &self.lu_cut_parameter_cache[order]
    }

    pub fn cut_esurface_for_order(&self, order: usize) -> &DualOrNot<F<T>> {
        &self.lu_cut_esurface_values[order]
    }

    pub fn non_dual_cut_esurface_value(&self) -> F<T> {
        match self.cut_esurface_for_order(0) {
            DualOrNot::NonDual(value) => value.clone(),
            DualOrNot::Dual(_) => {
                unreachable!("representative LU cut e-surface value must be non-dual")
            }
        }
    }

    pub fn new(unrescaled_sample: MomentumSample<T>) -> Self {
        Self {
            unrescaled_sample,
            dualized_momentum_sample_cache: Vec::new(),
            lu_cut_parameter_cache: Vec::new(),
            lu_cut_esurface_values: Vec::new(),
        }
    }
}

impl<T: FloatLike> Default for LUCTKinematicPoint<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl LUCounterTerm {
    pub(crate) fn cut_is_active(&self, cut_id: RaisedCutId) -> bool {
        self.active_cuts[cut_id]
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        frozen_mode: &FrozenCompilationMode,
    ) -> color_eyre::Result<()> {
        for (cut_id, evaluators) in self.evaluators.iter_mut_enumerated() {
            evaluators.compile(path.as_ref(), cut_id, frozen_mode)?;
        }
        Ok(())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> color_eyre::Result<()>,
    ) -> color_eyre::Result<()> {
        for evaluators in self.evaluators.iter_mut() {
            evaluators.for_each_generic_evaluator_mut(&mut f)?;
        }
        Ok(())
    }

    fn ensure_active_cut(&self, cut_id: RaisedCutId) -> Result<()> {
        if self.cut_is_active(cut_id) {
            return Ok(());
        }

        Err(eyre!(
            "Raised cut {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            cut_id.0
        ))
    }

    fn ensure_active_left_threshold(
        &self,
        cut_id: RaisedCutId,
        threshold_id: LeftThresholdId,
    ) -> Result<()> {
        if self.active_left_thresholds[cut_id][threshold_id] {
            return Ok(());
        }

        Err(eyre!(
            "Left threshold evaluator {} for raised cut {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            threshold_id.0,
            cut_id.0
        ))
    }

    fn ensure_active_right_threshold(
        &self,
        cut_id: RaisedCutId,
        threshold_id: RightThresholdId,
    ) -> Result<()> {
        if self.active_right_thresholds[cut_id][threshold_id] {
            return Ok(());
        }

        Err(eyre!(
            "Right threshold evaluator {} for raised cut {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            threshold_id.0,
            cut_id.0
        ))
    }

    fn ensure_active_iterated_threshold(
        &self,
        cut_id: RaisedCutId,
        iterated_index: (LeftThresholdId, RightThresholdId),
    ) -> Result<()> {
        if self.active_iterated_thresholds[cut_id][iterated_index] {
            return Ok(());
        }

        Err(eyre!(
            "Iterated threshold evaluator ({}, {}) for raised cut {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            iterated_index.0.0,
            iterated_index.1.0,
            cut_id.0
        ))
    }

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
    ) -> Result<Complex<F<T>>> {
        self.ensure_active_cut(cut_id)?;
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
            return Ok(Complex::new_re(F::from_f64(f64::NAN)));
        };

        let right_overlap = if let Ok(right_overlap) = overlap_subspace::find_maximal_overlap(
            &right_overlap_input,
            &right_existing_esurfaces,
            &sample_right_transformed_f64,
            &external_moms_f64,
        ) {
            right_overlap
        } else {
            return Ok(Complex::new_re(F::from_f64(f64::NAN)));
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

        let left_overlap_solutions = left_overlap_builders
            .iter()
            .map(|overlap_builder| {
                overlap_builder
                    .overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|esurface| {
                        overlap_builder
                            .new_esurface_builder(*esurface)
                            .solve_rstar(&mut self.rstar_dependence_calculator[cut_id])
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

        let right_overlap_solutions = right_overlap_builders
            .iter()
            .map(|overlap_builder| {
                overlap_builder
                    .overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|esurface| {
                        overlap_builder
                            .new_esurface_builder(*esurface)
                            .solve_rstar(&mut self.rstar_dependence_calculator[cut_id])
                    })
                    .collect_vec()
            })
            .collect_vec();

        let zero = kinematic_point.representative_sample().zero();
        let mut total_result = Complex::new_re(zero.clone());

        for order in 0..kinematic_point.lu_cut_parameter_cache.len() {
            let mut left_evaluations = zero_dual_or_not_complex(order, &zero);
            let lu_cut_params = kinematic_point.cut_params_for_order(order).clone();

            for solutions_group in left_overlap_solutions.iter() {
                for solution in solutions_group {
                    let sample = solution.rstar_sample_for_order(order);
                    debug!("left threshold parameters");

                    let left_threshold_params: ThresholdParams<T> =
                        sample.extract_threshold_parameters(true);
                    let inverse_transformed_sample = sample.get_inverse_transformed_sample();
                    let left_threshold_id = LeftThresholdId::from(sample.get_esurface_id().0);
                self.ensure_active_left_threshold(cut_id, left_threshold_id)?;

                    let params = T::get_parameters(
                        param_builder,
                        (false, false),
                        graph,
                        &inverse_transformed_sample,
                        settings.kinematics.externals.get_helicities(),
                        &settings.additional_params(),
                        Some(&left_threshold_params),
                        None,
                        Some(&lu_cut_params),
                    );

                    let result_of_this_ct = self.evaluators[cut_id].left_thresholds_evaluator
                        [left_threshold_id][order]
                        .evaluate(
                            params,
                            orientations,
                            settings,
                            evaluation_meta_data,
                            record_primary_timing,
                        )
                        .unwrap()
                        .pop()
                        .unwrap();

                    left_evaluations += multiply_dual_or_not_complex(
                        result_of_this_ct,
                        &sample.value_of_multi_channeling_factor,
                    );
                }
            }

            let mut right_evaluations = zero_dual_or_not_complex(order, &zero);

            for solutions_group in right_overlap_solutions.iter() {
                for solution in solutions_group {
                    let sample = solution.rstar_sample_for_order(order);
                    debug!("right threshold parameters");
                    let right_threshold_params: ThresholdParams<T> =
                        sample.extract_threshold_parameters(true);
                    let inverse_transformed_sample = sample.get_inverse_transformed_sample();
                    let right_threshold_id = RightThresholdId::from(sample.get_esurface_id().0);
                self.ensure_active_right_threshold(cut_id, right_threshold_id)?;

                    let params = T::get_parameters(
                        param_builder,
                        (false, false),
                        graph,
                        &inverse_transformed_sample,
                        settings.kinematics.externals.get_helicities(),
                        &settings.additional_params(),
                        None,
                        Some(&right_threshold_params),
                        Some(&lu_cut_params),
                    );

                    let result_of_this_ct = self.evaluators[cut_id].right_thresholds_evaluator
                        [right_threshold_id][order]
                        .evaluate(
                            params,
                            orientations,
                            settings,
                            evaluation_meta_data,
                            record_primary_timing,
                        )
                        .unwrap()
                        .pop()
                        .unwrap();

                    right_evaluations += multiply_dual_or_not_complex(
                        result_of_this_ct,
                        &sample.value_of_multi_channeling_factor,
                    );
                }
            }

            let mut cartesian_product_result = zero_dual_or_not_complex(order, &zero);

            for left_solutions_group in left_overlap_solutions.iter() {
                for left_solution in left_solutions_group {
                    let sample_left = left_solution.rstar_sample_for_order(order);
                    for right_solutions_group in right_overlap_solutions.iter() {
                        for right_solution in right_solutions_group {
                            let sample_right = right_solution.rstar_sample_for_order(order);

                            let left_threshold_params: ThresholdParams<T> =
                                sample_left.extract_threshold_parameters(false);
                            let right_threshold_params: ThresholdParams<T> =
                                sample_right.extract_threshold_parameters(false);
                            let multi_channeling_factor = multiply_dual_or_not_complex(
                                sample_left.value_of_multi_channeling_factor.clone(),
                                &sample_right.value_of_multi_channeling_factor,
                            );
                            let iterated_index = (
                                LeftThresholdId::from(sample_left.get_esurface_id().0),
                                RightThresholdId::from(sample_right.get_esurface_id().0),
                            );
                            self.ensure_active_iterated_threshold(cut_id, iterated_index)?;
            let inverse_transformed_momentum_sample =
                                merge_and_inverse_transform(&sample_left, &sample_right);

                            let params = T::get_parameters(
                                param_builder,
                                (false, false),
                                graph,
                                &inverse_transformed_momentum_sample,
                                settings.kinematics.externals.get_helicities(),
                                &settings.additional_params(),
                                Some(&left_threshold_params),
                                Some(&right_threshold_params),
                                Some(&lu_cut_params),
                            );

                            let result_of_this_ct = self.evaluators[cut_id].iterated_evaluator
                                [iterated_index][order]
                                .evaluate(
                                    params,
                                    orientations,
                                    settings,
                                    evaluation_meta_data,
                                    record_primary_timing,
                                )
                                .unwrap()
                                .pop()
                                .unwrap();

                            cartesian_product_result += multiply_dual_or_not_complex(
                                result_of_this_ct,
                                &multi_channeling_factor,
                            );
                        }
                    }
                }
            }

            let mut pass_one_result = left_evaluations.clone();
            pass_one_result += right_evaluations;
            pass_one_result += negate_dual_or_not_complex(cartesian_product_result);
            let pass_one_result = negate_dual_or_not_complex(pass_one_result);

            let mut params_for_pass_two = vec![];
            match pass_one_result {
                DualOrNot::Dual(dual_result) => {
                    params_for_pass_two
                        .extend_from_slice(&extract_t_derivatives_complex(dual_result));
                }
                DualOrNot::NonDual(non_dual_result) => {
                    params_for_pass_two.push(non_dual_result);
                }
            }

            match kinematic_point.cut_esurface_for_order(order).clone() {
                DualOrNot::Dual(dual_e_surface) => {
                    extract_t_derivatives(dual_e_surface)[1..]
                        .iter()
                        .for_each(|value| params_for_pass_two.push(Complex::new_re(value.clone())));
                }
                DualOrNot::NonDual(non_dual_e_surface) => {
                    params_for_pass_two.push(Complex::new_re(non_dual_e_surface));
                }
            }

            let pass_two_result = evaluate_evaluator_single(
                &mut self.evaluators[cut_id].pass_two_evaluator[order],
                &params_for_pass_two,
                evaluation_meta_data,
                record_primary_timing,
            );

            total_result += pass_two_result;
        }

        total_result
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
    fn solve_rstar(
        self,
        rstar_t_dependence_evaluator: &mut RstarTDependenceEvaluator,
    ) -> RstarSolution<'a, T> {
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

        let t_dependent_solution = if rstar_t_dependence_evaluator.supports_t_derivatives() {
            let t_star = match &self
                .overlap_builder
                .counterterm_builder
                .transformed_kinematic_point
                .non_dual_cut_params()
                .tstar
            {
                DualOrNot::NonDual(t_star) => t_star.clone(),
                DualOrNot::Dual(_) => {
                    unreachable!("representative LU cut parameters must stay non-dual")
                }
            };

            Some(
                rstar_t_dependence_evaluator.evaluate(
                    &t_star,
                    &solution.solution,
                    &self.overlap_builder.rotated_center,
                    subspace,
                    self.overlap_builder
                        .counterterm_builder
                        .transformed_kinematic_point
                        .unrescaled_sample(),
                    masses,
                    self.esurface,
                    &self
                        .overlap_builder
                        .counterterm_builder
                        .graph
                        .loop_momentum_basis,
                    lmbs,
                ),
            )
        } else {
            None
        };

        RstarSolution {
            esurface_ct_builder: self,
            solution,
            t_dependent_solution,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::LUCounterTerm;
    use crate::processes::RaisedCutId;
    use typed_index_collections::{TiVec, ti_vec};

    #[test]
    fn inactive_cut_guard_reports_runtime_access() {
        let counterterm = LUCounterTerm {
            evaluators: TiVec::new(),
            thresholds: TiVec::new(),
            subspaces: TiVec::new(),
            active_cuts: ti_vec![false],
            active_left_thresholds: ti_vec![TiVec::new()],
            active_right_thresholds: ti_vec![TiVec::new()],
            active_iterated_thresholds: TiVec::new(),
        };

        let error = counterterm.ensure_active_cut(RaisedCutId(0)).unwrap_err();
        assert!(error.to_string().contains("generation marked it inactive"));
    }
}

#[cfg(test)]
mod tests {
    use super::LUCounterTerm;
    use crate::processes::RaisedCutId;
    use typed_index_collections::{TiVec, ti_vec};

    #[test]
    fn inactive_cut_guard_reports_runtime_access() {
        let counterterm = LUCounterTerm {
            evaluators: TiVec::new(),
            thresholds: TiVec::new(),
            subspaces: TiVec::new(),
            active_cuts: ti_vec![false],
            active_left_thresholds: ti_vec![TiVec::new()],
            active_right_thresholds: ti_vec![TiVec::new()],
            active_iterated_thresholds: TiVec::new(),
        };

        let error = counterterm.ensure_active_cut(RaisedCutId(0)).unwrap_err();
        assert!(error.to_string().contains("generation marked it inactive"));
    }
}

struct RstarSolution<'a, T: FloatLike> {
    esurface_ct_builder: EsurfaceCTBuilder<'a, T>,
    solution: NewtonIterationResult<T>,
    t_dependent_solution: Option<HyperDual<F<T>>>,
}

struct DualRstarGeometry<T: FloatLike> {
    radius: HyperDual<F<T>>,
    radius_star: HyperDual<F<T>>,
    esurface_derivative: HyperDual<F<T>>,
    rstar_loop_momenta: LoopMomenta<HyperDual<F<T>>>,
    external_moms: ExternalFourMomenta<HyperDual<F<T>>>,
}

impl<'a, T: FloatLike> RstarSolution<'a, T> {
    fn max_supported_order(&self) -> usize {
        self.t_dependent_solution
            .as_ref()
            .map(|dual| dual.values.len().saturating_sub(1))
            .unwrap_or(0)
    }

    fn base_rstar_loop_momenta(&self) -> LoopMomenta<F<T>> {
        let subspace = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .subspace;

        &self
            .esurface_ct_builder
            .overlap_builder
            .unit_shifted_momenta
            .rescale(&self.solution.solution, subspace.as_subspace_simple())
            + &self.esurface_ct_builder.overlap_builder.rotated_center
    }

    fn truncated_rstar_solution(&self, order: usize) -> HyperDual<F<T>> {
        let dual_solution = self
            .t_dependent_solution
            .as_ref()
            .expect("higher-order LU threshold evaluation requires cached r_star(t)");
        HyperDual::from_values(
            shape_for_t_derivatives(order),
            dual_solution.values[..=order].to_vec(),
        )
    }

    fn dual_geometry_for_order(&self, order: usize) -> DualRstarGeometry<T> {
        let source_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .transformed_kinematic_point
            .sample_for_order(order);
        let dual_loop_momenta = source_sample
            .sample
            .dual_loop_moms
            .as_ref()
            .expect("higher-order LU threshold evaluation requires dual loop momenta")
            .clone();

        let radius_star = self.truncated_rstar_solution(order);
        let center = dualize_loop_momenta(
            &radius_star,
            &self.esurface_ct_builder.overlap_builder.rotated_center,
        );
        let external_moms = dualize_external_momenta(&radius_star, source_sample.external_moms());

        let shifted_momenta = dual_loop_momenta
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum.clone() - center.clone())
            .collect::<LoopMomenta<_>>();

        let radius = dual_shifted_radius(
            &shifted_momenta,
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .subspace,
        );
        let inverse_radius = new_constant(&radius, &radius.values[0].one()) / radius.clone();
        let unit_shifted_momenta = shifted_momenta.rescale(
            &inverse_radius,
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .subspace
                .as_subspace_simple(),
        );

        let rstar_loop_momenta = unit_shifted_momenta
            .rescale(
                &radius_star,
                self.esurface_ct_builder
                    .overlap_builder
                    .counterterm_builder
                    .subspace
                    .as_subspace_simple(),
            )
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum.clone() + center.clone())
            .collect();

        let (_, esurface_derivative) = compute_self_and_r_derivative_subspace_dual(
            self.esurface_ct_builder.esurface,
            &radius_star,
            &unit_shifted_momenta,
            &center,
            &external_moms,
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .real_mass_vector,
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .subspace,
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .all_lmbs,
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .graph,
        );

        DualRstarGeometry {
            radius,
            radius_star,
            esurface_derivative,
            rstar_loop_momenta,
            external_moms,
        }
    }

    fn non_dual_multichanneling_factor(&self, rstar_sample: &MomentumSample<T>) -> Complex<F<T>> {
        let subspace = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .subspace;

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
                    .fold(rstar_sample.one(), |acc, value| acc * value)
            })
            .fold(rstar_sample.zero(), |acc, value| acc + value);

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
            .fold(rstar_sample.one(), |acc, value| acc * value);

        Complex::new_re(multichanneling_numerator / multi_channeling_denominator)
    }

    fn dual_multichanneling_factor(&self, geometry: &DualRstarGeometry<T>) -> HyperDual<F<T>> {
        let subspace = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .subspace;
        let lmb = subspace.get_lmb(
            self.esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .all_lmbs,
        );
        let zero = new_constant(&geometry.radius, &geometry.radius.values[0].zero());
        let one = new_constant(&geometry.radius, &geometry.radius.values[0].one());

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
                        let esurface_value = esurface.compute_from_dual_momenta(
                            lmb,
                            self.esurface_ct_builder
                                .overlap_builder
                                .counterterm_builder
                                .real_mass_vector,
                            &geometry.rstar_loop_momenta,
                            &geometry.external_moms,
                        );

                        esurface_value.clone() * esurface_value
                    })
                    .fold(one.clone(), |acc, value| acc * value)
            })
            .fold(zero.clone(), |acc, value| acc + value);

        let multi_channeling_numerator = self
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
                let esurface_value = esurface.compute_from_dual_momenta(
                    lmb,
                    self.esurface_ct_builder
                        .overlap_builder
                        .counterterm_builder
                        .real_mass_vector,
                    &geometry.rstar_loop_momenta,
                    &geometry.external_moms,
                );

                esurface_value.clone() * esurface_value
            })
            .fold(one, |acc, value| acc * value);

        multi_channeling_numerator / multi_channeling_denominator
    }

    fn rstar_sample_for_order<'solution>(
        &'solution self,
        order: usize,
    ) -> RstarSample<'solution, 'a, T> {
        let base_rstar_loop_momenta = self.base_rstar_loop_momenta();

        if order == 0 {
            let mut rstar_sample = self
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .transformed_kinematic_point
                .representative_sample()
                .clone();
            rstar_sample.sample.loop_moms = base_rstar_loop_momenta;
            rstar_sample.sample.dual_loop_moms = None;

            let radius = self.esurface_ct_builder.overlap_builder.radius.clone();
            let radius_star = self.solution.solution.clone();
            let esurface_derivative = self.solution.derivative_at_solution.clone();
            let e_cm = &self
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .e_cm;
            let uv_localisation_settings = &self
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .settings
                .subtraction
                .local_ct_settings
                .uv_localisation;

            let threshold_params = ThresholdParams {
                radius: DualOrNot::NonDual(radius.clone()),
                radius_star: DualOrNot::NonDual(radius_star.clone()),
                esurface_derivative: DualOrNot::NonDual(esurface_derivative.clone()),
                uv_damp_plus: DualOrNot::NonDual(evaluate_uv_damper(
                    &radius,
                    &radius_star,
                    e_cm,
                    uv_localisation_settings,
                )),
                uv_damp_minus: DualOrNot::NonDual(evaluate_uv_damper(
                    &-radius.clone(),
                    &radius_star,
                    e_cm,
                    uv_localisation_settings,
                )),
                h_function: DualOrNot::NonDual(evaluate_integrated_ct_normalisation(
                    &radius,
                    &radius_star,
                    e_cm,
                    &self
                        .esurface_ct_builder
                        .overlap_builder
                        .counterterm_builder
                        .settings
                        .subtraction
                        .integrated_ct_settings,
                )),
            };

            let value_of_multi_channeling_factor =
                DualOrNot::NonDual(self.non_dual_multichanneling_factor(&rstar_sample));

            return RstarSample {
                rstar_solution: self,
                rstar_sample,
                threshold_params,
                value_of_multi_channeling_factor,
            };
        }

        assert!(
            order <= self.max_supported_order(),
            "requested LU threshold derivative order {} but only {} orders are cached",
            order,
            self.max_supported_order()
        );

        let geometry = self.dual_geometry_for_order(order);
        let mut rstar_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .transformed_kinematic_point
            .sample_for_order(order)
            .clone();
        rstar_sample.sample.loop_moms = base_rstar_loop_momenta;
        rstar_sample.sample.dual_loop_moms = Some(geometry.rstar_loop_momenta.clone());

        let e_cm = &self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .e_cm;
        let uv_localisation_settings = &self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .settings
            .subtraction
            .local_ct_settings
            .uv_localisation;

        let threshold_params = ThresholdParams {
            radius: DualOrNot::Dual(geometry.radius.clone()),
            radius_star: DualOrNot::Dual(geometry.radius_star.clone()),
            esurface_derivative: DualOrNot::Dual(geometry.esurface_derivative.clone()),
            uv_damp_plus: DualOrNot::Dual(evaluate_uv_damper_dual(
                &geometry.radius,
                &geometry.radius_star,
                e_cm,
                uv_localisation_settings,
            )),
            uv_damp_minus: DualOrNot::Dual(evaluate_uv_damper_dual(
                &-geometry.radius.clone(),
                &geometry.radius_star,
                e_cm,
                uv_localisation_settings,
            )),
            h_function: DualOrNot::Dual(evaluate_integrated_ct_normalisation_dual(
                &geometry.radius,
                &geometry.radius_star,
                e_cm,
                &self
                    .esurface_ct_builder
                    .overlap_builder
                    .counterterm_builder
                    .settings
                    .subtraction
                    .integrated_ct_settings,
            )),
        };

        let value_of_multi_channeling_factor = DualOrNot::Dual(real_dual_to_complex(
            self.dual_multichanneling_factor(&geometry),
        ));

        RstarSample {
            rstar_solution: self,
            rstar_sample,
            threshold_params,
            value_of_multi_channeling_factor,
        }
    }
}

struct RstarSample<'solution, 'a, T: FloatLike> {
    rstar_solution: &'solution RstarSolution<'a, T>,
    rstar_sample: MomentumSample<T>,
    threshold_params: ThresholdParams<T>,
    value_of_multi_channeling_factor: DualOrNot<Complex<F<T>>>,
}

impl<'solution, 'a, T: FloatLike> RstarSample<'solution, 'a, T> {
    fn extract_threshold_parameters(&self, is_first_call: bool) -> ThresholdParams<T> {
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
            debug!("radius: {}", self.threshold_params.radius);
            debug!("radius_star: {}", self.threshold_params.radius_star);
            debug!(
                "esurface_derivative: {}",
                self.threshold_params.esurface_derivative
            );
            debug!("uv_damp_plus: {}", self.threshold_params.uv_damp_plus);
            debug!("uv_damp_minus: {}", self.threshold_params.uv_damp_minus);
            debug!("h_function: {}", self.threshold_params.h_function);
            debug!(
                "value of multi-channeling factor: {}",
                self.value_of_multi_channeling_factor
            );
        }

        self.threshold_params.clone()
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
    left_sample: &RstarSample<'_, '_, T>,
    right_sample: &RstarSample<'_, '_, T>,
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

    match (
        merged_sample.sample.dual_loop_moms.as_mut(),
        right_sample.rstar_sample.sample.dual_loop_moms.as_ref(),
    ) {
        (Some(merged_dual_loop_moms), Some(right_dual_loop_moms)) => {
            for lmb_index in right_subspace.iter_lmb_indices() {
                merged_dual_loop_moms[lmb_index] = right_dual_loop_moms[lmb_index].clone();
            }
        }
        (None, None) => {}
        _ => {
            unreachable!("iterated LU samples must either both carry dual loop momenta or neither")
        }
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
