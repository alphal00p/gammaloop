use core::f64;
use std::{collections::BTreeMap, path::Path};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use smallvec::{SmallVec, smallvec};
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
use symbolica::domains::{
    dual::{DualNumberStructure, HyperDual},
    float::{Real, RealLike},
};
use tracing::{debug, warn};
use typed_index_collections::TiVec;

use crate::{
    GammaLoopContext,
    cff::{
        CutCFFIndex,
        esurface::{
            Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId,
            esurface_value_is_strictly_inside,
        },
        expression::OrientationID,
    },
    graph::{Graph, LmbIndex, LoopMomentumBasis},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            GenericEvaluator, ParamBuilder, ThresholdParams,
            evaluators::{
                EvaluatorStack, SingleOrAllOrientations, evaluate_evaluator,
                evaluate_evaluator_single,
            },
            param_builder::LUParams,
            threshold_multiplier::{
                ThresholdMultiplierEvaluatorCollection, ThresholdMultiplierInputWorkspace,
            },
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
    observables::{
        GenericThresholdCountertermComponentWeight, ThresholdCountertermComponentOccurrence,
    },
    processes::{
        CutGroupId, EvaluatorBuildTimings, IteratedCtCollection, IteratedThresholdPieces,
        LUCounterTermData, LUThresholdHelperOutputs, LeftThresholdId, RightThresholdId,
        SingleThresholdPieces, ThresholdCountertermComponentKind,
        ThresholdCountertermMetadataRegistry, ThresholdCountertermSide,
        ThresholdCountertermVariantId, build_derivative_structure,
    },
    settings::{GlobalSettings, RuntimeSettings, global::FrozenCompilationMode},
    subtraction::{
        RstarTDependenceEvaluator, RstarTDependenceInput, evaluate_integrated_ct_normalisation,
        evaluate_integrated_ct_normalisation_dual, evaluate_uv_damper, evaluate_uv_damper_dual,
        overlap_subspace::{self, OverlapGroup, OverlapInput, OverlapStructure},
    },
    utils::{
        F, FloatLike,
        hyperdual_utils::{
            DualOrNot, dualize_dual_t_to_dual_r_t, extract_coefficient_t_duals,
            extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            shape_from_cut_cff_index, simple_n_deriv_shape, variable_indices_from_cut_cff_index,
        },
        newton_solver::{NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity},
    },
};

fn zero_dual_or_not_complex<T: FloatLike>(order: usize, zero: &F<T>) -> DualOrNot<Complex<F<T>>> {
    if order == 0 {
        DualOrNot::NonDual(Complex::new_re(zero.clone()))
    } else {
        DualOrNot::Dual(HyperDual::new(simple_n_deriv_shape(order)))
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

struct ThresholdHelperEvaluation<T, P> {
    legacy: Option<T>,
    pieces: Option<P>,
}

type ThresholdHelperValue<T> = DualOrNot<Complex<F<T>>>;
type SingleThresholdHelperEvaluation<T> = ThresholdHelperEvaluation<
    ThresholdHelperValue<T>,
    SingleThresholdPieces<ThresholdHelperValue<T>>,
>;
type IteratedThresholdHelperEvaluation<T> = ThresholdHelperEvaluation<
    ThresholdHelperValue<T>,
    IteratedThresholdPieces<ThresholdHelperValue<T>>,
>;

struct ThresholdHelperEvaluationContext<'a, T: FloatLike> {
    helper: &'a mut GenericEvaluator,
    cut_cff_index: &'a CutCFFIndex,
    threshold_result: &'a ThresholdHelperValue<T>,
    outputs: LUThresholdHelperOutputs,
    evaluation_meta_data: &'a mut EvaluationMetaData,
    record_primary_timing: bool,
}

impl<T, P> ThresholdHelperEvaluation<T, P> {
    fn into_legacy(self) -> T {
        self.legacy
            .expect("legacy LU threshold helper output is unavailable")
    }

    fn into_pieces(self) -> P {
        self.pieces
            .expect("split LU threshold helper outputs are unavailable")
    }
}

#[derive(Clone, Debug)]
struct SingleMultiplierCoefficients<T: FloatLike> {
    local: F<T>,
    integrated: F<T>,
}

#[derive(Clone, Copy, Debug)]
enum SingleMultiplierComponent {
    Local,
    Integrated,
}

fn single_multiplier_context<'a, S>(
    component: SingleMultiplierComponent,
    base: &'a S,
    root: &'a S,
) -> (&'a S, &'a S) {
    // The absent side's threshold-rescaling map is the identity: its coordinates stay at `base`,
    // while the original integrand factor on that side remains part of the expanded term.
    match component {
        SingleMultiplierComponent::Local => (base, root),
        SingleMultiplierComponent::Integrated => (root, root),
    }
}

impl<T: FloatLike> SingleMultiplierCoefficients<T> {
    fn all_zero(&self) -> bool {
        self.local.is_zero() && self.integrated.is_zero()
    }
}

#[derive(Clone, Debug)]
struct IteratedMultiplierCoefficients<T: FloatLike> {
    local_local: PairMultiplierCoefficient<T>,
    local_integrated: PairMultiplierCoefficient<T>,
    integrated_local: PairMultiplierCoefficient<T>,
    integrated_integrated: PairMultiplierCoefficient<T>,
}

#[derive(Clone, Debug)]
struct PairMultiplierCoefficient<T: FloatLike> {
    left: F<T>,
    right: F<T>,
}

impl<T: FloatLike> PairMultiplierCoefficient<T> {
    fn effective(&self) -> F<T> {
        &self.left * &self.right
    }

    fn is_zero(&self) -> bool {
        self.left.is_zero() || self.right.is_zero()
    }
}

#[derive(Clone, Copy, Debug)]
enum IteratedMultiplierComponent {
    LocalLocal,
    LocalIntegrated,
    IntegratedLocal,
    IntegratedIntegrated,
}

fn iterated_multiplier_context<'a, S>(
    component: IteratedMultiplierComponent,
    base: &'a S,
    left_root: &'a S,
    right_root: &'a S,
    merged_root: &'a S,
) -> (&'a S, &'a S) {
    match component {
        IteratedMultiplierComponent::LocalLocal => (base, merged_root),
        IteratedMultiplierComponent::LocalIntegrated => (right_root, merged_root),
        IteratedMultiplierComponent::IntegratedLocal => (left_root, merged_root),
        IteratedMultiplierComponent::IntegratedIntegrated => (merged_root, merged_root),
    }
}

impl<T: FloatLike> IteratedMultiplierCoefficients<T> {
    fn all_zero(&self) -> bool {
        self.local_local.is_zero()
            && self.local_integrated.is_zero()
            && self.integrated_local.is_zero()
            && self.integrated_integrated.is_zero()
    }
}

#[allow(clippy::too_many_arguments)]
fn evaluate_single_multiplier_context<T: FloatLike>(
    collection: &mut ThresholdMultiplierEvaluatorCollection,
    workspace: &mut ThresholdMultiplierInputWorkspace<T>,
    variant_id: ThresholdCountertermVariantId,
    generation_lmb: &LoopMomentumBasis,
    masses: &EdgeVec<F<T>>,
    model_values: &[Complex<F<f64>>],
    additional_values: &[F<T>],
    effective: &MomentumSample<T>,
    star: &MomentumSample<T>,
    evaluation_meta_data: &mut EvaluationMetaData,
    record_primary_timing: bool,
) -> Result<F<T>> {
    let Some(evaluator_id) = collection.evaluator_id_for_variant(variant_id)? else {
        return Ok(effective.sample.one());
    };
    let values = workspace.bind(
        collection.layout(),
        generation_lmb,
        masses,
        model_values,
        additional_values,
        effective,
        star,
    )?;
    collection.evaluators_mut()[evaluator_id.0].evaluate(
        values,
        evaluation_meta_data,
        record_primary_timing,
    )
}

#[allow(clippy::too_many_arguments)]
fn evaluate_pair_multiplier_context<T: FloatLike>(
    collection: &mut ThresholdMultiplierEvaluatorCollection,
    workspace: &mut ThresholdMultiplierInputWorkspace<T>,
    left_variant_id: ThresholdCountertermVariantId,
    right_variant_id: ThresholdCountertermVariantId,
    generation_lmb: &LoopMomentumBasis,
    masses: &EdgeVec<F<T>>,
    model_values: &[Complex<F<f64>>],
    additional_values: &[F<T>],
    effective: &MomentumSample<T>,
    star: &MomentumSample<T>,
    evaluation_meta_data: &mut EvaluationMetaData,
    record_primary_timing: bool,
) -> Result<(F<T>, F<T>)> {
    let left_evaluator_id = collection.evaluator_id_for_variant(left_variant_id)?;
    let right_evaluator_id = collection.evaluator_id_for_variant(right_variant_id)?;
    if left_evaluator_id.is_none() && right_evaluator_id.is_none() {
        return Ok((effective.sample.one(), effective.sample.one()));
    }

    let values = workspace.bind(
        collection.layout(),
        generation_lmb,
        masses,
        model_values,
        additional_values,
        effective,
        star,
    )?;
    if let (Some(left_evaluator_id), Some(right_evaluator_id)) =
        (left_evaluator_id, right_evaluator_id)
        && left_evaluator_id == right_evaluator_id
    {
        // Expression interning gives equal functions one evaluator ID. Both sides receive the
        // same bound input slice, so evaluate the shared function once and reuse its value.
        let value = collection.evaluators_mut()[left_evaluator_id.0].evaluate(
            values,
            evaluation_meta_data,
            record_primary_timing,
        )?;
        return Ok((value.clone(), value));
    }
    let left = if let Some(evaluator_id) = left_evaluator_id {
        collection.evaluators_mut()[evaluator_id.0].evaluate(
            values,
            evaluation_meta_data,
            record_primary_timing,
        )?
    } else {
        effective.sample.one()
    };
    // Deliberately evaluate the right factor even when the left factor is exactly zero. This
    // keeps finite/real validation symmetric and makes one global context authoritative.
    let right = if let Some(evaluator_id) = right_evaluator_id {
        collection.evaluators_mut()[evaluator_id.0].evaluate(
            values,
            evaluation_meta_data,
            record_primary_timing,
        )?
    } else {
        effective.sample.one()
    };
    Ok((left, right))
}

fn weight_single_multiplier_pieces<T: FloatLike>(
    pieces: SingleThresholdPieces<DualOrNot<Complex<F<T>>>>,
    coefficients: &SingleMultiplierCoefficients<T>,
    order: usize,
    zero: &F<T>,
) -> DualOrNot<Complex<F<T>>> {
    let mut result = zero_dual_or_not_complex(order, zero);
    if !coefficients.local.is_zero() {
        result += multiply_dual_or_not_complex(
            pieces.local,
            &DualOrNot::NonDual(Complex::new_re(coefficients.local.clone())),
        );
    }
    if !coefficients.integrated.is_zero() {
        result += multiply_dual_or_not_complex(
            pieces.integrated,
            &DualOrNot::NonDual(Complex::new_re(coefficients.integrated.clone())),
        );
    }
    result
}

fn weight_iterated_multiplier_pieces<T: FloatLike>(
    pieces: IteratedThresholdPieces<DualOrNot<Complex<F<T>>>>,
    coefficients: &IteratedMultiplierCoefficients<T>,
    order: usize,
    zero: &F<T>,
) -> DualOrNot<Complex<F<T>>> {
    let mut result = zero_dual_or_not_complex(order, zero);
    if !coefficients.local_local.is_zero() {
        result += multiply_dual_or_not_complex(
            pieces.local_local,
            &DualOrNot::NonDual(Complex::new_re(coefficients.local_local.effective())),
        );
    }
    if !coefficients.local_integrated.is_zero() {
        result += multiply_dual_or_not_complex(
            pieces.local_integrated,
            &DualOrNot::NonDual(Complex::new_re(coefficients.local_integrated.effective())),
        );
    }
    if !coefficients.integrated_local.is_zero() {
        result += multiply_dual_or_not_complex(
            pieces.integrated_local,
            &DualOrNot::NonDual(Complex::new_re(coefficients.integrated_local.effective())),
        );
    }
    if !coefficients.integrated_integrated.is_zero() {
        result += multiply_dual_or_not_complex(
            pieces.integrated_integrated,
            &DualOrNot::NonDual(Complex::new_re(
                coefficients.integrated_integrated.effective(),
            )),
        );
    }
    result
}

fn plain_t_dual_or_scalar_complex<T: FloatLike>(
    value: &DualOrNot<Complex<F<T>>>,
    t_variable: Option<usize>,
) -> DualOrNot<Complex<F<T>>> {
    match value {
        DualOrNot::Dual(dual) => {
            if let Some(t_variable) = t_variable {
                DualOrNot::Dual(extract_zero_threshold_coefficient_t_dual(dual, t_variable))
            } else {
                DualOrNot::NonDual(dual.values[0].clone())
            }
        }
        DualOrNot::NonDual(value) => DualOrNot::NonDual(value.clone()),
    }
}

fn extract_zero_threshold_coefficient_t_dual<
    T: Clone + spenso::algebra::algebraic_traits::RefZero + Default,
>(
    dual: &HyperDual<T>,
    t_variable: usize,
) -> HyperDual<T> {
    let (coefficient_orders, coefficient_duals) = extract_coefficient_t_duals(dual, t_variable);
    let zero_index = coefficient_orders
        .iter()
        .position(|orders| orders.iter().all(|order| *order == 0))
        .expect("Could not find zero-threshold coefficient in mixed dual");

    coefficient_duals[zero_index].clone()
}

fn embed_t_dual_in_target_shape<T: FloatLike>(
    t_dual: &HyperDual<F<T>>,
    target_shape: &HyperDual<F<T>>,
    t_variable: Option<usize>,
) -> HyperDual<F<T>> {
    match t_variable {
        Some(t_variable) => {
            dualize_dual_t_to_dual_r_t(t_dual.clone(), target_shape.clone(), t_variable)
        }
        None => new_constant(target_shape, &t_dual.values[0]),
    }
}

fn activate_threshold_variable_in_target_shape<T: FloatLike>(
    dual: &mut HyperDual<F<T>>,
    variable: Option<usize>,
) {
    let Some(variable) = variable else {
        return;
    };

    let n_variables = dual.get_shape()[0].len();
    let mut derivative_shape = vec![0; n_variables];
    derivative_shape[variable] = 1;

    let derivative_index = dual
        .get_shape()
        .iter()
        .position(|shape| shape == &derivative_shape)
        .expect("Could not find threshold-variable derivative slot in target shape");

    dual.values[derivative_index] = dual.values[0].one();
}

fn embed_real_dual_or_not_in_target_shape<T: FloatLike>(
    value: &DualOrNot<F<T>>,
    target_shape: &Option<HyperDual<F<T>>>,
    t_variable: Option<usize>,
) -> DualOrNot<F<T>> {
    match (value, target_shape) {
        (DualOrNot::Dual(dual), Some(target_shape)) => {
            DualOrNot::Dual(embed_t_dual_in_target_shape(dual, target_shape, t_variable))
        }
        (DualOrNot::NonDual(value), Some(target_shape)) => {
            DualOrNot::Dual(new_constant(target_shape, value))
        }
        _ => value.clone(),
    }
}

fn embedded_dual_loop_momenta_for_cut_cff_index<T: FloatLike>(
    sample: &MomentumSample<T>,
    target_shape: &Option<HyperDual<F<T>>>,
    t_variable: Option<usize>,
) -> Option<LoopMomenta<HyperDual<F<T>>>> {
    match target_shape {
        Some(target_shape) => Some(match sample.sample.dual_loop_moms.as_ref() {
            Some(dual_loop_moms) => dual_loop_moms
                .iter()
                .map(|momentum| {
                    momentum.map_ref(&|component| {
                        embed_t_dual_in_target_shape(component, target_shape, t_variable)
                    })
                })
                .collect(),
            None => sample
                .loop_moms()
                .iter()
                .map(|momentum| {
                    momentum.map_ref(&|component| new_constant(target_shape, component))
                })
                .collect(),
        }),
        None => sample.sample.dual_loop_moms.clone(),
    }
}

fn extend_plain_helper_real_params<T: FloatLike>(
    params: &mut Vec<Complex<F<T>>>,
    value: &DualOrNot<F<T>>,
    t_variable: Option<usize>,
) {
    match value {
        DualOrNot::Dual(dual) => {
            if let Some(t_variable) = t_variable {
                params.extend(
                    extract_t_derivatives(extract_zero_threshold_coefficient_t_dual(
                        dual, t_variable,
                    ))
                    .into_iter()
                    .map(Complex::new_re),
                );
            } else {
                params.push(Complex::new_re(dual.values[0].clone()));
            }
        }
        DualOrNot::NonDual(value) => params.push(Complex::new_re(value.clone())),
    }
}

fn extend_extracted_f_params<T: FloatLike>(
    params: &mut Vec<Complex<F<T>>>,
    value: &DualOrNot<Complex<F<T>>>,
    t_variable: Option<usize>,
) {
    match value {
        DualOrNot::Dual(dual) => {
            if let Some(t_variable) = t_variable {
                let (_, coefficient_duals) = extract_coefficient_t_duals(dual, t_variable);
                for coefficient_dual in coefficient_duals {
                    params.extend(extract_t_derivatives_complex(coefficient_dual));
                }
            } else {
                params.extend(dual.values.iter().cloned());
            }
        }
        DualOrNot::NonDual(value) => params.push(value.clone()),
    }
}

fn extend_extracted_eta_params<T: FloatLike>(
    params: &mut Vec<Complex<F<T>>>,
    value: &DualOrNot<F<T>>,
    lu_variable: Option<usize>,
    threshold_variable: Option<usize>,
) {
    match value {
        DualOrNot::Dual(dual) => {
            if let Some(lu_variable) = lu_variable {
                let (coefficient_orders, coefficient_duals) =
                    extract_coefficient_t_duals(dual, lu_variable);
                let projected_threshold_variable = threshold_variable.map(|variable| {
                    if variable > lu_variable {
                        variable - 1
                    } else {
                        variable
                    }
                });
                for (orders, coefficient_dual) in
                    coefficient_orders.into_iter().zip(coefficient_duals)
                {
                    if orders.iter().enumerate().any(|(variable, &order)| {
                        Some(variable) != projected_threshold_variable && order != 0
                    }) {
                        continue;
                    }
                    params.extend(
                        extract_t_derivatives(coefficient_dual)
                            .into_iter()
                            .map(Complex::new_re),
                    );
                }
            } else {
                params.extend(
                    dual.get_shape()
                        .iter()
                        .zip(&dual.values)
                        .filter(|(orders, _)| {
                            !orders.iter().enumerate().any(|(variable, &order)| {
                                Some(variable) != threshold_variable && order != 0
                            })
                        })
                        .map(|(_, value)| Complex::new_re(value.clone())),
                );
            }
        }
        DualOrNot::NonDual(value) => params.push(Complex::new_re(value.clone())),
    }
}

fn evaluate_threshold_helper_single<
    T: FloatLike + crate::integrands::process::evaluators::GenericEvaluatorFloat,
>(
    context: ThresholdHelperEvaluationContext<'_, T>,
    threshold_params: &ThresholdParams<T>,
) -> SingleThresholdHelperEvaluation<T> {
    let ThresholdHelperEvaluationContext {
        helper,
        cut_cff_index,
        threshold_result,
        outputs,
        evaluation_meta_data,
        record_primary_timing,
    } = context;
    let variable_indices = variable_indices_from_cut_cff_index(cut_cff_index);
    let t_variable = variable_indices.lu_cut;
    let mut helper_params = Vec::new();
    extend_extracted_f_params(&mut helper_params, threshold_result, t_variable);
    extend_extracted_eta_params(
        &mut helper_params,
        &threshold_params.esurface_derivative,
        t_variable,
        variable_indices
            .left_threshold
            .or(variable_indices.right_threshold),
    );
    extend_plain_helper_real_params(&mut helper_params, &threshold_params.radius, t_variable);
    extend_plain_helper_real_params(
        &mut helper_params,
        &threshold_params.radius_star,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &threshold_params.uv_damp_plus,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &threshold_params.uv_damp_minus,
        t_variable,
    );
    extend_plain_helper_real_params(&mut helper_params, &threshold_params.h_function, t_variable);

    debug!(
        "LU threshold helper input (single): cut_cff_index={:?}, t_variable={:?}, param_count={}",
        cut_cff_index,
        t_variable,
        helper_params.len()
    );
    for (idx, value) in helper_params.iter().enumerate() {
        debug!(
            "LU threshold helper param (single): cut_cff_index={:?}, index={}, value={}",
            cut_cff_index, idx, value
        );
    }

    let mut result = evaluate_evaluator(
        helper,
        &helper_params,
        evaluation_meta_data,
        record_primary_timing,
    );

    match outputs {
        LUThresholdHelperOutputs::Legacy => {
            assert_eq!(
                result.len(),
                1,
                "legacy single LU threshold helper must return one combined output",
            );
            ThresholdHelperEvaluation {
                legacy: result.pop(),
                pieces: None,
            }
        }
        LUThresholdHelperOutputs::Pieces => {
            assert_eq!(
                result.len(),
                2,
                "single LU threshold helper must return local and integrated pieces",
            );
            let integrated = result.pop().unwrap();
            let local = result.pop().unwrap();
            ThresholdHelperEvaluation {
                legacy: None,
                pieces: Some(SingleThresholdPieces { local, integrated }),
            }
        }
        LUThresholdHelperOutputs::LegacyAndPieces => {
            assert_eq!(
                result.len(),
                3,
                "recording single LU threshold helper must return legacy, local, and integrated outputs",
            );
            let integrated = result.pop().unwrap();
            let local = result.pop().unwrap();
            let legacy = result.pop().unwrap();
            ThresholdHelperEvaluation {
                legacy: Some(legacy),
                pieces: Some(SingleThresholdPieces { local, integrated }),
            }
        }
    }
}

fn evaluate_threshold_helper_iterated<
    T: FloatLike + crate::integrands::process::evaluators::GenericEvaluatorFloat,
>(
    context: ThresholdHelperEvaluationContext<'_, T>,
    left_threshold_params: &ThresholdParams<T>,
    right_threshold_params: &ThresholdParams<T>,
) -> IteratedThresholdHelperEvaluation<T> {
    let ThresholdHelperEvaluationContext {
        helper,
        cut_cff_index,
        threshold_result,
        outputs,
        evaluation_meta_data,
        record_primary_timing,
    } = context;
    let variable_indices = variable_indices_from_cut_cff_index(cut_cff_index);
    let t_variable = variable_indices.lu_cut;
    let mut helper_params = Vec::new();
    extend_extracted_f_params(&mut helper_params, threshold_result, t_variable);
    extend_extracted_eta_params(
        &mut helper_params,
        &left_threshold_params.esurface_derivative,
        t_variable,
        variable_indices.left_threshold,
    );
    extend_extracted_eta_params(
        &mut helper_params,
        &right_threshold_params.esurface_derivative,
        t_variable,
        variable_indices.right_threshold,
    );

    extend_plain_helper_real_params(
        &mut helper_params,
        &left_threshold_params.radius,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &left_threshold_params.radius_star,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &left_threshold_params.uv_damp_plus,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &left_threshold_params.uv_damp_minus,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &left_threshold_params.h_function,
        t_variable,
    );

    extend_plain_helper_real_params(
        &mut helper_params,
        &right_threshold_params.radius,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &right_threshold_params.radius_star,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &right_threshold_params.uv_damp_plus,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &right_threshold_params.uv_damp_minus,
        t_variable,
    );
    extend_plain_helper_real_params(
        &mut helper_params,
        &right_threshold_params.h_function,
        t_variable,
    );

    debug!(
        "LU threshold helper input (iterated): cut_cff_index={:?}, t_variable={:?}, param_count={}",
        cut_cff_index,
        t_variable,
        helper_params.len()
    );
    for (idx, value) in helper_params.iter().enumerate() {
        debug!(
            "LU threshold helper param (iterated): cut_cff_index={:?}, index={}, value={}",
            cut_cff_index, idx, value
        );
    }

    let mut result = evaluate_evaluator(
        helper,
        &helper_params,
        evaluation_meta_data,
        record_primary_timing,
    );

    match outputs {
        LUThresholdHelperOutputs::Legacy => {
            assert_eq!(
                result.len(),
                1,
                "legacy iterated LU threshold helper must return one combined output",
            );
            ThresholdHelperEvaluation {
                legacy: result.pop(),
                pieces: None,
            }
        }
        LUThresholdHelperOutputs::Pieces | LUThresholdHelperOutputs::LegacyAndPieces => {
            let expected = if outputs == LUThresholdHelperOutputs::Pieces {
                4
            } else {
                5
            };
            assert_eq!(
                result.len(),
                expected,
                "iterated LU threshold helper must return the configured legacy and/or LL, LI, IL, and II outputs",
            );
            let integrated_integrated = result.pop().unwrap();
            let integrated_local = result.pop().unwrap();
            let local_integrated = result.pop().unwrap();
            let local_local = result.pop().unwrap();
            let legacy = (outputs == LUThresholdHelperOutputs::LegacyAndPieces)
                .then(|| result.pop().unwrap());
            ThresholdHelperEvaluation {
                legacy,
                pieces: Some(IteratedThresholdPieces {
                    local_local,
                    local_integrated,
                    integrated_local,
                    integrated_integrated,
                }),
            }
        }
    }
}

fn real_dual_to_complex<T: FloatLike>(dual: HyperDual<F<T>>) -> HyperDual<Complex<F<T>>> {
    let shape = dual.get_shape().iter().map(|v| v.to_vec()).collect();
    let values = dual.values.into_iter().map(Complex::new_re).collect_vec();
    HyperDual::from_values(shape, values)
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
pub struct LUThresholdHelperEvaluators {
    pub left_thresholds: TiVec<LeftThresholdId, BTreeMap<CutCFFIndex, GenericEvaluator>>,
    pub right_thresholds: TiVec<RightThresholdId, BTreeMap<CutCFFIndex, GenericEvaluator>>,
    pub iterated: IteratedCtCollection<BTreeMap<CutCFFIndex, GenericEvaluator>>,
}

impl LUThresholdHelperEvaluators {
    fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> color_eyre::Result<()>,
    ) -> color_eyre::Result<()> {
        for evaluators in self.left_thresholds.iter_mut() {
            for evaluator in evaluators.values_mut() {
                f(evaluator)?;
            }
        }
        for evaluators in self.right_thresholds.iter_mut() {
            for evaluator in evaluators.values_mut() {
                f(evaluator)?;
            }
        }
        for evaluators in self.iterated.iter_mut() {
            for evaluator in evaluators.values_mut() {
                f(evaluator)?;
            }
        }

        Ok(())
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LUCounterTermEvaluators {
    pub left_thresholds_evaluator: TiVec<LeftThresholdId, BTreeMap<CutCFFIndex, EvaluatorStack>>,
    pub right_thresholds_evaluator: TiVec<RightThresholdId, BTreeMap<CutCFFIndex, EvaluatorStack>>,
    pub iterated_evaluator: IteratedCtCollection<BTreeMap<CutCFFIndex, EvaluatorStack>>,
    pub threshold_helpers: LUThresholdHelperEvaluators,
    pub threshold_multipliers: Option<ThresholdMultiplierEvaluatorCollection>,
    pub residue_from_e_surface_evaluators: Vec<GenericEvaluator>,
}

impl LUCounterTermEvaluators {
    pub(crate) fn generic_compileable_evaluator_count(&self) -> usize {
        let left = self
            .left_thresholds_evaluator
            .iter()
            .flat_map(|evaluators| evaluators.values())
            .map(EvaluatorStack::generic_evaluator_count)
            .sum::<usize>();
        let right = self
            .right_thresholds_evaluator
            .iter()
            .flat_map(|evaluators| evaluators.values())
            .map(EvaluatorStack::generic_evaluator_count)
            .sum::<usize>();
        let iterated = self
            .iterated_evaluator
            .iter()
            .flat_map(|evaluators| evaluators.values())
            .map(EvaluatorStack::generic_evaluator_count)
            .sum::<usize>();

        // Ignore eager-only helper and multiplier evaluators, as well as pass-two evaluators.
        left + right + iterated
    }

    fn evaluate_residue_from_esurface<T: FloatLike>(
        &mut self,
        cut_group_id: CutGroupId,
        order: usize,
        pass_one_result: DualOrNot<Complex<F<T>>>,
        cut_esurface: &DualOrNot<F<T>>,
        evaluation_meta_data: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Complex<F<T>> {
        let mut params_for_pass_two = vec![];
        match pass_one_result {
            DualOrNot::Dual(dual_result) => {
                params_for_pass_two.extend_from_slice(&extract_t_derivatives_complex(dual_result));
            }
            DualOrNot::NonDual(non_dual_result) => {
                params_for_pass_two.push(non_dual_result);
            }
        }

        match cut_esurface {
            DualOrNot::Dual(dual_e_surface) => {
                extract_t_derivatives(dual_e_surface.clone())[1..]
                    .iter()
                    .for_each(|value| params_for_pass_two.push(Complex::new_re(value.clone())));
            }
            DualOrNot::NonDual(non_dual_e_surface) => {
                debug!("non dual esurface: {}", non_dual_e_surface);
                params_for_pass_two.push(Complex::new_re(non_dual_e_surface.clone()));
            }
        }

        debug!(
            "LU pass-two evaluator input: cut_group_id={}, lu_order={}, param_count={}",
            cut_group_id.0,
            order + 1,
            params_for_pass_two.len()
        );
        for (idx, value) in params_for_pass_two.iter().enumerate() {
            debug!(
                "LU pass-two evaluator param: cut_group_id={}, lu_order={}, index={}, value={}",
                cut_group_id.0,
                order + 1,
                idx,
                value
            );
        }

        evaluate_evaluator_single(
            &mut self.residue_from_e_surface_evaluators[order],
            &params_for_pass_two,
            evaluation_meta_data,
            record_primary_timing,
        )
    }

    pub fn from_atoms(
        counterterm_data: &LUCounterTermData,
        max_cut_order: usize,
        threshold_helpers: LUThresholdHelperEvaluators,
        threshold_multipliers: Option<ThresholdMultiplierEvaluatorCollection>,
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
                    .map(|(i, atom)| {
                        let dual_shape = shape_from_cut_cff_index(i);

                        let (evaluator, evaluator_timings) = EvaluatorStack::new_with_timings(
                            std::slice::from_ref(atom),
                            param_builder,
                            &orientations.raw,
                            dual_shape,
                            &settings.generation.evaluator,
                        )
                        .unwrap();
                        timings += evaluator_timings;
                        (*i, evaluator)
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
                    .map(|(i, atom)| {
                        let dual_shape = shape_from_cut_cff_index(i);

                        let (evaluator, evaluator_timings) = EvaluatorStack::new_with_timings(
                            std::slice::from_ref(atom),
                            param_builder,
                            &orientations.raw,
                            dual_shape,
                            &settings.generation.evaluator,
                        )
                        .unwrap();
                        timings += evaluator_timings;
                        (*i, evaluator)
                    })
                    .collect()
            })
            .collect();

        let iterated_timings = std::cell::Cell::new(EvaluatorBuildTimings::default());
        let iterated_evaluator = counterterm_data.iterated.map_ref(|parametric_integrands| {
            parametric_integrands
                .integrands
                .iter()
                .map(|(i, atom)| {
                    let dual_shape = shape_from_cut_cff_index(i);

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
                    (*i, evaluator)
                })
                .collect()
        });
        timings += iterated_timings.get();

        let symbolica_started = std::time::Instant::now();
        let pass_two_evaluator = (1..=max_cut_order)
            .map(|order| {
                build_derivative_structure(order as u8, -1, &settings.generation.evaluator)
            })
            .collect();

        timings.symbolica_time += symbolica_started.elapsed();

        if let Some(multipliers) = &threshold_multipliers {
            multipliers
                .validate(
                    counterterm_data.left_thresholds.len(),
                    counterterm_data.right_thresholds.len(),
                )
                .expect("invalid threshold-multiplier registry generated for LU counterterms");
        }

        (
            LUCounterTermEvaluators {
                left_thresholds_evaluator,
                right_thresholds_evaluator,
                iterated_evaluator,
                threshold_helpers,
                threshold_multipliers,
                residue_from_e_surface_evaluators: pass_two_evaluator,
            },
            timings,
        )
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        cut_group_id: CutGroupId,
        frozen_mode: &FrozenCompilationMode,
    ) -> color_eyre::Result<()> {
        for (threshold_id, evaluators) in self.left_thresholds_evaluator.iter_mut_enumerated() {
            for (index, evaluator) in evaluators.iter_mut() {
                let name = format!(
                    "cut_group_{}_left_threshold_{}_index_{}",
                    cut_group_id.0, threshold_id.0, index
                );
                evaluator.compile(&name, path.as_ref(), frozen_mode)?;
            }
        }

        for (threshold_id, evaluators) in self.right_thresholds_evaluator.iter_mut_enumerated() {
            for (index, evaluator) in evaluators.iter_mut() {
                let name = format!(
                    "cut_group_{}_right_threshold_{}_index_{}",
                    cut_group_id.0, threshold_id.0, index
                );
                evaluator.compile(&name, path.as_ref(), frozen_mode)?;
            }
        }

        for (iterated_index, evaluators) in self.iterated_evaluator.iter_mut().enumerate() {
            for (index, evaluator) in evaluators.iter_mut() {
                let name = format!(
                    "cut_group_{}_iterated_{}_index_{}",
                    cut_group_id.0, iterated_index, index
                );
                evaluator.compile(&name, path.as_ref(), frozen_mode)?;
            }
        }

        for (order, pass_to_evaluator) in self
            .residue_from_e_surface_evaluators
            .iter_mut()
            .enumerate()
        {
            let name = format!("cut_group_{}_pass_two_{}", cut_group_id.0, order);
            pass_to_evaluator.compile_external(
                path.as_ref().join(&name).with_extension("cpp"),
                &name,
                path.as_ref().join(&name).with_extension("so"),
                frozen_mode,
            )?;
        }

        Ok(())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> color_eyre::Result<()>,
    ) -> color_eyre::Result<()> {
        for evaluators in self.left_thresholds_evaluator.iter_mut() {
            for evaluator in evaluators.values_mut() {
                evaluator.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        for evaluators in self.right_thresholds_evaluator.iter_mut() {
            for evaluator in evaluators.values_mut() {
                evaluator.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        for evaluators in self.iterated_evaluator.iter_mut() {
            for evaluator in evaluators.values_mut() {
                evaluator.for_each_generic_evaluator_mut(&mut f)?;
            }
        }
        self.threshold_helpers
            .for_each_generic_evaluator_mut(&mut f)?;
        if let Some(multipliers) = &mut self.threshold_multipliers {
            multipliers.for_each_generic_evaluator_mut(&mut f)?;
        }
        for pass_to_evaluator in self.residue_from_e_surface_evaluators.iter_mut() {
            f(pass_to_evaluator)?;
        }

        Ok(())
    }
}

type CutThresholds = (
    TiVec<LeftThresholdId, Esurface>,
    TiVec<RightThresholdId, Esurface>,
);

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct LUVariantSubspaces {
    pub left_variant_ids: TiVec<LeftThresholdId, ThresholdCountertermVariantId>,
    pub right_variant_ids: TiVec<RightThresholdId, ThresholdCountertermVariantId>,
    pub left: TiVec<LeftThresholdId, SubspaceData>,
    pub right: TiVec<RightThresholdId, SubspaceData>,
    pub left_common: Option<SubspaceData>,
    pub right_common: Option<SubspaceData>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct LUCounterTerm {
    pub evaluators: TiVec<CutGroupId, LUCounterTermEvaluators>,
    pub thresholds: TiVec<CutGroupId, CutThresholds>,
    pub subspaces: TiVec<CutGroupId, (SubspaceData, SubspaceData)>,
    /// Present only when directives require per-variant projected E-surface instances.
    pub variant_subspaces: Option<TiVec<CutGroupId, LUVariantSubspaces>>,
    /// Stable graph-local static metadata. The absent/no-directive path retains no allocation.
    pub metadata_registry: Option<ThresholdCountertermMetadataRegistry>,
    pub active_cut_groups: TiVec<CutGroupId, bool>,
    pub active_left_thresholds: TiVec<CutGroupId, TiVec<LeftThresholdId, bool>>,
    pub active_right_thresholds: TiVec<CutGroupId, TiVec<RightThresholdId, bool>>,
    pub active_iterated_thresholds: TiVec<CutGroupId, IteratedCtCollection<bool>>,
    pub rstar_dependence_calculator: TiVec<CutGroupId, RstarTDependenceEvaluator>,
}

pub(crate) struct LUCountertermEvaluation<T: FloatLike> {
    pub total: Complex<F<T>>,
    pub components: Option<Vec<GenericThresholdCountertermComponentWeight<T>>>,
}

struct PendingLUCountertermComponent<T: FloatLike> {
    component_id: usize,
    occurrence: ThresholdCountertermComponentOccurrence,
    multiplier_values: SmallVec<[F<T>; 2]>,
    effective_multiplier: F<T>,
    /// `None` is reserved for an exact-zero multiplier and avoids pass-one/helper/pass-two work.
    pass_one: Option<DualOrNot<Complex<F<T>>>>,
}

impl<T: FloatLike> PendingLUCountertermComponent<T> {
    fn finish(
        self,
        zero: &F<T>,
        evaluate_pass_two: impl FnOnce(DualOrNot<Complex<F<T>>>) -> Complex<F<T>>,
    ) -> (Complex<F<T>>, GenericThresholdCountertermComponentWeight<T>) {
        let Some(pass_one) = self.pass_one else {
            let weighted = Complex::new_re(zero.clone());
            return (
                weighted.clone(),
                GenericThresholdCountertermComponentWeight {
                    component_id: self.component_id,
                    occurrence: self.occurrence,
                    multiplier_values: self.multiplier_values,
                    effective_multiplier: self.effective_multiplier,
                    bare: None,
                    weighted,
                    evaluation_skipped: true,
                },
            );
        };

        let bare = evaluate_pass_two(pass_one);
        let weighted = &bare * &Complex::new_re(self.effective_multiplier.clone());
        (
            weighted.clone(),
            GenericThresholdCountertermComponentWeight {
                component_id: self.component_id,
                occurrence: self.occurrence,
                multiplier_values: self.multiplier_values,
                effective_multiplier: self.effective_multiplier,
                bare: Some(bare),
                weighted,
                evaluation_skipped: false,
            },
        )
    }
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

    pub fn cut_params_for_cut_cff_index(&self, cut_cff_index: &CutCFFIndex) -> LUParams<T> {
        let lu_order = cut_cff_index
            .lu_cut_order
            .expect("LU cut parameter lookup requires lu_cut_order")
            - 1;
        let base_params = self.cut_params_for_order(lu_order).clone();
        let target_shape = shape_from_cut_cff_index(cut_cff_index).map(HyperDual::new);
        let t_variable = variable_indices_from_cut_cff_index(cut_cff_index).lu_cut;

        LUParams {
            tstar: embed_real_dual_or_not_in_target_shape(
                &base_params.tstar,
                &target_shape,
                t_variable,
            ),
            h_function: embed_real_dual_or_not_in_target_shape(
                &base_params.h_function,
                &target_shape,
                t_variable,
            ),
        }
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

impl LUCounterTerm {
    fn radial_root_identity(
        graph_name: &str,
        cut_group_id: CutGroupId,
        side: &str,
        overlap_group: usize,
        esurface_id: EsurfaceID,
        probe_rotation: &Rotation,
    ) -> RadialRootIdentity {
        RadialRootIdentity::new(format!(
            "LU graph '{graph_name}' cut group {} {side} overlap group {overlap_group} E-surface {} probe rotation {}",
            cut_group_id.0, esurface_id.0, probe_rotation.method,
        ))
    }

    pub(crate) fn cut_group_is_active(&self, cut_group_id: CutGroupId) -> bool {
        self.active_cut_groups[cut_group_id]
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        frozen_mode: &FrozenCompilationMode,
    ) -> color_eyre::Result<()> {
        for (cut_group_id, evaluators) in self.evaluators.iter_mut_enumerated() {
            evaluators.compile(path.as_ref(), cut_group_id, frozen_mode)?;
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

    fn ensure_active_cut_group(&self, cut_group_id: CutGroupId) -> Result<()> {
        if self.cut_group_is_active(cut_group_id) {
            return Ok(());
        }

        Err(eyre!(
            "Cut group {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            cut_group_id.0
        ))
    }

    fn ensure_active_left_threshold(
        &self,
        cut_group_id: CutGroupId,
        threshold_id: LeftThresholdId,
    ) -> Result<()> {
        if self.active_left_thresholds[cut_group_id][threshold_id] {
            return Ok(());
        }

        Err(eyre!(
            "Left threshold evaluator {} for cut group {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            threshold_id.0,
            cut_group_id.0
        ))
    }

    fn ensure_active_right_threshold(
        &self,
        cut_group_id: CutGroupId,
        threshold_id: RightThresholdId,
    ) -> Result<()> {
        if self.active_right_thresholds[cut_group_id][threshold_id] {
            return Ok(());
        }

        Err(eyre!(
            "Right threshold evaluator {} for cut group {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            threshold_id.0,
            cut_group_id.0
        ))
    }

    fn ensure_active_iterated_threshold(
        &self,
        cut_group_id: CutGroupId,
        iterated_index: (LeftThresholdId, RightThresholdId),
    ) -> Result<()> {
        if self.active_iterated_thresholds[cut_group_id][iterated_index] {
            return Ok(());
        }

        Err(eyre!(
            "Iterated threshold evaluator ({}, {}) for cut group {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            iterated_index.0.0,
            iterated_index.1.0,
            cut_group_id.0
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate<T: FloatLike>(
        &mut self,
        kinematic_point: &LUCTKinematicPoint<T>,
        cut_group_id: CutGroupId,
        reversed_edges: &[EdgeIndex],
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        probe_rotation: &Rotation,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder<f64>,
        orientations: SingleOrAllOrientations<'_, OrientationID>,
        evaluation_meta_data: &mut EvaluationMetaData,
        record_primary_timing: bool,
        record_components: bool,
    ) -> Result<LUCountertermEvaluation<T>> {
        self.ensure_active_cut_group(cut_group_id)?;

        // Metadata alone does not request runtime allocations. Detailed storage is enabled only
        // for an accepted event whose caller explicitly asked for additional weights.
        let record_components = record_components && self.metadata_registry.is_some();
        let additional_params = settings.additional_params::<T>();

        let (left_thresholds_typed, right_thresholds_typed) = &self.thresholds[cut_group_id];

        let left_thresholds = TiVec::from_ref(&left_thresholds_typed.raw);
        let right_thresholds = TiVec::from_ref(&right_thresholds_typed.raw);

        if left_thresholds.is_empty() && right_thresholds.is_empty() {
            return Ok(LUCountertermEvaluation {
                total: Complex::new_re(kinematic_point.representative_sample().zero()),
                components: record_components.then(Vec::new),
            });
        }

        let (legacy_left_subspace, legacy_right_subspace) = &self.subspaces[cut_group_id];
        let variant_subspaces = self
            .variant_subspaces
            .as_ref()
            .map(|subspaces| &subspaces[cut_group_id]);
        let threshold_helper_outputs = if variant_subspaces.is_some() {
            LUThresholdHelperOutputs::Pieces
        } else if self.metadata_registry.is_some() {
            LUThresholdHelperOutputs::LegacyAndPieces
        } else {
            LUThresholdHelperOutputs::Legacy
        };
        let left_subspace = variant_subspaces
            .and_then(|subspaces| subspaces.left_common.as_ref())
            .unwrap_or(legacy_left_subspace);
        let right_subspace = variant_subspaces
            .and_then(|subspaces| subspaces.right_common.as_ref())
            .unwrap_or(legacy_right_subspace);
        let left_threshold_subspaces =
            variant_subspaces.map(|subspaces| subspaces.left.raw.as_slice());
        let right_threshold_subspaces =
            variant_subspaces.map(|subspaces| subspaces.right.raw.as_slice());
        let detailed_variant_ids = if record_components {
            let registry = self
                .metadata_registry
                .as_ref()
                .expect("record_components requires threshold metadata");
            let (left, right) = if let Some(subspaces) = variant_subspaces {
                (
                    subspaces.left_variant_ids.iter().copied().collect_vec(),
                    subspaces.right_variant_ids.iter().copied().collect_vec(),
                )
            } else {
                (
                    registry.lu_variant_ids(cut_group_id, ThresholdCountertermSide::Left)?,
                    registry.lu_variant_ids(cut_group_id, ThresholdCountertermSide::Right)?,
                )
            };
            if left.len() != left_thresholds.len() || right.len() != right_thresholds.len() {
                return Err(eyre!(
                    "graph '{}' cut group {} threshold metadata maps {} left/{} right variants onto {} left/{} right runtime thresholds",
                    graph.name,
                    cut_group_id.0,
                    left.len(),
                    right.len(),
                    left_thresholds.len(),
                    right_thresholds.len(),
                ));
            }
            Some((left, right))
        } else {
            None
        };
        let (sample_left_transformed, sample_right_transformed) = (
            kinematic_point
                .lmb_transform(&graph.loop_momentum_basis, left_subspace.get_lmb(all_lmbs)),
            kinematic_point
                .lmb_transform(&graph.loop_momentum_basis, right_subspace.get_lmb(all_lmbs)),
        );

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
        let esurface_existence_threshold =
            F::from_f64(settings.subtraction.esurface_existence_threshold);

        let left_existing_esurfaces = self.thresholds[cut_group_id]
            .0
            .iter_enumerated()
            .filter_map(|(left_id, esurface)| {
                let threshold_subspace = left_threshold_subspaces
                    .map_or(left_subspace, |subspaces| &subspaces[left_id.0]);
                let classification = esurface.classify_existence_subspace(
                    sample_left_transformed.representative_sample().loop_moms(),
                    kinematic_point.representative_sample().external_moms(),
                    threshold_subspace,
                    all_lmbs,
                    graph,
                    masses,
                    reversed_edges,
                    &e_cm,
                    &esurface_existence_threshold,
                );
                debug!(
                    graph = %graph.name,
                    side = "left",
                    esurface_id = left_id.0,
                    status = classification.label(),
                    normalized_margin = ?classification.normalized_margin(),
                    non_existing_reason = ?classification.non_existing_reason(),
                    "classified LU threshold surface"
                );
                if classification.is_existing() {
                    Some(EsurfaceID::from(left_id.0))
                } else {
                    None
                }
            })
            .collect::<TiVec<ExistingEsurfaceId, _>>();

        let right_existing_esurfaces = self.thresholds[cut_group_id]
            .1
            .iter_enumerated()
            .filter_map(|(right_id, esurface)| {
                let threshold_subspace = right_threshold_subspaces
                    .map_or(right_subspace, |subspaces| &subspaces[right_id.0]);
                let classification = esurface.classify_existence_subspace(
                    sample_right_transformed.representative_sample().loop_moms(),
                    kinematic_point.representative_sample().external_moms(),
                    threshold_subspace,
                    all_lmbs,
                    graph,
                    masses,
                    reversed_edges,
                    &e_cm,
                    &esurface_existence_threshold,
                );
                debug!(
                    graph = %graph.name,
                    side = "right",
                    esurface_id = right_id.0,
                    status = classification.label(),
                    normalized_margin = ?classification.normalized_margin(),
                    non_existing_reason = ?classification.non_existing_reason(),
                    "classified LU threshold surface"
                );
                if classification.is_existing() {
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

        if left_existing_esurfaces.is_empty() && right_existing_esurfaces.is_empty() {
            return Ok(LUCountertermEvaluation {
                total: Complex::new_re(kinematic_point.representative_sample().zero()),
                components: record_components.then(Vec::new),
            });
        }

        let left_overlap_input = OverlapInput {
            graph,
            subspace: left_subspace,
            settings,
            edge_masses: masses_f64.clone(),
            lmbs: all_lmbs,
            thresholds: left_thresholds,
            threshold_subspaces: left_threshold_subspaces,
        };

        let right_overlap_input = OverlapInput {
            graph,
            subspace: right_subspace,
            settings,
            edge_masses: masses_f64,
            lmbs: all_lmbs,
            thresholds: right_thresholds,
            threshold_subspaces: right_threshold_subspaces,
        };

        let left_overlap = match overlap_subspace::find_maximal_overlap(
            &left_overlap_input,
            &left_existing_esurfaces,
            &sample_left_transformed_f64,
            &external_moms_f64,
            probe_rotation,
        ) {
            Ok(left_overlap) => left_overlap,
            Err(error) => {
                evaluation_meta_data.record_threshold_counterterm_error(format!(
                    "LU graph '{}' cut group {} left subspace failed overlap-center construction in probe rotation {}: {}",
                    graph.name,
                    cut_group_id.0,
                    probe_rotation.method,
                    error,
                ));
                return Ok(LUCountertermEvaluation {
                    total: Complex::new_re(F::from_f64(f64::NAN)),
                    components: None,
                });
            }
        };

        let right_overlap = match overlap_subspace::find_maximal_overlap(
            &right_overlap_input,
            &right_existing_esurfaces,
            &sample_right_transformed_f64,
            &external_moms_f64,
            probe_rotation,
        ) {
            Ok(right_overlap) => right_overlap,
            Err(error) => {
                evaluation_meta_data.record_threshold_counterterm_error(format!(
                    "LU graph '{}' cut group {} right subspace failed overlap-center construction in probe rotation {}: {}",
                    graph.name,
                    cut_group_id.0,
                    probe_rotation.method,
                    error,
                ));
                return Ok(LUCountertermEvaluation {
                    total: Complex::new_re(F::from_f64(f64::NAN)),
                    components: None,
                });
            }
        };

        debug!("left overlap structure: {}", left_overlap);

        debug!("right overlap structure: {}", right_overlap);

        let left_counterterm_builder = CounterTermBuilder::new(
            graph,
            settings,
            left_overlap_input.thresholds,
            sample_left_transformed,
            &left_overlap,
            masses,
            all_lmbs,
            left_subspace,
            left_threshold_subspaces,
            probe_rotation,
        );

        let left_overlap_builders = if left_threshold_subspaces.is_some() {
            SideOverlapBuilders::Projected(
                left_overlap
                    .overlap_groups
                    .iter()
                    .map(|overlap_group| {
                        overlap_group
                            .existing_esurfaces
                            .iter()
                            .map(|existing_esurface_id| {
                                let esurface_id =
                                    left_overlap.existing_esurfaces[*existing_esurface_id];
                                left_counterterm_builder
                                    .new_overlap_builder(overlap_group, Some(esurface_id))
                            })
                            .collect_vec()
                    })
                    .collect_vec(),
            )
        } else {
            SideOverlapBuilders::Homogeneous(
                left_overlap
                    .overlap_groups
                    .iter()
                    .map(|overlap_group| {
                        left_counterterm_builder.new_overlap_builder(overlap_group, None)
                    })
                    .collect_vec(),
            )
        };

        let mut radial_root_failed = false;
        let mut left_overlap_solutions = Vec::with_capacity(left_overlap.overlap_groups.len());
        match &left_overlap_builders {
            SideOverlapBuilders::Homogeneous(builders) => {
                for (overlap_group_index, overlap_builder) in builders.iter().enumerate() {
                    let mut group_solutions =
                        Vec::with_capacity(overlap_builder.overlap_group.existing_esurfaces.len());
                    for &existing_esurface_id in &overlap_builder.overlap_group.existing_esurfaces {
                        let esurface_id = left_overlap.existing_esurfaces[existing_esurface_id];
                        let radial_root_identity = Self::radial_root_identity(
                            &graph.name,
                            cut_group_id,
                            "left",
                            overlap_group_index,
                            esurface_id,
                            probe_rotation,
                        );
                        let Some(solution) = overlap_builder
                            .new_esurface_builder(existing_esurface_id)
                            .solve_rstar(
                                &mut self.rstar_dependence_calculator[cut_group_id],
                                &radial_root_identity,
                                &mut evaluation_meta_data.radial_root_diagnostics,
                            )
                        else {
                            evaluation_meta_data.record_threshold_counterterm_error(format!(
                                "LU graph '{}' cut group {} left overlap group {} E-surface {} failed center or radial-root validation in probe rotation {}",
                                graph.name,
                                cut_group_id.0,
                                overlap_group_index,
                                esurface_id.0,
                                probe_rotation.method,
                            ));
                            radial_root_failed = true;
                            continue;
                        };
                        group_solutions.push(solution);
                    }
                    left_overlap_solutions.push(group_solutions);
                }
            }
            SideOverlapBuilders::Projected(builders) => {
                for (overlap_group_index, (overlap_group, group_builders)) in
                    left_overlap.overlap_groups.iter().zip(builders).enumerate()
                {
                    let mut group_solutions = Vec::with_capacity(group_builders.len());
                    for (&existing_esurface_id, overlap_builder) in
                        overlap_group.existing_esurfaces.iter().zip(group_builders)
                    {
                        let esurface_id = left_overlap.existing_esurfaces[existing_esurface_id];
                        let radial_root_identity = Self::radial_root_identity(
                            &graph.name,
                            cut_group_id,
                            "left",
                            overlap_group_index,
                            esurface_id,
                            probe_rotation,
                        );
                        let Some(solution) = overlap_builder
                            .new_esurface_builder(existing_esurface_id)
                            .solve_rstar(
                                &mut self.rstar_dependence_calculator[cut_group_id],
                                &radial_root_identity,
                                &mut evaluation_meta_data.radial_root_diagnostics,
                            )
                        else {
                            evaluation_meta_data.record_threshold_counterterm_error(format!(
                                "LU graph '{}' cut group {} left projected E-surface instance {} in overlap group {} failed center or radial-root validation in probe rotation {}",
                                graph.name,
                                cut_group_id.0,
                                esurface_id.0,
                                overlap_group_index,
                                probe_rotation.method,
                            ));
                            radial_root_failed = true;
                            continue;
                        };
                        group_solutions.push(solution);
                    }
                    left_overlap_solutions.push(group_solutions);
                }
            }
        }

        let right_counterterm_builder = CounterTermBuilder::new(
            graph,
            settings,
            right_overlap_input.thresholds,
            sample_right_transformed,
            &right_overlap,
            masses,
            all_lmbs,
            right_subspace,
            right_threshold_subspaces,
            probe_rotation,
        );

        let right_overlap_builders = if right_threshold_subspaces.is_some() {
            SideOverlapBuilders::Projected(
                right_overlap
                    .overlap_groups
                    .iter()
                    .map(|overlap_group| {
                        overlap_group
                            .existing_esurfaces
                            .iter()
                            .map(|existing_esurface_id| {
                                let esurface_id =
                                    right_overlap.existing_esurfaces[*existing_esurface_id];
                                right_counterterm_builder
                                    .new_overlap_builder(overlap_group, Some(esurface_id))
                            })
                            .collect_vec()
                    })
                    .collect_vec(),
            )
        } else {
            SideOverlapBuilders::Homogeneous(
                right_overlap
                    .overlap_groups
                    .iter()
                    .map(|overlap_group| {
                        right_counterterm_builder.new_overlap_builder(overlap_group, None)
                    })
                    .collect_vec(),
            )
        };

        let mut right_overlap_solutions = Vec::with_capacity(right_overlap.overlap_groups.len());
        match &right_overlap_builders {
            SideOverlapBuilders::Homogeneous(builders) => {
                for (overlap_group_index, overlap_builder) in builders.iter().enumerate() {
                    let mut group_solutions =
                        Vec::with_capacity(overlap_builder.overlap_group.existing_esurfaces.len());
                    for &existing_esurface_id in &overlap_builder.overlap_group.existing_esurfaces {
                        let esurface_id = right_overlap.existing_esurfaces[existing_esurface_id];
                        let radial_root_identity = Self::radial_root_identity(
                            &graph.name,
                            cut_group_id,
                            "right",
                            overlap_group_index,
                            esurface_id,
                            probe_rotation,
                        );
                        let Some(solution) = overlap_builder
                            .new_esurface_builder(existing_esurface_id)
                            .solve_rstar(
                                &mut self.rstar_dependence_calculator[cut_group_id],
                                &radial_root_identity,
                                &mut evaluation_meta_data.radial_root_diagnostics,
                            )
                        else {
                            evaluation_meta_data.record_threshold_counterterm_error(format!(
                                "LU graph '{}' cut group {} right overlap group {} E-surface {} failed center or radial-root validation in probe rotation {}",
                                graph.name,
                                cut_group_id.0,
                                overlap_group_index,
                                esurface_id.0,
                                probe_rotation.method,
                            ));
                            radial_root_failed = true;
                            continue;
                        };
                        group_solutions.push(solution);
                    }
                    right_overlap_solutions.push(group_solutions);
                }
            }
            SideOverlapBuilders::Projected(builders) => {
                for (overlap_group_index, (overlap_group, group_builders)) in right_overlap
                    .overlap_groups
                    .iter()
                    .zip(builders)
                    .enumerate()
                {
                    let mut group_solutions = Vec::with_capacity(group_builders.len());
                    for (&existing_esurface_id, overlap_builder) in
                        overlap_group.existing_esurfaces.iter().zip(group_builders)
                    {
                        let esurface_id = right_overlap.existing_esurfaces[existing_esurface_id];
                        let radial_root_identity = Self::radial_root_identity(
                            &graph.name,
                            cut_group_id,
                            "right",
                            overlap_group_index,
                            esurface_id,
                            probe_rotation,
                        );
                        let Some(solution) = overlap_builder
                            .new_esurface_builder(existing_esurface_id)
                            .solve_rstar(
                                &mut self.rstar_dependence_calculator[cut_group_id],
                                &radial_root_identity,
                                &mut evaluation_meta_data.radial_root_diagnostics,
                            )
                        else {
                            evaluation_meta_data.record_threshold_counterterm_error(format!(
                                "LU graph '{}' cut group {} right projected E-surface instance {} in overlap group {} failed center or radial-root validation in probe rotation {}",
                                graph.name,
                                cut_group_id.0,
                                esurface_id.0,
                                overlap_group_index,
                                probe_rotation.method,
                            ));
                            radial_root_failed = true;
                            continue;
                        };
                        group_solutions.push(solution);
                    }
                    right_overlap_solutions.push(group_solutions);
                }
            }
        }

        if radial_root_failed {
            return Ok(LUCountertermEvaluation {
                total: Complex::new_re(F::from_f64(f64::NAN)),
                components: None,
            });
        }

        let zero = kinematic_point.representative_sample().zero();
        let mut total_result = Complex::new_re(zero.clone());
        let mut completed_components = record_components.then(Vec::new);
        let mut threshold_multiplier_workspace = self.evaluators[cut_group_id]
            .threshold_multipliers
            .as_ref()
            .map(|collection| {
                ThresholdMultiplierInputWorkspace::new(collection.layout(), zero.clone())
            });
        if threshold_multiplier_workspace.is_some() && variant_subspaces.is_none() {
            return Err(eyre!(
                "graph '{}' cut group {} has threshold multipliers but no variant/subspace registry",
                graph.name,
                cut_group_id.0,
            ));
        }

        for order in 0..kinematic_point.lu_cut_parameter_cache.len() {
            let mut pending_components = record_components.then(Vec::new);
            let mut left_evaluations = zero_dual_or_not_complex(order, &zero);

            for (overlap_group_index, solutions_group) in left_overlap_solutions.iter().enumerate()
            {
                for solution in solutions_group {
                    let representative_sample = solution.rstar_sample_for_order(order);
                    let left_threshold_id =
                        LeftThresholdId::from(representative_sample.get_esurface_id().0);
                    self.ensure_active_left_threshold(cut_group_id, left_threshold_id)?;

                    let matching_cut_indices = self.evaluators[cut_group_id]
                        .left_thresholds_evaluator[left_threshold_id]
                        .keys()
                        .filter(|cut_cff_index| {
                            cut_cff_index.lu_cut_order == Some(order + 1)
                                && cut_cff_index.right_threshold_order.is_none()
                        })
                        .copied()
                        .collect_vec();

                    for cut_cff_index in matching_cut_indices {
                        let variable_indices = variable_indices_from_cut_cff_index(&cut_cff_index);
                        let lu_cut_params =
                            kinematic_point.cut_params_for_cut_cff_index(&cut_cff_index);
                        let sample = solution.rstar_sample_for_cut_cff_index(
                            &cut_cff_index,
                            variable_indices.left_threshold,
                        );
                        debug!("left threshold parameters");

                        let left_threshold_params: ThresholdParams<T> =
                            sample.extract_threshold_parameters(true);
                        debug!(
                            "LU left evaluator input: cut_group_id={}, left_threshold_id={}, cut_cff_index={:?}, lu_order={}, r={}, rstar={}",
                            cut_group_id.0,
                            left_threshold_id.0,
                            cut_cff_index,
                            order + 1,
                            left_threshold_params.radius,
                            left_threshold_params.radius_star,
                        );
                        let inverse_transformed_sample = sample.get_inverse_transformed_sample();
                        let detailed_variant_id = detailed_variant_ids
                            .as_ref()
                            .map(|(left, _)| left[left_threshold_id.0]);
                        let multiplier_coefficients = if let Some(workspace) =
                            threshold_multiplier_workspace.as_mut()
                        {
                            let variant_id = variant_subspaces
                                .expect("multiplier registry requires variant subspaces")
                                .left_variant_ids[left_threshold_id];
                            let collection = self.evaluators[cut_group_id]
                                .threshold_multipliers
                                .as_mut()
                                .expect("multiplier workspace requires evaluator collection");
                            let (effective, star) = single_multiplier_context(
                                SingleMultiplierComponent::Local,
                                kinematic_point.representative_sample(),
                                &inverse_transformed_sample,
                            );
                            let local = evaluate_single_multiplier_context(
                                collection,
                                workspace,
                                variant_id,
                                &graph.loop_momentum_basis,
                                masses,
                                param_builder.model_values(),
                                &additional_params,
                                effective,
                                star,
                                evaluation_meta_data,
                                record_primary_timing,
                            )
                            .with_context(|| {
                                format!(
                                    "Failed to evaluate local multiplier for graph '{}' cut group {} left variant {}",
                                    graph.name, cut_group_id.0, variant_id.0,
                                )
                            })?;
                            let (effective, star) = single_multiplier_context(
                                SingleMultiplierComponent::Integrated,
                                kinematic_point.representative_sample(),
                                &inverse_transformed_sample,
                            );
                            let integrated = evaluate_single_multiplier_context(
                                collection,
                                workspace,
                                variant_id,
                                &graph.loop_momentum_basis,
                                masses,
                                param_builder.model_values(),
                                &additional_params,
                                effective,
                                star,
                                evaluation_meta_data,
                                record_primary_timing,
                            )
                            .with_context(|| {
                                format!(
                                    "Failed to evaluate integrated multiplier for graph '{}' cut group {} left variant {}",
                                    graph.name, cut_group_id.0, variant_id.0,
                                )
                            })?;
                            let coefficients = SingleMultiplierCoefficients { local, integrated };
                            Some(coefficients)
                        } else if record_components {
                            let one = kinematic_point.representative_sample().one();
                            Some(SingleMultiplierCoefficients {
                                local: one.clone(),
                                integrated: one,
                            })
                        } else {
                            None
                        };

                        if multiplier_coefficients
                            .as_ref()
                            .is_some_and(SingleMultiplierCoefficients::all_zero)
                        {
                            if let Some(pending) = pending_components.as_mut() {
                                let coefficients = multiplier_coefficients.as_ref().unwrap();
                                let variant_id = detailed_variant_id.expect(
                                    "detailed zero multiplier requires a stable left variant ID",
                                );
                                let registry = self.metadata_registry.as_ref().unwrap();
                                for (kind, coefficient) in [
                                    (
                                        ThresholdCountertermComponentKind::Local,
                                        &coefficients.local,
                                    ),
                                    (
                                        ThresholdCountertermComponentKind::Integrated,
                                        &coefficients.integrated,
                                    ),
                                ] {
                                    pending.push(PendingLUCountertermComponent {
                                        component_id: registry.component_id(
                                            Some(cut_group_id),
                                            kind,
                                            &[variant_id],
                                        )?,
                                        occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                            overlap_groups: SmallVec::from_slice(&[
                                                overlap_group_index,
                                            ]),
                                            left_threshold_order: cut_cff_index.left_threshold_order,
                                            right_threshold_order: None,
                                            lu_cut_order: cut_cff_index.lu_cut_order,
                                        },
                                        multiplier_values: smallvec![coefficient.clone()],
                                        effective_multiplier: coefficient.clone(),
                                        pass_one: None,
                                    });
                                }
                            }
                            continue;
                        }

                        let params = T::get_parameters(
                            param_builder,
                            (false, false),
                            graph,
                            &inverse_transformed_sample,
                            settings.kinematics.externals.get_helicities(),
                            &additional_params,
                            Some(&left_threshold_params),
                            None,
                            Some(&lu_cut_params),
                        );

                        let result_of_this_ct = self.evaluators[cut_group_id]
                            .left_thresholds_evaluator[left_threshold_id]
                            .get_mut(&cut_cff_index)
                            .unwrap()
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

                        debug!(
                            "result of left threshold evaluator {:?}: {}",
                            cut_cff_index, result_of_this_ct
                        );

                        let helper_evaluation = evaluate_threshold_helper_single(
                            ThresholdHelperEvaluationContext {
                                helper: self.evaluators[cut_group_id]
                                    .threshold_helpers
                                    .left_thresholds[left_threshold_id]
                                    .get_mut(&cut_cff_index)
                                    .unwrap(),
                                cut_cff_index: &cut_cff_index,
                                threshold_result: &result_of_this_ct,
                                outputs: threshold_helper_outputs,
                                evaluation_meta_data,
                                record_primary_timing,
                            },
                            &left_threshold_params,
                        );
                        let multi_channeling_factor = plain_t_dual_or_scalar_complex(
                            &sample.value_of_multi_channeling_factor,
                            variable_indices.lu_cut,
                        );
                        if let Some(pending) = pending_components.as_mut() {
                            let helper_completed_pieces = helper_evaluation.into_pieces();
                            let coefficients = multiplier_coefficients
                                .as_ref()
                                .expect("detailed components always materialize identity factors");
                            let variant_id = detailed_variant_id
                                .expect("detailed left component requires a stable variant ID");
                            let registry = self.metadata_registry.as_ref().unwrap();
                            for (kind, piece, coefficient) in [
                                (
                                    ThresholdCountertermComponentKind::Local,
                                    helper_completed_pieces.local,
                                    coefficients.local.clone(),
                                ),
                                (
                                    ThresholdCountertermComponentKind::Integrated,
                                    helper_completed_pieces.integrated,
                                    coefficients.integrated.clone(),
                                ),
                            ] {
                                let pass_one = (!coefficient.is_zero()).then(|| {
                                    negate_dual_or_not_complex(multiply_dual_or_not_complex(
                                        piece,
                                        &multi_channeling_factor,
                                    ))
                                });
                                pending.push(PendingLUCountertermComponent {
                                    component_id: registry.component_id(
                                        Some(cut_group_id),
                                        kind,
                                        &[variant_id],
                                    )?,
                                    occurrence:
                                        ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                            overlap_groups: SmallVec::from_slice(&[
                                                overlap_group_index,
                                            ]),
                                            left_threshold_order: cut_cff_index
                                                .left_threshold_order,
                                            right_threshold_order: None,
                                            lu_cut_order: cut_cff_index.lu_cut_order,
                                        },
                                    multiplier_values: smallvec![coefficient.clone()],
                                    effective_multiplier: coefficient,
                                    pass_one,
                                });
                            }
                        } else {
                            let helper_completed_result = if let Some(coefficients) =
                                &multiplier_coefficients
                            {
                                weight_single_multiplier_pieces(
                                    helper_evaluation.into_pieces(),
                                    coefficients,
                                    order,
                                    &zero,
                                )
                            } else if threshold_helper_outputs == LUThresholdHelperOutputs::Pieces {
                                let helper_completed_pieces = helper_evaluation.into_pieces();
                                let mut result = helper_completed_pieces.local;
                                result += helper_completed_pieces.integrated;
                                result
                            } else {
                                helper_evaluation.into_legacy()
                            };
                            left_evaluations += multiply_dual_or_not_complex(
                                helper_completed_result,
                                &multi_channeling_factor,
                            );
                        }
                    }
                }
            }

            let mut right_evaluations = zero_dual_or_not_complex(order, &zero);

            for (overlap_group_index, solutions_group) in right_overlap_solutions.iter().enumerate()
            {
                for solution in solutions_group {
                    let representative_sample = solution.rstar_sample_for_order(order);
                    let right_threshold_id =
                        RightThresholdId::from(representative_sample.get_esurface_id().0);
                    self.ensure_active_right_threshold(cut_group_id, right_threshold_id)?;

                    let matching_cut_indices = self.evaluators[cut_group_id]
                        .right_thresholds_evaluator[right_threshold_id]
                        .keys()
                        .filter(|cut_cff_index| {
                            cut_cff_index.lu_cut_order == Some(order + 1)
                                && cut_cff_index.left_threshold_order.is_none()
                        })
                        .copied()
                        .collect_vec();

                    for cut_cff_index in matching_cut_indices {
                        let variable_indices = variable_indices_from_cut_cff_index(&cut_cff_index);
                        let lu_cut_params =
                            kinematic_point.cut_params_for_cut_cff_index(&cut_cff_index);
                        let sample = solution.rstar_sample_for_cut_cff_index(
                            &cut_cff_index,
                            variable_indices.right_threshold,
                        );
                        debug!("right threshold parameters");
                        let right_threshold_params: ThresholdParams<T> =
                            sample.extract_threshold_parameters(true);
                        debug!(
                            "LU right evaluator input: cut_group_id={}, right_threshold_id={}, cut_cff_index={:?}, lu_order={}, r={}, rstar={}",
                            cut_group_id.0,
                            right_threshold_id.0,
                            cut_cff_index,
                            order + 1,
                            right_threshold_params.radius,
                            right_threshold_params.radius_star,
                        );
                        let inverse_transformed_sample = sample.get_inverse_transformed_sample();
                        let detailed_variant_id = detailed_variant_ids
                            .as_ref()
                            .map(|(_, right)| right[right_threshold_id.0]);
                        let multiplier_coefficients = if let Some(workspace) =
                            threshold_multiplier_workspace.as_mut()
                        {
                            let variant_id = variant_subspaces
                                .expect("multiplier registry requires variant subspaces")
                                .right_variant_ids[right_threshold_id];
                            let collection = self.evaluators[cut_group_id]
                                .threshold_multipliers
                                .as_mut()
                                .expect("multiplier workspace requires evaluator collection");
                            let (effective, star) = single_multiplier_context(
                                SingleMultiplierComponent::Local,
                                kinematic_point.representative_sample(),
                                &inverse_transformed_sample,
                            );
                            let local = evaluate_single_multiplier_context(
                                collection,
                                workspace,
                                variant_id,
                                &graph.loop_momentum_basis,
                                masses,
                                param_builder.model_values(),
                                &additional_params,
                                effective,
                                star,
                                evaluation_meta_data,
                                record_primary_timing,
                            )
                            .with_context(|| {
                                format!(
                                    "Failed to evaluate local multiplier for graph '{}' cut group {} right variant {}",
                                    graph.name, cut_group_id.0, variant_id.0,
                                )
                            })?;
                            let (effective, star) = single_multiplier_context(
                                SingleMultiplierComponent::Integrated,
                                kinematic_point.representative_sample(),
                                &inverse_transformed_sample,
                            );
                            let integrated = evaluate_single_multiplier_context(
                                collection,
                                workspace,
                                variant_id,
                                &graph.loop_momentum_basis,
                                masses,
                                param_builder.model_values(),
                                &additional_params,
                                effective,
                                star,
                                evaluation_meta_data,
                                record_primary_timing,
                            )
                            .with_context(|| {
                                format!(
                                    "Failed to evaluate integrated multiplier for graph '{}' cut group {} right variant {}",
                                    graph.name, cut_group_id.0, variant_id.0,
                                )
                            })?;
                            let coefficients = SingleMultiplierCoefficients { local, integrated };
                            Some(coefficients)
                        } else if record_components {
                            let one = kinematic_point.representative_sample().one();
                            Some(SingleMultiplierCoefficients {
                                local: one.clone(),
                                integrated: one,
                            })
                        } else {
                            None
                        };

                        if multiplier_coefficients
                            .as_ref()
                            .is_some_and(SingleMultiplierCoefficients::all_zero)
                        {
                            if let Some(pending) = pending_components.as_mut() {
                                let coefficients = multiplier_coefficients.as_ref().unwrap();
                                let variant_id = detailed_variant_id.expect(
                                    "detailed zero multiplier requires a stable right variant ID",
                                );
                                let registry = self.metadata_registry.as_ref().unwrap();
                                for (kind, coefficient) in [
                                    (
                                        ThresholdCountertermComponentKind::Local,
                                        &coefficients.local,
                                    ),
                                    (
                                        ThresholdCountertermComponentKind::Integrated,
                                        &coefficients.integrated,
                                    ),
                                ] {
                                    pending.push(PendingLUCountertermComponent {
                                        component_id: registry.component_id(
                                            Some(cut_group_id),
                                            kind,
                                            &[variant_id],
                                        )?,
                                        occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                            overlap_groups: SmallVec::from_slice(&[
                                                overlap_group_index,
                                            ]),
                                            left_threshold_order: None,
                                            right_threshold_order: cut_cff_index.right_threshold_order,
                                            lu_cut_order: cut_cff_index.lu_cut_order,
                                        },
                                        multiplier_values: smallvec![coefficient.clone()],
                                        effective_multiplier: coefficient.clone(),
                                        pass_one: None,
                                    });
                                }
                            }
                            continue;
                        }

                        let params = T::get_parameters(
                            param_builder,
                            (false, false),
                            graph,
                            &inverse_transformed_sample,
                            settings.kinematics.externals.get_helicities(),
                            &additional_params,
                            None,
                            Some(&right_threshold_params),
                            Some(&lu_cut_params),
                        );

                        let result_of_this_ct = self.evaluators[cut_group_id]
                            .right_thresholds_evaluator[right_threshold_id]
                            .get_mut(&cut_cff_index)
                            .unwrap()
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

                        debug!(
                            "result of right threshold evaluator {:?}: {}",
                            cut_cff_index, result_of_this_ct
                        );

                        let helper_evaluation = evaluate_threshold_helper_single(
                            ThresholdHelperEvaluationContext {
                                helper: self.evaluators[cut_group_id]
                                    .threshold_helpers
                                    .right_thresholds[right_threshold_id]
                                    .get_mut(&cut_cff_index)
                                    .unwrap(),
                                cut_cff_index: &cut_cff_index,
                                threshold_result: &result_of_this_ct,
                                outputs: threshold_helper_outputs,
                                evaluation_meta_data,
                                record_primary_timing,
                            },
                            &right_threshold_params,
                        );
                        let multi_channeling_factor = plain_t_dual_or_scalar_complex(
                            &sample.value_of_multi_channeling_factor,
                            variable_indices.lu_cut,
                        );
                        if let Some(pending) = pending_components.as_mut() {
                            let helper_completed_pieces = helper_evaluation.into_pieces();
                            let coefficients = multiplier_coefficients
                                .as_ref()
                                .expect("detailed components always materialize identity factors");
                            let variant_id = detailed_variant_id
                                .expect("detailed right component requires a stable variant ID");
                            let registry = self.metadata_registry.as_ref().unwrap();
                            for (kind, piece, coefficient) in [
                                (
                                    ThresholdCountertermComponentKind::Local,
                                    helper_completed_pieces.local,
                                    coefficients.local.clone(),
                                ),
                                (
                                    ThresholdCountertermComponentKind::Integrated,
                                    helper_completed_pieces.integrated,
                                    coefficients.integrated.clone(),
                                ),
                            ] {
                                let pass_one = (!coefficient.is_zero()).then(|| {
                                    negate_dual_or_not_complex(multiply_dual_or_not_complex(
                                        piece,
                                        &multi_channeling_factor,
                                    ))
                                });
                                pending.push(PendingLUCountertermComponent {
                                    component_id: registry.component_id(
                                        Some(cut_group_id),
                                        kind,
                                        &[variant_id],
                                    )?,
                                    occurrence:
                                        ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                            overlap_groups: SmallVec::from_slice(&[
                                                overlap_group_index,
                                            ]),
                                            left_threshold_order: None,
                                            right_threshold_order: cut_cff_index
                                                .right_threshold_order,
                                            lu_cut_order: cut_cff_index.lu_cut_order,
                                        },
                                    multiplier_values: smallvec![coefficient.clone()],
                                    effective_multiplier: coefficient,
                                    pass_one,
                                });
                            }
                        } else {
                            let helper_completed_result = if let Some(coefficients) =
                                &multiplier_coefficients
                            {
                                weight_single_multiplier_pieces(
                                    helper_evaluation.into_pieces(),
                                    coefficients,
                                    order,
                                    &zero,
                                )
                            } else if threshold_helper_outputs == LUThresholdHelperOutputs::Pieces {
                                let helper_completed_pieces = helper_evaluation.into_pieces();
                                let mut result = helper_completed_pieces.local;
                                result += helper_completed_pieces.integrated;
                                result
                            } else {
                                helper_evaluation.into_legacy()
                            };
                            right_evaluations += multiply_dual_or_not_complex(
                                helper_completed_result,
                                &multi_channeling_factor,
                            );
                        }
                    }
                }
            }

            let mut cartesian_product_result = zero_dual_or_not_complex(order, &zero);

            for (left_overlap_group_index, left_solutions_group) in
                left_overlap_solutions.iter().enumerate()
            {
                for left_solution in left_solutions_group {
                    let left_sample_representative = left_solution.rstar_sample_for_order(order);
                    for (right_overlap_group_index, right_solutions_group) in
                        right_overlap_solutions.iter().enumerate()
                    {
                        for right_solution in right_solutions_group {
                            let right_sample_representative =
                                right_solution.rstar_sample_for_order(order);
                            let iterated_index = (
                                LeftThresholdId::from(
                                    left_sample_representative.get_esurface_id().0,
                                ),
                                RightThresholdId::from(
                                    right_sample_representative.get_esurface_id().0,
                                ),
                            );
                            self.ensure_active_iterated_threshold(cut_group_id, iterated_index)?;

                            let matching_cut_indices = self.evaluators[cut_group_id]
                                .iterated_evaluator[iterated_index]
                                .keys()
                                .filter(|cut_cff_index| {
                                    cut_cff_index.lu_cut_order == Some(order + 1)
                                })
                                .copied()
                                .collect_vec();

                            for cut_cff_index in matching_cut_indices {
                                let variable_indices =
                                    variable_indices_from_cut_cff_index(&cut_cff_index);
                                let lu_cut_params =
                                    kinematic_point.cut_params_for_cut_cff_index(&cut_cff_index);
                                let sample_left = left_solution.rstar_sample_for_cut_cff_index(
                                    &cut_cff_index,
                                    variable_indices.left_threshold,
                                );
                                let sample_right = right_solution.rstar_sample_for_cut_cff_index(
                                    &cut_cff_index,
                                    variable_indices.right_threshold,
                                );

                                let left_threshold_params: ThresholdParams<T> =
                                    sample_left.extract_threshold_parameters(false);
                                let right_threshold_params: ThresholdParams<T> =
                                    sample_right.extract_threshold_parameters(false);
                                debug!(
                                    "LU iterated evaluator input: cut_group_id={}, left_threshold_id={}, right_threshold_id={}, cut_cff_index={:?}, lu_order={}, left_r={}, left_rstar={}, right_r={}, right_rstar={}",
                                    cut_group_id.0,
                                    iterated_index.0.0,
                                    iterated_index.1.0,
                                    cut_cff_index,
                                    order + 1,
                                    left_threshold_params.radius,
                                    left_threshold_params.radius_star,
                                    right_threshold_params.radius,
                                    right_threshold_params.radius_star,
                                );
                                let inverse_transformed_momentum_sample =
                                    merge_and_inverse_transform(&sample_left, &sample_right);
                                let detailed_pair_variant_ids =
                                    detailed_variant_ids.as_ref().map(|(left, right)| {
                                        (left[iterated_index.0.0], right[iterated_index.1.0])
                                    });
                                let multiplier_coefficients = if let Some(workspace) =
                                    threshold_multiplier_workspace.as_mut()
                                {
                                    let variant_subspaces = variant_subspaces
                                        .expect("multiplier registry requires variant subspaces");
                                    let left_variant_id =
                                        variant_subspaces.left_variant_ids[iterated_index.0];
                                    let right_variant_id =
                                        variant_subspaces.right_variant_ids[iterated_index.1];
                                    let left_root = sample_left.get_inverse_transformed_sample();
                                    let right_root = sample_right.get_inverse_transformed_sample();
                                    let collection = self.evaluators[cut_group_id]
                                        .threshold_multipliers
                                        .as_mut()
                                        .expect(
                                            "multiplier workspace requires evaluator collection",
                                        );

                                    let (effective, star) = iterated_multiplier_context(
                                        IteratedMultiplierComponent::LocalLocal,
                                        kinematic_point.representative_sample(),
                                        &left_root,
                                        &right_root,
                                        &inverse_transformed_momentum_sample,
                                    );
                                    let (left, right) = evaluate_pair_multiplier_context(
                                        collection,
                                        workspace,
                                        left_variant_id,
                                        right_variant_id,
                                        &graph.loop_momentum_basis,
                                        masses,
                                        param_builder.model_values(),
                                        &additional_params,
                                        effective,
                                        star,
                                        evaluation_meta_data,
                                        record_primary_timing,
                                    )
                                    .with_context(|| {
                                        format!(
                                            "Failed to evaluate LL multipliers for graph '{}' cut group {} variants ({}, {})",
                                            graph.name,
                                            cut_group_id.0,
                                            left_variant_id.0,
                                            right_variant_id.0,
                                        )
                                    })?;
                                    let local_local = PairMultiplierCoefficient { left, right };

                                    let (effective, star) = iterated_multiplier_context(
                                        IteratedMultiplierComponent::IntegratedLocal,
                                        kinematic_point.representative_sample(),
                                        &left_root,
                                        &right_root,
                                        &inverse_transformed_momentum_sample,
                                    );
                                    let (left, right) = evaluate_pair_multiplier_context(
                                        collection,
                                        workspace,
                                        left_variant_id,
                                        right_variant_id,
                                        &graph.loop_momentum_basis,
                                        masses,
                                        param_builder.model_values(),
                                        &additional_params,
                                        effective,
                                        star,
                                        evaluation_meta_data,
                                        record_primary_timing,
                                    )
                                    .with_context(|| {
                                        format!(
                                            "Failed to evaluate IL multipliers for graph '{}' cut group {} variants ({}, {})",
                                            graph.name,
                                            cut_group_id.0,
                                            left_variant_id.0,
                                            right_variant_id.0,
                                        )
                                    })?;
                                    let integrated_local =
                                        PairMultiplierCoefficient { left, right };

                                    let (effective, star) = iterated_multiplier_context(
                                        IteratedMultiplierComponent::LocalIntegrated,
                                        kinematic_point.representative_sample(),
                                        &left_root,
                                        &right_root,
                                        &inverse_transformed_momentum_sample,
                                    );
                                    let (left, right) = evaluate_pair_multiplier_context(
                                        collection,
                                        workspace,
                                        left_variant_id,
                                        right_variant_id,
                                        &graph.loop_momentum_basis,
                                        masses,
                                        param_builder.model_values(),
                                        &additional_params,
                                        effective,
                                        star,
                                        evaluation_meta_data,
                                        record_primary_timing,
                                    )
                                    .with_context(|| {
                                        format!(
                                            "Failed to evaluate LI multipliers for graph '{}' cut group {} variants ({}, {})",
                                            graph.name,
                                            cut_group_id.0,
                                            left_variant_id.0,
                                            right_variant_id.0,
                                        )
                                    })?;
                                    let local_integrated =
                                        PairMultiplierCoefficient { left, right };

                                    let (effective, star) = iterated_multiplier_context(
                                        IteratedMultiplierComponent::IntegratedIntegrated,
                                        kinematic_point.representative_sample(),
                                        &left_root,
                                        &right_root,
                                        &inverse_transformed_momentum_sample,
                                    );
                                    let (left, right) = evaluate_pair_multiplier_context(
                                        collection,
                                        workspace,
                                        left_variant_id,
                                        right_variant_id,
                                        &graph.loop_momentum_basis,
                                        masses,
                                        param_builder.model_values(),
                                        &additional_params,
                                        effective,
                                        star,
                                        evaluation_meta_data,
                                        record_primary_timing,
                                    )
                                    .with_context(|| {
                                        format!(
                                            "Failed to evaluate II multipliers for graph '{}' cut group {} variants ({}, {})",
                                            graph.name,
                                            cut_group_id.0,
                                            left_variant_id.0,
                                            right_variant_id.0,
                                        )
                                    })?;
                                    let integrated_integrated =
                                        PairMultiplierCoefficient { left, right };

                                    let coefficients = IteratedMultiplierCoefficients {
                                        local_local,
                                        local_integrated,
                                        integrated_local,
                                        integrated_integrated,
                                    };
                                    Some(coefficients)
                                } else if record_components {
                                    let one = kinematic_point.representative_sample().one();
                                    let identity = || PairMultiplierCoefficient {
                                        left: one.clone(),
                                        right: one.clone(),
                                    };
                                    Some(IteratedMultiplierCoefficients {
                                        local_local: identity(),
                                        local_integrated: identity(),
                                        integrated_local: identity(),
                                        integrated_integrated: identity(),
                                    })
                                } else {
                                    None
                                };

                                if multiplier_coefficients
                                    .as_ref()
                                    .is_some_and(IteratedMultiplierCoefficients::all_zero)
                                {
                                    if let Some(pending) = pending_components.as_mut() {
                                        let coefficients =
                                            multiplier_coefficients.as_ref().unwrap();
                                        let (left_variant_id, right_variant_id) =
                                            detailed_pair_variant_ids.expect(
                                                "detailed zero multipliers require stable pair variant IDs",
                                            );
                                        let registry = self.metadata_registry.as_ref().unwrap();
                                        for (kind, coefficient) in [
                                            (
                                                ThresholdCountertermComponentKind::LocalLocal,
                                                &coefficients.local_local,
                                            ),
                                            (
                                                ThresholdCountertermComponentKind::LocalIntegrated,
                                                &coefficients.local_integrated,
                                            ),
                                            (
                                                ThresholdCountertermComponentKind::IntegratedLocal,
                                                &coefficients.integrated_local,
                                            ),
                                            (
                                                ThresholdCountertermComponentKind::IntegratedIntegrated,
                                                &coefficients.integrated_integrated,
                                            ),
                                        ] {
                                            pending.push(PendingLUCountertermComponent {
                                                component_id: registry.component_id(
                                                    Some(cut_group_id),
                                                    kind,
                                                    &[left_variant_id, right_variant_id],
                                                )?,
                                                occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                                    overlap_groups: SmallVec::from_slice(&[
                                                        left_overlap_group_index,
                                                        right_overlap_group_index,
                                                    ]),
                                                    left_threshold_order: cut_cff_index.left_threshold_order,
                                                    right_threshold_order: cut_cff_index.right_threshold_order,
                                                    lu_cut_order: cut_cff_index.lu_cut_order,
                                                },
                                                multiplier_values: smallvec![
                                                    coefficient.left.clone(),
                                                    coefficient.right.clone(),
                                                ],
                                                effective_multiplier: coefficient.effective(),
                                                pass_one: None,
                                            });
                                        }
                                    }
                                    continue;
                                }

                                let multi_channeling_factor = plain_t_dual_or_scalar_complex(
                                    &multiply_dual_or_not_complex(
                                        sample_left.value_of_multi_channeling_factor.clone(),
                                        &sample_right.value_of_multi_channeling_factor,
                                    ),
                                    variable_indices.lu_cut,
                                );

                                let params = T::get_parameters(
                                    param_builder,
                                    (false, false),
                                    graph,
                                    &inverse_transformed_momentum_sample,
                                    settings.kinematics.externals.get_helicities(),
                                    &additional_params,
                                    Some(&left_threshold_params),
                                    Some(&right_threshold_params),
                                    Some(&lu_cut_params),
                                );

                                let result_of_this_ct = self.evaluators[cut_group_id]
                                    .iterated_evaluator[iterated_index]
                                    .get_mut(&cut_cff_index)
                                    .unwrap()
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

                                let helper_evaluation = evaluate_threshold_helper_iterated(
                                    ThresholdHelperEvaluationContext {
                                        helper: self.evaluators[cut_group_id]
                                            .threshold_helpers
                                            .iterated[iterated_index]
                                            .get_mut(&cut_cff_index)
                                            .unwrap(),
                                        cut_cff_index: &cut_cff_index,
                                        threshold_result: &result_of_this_ct,
                                        outputs: threshold_helper_outputs,
                                        evaluation_meta_data,
                                        record_primary_timing,
                                    },
                                    &left_threshold_params,
                                    &right_threshold_params,
                                );
                                if let Some(pending) = pending_components.as_mut() {
                                    let helper_completed_pieces = helper_evaluation.into_pieces();
                                    let coefficients = multiplier_coefficients.as_ref().expect(
                                        "detailed components always materialize identity factors",
                                    );
                                    let (left_variant_id, right_variant_id) =
                                        detailed_pair_variant_ids.expect(
                                            "detailed pair component requires stable variant IDs",
                                        );
                                    let registry = self.metadata_registry.as_ref().unwrap();
                                    for (kind, piece, coefficient) in [
                                        (
                                            ThresholdCountertermComponentKind::LocalLocal,
                                            helper_completed_pieces.local_local,
                                            &coefficients.local_local,
                                        ),
                                        (
                                            ThresholdCountertermComponentKind::LocalIntegrated,
                                            helper_completed_pieces.local_integrated,
                                            &coefficients.local_integrated,
                                        ),
                                        (
                                            ThresholdCountertermComponentKind::IntegratedLocal,
                                            helper_completed_pieces.integrated_local,
                                            &coefficients.integrated_local,
                                        ),
                                        (
                                            ThresholdCountertermComponentKind::IntegratedIntegrated,
                                            helper_completed_pieces.integrated_integrated,
                                            &coefficients.integrated_integrated,
                                        ),
                                    ] {
                                        let pass_one = (!coefficient.is_zero()).then(|| {
                                            multiply_dual_or_not_complex(
                                                piece,
                                                &multi_channeling_factor,
                                            )
                                        });
                                        pending.push(PendingLUCountertermComponent {
                                            component_id: registry.component_id(
                                                Some(cut_group_id),
                                                kind,
                                                &[left_variant_id, right_variant_id],
                                            )?,
                                            occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                                                overlap_groups: SmallVec::from_slice(&[
                                                    left_overlap_group_index,
                                                    right_overlap_group_index,
                                                ]),
                                                left_threshold_order: cut_cff_index.left_threshold_order,
                                                right_threshold_order: cut_cff_index.right_threshold_order,
                                                lu_cut_order: cut_cff_index.lu_cut_order,
                                            },
                                            multiplier_values: smallvec![
                                                coefficient.left.clone(),
                                                coefficient.right.clone(),
                                            ],
                                            effective_multiplier: coefficient.effective(),
                                            pass_one,
                                        });
                                    }
                                } else {
                                    let helper_completed_result =
                                        if let Some(coefficients) = &multiplier_coefficients {
                                            weight_iterated_multiplier_pieces(
                                                helper_evaluation.into_pieces(),
                                                coefficients,
                                                order,
                                                &zero,
                                            )
                                        } else if threshold_helper_outputs
                                            == LUThresholdHelperOutputs::Pieces
                                        {
                                            let helper_completed_pieces =
                                                helper_evaluation.into_pieces();
                                            let mut result = helper_completed_pieces.local_local;
                                            result += helper_completed_pieces.local_integrated;
                                            result += helper_completed_pieces.integrated_local;
                                            result += helper_completed_pieces.integrated_integrated;
                                            result
                                        } else {
                                            helper_evaluation.into_legacy()
                                        };

                                    cartesian_product_result += multiply_dual_or_not_complex(
                                        helper_completed_result,
                                        &multi_channeling_factor,
                                    );
                                }
                            }
                        }
                    }
                }
            }

            let cut_esurface = kinematic_point.cut_esurface_for_order(order);
            if let Some(pending_components) = pending_components {
                let completed_components = completed_components
                    .as_mut()
                    .expect("detailed pending components require completed storage");
                for component in pending_components {
                    let (weighted, completed) = component.finish(&zero, |pass_one| {
                        self.evaluators[cut_group_id].evaluate_residue_from_esurface(
                            cut_group_id,
                            order,
                            pass_one,
                            cut_esurface,
                            evaluation_meta_data,
                            record_primary_timing,
                        )
                    });
                    total_result += weighted;
                    completed_components.push(completed);
                }
            } else {
                let mut pass_one_result = left_evaluations;
                pass_one_result += right_evaluations;
                pass_one_result += negate_dual_or_not_complex(cartesian_product_result);
                let pass_one_result = negate_dual_or_not_complex(pass_one_result);
                total_result += self.evaluators[cut_group_id].evaluate_residue_from_esurface(
                    cut_group_id,
                    order,
                    pass_one_result,
                    cut_esurface,
                    evaluation_meta_data,
                    record_primary_timing,
                );
            }
        }

        Ok(LUCountertermEvaluation {
            total: total_result,
            components: completed_components,
        })
    }
}

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    real_mass_vector: &'a EdgeVec<F<T>>,
    e_cm: F<T>,
    graph: &'a Graph,
    subspace: &'a SubspaceData,
    threshold_subspaces: Option<&'a [SubspaceData]>,
    all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
    settings: &'a RuntimeSettings,
    esurface_collection: &'a EsurfaceCollection,
    transformed_kinematic_point: LUCTKinematicPoint<T>,
    probe_rotation: &'a Rotation,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        graph: &'a Graph,
        settings: &'a RuntimeSettings,
        esurface_collection: &'a EsurfaceCollection,
        transformed_kinematic_point: LUCTKinematicPoint<T>,
        overlap_structure: &'a OverlapStructure,
        masses: &'a EdgeVec<F<T>>,
        all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
        subspace: &'a SubspaceData,
        threshold_subspaces: Option<&'a [SubspaceData]>,
        probe_rotation: &'a Rotation,
    ) -> Self {
        let e_cm = F::from_f64(settings.kinematics.e_cm);

        Self {
            real_mass_vector: masses,
            e_cm,
            graph,
            settings,
            esurface_collection,
            overlap_structure,
            transformed_kinematic_point,
            all_lmbs,
            subspace,
            threshold_subspaces,
            probe_rotation,
        }
    }

    fn threshold_subspace(&self, esurface_id: EsurfaceID) -> &SubspaceData {
        self.threshold_subspaces
            .map_or(self.subspace, |subspaces| &subspaces[esurface_id.0])
    }

    fn new_overlap_builder(
        &'a self,
        overlap_group: &'a OverlapGroup,
        selected_esurface: Option<EsurfaceID>,
    ) -> OverlapBuilder<'a, T> {
        let subspace = selected_esurface
            .map(|esurface_id| self.threshold_subspace(esurface_id))
            .unwrap_or(self.subspace);

        let center = overlap_group.center.cast();

        let shifted_loop_momenta = self
            .transformed_kinematic_point
            .representative_sample()
            .loop_moms()
            - &center;

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
            subspace,
            center,
            unit_shifted_momenta,
            radius,
        }
    }
}

struct OverlapBuilder<'a, T: FloatLike> {
    counterterm_builder: &'a CounterTermBuilder<'a, T>,
    overlap_group: &'a OverlapGroup,
    subspace: &'a SubspaceData,
    /// Solver-derived centers already belong to the current probe and cut-side LMB frame.
    center: LoopMomenta<F<T>>,
    unit_shifted_momenta: LoopMomenta<F<T>>,
    radius: F<T>,
}

enum SideOverlapBuilders<'a, T: FloatLike> {
    Homogeneous(Vec<OverlapBuilder<'a, T>>),
    /// One builder per projected E-surface instance, grouped like the overlap result.
    Projected(Vec<Vec<OverlapBuilder<'a, T>>>),
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
        radial_root_identity: &RadialRootIdentity,
        radial_root_diagnostics: &mut RadialRootDiagnostics,
    ) -> Option<RstarSolution<'a, T>> {
        let subspace = self.overlap_builder.subspace;
        let graph = self.overlap_builder.counterterm_builder.graph;
        let lmbs = self.overlap_builder.counterterm_builder.all_lmbs;
        let masses = self.overlap_builder.counterterm_builder.real_mass_vector;

        debug!("subspace: {:?}", subspace);

        let representative_sample = self
            .overlap_builder
            .counterterm_builder
            .transformed_kinematic_point
            .representative_sample();
        let mut center_with_fixed_complement = representative_sample.loop_moms().clone();
        for loop_index in subspace.iter_lmb_indices() {
            center_with_fixed_complement[loop_index] =
                self.overlap_builder.center[loop_index].clone();
        }

        let center_surface_values = self
            .overlap_builder
            .overlap_group
            .existing_esurfaces
            .iter()
            .map(|&existing_esurface_id| {
                let esurface_id = self
                    .overlap_builder
                    .counterterm_builder
                    .overlap_structure
                    .existing_esurfaces[existing_esurface_id];
                let esurface =
                    &self.overlap_builder.counterterm_builder.esurface_collection[esurface_id];
                let surface_subspace = self
                    .overlap_builder
                    .counterterm_builder
                    .threshold_subspace(esurface_id);
                let mut surface_center_with_fixed_complement =
                    representative_sample.loop_moms().clone();
                for loop_index in surface_subspace.iter_lmb_indices() {
                    surface_center_with_fixed_complement[loop_index] =
                        self.overlap_builder.center[loop_index].clone();
                }
                let value = esurface.compute_from_momenta(
                    surface_subspace.get_lmb(lmbs),
                    masses,
                    &surface_center_with_fixed_complement,
                    representative_sample.external_moms(),
                );
                let is_valid = esurface_value_is_strictly_inside(
                    &value,
                    &self.overlap_builder.counterterm_builder.e_cm,
                );
                (esurface_id, value, is_valid)
            })
            .collect_vec();
        let all_center_values_valid = center_surface_values
            .iter()
            .all(|(_, _, is_valid)| *is_valid);

        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #center;
            stage = "lu_threshold_center_values",
            graph = %graph.name,
            selected_esurface_id = self.esurface_id.0,
            rotation_id = %self.overlap_builder.counterterm_builder.probe_rotation.method,
            center_provenance = "current_probe_cut_lmb_frame",
            subspace_loop_indices = ?subspace.iter_lmb_indices().collect_vec(),
            file.active_center = %format!("{}", self.overlap_builder.center),
            file.center_with_fixed_complement = %format!("{}", center_with_fixed_complement),
            file.surface_values = %center_surface_values
                .iter()
                .map(|(esurface_id, value, is_valid)| format!(
                    "local={} value={:+16e} inside={}",
                    esurface_id.0,
                    value,
                    is_valid,
                ))
                .join("\n"),
            "LU threshold center values"
        );

        if !all_center_values_valid {
            warn!(
                graph = %graph.name,
                selected_esurface_id = self.esurface_id.0,
                rotation_id = %self.overlap_builder.counterterm_builder.probe_rotation.method,
                center_provenance = "current_probe_cut_lmb_frame",
                center = %center_with_fixed_complement,
                surface_values = %center_surface_values
                    .iter()
                    .map(|(esurface_id, value, is_valid)| format!(
                        "local={} value={:+16e} inside={}",
                        esurface_id.0,
                        value,
                        is_valid,
                    ))
                    .join("; "),
                "refusing to evaluate an LU threshold counterterm with an invalid probe-frame overlap center"
            );
            return None;
        }

        let (raw_radius_guess, _) = self.esurface.get_radius_guess_subspace(
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
                &self.overlap_builder.center,
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

        let zero = raw_radius_guess.zero();
        let mut radius_guess = raw_radius_guess.clone();
        if radius_guess.is_nan() || radius_guess.is_infinite() || radius_guess <= zero {
            radius_guess = self.overlap_builder.counterterm_builder.e_cm.clone();
        }

        debug!("initial radius guess: {:?}", radius_guess);

        // Some residual is expected when Newton stagnates at the representable value nearest the
        // root. Direct convergence uses the active precision; the shared diagnostics may also
        // certify a residual-limited higher-precision result against the preceding precision.
        let tolerance = F::from_f64(
            self.overlap_builder
                .counterterm_builder
                .settings
                .subtraction
                .radial_root_residual_tolerance,
        );
        let solution = match radial_root_diagnostics.solve(
            radial_root_identity,
            &zero,
            &radius_guess,
            function,
            &tolerance,
            MAX_ITERATIONS,
            64,
            &self.overlap_builder.counterterm_builder.e_cm,
        ) {
            Ok(solution) => solution,
            Err(error) => {
                warn!(
                    graph = %graph.name,
                    esurface_id = self.esurface_id.0,
                    rotation_id = %self.overlap_builder.counterterm_builder.probe_rotation.method,
                    center_provenance = "current_probe_cut_lmb_frame",
                    raw_radius_guess = %raw_radius_guess,
                    radius_guess = %radius_guess,
                    error = %error,
                    "refusing to evaluate a threshold counterterm with an invalid radial solution"
                );
                return None;
            }
        };

        debug!("r* solution: {:?}", solution);

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
                rstar_t_dependence_evaluator.evaluate(RstarTDependenceInput {
                    t_star: &t_star,
                    radius_star: &solution.solution,
                    overlap_center: &self.overlap_builder.center,
                    subspace,
                    unrescaled_momentum_sample: self
                        .overlap_builder
                        .counterterm_builder
                        .transformed_kinematic_point
                        .unrescaled_sample(),
                    masses,
                    threshold_esurface: self.esurface,
                    lmb: &self
                        .overlap_builder
                        .counterterm_builder
                        .graph
                        .loop_momentum_basis,
                    all_lmbs: lmbs,
                }),
            )
        } else {
            None
        };

        Some(RstarSolution {
            esurface_ct_builder: self,
            solution,
            t_dependent_solution,
        })
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
    fn subspace(&self) -> &SubspaceData {
        self.esurface_ct_builder.overlap_builder.subspace
    }

    fn max_supported_order(&self) -> usize {
        self.t_dependent_solution
            .as_ref()
            .map(|dual| dual.values.len().saturating_sub(1))
            .unwrap_or(0)
    }

    fn base_rstar_loop_momenta(&self) -> LoopMomenta<F<T>> {
        let subspace = self.subspace();

        &self
            .esurface_ct_builder
            .overlap_builder
            .unit_shifted_momenta
            .rescale(&self.solution.solution, subspace.as_subspace_simple())
            + &self.esurface_ct_builder.overlap_builder.center
    }

    fn truncated_rstar_solution(&self, order: usize) -> HyperDual<F<T>> {
        let dual_solution = self
            .t_dependent_solution
            .as_ref()
            .expect("higher-order LU threshold evaluation requires cached r_star(t)");
        HyperDual::from_values(
            simple_n_deriv_shape(order),
            dual_solution.values[..=order].to_vec(),
        )
    }

    fn embedded_truncated_rstar_solution(&self, cut_cff_index: &CutCFFIndex) -> HyperDual<F<T>> {
        let lu_order = cut_cff_index
            .lu_cut_order
            .expect("mixed LU geometry requires lu_cut_order")
            - 1;
        let target_shape = shape_from_cut_cff_index(cut_cff_index)
            .map(HyperDual::new)
            .expect("mixed LU geometry requires a dual shape");
        let t_variable = variable_indices_from_cut_cff_index(cut_cff_index).lu_cut;
        let t_dual = if lu_order == 0 {
            HyperDual::from_values(
                simple_n_deriv_shape(0),
                vec![self.solution.solution.clone()],
            )
        } else {
            self.truncated_rstar_solution(lu_order)
        };

        embed_t_dual_in_target_shape(&t_dual, &target_shape, t_variable)
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
            &self.esurface_ct_builder.overlap_builder.center,
        );
        let external_moms = dualize_external_momenta(&radius_star, source_sample.external_moms());

        let shifted_momenta = dual_loop_momenta
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum.clone() - center.clone())
            .collect::<LoopMomenta<_>>();

        let radius = dual_shifted_radius(&shifted_momenta, self.subspace());
        let inverse_radius = new_constant(&radius, &radius.values[0].one()) / radius.clone();
        let unit_shifted_momenta =
            shifted_momenta.rescale(&inverse_radius, self.subspace().as_subspace_simple());

        let rstar_loop_momenta = unit_shifted_momenta
            .rescale(&radius_star, self.subspace().as_subspace_simple())
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
            self.subspace(),
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

    fn dual_geometry_for_cut_cff_index(
        &self,
        cut_cff_index: &CutCFFIndex,
        threshold_variable: Option<usize>,
    ) -> DualRstarGeometry<T> {
        let lu_order = cut_cff_index
            .lu_cut_order
            .expect("mixed LU geometry requires lu_cut_order")
            - 1;
        let source_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .transformed_kinematic_point
            .sample_for_order(lu_order);
        let target_shape = shape_from_cut_cff_index(cut_cff_index)
            .map(HyperDual::new)
            .expect("mixed LU geometry requires a dual shape");
        let t_variable = variable_indices_from_cut_cff_index(cut_cff_index).lu_cut;
        let dual_loop_momenta = embedded_dual_loop_momenta_for_cut_cff_index(
            source_sample,
            &Some(target_shape.clone()),
            t_variable,
        )
        .expect("mixed LU geometry requires dual loop momenta");

        let mut radius_star = self.embedded_truncated_rstar_solution(cut_cff_index);
        activate_threshold_variable_in_target_shape(&mut radius_star, threshold_variable);
        let center = dualize_loop_momenta(
            &radius_star,
            &self.esurface_ct_builder.overlap_builder.center,
        );
        let external_moms = dualize_external_momenta(&radius_star, source_sample.external_moms());

        let shifted_momenta = dual_loop_momenta
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum.clone() - center.clone())
            .collect::<LoopMomenta<_>>();

        let radius = dual_shifted_radius(&shifted_momenta, self.subspace());
        let inverse_radius = new_constant(&radius, &radius.values[0].one()) / radius.clone();
        let unit_shifted_momenta =
            shifted_momenta.rescale(&inverse_radius, self.subspace().as_subspace_simple());

        let rstar_loop_momenta = unit_shifted_momenta
            .rescale(&radius_star, self.subspace().as_subspace_simple())
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
            self.subspace(),
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
        let subspace = self.subspace();

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
        let subspace = self.subspace();
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

    fn rstar_sample_for_cut_cff_index<'solution>(
        &'solution self,
        cut_cff_index: &CutCFFIndex,
        threshold_variable: Option<usize>,
    ) -> RstarSample<'solution, 'a, T> {
        let lu_order = cut_cff_index
            .lu_cut_order
            .expect("LU threshold sampling requires lu_cut_order")
            - 1;

        if shape_from_cut_cff_index(cut_cff_index).is_none() {
            return self.rstar_sample_for_order(lu_order);
        }

        let base_rstar_loop_momenta = self.base_rstar_loop_momenta();
        let geometry = self.dual_geometry_for_cut_cff_index(cut_cff_index, threshold_variable);
        let mut rstar_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .transformed_kinematic_point
            .sample_for_order(lu_order)
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
        let subspace = self.rstar_solution.subspace();
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
    let left_subspace = left_sample.rstar_solution.subspace();

    let right_subspace = right_sample.rstar_solution.subspace();

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

#[cfg(test)]
mod tests {
    use super::{
        IteratedMultiplierCoefficients, IteratedMultiplierComponent, LUCounterTerm,
        PairMultiplierCoefficient, PendingLUCountertermComponent, SingleMultiplierCoefficients,
        SingleMultiplierComponent, extend_extracted_eta_params,
        extract_zero_threshold_coefficient_t_dual, iterated_multiplier_context,
        single_multiplier_context, weight_iterated_multiplier_pieces,
        weight_single_multiplier_pieces,
    };
    use crate::cff::CutCFFIndex;
    use crate::observables::events::ThresholdCountertermComponentOccurrence;
    use crate::processes::{CutGroupId, IteratedThresholdPieces, SingleThresholdPieces};
    use crate::utils::{
        ArbPrec, F, FloatLike, f128,
        hyperdual_utils::{DualOrNot, new_from_values, shape_from_cut_cff_index},
    };
    use spenso::algebra::complex::Complex;
    use symbolica::domains::dual::HyperDual;
    use typed_index_collections::{TiVec, ti_vec};

    #[test]
    fn inactive_cut_group_guard_reports_runtime_access() {
        let counterterm = LUCounterTerm {
            evaluators: TiVec::new(),
            thresholds: TiVec::new(),
            subspaces: TiVec::new(),
            variant_subspaces: None,
            metadata_registry: None,
            rstar_dependence_calculator: TiVec::new(),
            active_cut_groups: ti_vec![false],
            active_left_thresholds: ti_vec![TiVec::new()],
            active_right_thresholds: ti_vec![TiVec::new()],
            active_iterated_thresholds: TiVec::new(),
        };

        let error = counterterm
            .ensure_active_cut_group(CutGroupId(0))
            .unwrap_err();
        assert!(error.to_string().contains("generation marked it inactive"));
    }

    #[test]
    fn zero_threshold_coefficient_lookup_is_explicit() {
        let shape = HyperDual::new(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: Some(2),
                lu_cut_order: Some(2),
            })
            .unwrap(),
        );
        let mixed_dual = new_from_values(
            &shape,
            &[
                F(10.0_f64),
                F(11.0_f64),
                F(20.0_f64),
                F(30.0_f64),
                F(40.0_f64),
                F(21.0_f64),
                F(31.0_f64),
                F(41.0_f64),
            ],
        );

        let zero_threshold = extract_zero_threshold_coefficient_t_dual(&mixed_dual, 0);

        assert_eq!(zero_threshold.values, vec![F(10.0_f64), F(11.0_f64)]);
    }

    #[test]
    fn iterated_eta_parameters_project_onto_their_own_threshold_axis() {
        let cut_cff_index = CutCFFIndex {
            left_threshold_order: Some(1),
            right_threshold_order: Some(2),
            lu_cut_order: Some(1),
        };
        let shape = HyperDual::new(shape_from_cut_cff_index(&cut_cff_index).unwrap());
        let mixed_eta = DualOrNot::Dual(new_from_values(&shape, &[F(10.0_f64), F(20.0_f64)]));

        let mut left_params = Vec::new();
        extend_extracted_eta_params(&mut left_params, &mixed_eta, None, None);
        assert_eq!(
            left_params,
            vec![Complex::new_re(F(10.0_f64))],
            "an order-one left eta must take the zero coefficient of the right threshold axis",
        );

        let mut right_params = Vec::new();
        extend_extracted_eta_params(&mut right_params, &mixed_eta, None, Some(0));
        assert_eq!(
            right_params,
            vec![Complex::new_re(F(10.0_f64)), Complex::new_re(F(20.0_f64)),],
            "the right eta must retain derivatives along its own threshold axis",
        );
    }

    #[test]
    fn multiplier_contexts_follow_expanded_term_semantics() {
        let base = "base";
        let left_root = "left_root";
        let right_root = "right_root";
        let merged_root = "merged_root";

        assert_eq!(
            single_multiplier_context(SingleMultiplierComponent::Local, &base, &left_root),
            (&base, &left_root),
        );
        assert_eq!(
            single_multiplier_context(SingleMultiplierComponent::Integrated, &base, &left_root),
            (&left_root, &left_root),
        );
        assert_eq!(
            iterated_multiplier_context(
                IteratedMultiplierComponent::LocalLocal,
                &base,
                &left_root,
                &right_root,
                &merged_root,
            ),
            (&base, &merged_root),
        );
        assert_eq!(
            iterated_multiplier_context(
                IteratedMultiplierComponent::IntegratedLocal,
                &base,
                &left_root,
                &right_root,
                &merged_root,
            ),
            (&left_root, &merged_root),
        );
        assert_eq!(
            iterated_multiplier_context(
                IteratedMultiplierComponent::LocalIntegrated,
                &base,
                &left_root,
                &right_root,
                &merged_root,
            ),
            (&right_root, &merged_root),
        );
        assert_eq!(
            iterated_multiplier_context(
                IteratedMultiplierComponent::IntegratedIntegrated,
                &base,
                &left_root,
                &right_root,
                &merged_root,
            ),
            (&merged_root, &merged_root),
        );

        // One original term, two single-sided (local/integrated) terms on each side, and four
        // paired terms exhaust the 3x3 expansion used by the GL638 two-sided regression.
        let left_only_terms = 2;
        let right_only_terms = 2;
        let paired_terms = 4;
        assert_eq!(1 + left_only_terms + right_only_terms + paired_terms, 9);
        assert_eq!((1 + left_only_terms) * (1 + right_only_terms), 9);
    }

    #[test]
    fn zero_multiplier_coefficients_skip_unused_pieces() {
        let zero = F(0.0_f64);
        let single = weight_single_multiplier_pieces(
            SingleThresholdPieces {
                local: DualOrNot::NonDual(Complex::new_re(F(f64::NAN))),
                integrated: DualOrNot::NonDual(Complex::new_re(F(3.0))),
            },
            &SingleMultiplierCoefficients {
                local: F(0.0),
                integrated: F(5.0),
            },
            0,
            &zero,
        )
        .unwrap_real();
        assert_eq!(single, Complex::new_re(F(15.0)));

        let iterated = weight_iterated_multiplier_pieces(
            IteratedThresholdPieces {
                local_local: DualOrNot::NonDual(Complex::new_re(F(1.0))),
                local_integrated: DualOrNot::NonDual(Complex::new_re(F(f64::NAN))),
                integrated_local: DualOrNot::NonDual(Complex::new_re(F(100.0))),
                integrated_integrated: DualOrNot::NonDual(Complex::new_re(F(1000.0))),
            },
            &IteratedMultiplierCoefficients {
                local_local: PairMultiplierCoefficient {
                    left: F(2.0),
                    right: F(1.0),
                },
                local_integrated: PairMultiplierCoefficient {
                    left: F(0.0),
                    right: F(1.0),
                },
                integrated_local: PairMultiplierCoefficient {
                    left: F(5.0),
                    right: F(1.0),
                },
                integrated_integrated: PairMultiplierCoefficient {
                    left: F(7.0),
                    right: F(1.0),
                },
            },
            0,
            &zero,
        )
        .unwrap_real();
        assert_eq!(iterated, Complex::new_re(F(7502.0)));

        assert!(
            SingleMultiplierCoefficients {
                local: zero,
                integrated: zero,
            }
            .all_zero()
        );
        assert!(
            IteratedMultiplierCoefficients {
                local_local: PairMultiplierCoefficient {
                    left: zero,
                    right: F(1.0),
                },
                local_integrated: PairMultiplierCoefficient {
                    left: zero,
                    right: F(1.0),
                },
                integrated_local: PairMultiplierCoefficient {
                    left: zero,
                    right: F(1.0),
                },
                integrated_integrated: PairMultiplierCoefficient {
                    left: zero,
                    right: F(1.0),
                },
            }
            .all_zero()
        );
    }

    #[test]
    fn detailed_component_finishing_skips_exact_zero_and_reconciles_weights() {
        let occurrence = ThresholdCountertermComponentOccurrence::LocalUnitarity {
            overlap_groups: smallvec::smallvec![2],
            left_threshold_order: Some(1),
            right_threshold_order: None,
            lu_cut_order: Some(2),
        };
        let nonzero = PendingLUCountertermComponent {
            component_id: 4,
            occurrence: occurrence.clone(),
            multiplier_values: smallvec::smallvec![F(2.0_f64)],
            effective_multiplier: F(2.0),
            // Single-sided signs are already carried by the pass-one component.
            pass_one: Some(DualOrNot::NonDual(Complex::new_re(F(-3.0)))),
        };
        let mut pass_two_calls = 0;
        let (weighted, completed) = nonzero.finish(&F(0.0), |pass_one| {
            pass_two_calls += 1;
            pass_one.unwrap_real() * Complex::new_re(F(4.0))
        });
        assert_eq!(pass_two_calls, 1);
        assert_eq!(completed.bare, Some(Complex::new_re(F(-12.0))));
        assert_eq!(weighted, Complex::new_re(F(-24.0)));
        assert_eq!(completed.weighted, weighted);

        let skipped = PendingLUCountertermComponent {
            component_id: 5,
            occurrence,
            multiplier_values: smallvec::smallvec![F(0.0_f64)],
            effective_multiplier: F(0.0),
            pass_one: None,
        };
        let (skipped_weight, skipped) = skipped.finish(&F(0.0), |_| {
            panic!("exact-zero component must not reach pass two")
        });
        assert!(skipped.evaluation_skipped);
        assert!(skipped.bare.is_none());
        assert_eq!(skipped_weight, Complex::new_re(F(0.0)));

        let decomposition = crate::observables::events::GenericThresholdCountertermEventInfo {
            original: Complex::new_re(F(30.0)),
            components: vec![completed, skipped],
        };
        assert_eq!(decomposition.total(), Complex::new_re(F(6.0)));
    }

    fn assert_detailed_component_precision<T: FloatLike>() {
        let value = |value| F(T::from_f64(value));
        let component = PendingLUCountertermComponent {
            component_id: 0,
            occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                overlap_groups: smallvec::smallvec![0],
                left_threshold_order: Some(1),
                right_threshold_order: None,
                lu_cut_order: Some(1),
            },
            multiplier_values: smallvec::smallvec![value(2.0)],
            effective_multiplier: value(2.0),
            pass_one: Some(DualOrNot::NonDual(Complex::new_re(value(-3.0)))),
        };
        let (_, completed) = component.finish(&value(0.0), |pass_one| {
            pass_one.unwrap_real() * Complex::new_re(value(4.0))
        });

        assert_eq!(completed.bare.unwrap().re.0.into_f64(), -12.0);
        assert_eq!(completed.weighted.re.0.into_f64(), -24.0);
        assert_eq!(completed.effective_multiplier.0.into_f64(), 2.0);
    }

    #[test]
    fn detailed_components_preserve_all_runtime_precision_lanes() {
        assert_detailed_component_precision::<f64>();
        assert_detailed_component_precision::<f128>();
        assert_detailed_component_precision::<ArbPrec>();
    }
}
