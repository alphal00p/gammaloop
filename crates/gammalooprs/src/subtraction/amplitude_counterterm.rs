use std::{collections::BTreeMap, path::Path, slice};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
use symbolica::atom::Atom;

use color_eyre::Result;
use tracing::{debug, instrument};
use typed_index_collections::{TiVec, ti_vec};

use crate::{
    GammaLoopContext,
    cff::{
        CutCFFIndex,
        esurface::{
            Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId, ExistingEsurfaces,
            GroupEsurfaceId, RaisedEsurfaceData, RaisedEsurfaceId,
        },
        expression::OrientationID,
    },
    graph::{FeynmanGraph, Graph, GraphGroupPosition},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            GenericEvaluator, ParamBuilder, ThresholdParams,
            evaluators::{
                EvaluatorStack, SingleOrAllOrientations, evaluate_evaluator,
                evaluate_evaluator_single,
            },
        },
    },
    model::Model,
    momentum::{
        Energy, FourMomentum, Rotation,
        sample::{LoopMomenta, MomentumSample},
    },
    processes::EvaluatorBuildTimings,
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        evaluate_integrated_ct_normalisation, evaluate_uv_damper,
        overlap::{OverlapGroup, OverlapStructure},
    },
    utils::{
        F, FloatLike,
        hyperdual_utils::{
            DualOrNot, extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            shape_from_cut_cff_index, simple_n_deriv_shape,
        },
        newton_solver::{NewtonIterationResult, newton_iteration_and_derivative},
    },
};
use symbolica::domains::dual::HyperDual;

const MAX_ITERATIONS: usize = 40;
const TOLERANCE: f64 = 1.0;

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

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermData {
    pub overlap: OverlapStructure,
    pub evaluators: TiVec<RaisedEsurfaceId, AmplitudeCountertermEvaluator>,
    pub helper_evaluators: Vec<GenericEvaluator>,
    // `generated_mask` tracks whether a threshold slot actually has symbolic content.
    // This is independent of any orientation filtering and reflects the underlying
    // threshold-generation logic itself.
    pub generated_mask: TiVec<RaisedEsurfaceId, bool>,
    // `active_mask` is an additional generation-time gate coming from the selected
    // orientation subset. A slot can be generated in principle but inactive for the
    // current generated evaluator set, in which case we keep its index but compile it
    // to a zero/dummy evaluator and hard-fail if runtime ever reaches it.
    pub active_mask: TiVec<RaisedEsurfaceId, bool>,
    pub raised_data: RaisedEsurfaceData,
    pub esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
    pub own_group_position: GraphGroupPosition,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermAtom {
    pub parametric: BTreeMap<CutCFFIndex, Atom>,
}

impl AmplitudeCountertermAtom {
    pub(crate) fn is_generated(&self) -> bool {
        !self.parametric.is_empty()
    }

    pub(crate) fn zero_like(&self) -> Self {
        Self {
            parametric: self
                .parametric
                .iter()
                .map(|(k, _)| (*k, Atom::Zero))
                .collect(),
        }
    }

    #[instrument(
           skip_all,
           fields(indicatif.pb_show = true, indicatif.pb_msg = "Building Threshold CT Evaluator"),
       )]
    pub(crate) fn to_evaluator_with_timings(
        &self,
        param_builder: &ParamBuilder,
        orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
        global_settings: &GlobalSettings,
    ) -> (AmplitudeCountertermEvaluator, EvaluatorBuildTimings) {
        let mut evaluator_stacks = BTreeMap::new();
        let mut timings = EvaluatorBuildTimings::default();

        for (index, integrand) in self.parametric.iter() {
            let dual_shape = shape_from_cut_cff_index(index);

            let (evaluator_stack, evaluator_timings) = EvaluatorStack::new_with_timings(
                slice::from_ref(integrand),
                param_builder,
                orientations.as_slice().as_ref(),
                dual_shape,
                &global_settings.generation.evaluator,
            )
            .unwrap();
            timings += evaluator_timings;
            evaluator_stacks.insert(*index, evaluator_stack);
        }

        (AmplitudeCountertermEvaluator { evaluator_stacks }, timings)
    }

    pub(crate) fn new() -> Self {
        Self {
            parametric: BTreeMap::new(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermEvaluator {
    pub evaluator_stacks: BTreeMap<CutCFFIndex, EvaluatorStack>,
}

impl AmplitudeCountertermEvaluator {
    pub(crate) fn generic_evaluator_count(&self) -> usize {
        self.evaluator_stacks
            .values()
            .map(EvaluatorStack::generic_evaluator_count)
            .sum()
    }
}

#[derive(Debug, Clone)]
pub struct AmplitudeLocalCountertermEvaluation<T: FloatLike> {
    pub esurface_id: RaisedEsurfaceId,
    pub overlap_group: usize,
    pub value: Complex<F<T>>,
}

#[derive(Debug, Clone)]
pub struct AmplitudeCountertermEvaluation<T: FloatLike> {
    pub total: Complex<F<T>>,
    pub local_counterterms: Vec<AmplitudeLocalCountertermEvaluation<T>>,
}

impl AmplitudeCountertermData {
    pub(crate) fn generic_evaluator_count(&self) -> usize {
        let stack_count = self
            .evaluators
            .iter()
            .map(AmplitudeCountertermEvaluator::generic_evaluator_count)
            .sum::<usize>();
        let overlap_count = self
            .overlap
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .prefactor_evaluator
                    .as_ref()
                    .map(|evaluators| evaluators.len())
                    .unwrap_or(0)
            })
            .sum::<usize>();
        stack_count + overlap_count + self.helper_evaluators.len()
    }

    pub fn new_empty(own_group_position: GraphGroupPosition) -> Self {
        Self {
            overlap: OverlapStructure::new_empty(),
            evaluators: TiVec::new(),
            helper_evaluators: vec![],
            generated_mask: TiVec::new(),
            active_mask: TiVec::new(),
            raised_data: RaisedEsurfaceData {
                raised_groups: TiVec::new(),
                pass_two_evaluator: None,
            },
            esurface_map: TiVec::new(),
            own_group_position,
        }
    }

    fn ensure_active_raised_esurface(&self, raised_esurface_id: RaisedEsurfaceId) -> Result<()> {
        if self.active_mask[raised_esurface_id] {
            return Ok(());
        }

        Err(color_eyre::eyre::eyre!(
            "Amplitude threshold evaluator {} was reached at runtime even though generation marked it inactive for the selected orientation subset",
            raised_esurface_id.0
        ))
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        _override_existing: bool,
        frozen_mode: &crate::settings::global::FrozenCompilationMode,
    ) -> Result<()> {
        for (i, e) in self.evaluators.iter_mut_enumerated() {
            for (cff_index, evaluator_stack) in e.evaluator_stacks.iter_mut() {
                let order_index = cff_index.left_threshold_order.unwrap() - 1;
                evaluator_stack.compile(
                    format!("esurface_{}_order_{}", i.0, order_index + 1),
                    path.as_ref(),
                    frozen_mode,
                )?;
            }
        }

        for (group_index, group) in self.overlap.overlap_groups.iter_mut().enumerate() {
            if let Some(prefactor_evaluators) = group.prefactor_evaluator.as_mut() {
                for (order_index, prefactor_evaluator) in
                    prefactor_evaluators.iter_mut().enumerate()
                {
                    prefactor_evaluator.borrow_mut().compile_external(
                        path.as_ref()
                            .join(format!(
                                "overlap_prefactor_{group_index}_order_{}",
                                order_index + 1
                            ))
                            .with_extension("cpp"),
                        format!("overlap_prefactor_{group_index}_order_{}", order_index + 1),
                        path.as_ref()
                            .join(format!(
                                "overlap_prefactor_{group_index}_order_{}",
                                order_index + 1
                            ))
                            .with_extension("so"),
                        frozen_mode,
                    )?;
                }
            }
        }

        for (order_index, evaluator) in self.helper_evaluators.iter_mut().enumerate() {
            evaluator.compile_external(
                path.as_ref()
                    .join(format!("threshold_helper_{order_index}"))
                    .with_extension("cpp"),
                format!("threshold_helper_{order_index}"),
                path.as_ref()
                    .join(format!("threshold_helper_{order_index}"))
                    .with_extension("so"),
                frozen_mode,
            )?;
        }
        Ok(())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        for evaluator in self.evaluators.iter_mut() {
            for (_cut_index, evaluator_stack) in &mut evaluator.evaluator_stacks {
                evaluator_stack.for_each_generic_evaluator_mut(&mut f)?;
            }
        }

        for group in &mut self.overlap.overlap_groups {
            if let Some(prefactor_evaluators) = group.prefactor_evaluator.as_mut() {
                for prefactor_evaluator in prefactor_evaluators.iter_mut() {
                    f(prefactor_evaluator.get_mut())?;
                }
            }
        }

        for evaluator in &mut self.helper_evaluators {
            f(evaluator)?;
        }

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        model: &Model,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder<f64>,
        orientation: SingleOrAllOrientations<'_, OrientationID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<AmplitudeCountertermEvaluation<T>> {
        debug!("start evaluate threshold counterterm");
        let existing_esurfaces = self
            .overlap
            .existing_esurfaces
            .iter()
            .map(|e| e.0)
            .collect::<Vec<_>>();
        debug!("subtracting esurfaces: {:?}", existing_esurfaces);
        debug!("overlap structure\n: {}", self.overlap);

        let counter_term_builder = CounterTermBuilder::new(
            graph,
            model,
            rotation,
            settings,
            esurfaces,
            momentum_sample,
            &self.overlap,
            &self.raised_data,
            self.own_group_position,
            &self.esurface_map,
        );

        let mut result = Complex::new_re(momentum_sample.zero());
        let mut local_counterterms = Vec::new();

        for (overlap_group, group) in self.overlap.overlap_groups.iter().enumerate() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            for existing_esurface_id in group.existing_esurfaces.iter() {
                let Some(esurface_builder) =
                    overlap_builder.new_esurface_builder(*existing_esurface_id)
                else {
                    continue;
                };

                let raised_esurface_id = esurface_builder.raised_esurface_id;
                self.ensure_active_raised_esurface(raised_esurface_id)?;
                let single_result = esurface_builder.solve_rstar().rstar_samples().evaluate(
                    param_builder,
                    orientation,
                    evaluation_metadata,
                    record_primary_timing,
                    &mut self.evaluators[raised_esurface_id],
                    &mut self.helper_evaluators,
                )?;

                if !single_result.is_zero() {
                    //    debug!(
                    //        "Param Builder for {}:\n{}",
                    //        existing_esurface_id, param_builder
                    //    );
                    debug!(
                        "Counterterm for esurface {}: {:+16e}",
                        existing_esurface_id, single_result
                    );
                }

                local_counterterms.push(AmplitudeLocalCountertermEvaluation {
                    esurface_id: raised_esurface_id,
                    overlap_group,
                    value: single_result.clone(),
                });
                result += single_result;
            }
        }
        Ok(AmplitudeCountertermEvaluation {
            total: result,
            local_counterterms,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn kinematics_for_approach<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        model: &Model,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
    ) -> Result<OverlapStructureWithKinematics<T>> {
        let counter_term_builder = CounterTermBuilder::new(
            graph,
            model,
            rotation,
            settings,
            esurfaces,
            momentum_sample,
            &self.overlap,
            &self.raised_data,
            self.own_group_position,
            &self.esurface_map,
        );

        let mut overlap_sturcture_with_kinematics = OverlapStructureWithKinematics {
            existing_esurfaces: self.overlap.existing_esurfaces.clone(),
            overlap_groups_with_kinematics: Vec::new(),
        };

        for group in self.overlap.overlap_groups.iter() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            let mut loop_momenta_at_esurface: TiVec<ExistingEsurfaceId, Option<MomentumSample<T>>> =
                ti_vec![];
            for existing_esurface_id in group.existing_esurfaces.iter() {
                let single_result = overlap_builder
                    .new_esurface_builder(*existing_esurface_id)
                    .map(|esurface_builder| -> Result<_> {
                        self.ensure_active_raised_esurface(esurface_builder.raised_esurface_id)?;
                        let rstar_sample = esurface_builder.solve_rstar().rstar_samples();
                        Result::Ok(rstar_sample.rstar_sample)
                    })
                    .transpose()?;

                loop_momenta_at_esurface.push(single_result);
            }

            overlap_sturcture_with_kinematics
                .overlap_groups_with_kinematics
                .push(OverlapGroupWithKinematics {
                    overlap_group: group.clone(),
                    loop_momenta_at_esurface,
                });
        }

        Ok(overlap_sturcture_with_kinematics)
    }
}

pub struct OverlapGroupWithKinematics<T: FloatLike> {
    pub overlap_group: OverlapGroup,
    pub loop_momenta_at_esurface: TiVec<ExistingEsurfaceId, Option<MomentumSample<T>>>,
}

pub struct OverlapStructureWithKinematics<T: FloatLike> {
    pub existing_esurfaces: ExistingEsurfaces,
    pub overlap_groups_with_kinematics: Vec<OverlapGroupWithKinematics<T>>,
}

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    raised_data: &'a RaisedEsurfaceData,
    own_group_position: GraphGroupPosition,
    esurface_map: &'a TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
    real_mass_vector: EdgeVec<F<T>>,
    e_cm: F<T>,
    graph: &'a Graph,
    rotation_for_overlap: &'a Rotation,
    settings: &'a RuntimeSettings,
    esurface_collection: &'a EsurfaceCollection,
    sample: &'a MomentumSample<T>,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        graph: &'a Graph,
        model: &'a Model,
        rotation_for_overlap: &'a Rotation,
        settings: &'a RuntimeSettings,
        esurface_collection: &'a EsurfaceCollection,
        sample: &'a MomentumSample<T>,
        overlap_structure: &'a OverlapStructure,
        raised_data: &'a RaisedEsurfaceData,
        own_group_position: GraphGroupPosition,
        esurface_map: &'a TiVec<
            GroupEsurfaceId,
            TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>,
        >,
    ) -> Self {
        let real_mass_vector = graph.get_real_mass_vector(model);
        let e_cm = F::from_f64(settings.kinematics.e_cm);

        Self {
            real_mass_vector,
            e_cm,
            graph,
            rotation_for_overlap,
            settings,
            esurface_collection,
            overlap_structure,
            raised_data,
            sample,
            own_group_position,
            esurface_map,
        }
    }

    fn new_overlap_builder(&'a self, overlap_group: &'a OverlapGroup) -> OverlapBuilder<'a, T> {
        let center = &overlap_group.center;

        let (unrotated_center, rotated_center) = (
            center.cast(),
            center.rotate(self.rotation_for_overlap).cast(),
        );

        let shifted_loop_momenta = self.sample.loop_moms() - &rotated_center;
        let radius = shifted_loop_momenta.hyper_radius_squared(None).sqrt();
        let unit_shifted_momenta = shifted_loop_momenta.rescale(&radius.inv(), None);

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
    ) -> Option<EsurfaceCTBuilder<'a, T>> {
        let group_esurface_id = self
            .counterterm_builder
            .overlap_structure
            .existing_esurfaces[existing_esurface_id];

        let raised_esurface_id = self.counterterm_builder.esurface_map[group_esurface_id]
            [self.counterterm_builder.own_group_position];

        raised_esurface_id.map(|raised_esurface_id| {
            let esurface_id = self.counterterm_builder.raised_data.raised_groups
                [raised_esurface_id]
                .esurface_ids[0];

            EsurfaceCTBuilder {
                overlap_builder: self,
                _existing_esurface_id: existing_esurface_id,
                esurface: &self.counterterm_builder.esurface_collection[esurface_id],
                esurface_id,
                raised_esurface_id,
            }
        })
    }
}

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    _existing_esurface_id: ExistingEsurfaceId,
    esurface: &'a Esurface,
    esurface_id: EsurfaceID,
    raised_esurface_id: RaisedEsurfaceId,
}

impl<'a, T: FloatLike> EsurfaceCTBuilder<'a, T> {
    fn solve_rstar(self) -> RstarSolution<'a, T> {
        let (radius_guess, _) = self.esurface.get_radius_guess(
            &self.overlap_builder.unit_shifted_momenta,
            self.overlap_builder
                .counterterm_builder
                .sample
                .external_moms(),
            &self
                .overlap_builder
                .counterterm_builder
                .graph
                .loop_momentum_basis,
        );

        let function = |r: &_| {
            self.esurface.compute_self_and_r_derivative(
                r,
                &self.overlap_builder.unit_shifted_momenta,
                &self.overlap_builder.rotated_center,
                self.overlap_builder
                    .counterterm_builder
                    .sample
                    .external_moms(),
                &self.overlap_builder.counterterm_builder.real_mass_vector,
                &self
                    .overlap_builder
                    .counterterm_builder
                    .graph
                    .loop_momentum_basis,
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
        let rstar_loop_momenta = &self
            .esurface_ct_builder
            .overlap_builder
            .unit_shifted_momenta
            .rescale(&self.solution.solution, None)
            + &self.esurface_ct_builder.overlap_builder.rotated_center;

        let mut rstar_sample = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .sample
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

impl<'a, T: FloatLike> RstarSample<'a, T> {
    fn evaluate<'b, 'c: 'b>(
        self,
        param_builder: &mut ParamBuilder<f64>,
        orientations: SingleOrAllOrientations<'a, OrientationID>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        ct_evaluator: &mut AmplitudeCountertermEvaluator,
        helper_evaluators: &mut [GenericEvaluator],
    ) -> Result<Complex<F<T>>> {
        let esurface_ct_builder = &self.rstar_solution.esurface_ct_builder;
        let ct_builder = esurface_ct_builder.overlap_builder.counterterm_builder;

        let esurface_id = esurface_ct_builder.esurface_id;

        let model_params = param_builder
            .model_values()
            .iter()
            .map(|c| Complex::new(F::from_ff64(c.re), F::from_ff64(c.im)))
            .collect::<Vec<_>>();

        let radius = self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .radius
            .clone();

        let radius_star = self.rstar_solution.solution.solution.clone();
        let e_cm = &ct_builder.e_cm;
        let settings = &ct_builder
            .settings
            .subtraction
            .local_ct_settings
            .uv_localisation;

        let integrated_settings = &ct_builder.settings.subtraction.integrated_ct_settings;

        let uv_damp_plus = evaluate_uv_damper(&radius, &radius_star, e_cm, settings);
        let uv_damp_minus = evaluate_uv_damper(&-&radius, &radius_star, e_cm, settings);

        debug!("uv_damp_plus: {:?}", uv_damp_plus);
        debug!("uv_damp_minus: {:?}", uv_damp_minus);
        debug!("radius: {:?}", radius);
        debug!("radius_star: {:?}", radius_star);

        let h_function =
            evaluate_integrated_ct_normalisation(&radius, &radius_star, e_cm, integrated_settings);

        let mut total_ct = Complex::new_re(self.rstar_sample.zero());

        for (cut_cff_index, evaluator_stack) in ct_evaluator.evaluator_stacks.iter_mut() {
            let order_index = cut_cff_index.left_threshold_order.unwrap() - 1;
            let (sample_for_order, threshold_params) = if order_index == 0 {
                debug!(
                    "rescaled loop momenta at rstar:\n{}",
                    self.rstar_sample.loop_moms()
                );
                (
                    self.rstar_sample.clone(),
                    ThresholdParams {
                        radius: DualOrNot::NonDual(radius.clone()),
                        radius_star: DualOrNot::NonDual(radius_star.clone()),
                        esurface_derivative: DualOrNot::NonDual(
                            self.rstar_solution.solution.derivative_at_solution.clone(),
                        ),
                        uv_damp_plus: DualOrNot::NonDual(uv_damp_plus.clone()),
                        uv_damp_minus: DualOrNot::NonDual(uv_damp_minus.clone()),
                        h_function: DualOrNot::NonDual(h_function.clone()),
                    },
                )
            } else {
                let dual_shape = HyperDual::<F<T>>::new(simple_n_deriv_shape(order_index));
                let dual_radius_star = dual_shape.variable(0, radius_star.clone());

                let dualized_center = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .rotated_center
                    .iter()
                    .map(|momentum| {
                        momentum.map_ref(&|value| new_constant(&dual_radius_star, value))
                    })
                    .collect::<LoopMomenta<_>>();
                let dual_loop_momenta = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .unit_shifted_momenta
                    .rescale_with_hyper_dual(&dual_radius_star, None)
                    .iter()
                    .zip(dualized_center.iter())
                    .map(|(momentum, center)| momentum.clone() + center.clone())
                    .collect::<LoopMomenta<_>>();
                debug!("rescaled loop momenta at rstar:\n{}", dual_loop_momenta);
                let mut sample_with_duals = self.rstar_sample.clone();
                sample_with_duals.sample.dual_loop_moms = Some(dual_loop_momenta);

                (
                    sample_with_duals,
                    ThresholdParams {
                        radius: DualOrNot::Dual(new_constant(&dual_radius_star, &radius)),
                        radius_star: DualOrNot::Dual(dual_radius_star.clone()),
                        esurface_derivative: DualOrNot::Dual(new_constant(
                            &dual_radius_star,
                            &self.rstar_solution.solution.derivative_at_solution.clone(),
                        )),
                        uv_damp_plus: DualOrNot::Dual(new_constant(
                            &dual_radius_star,
                            &uv_damp_plus,
                        )),
                        uv_damp_minus: DualOrNot::Dual(new_constant(
                            &dual_radius_star,
                            &uv_damp_minus,
                        )),
                        h_function: DualOrNot::Dual(new_constant(&dual_radius_star, &h_function)),
                    },
                )
            };

            let esurface_derivatives = if order_index == 0 {
                DualOrNot::NonDual(self.rstar_solution.solution.derivative_at_solution.clone())
            } else {
                let dual_shape_for_esurface =
                    HyperDual::<F<T>>::new(simple_n_deriv_shape(order_index + 1));
                let dual_radius_star_for_esurface =
                    dual_shape_for_esurface.variable(0, radius_star.clone());
                let dualized_center_for_esurface = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .rotated_center
                    .iter()
                    .map(|momentum| {
                        momentum
                            .map_ref(&|value| new_constant(&dual_radius_star_for_esurface, value))
                    })
                    .collect::<LoopMomenta<_>>();
                let dual_loop_momenta_for_esurface = self
                    .rstar_solution
                    .esurface_ct_builder
                    .overlap_builder
                    .unit_shifted_momenta
                    .rescale_with_hyper_dual(&dual_radius_star_for_esurface, None)
                    .iter()
                    .zip(dualized_center_for_esurface.iter())
                    .map(|(momentum, center)| momentum.clone() + center.clone())
                    .collect::<LoopMomenta<_>>();

                let dualized_externals = self
                    .rstar_sample
                    .external_moms()
                    .iter()
                    .map(|momentum| FourMomentum {
                        temporal: Energy {
                            value: new_constant(
                                &dual_radius_star_for_esurface,
                                &momentum.temporal.value,
                            ),
                        },
                        spatial: momentum
                            .spatial
                            .map_ref(&|value| new_constant(&dual_radius_star_for_esurface, value)),
                    })
                    .collect();

                DualOrNot::Dual(
                    self.rstar_solution
                        .esurface_ct_builder
                        .esurface
                        .compute_from_dual_momenta(
                            &ct_builder.graph.loop_momentum_basis,
                            &ct_builder.real_mass_vector,
                            &dual_loop_momenta_for_esurface,
                            &dualized_externals,
                        ),
                )
            };

            let params = T::get_parameters(
                param_builder,
                (false, false),
                ct_builder.graph,
                &sample_for_order,
                ct_builder.settings.kinematics.externals.get_helicities(),
                &ct_builder.settings.additional_params(),
                Some(&threshold_params),
                None,
                None,
            );

            let pass_one_result = evaluator_stack
                .evaluate(
                    params,
                    orientations,
                    ct_builder.settings,
                    evaluation_metadata,
                    record_primary_timing,
                )
                .expect("Amplitude counterterm evaluator stack failed")
                .pop()
                .unwrap();

            let prefactor = self.evaluate_multichanneling_prefactor(
                &sample_for_order,
                &model_params,
                evaluation_metadata,
                record_primary_timing,
                order_index,
            );

            let pass_one_result = multiply_dual_or_not_complex(pass_one_result, &prefactor);

            debug!(
                "Pass one result for esurface {} order {}: {}",
                esurface_id.0,
                order_index + 1,
                pass_one_result
            );

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

            match esurface_derivatives {
                DualOrNot::Dual(dual_e_surface) => {
                    extract_t_derivatives(dual_e_surface)[1..]
                        .iter()
                        .for_each(|value| params_for_pass_two.push(Complex::new_re(value.clone())));
                }
                DualOrNot::NonDual(non_dual_e_surface) => {
                    params_for_pass_two.push(Complex::new_re(non_dual_e_surface));
                }
            }

            params_for_pass_two.push(Complex::new_re(radius.clone()));
            params_for_pass_two.push(Complex::new_re(radius_star.clone()));
            params_for_pass_two.push(Complex::new_re(uv_damp_plus.clone()));
            params_for_pass_two.push(Complex::new_re(uv_damp_minus.clone()));
            params_for_pass_two.push(Complex::new_re(h_function.clone()));

            let pass_two_result = evaluate_evaluator_single(
                &mut helper_evaluators[order_index],
                &params_for_pass_two,
                evaluation_metadata,
                record_primary_timing,
            );

            debug!(
                "Pass two result for esurface {} order {}: {}",
                esurface_id.0,
                order_index + 1,
                pass_two_result
            );

            total_ct += pass_two_result;
        }

        debug!(
            ct_eval = format!("{:+16e}", total_ct),
            "esurface {}", esurface_id.0
        );

        Ok(total_ct)
    }

    fn evaluate_multichanneling_prefactor(
        &self,
        momentum_sample: &MomentumSample<T>,
        model_params: &[Complex<F<T>>],
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        order_index: usize,
    ) -> DualOrNot<Complex<F<T>>> {
        let overlap_builder = self.rstar_solution.esurface_ct_builder.overlap_builder;
        let overlap = overlap_builder.counterterm_builder.overlap_structure;

        if overlap.overlap_groups.len() < 2 {
            return DualOrNot::NonDual(Complex::new_re(momentum_sample.one()));
        }

        let multiplicative_offset = momentum_sample
            .sample
            .dual_loop_moms
            .as_ref()
            .map(|dual_loop_moms| dual_loop_moms.first().unwrap().px.values.len())
            .unwrap_or(1);
        let zero = Complex::new_re(momentum_sample.zero());
        let mut params = if let Some(dual_loop_moms) = &momentum_sample.sample.dual_loop_moms {
            dual_loop_moms
                .iter()
                .flat_map(|mom| {
                    [
                        mom.px.values.clone(),
                        mom.py.values.clone(),
                        mom.pz.values.clone(),
                    ]
                    .into_iter()
                    .flatten()
                    .map(Complex::new_re)
                    .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        } else {
            momentum_sample
                .loop_moms()
                .iter()
                .flat_map(|momentum| {
                    [
                        momentum.px.clone(),
                        momentum.py.clone(),
                        momentum.pz.clone(),
                    ]
                    .into_iter()
                    .map(Complex::new_re)
                    .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };

        params.extend(momentum_sample.external_moms().iter().flat_map(|momentum| {
            [
                momentum.temporal.value.clone(),
                momentum.spatial.px.clone(),
                momentum.spatial.py.clone(),
                momentum.spatial.pz.clone(),
            ]
            .into_iter()
            .flat_map(|value| {
                std::iter::once(Complex::new_re(value))
                    .chain((1..multiplicative_offset).map(|_| zero.clone()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
        }));
        params.extend(model_params.iter().cloned().flat_map(|value| {
            std::iter::once(value)
                .chain((1..multiplicative_offset).map(|_| zero.clone()))
                .collect::<Vec<_>>()
        }));

        let evaluator = overlap_builder
            .overlap_group
            .prefactor_evaluator
            .as_ref()
            .unwrap()
            .get(order_index)
            .expect("missing overlap prefactor evaluator for amplitude threshold order");

        evaluate_evaluator(
            &mut evaluator.borrow_mut(),
            &params,
            evaluation_metadata,
            record_primary_timing,
        )
        .pop()
        .expect("overlap prefactor evaluator should return exactly one value")
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{AmplitudeCountertermAtom, AmplitudeCountertermData};
    use crate::{
        cff::{CutCFFIndex, esurface::RaisedEsurfaceId},
        graph::GraphGroupPosition,
    };
    use symbolica::{atom::Atom, symbol};
    use typed_index_collections::ti_vec;

    #[test]
    fn empty_amplitude_counterterm_atom_is_not_generated() {
        let atom = AmplitudeCountertermAtom {
            parametric: BTreeMap::new(),
        };

        assert!(!atom.is_generated());
    }

    #[test]
    fn non_empty_amplitude_counterterm_atom_is_generated() {
        let atom = AmplitudeCountertermAtom {
            parametric: BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::var(symbol!("x")))]),
        };

        assert!(atom.is_generated());
    }

    #[test]
    fn zero_like_amplitude_counterterm_atom_is_zero() {
        let atom = AmplitudeCountertermAtom {
            parametric: BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::var(symbol!("x")))]),
        };

        let zeroed = atom.zero_like();

        assert_eq!(
            zeroed.parametric,
            BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::Zero)])
        );
    }

    #[test]
    fn inactive_amplitude_esurface_guard_reports_runtime_access() {
        let mut data = AmplitudeCountertermData::new_empty(GraphGroupPosition(0));
        data.active_mask = ti_vec![false];

        let error = data
            .ensure_active_raised_esurface(RaisedEsurfaceId(0))
            .unwrap_err();
        assert!(error.to_string().contains("generation marked it inactive"));
    }
}
