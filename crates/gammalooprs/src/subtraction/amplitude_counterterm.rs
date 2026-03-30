use std::{cell::RefCell, path::Path};

use bincode_trait_derive::{Decode, Encode};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
use symbolica::atom::Atom;

use color_eyre::Result;
use tracing::{debug, instrument};
use typed_index_collections::TiVec;

use crate::{
    GammaLoopContext,
    cff::{
        esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId, GroupEsurfaceId},
        expression::OrientationID,
    },
    graph::{FeynmanGraph, Graph, GraphGroupPosition},
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            ParamBuilder, ThresholdParams,
            evaluators::{EvaluatorStack, SingleOrAllOrientations, evaluate_evaluator_single},
        },
    },
    model::Model,
    momentum::{
        Rotation,
        sample::{LoopMomenta, MomentumSample},
    },
    settings::{GlobalSettings, RuntimeSettings},
    subtraction::{
        evaluate_integrated_ct_normalisation, evaluate_uv_damper,
        overlap::{OverlapGroup, OverlapStructure},
    },
    utils::{
        F, FloatLike,
        hyperdual_utils::DualOrNot,
        newton_solver::{NewtonIterationResult, newton_iteration_and_derivative},
    },
};

const MAX_ITERATIONS: usize = 40;
const TOLERANCE: f64 = 1.0;

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermData {
    pub overlap: OverlapStructure,
    pub evaluators: TiVec<EsurfaceID, AmplitudeCountertermEvaluator>,
    pub generated_mask: TiVec<EsurfaceID, bool>,
    pub esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
    pub own_group_position: GraphGroupPosition,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermAtom {
    pub parametric_local: Atom,
    pub parametric_integrated: Atom,
}

impl AmplitudeCountertermAtom {
    pub(crate) fn is_generated(&self) -> bool {
        self.parametric_local != Atom::new() || self.parametric_integrated != Atom::new()
    }

    #[instrument(
           skip_all,
           fields(indicatif.pb_show = true, indicatif.pb_msg = "Building Threshold CT Evaluator"),
       )]
    pub(crate) fn to_evaluator(
        &self,
        param_builder: &ParamBuilder,
        orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
        global_settings: &GlobalSettings,
    ) -> AmplitudeCountertermEvaluator {
        AmplitudeCountertermEvaluator {
            evaluator_stack: RefCell::new(
                EvaluatorStack::new(
                    &[&self.parametric_local, &self.parametric_integrated],
                    param_builder,
                    orientations.as_slice().as_ref(),
                    None,
                    &global_settings.generation.evaluator,
                )
                .unwrap(),
            ),
        }
    }

    pub(crate) fn new() -> Self {
        Self {
            parametric_local: Atom::new(),
            parametric_integrated: Atom::new(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermEvaluator {
    pub evaluator_stack: RefCell<EvaluatorStack>,
}

#[derive(Debug, Clone)]
pub struct AmplitudeCountertermEvaluation<T: FloatLike> {
    pub total: Complex<F<T>>,
    pub local_counterterms: Vec<Complex<F<T>>>,
}

impl AmplitudeCountertermData {
    pub fn new_empty(own_group_position: GraphGroupPosition) -> Self {
        Self {
            overlap: OverlapStructure::new_empty(),
            evaluators: TiVec::new(),
            generated_mask: TiVec::new(),
            esurface_map: TiVec::new(),
            own_group_position,
        }
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        _override_existing: bool,
        settings: &GlobalSettings,
    ) -> Result<()> {
        for (i, e) in self.evaluators.iter_mut_enumerated() {
            let mut evaluator_stack = e.evaluator_stack.borrow_mut();
            evaluator_stack.compile(format!("esurface_{}", i.0), path.as_ref(), settings)?;
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
    ) -> AmplitudeCountertermEvaluation<T> {
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
            &self.evaluators,
            self.own_group_position,
            &self.esurface_map,
        );

        let mut result = Complex::new_re(momentum_sample.zero());
        let mut local_counterterms = Vec::new();

        for group in self.overlap.overlap_groups.iter() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            for existing_esurface_id in group.existing_esurfaces.iter() {
                let single_result = overlap_builder
                    .new_esurface_builder(*existing_esurface_id)
                    .map(|esurface_builder| {
                        esurface_builder.solve_rstar().rstar_samples().evaluate(
                            param_builder,
                            orientation,
                            evaluation_metadata,
                            record_primary_timing,
                        )
                    })
                    .unwrap_or_else(|| Complex::new_re(momentum_sample.zero()));

                if !single_result.is_zero() {
                    debug!(
                        "Param Builder for {}:\n{}",
                        existing_esurface_id, param_builder
                    );
                    debug!(
                        "Counterterm for esurface {}: {:+16e}",
                        existing_esurface_id, single_result
                    );
                }

                local_counterterms.push(single_result.clone());
                result += single_result;
            }
        }
        AmplitudeCountertermEvaluation {
            total: result,
            local_counterterms,
        }
    }
}

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    evaluators: &'a TiVec<EsurfaceID, AmplitudeCountertermEvaluator>,
    own_group_position: GraphGroupPosition,
    esurface_map: &'a TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
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
        evaluators: &'a TiVec<EsurfaceID, AmplitudeCountertermEvaluator>,
        own_group_position: GraphGroupPosition,
        esurface_map: &'a TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
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
            sample,
            evaluators,
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

        let esurface_id = self.counterterm_builder.esurface_map[group_esurface_id]
            [self.counterterm_builder.own_group_position];

        esurface_id.map(|esurface_id| EsurfaceCTBuilder {
            overlap_builder: self,
            _existing_esurface_id: existing_esurface_id,
            esurface: &self.counterterm_builder.esurface_collection[esurface_id],
            esurface_id,
        })
    }
}

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    _existing_esurface_id: ExistingEsurfaceId,
    esurface: &'a Esurface,
    esurface_id: EsurfaceID,
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
    ) -> Complex<F<T>> {
        let esurface_ct_builder = &self.rstar_solution.esurface_ct_builder;
        let ct_builder = esurface_ct_builder.overlap_builder.counterterm_builder;

        let esurface_id = esurface_ct_builder.esurface_id;

        let model_params = param_builder
            .model_values()
            .iter()
            .map(|c| Complex::new(F::from_ff64(c.re), F::from_ff64(c.im)))
            .collect::<Vec<_>>();

        let prefactor = self.evaluate_multichanneling_prefactor(
            &self.rstar_sample,
            model_params,
            evaluation_metadata,
            record_primary_timing,
        );

        let radius = self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .radius
            .clone();

        let radius_star = self.rstar_solution.solution.solution;
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

        let threshold_params = ThresholdParams {
            radius,
            radius_star,
            esurface_derivative: self.rstar_solution.solution.derivative_at_solution,

            uv_damp_plus,
            uv_damp_minus,
            h_function,
        };

        let params = T::get_parameters(
            param_builder,
            (false, false),
            ct_builder.graph,
            &self.rstar_sample,
            ct_builder.settings.kinematics.externals.get_helicities(),
            &ct_builder.settings.additional_params(),
            Some(&threshold_params),
            None,
            None,
        );

        let results = ct_builder.evaluators[esurface_id]
            .evaluator_stack
            .borrow_mut()
            .evaluate(
                params,
                orientations,
                ct_builder.settings,
                evaluation_metadata,
                record_primary_timing,
            )
            .expect("Amplitude counterterm evaluator stack failed")
            .into_iter()
            .map(DualOrNot::unwrap_real)
            .collect_vec();

        let final_result = if results.len() >= 2 {
            debug!(
                "results\nlocal ct:      {:+16e}\nintegrated ct: {:+16e}\nprefactor:     {:+16e}",
                results[0], results[1], prefactor
            );
            (&results[0] + &results[1]) * prefactor
        } else {
            let mut total_ct = Complex::new_re(self.rstar_sample.zero());
            for value in results.iter() {
                total_ct += value;
            }
            total_ct * prefactor
        };

        debug!(
            ct_eval = format!("{:+16e}", final_result),
            "esurface {}", esurface_id.0
        );

        final_result
    }

    fn evaluate_multichanneling_prefactor(
        &self,
        momentum_sample: &MomentumSample<T>,
        model_params: Vec<Complex<F<T>>>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Complex<F<T>> {
        let overlap_builder = self.rstar_solution.esurface_ct_builder.overlap_builder;
        let overlap = overlap_builder.counterterm_builder.overlap_structure;

        if overlap.overlap_groups.len() < 2 {
            return Complex::new_re(momentum_sample.one());
        }

        let params = momentum_sample
            .loop_moms()
            .iter()
            .flat_map(|x| x.into_iter().cloned().map(Complex::new_re))
            .chain(
                momentum_sample
                    .external_moms()
                    .iter()
                    .flat_map(|x| x.into_iter().cloned().map(Complex::new_re)),
            )
            .chain(model_params)
            .collect::<Vec<_>>();

        let evaluator = overlap_builder
            .overlap_group
            .prefactor_evaluator
            .as_ref()
            .unwrap();

        evaluate_evaluator_single(
            &mut evaluator.borrow_mut(),
            &params,
            evaluation_metadata,
            record_primary_timing,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::AmplitudeCountertermAtom;
    use symbolica::{atom::Atom, symbol};

    #[test]
    fn empty_amplitude_counterterm_atom_is_not_generated() {
        let atom = AmplitudeCountertermAtom {
            parametric_local: Atom::new(),
            parametric_integrated: Atom::new(),
        };

        assert!(!atom.is_generated());
    }

    #[test]
    fn non_empty_amplitude_counterterm_atom_is_generated() {
        let atom = AmplitudeCountertermAtom {
            parametric_local: Atom::var(symbol!("x")),
            parametric_integrated: Atom::new(),
        };

        assert!(atom.is_generated());
    }
}
