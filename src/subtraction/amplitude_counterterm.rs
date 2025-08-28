use std::path::Path;

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::involution::EdgeVec;
use log::debug;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::Atom,
    domains::float::{NumericalFloatLike, Real},
    evaluate::OptimizationSettings,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        esurface::{
            self, Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId, GroupEsurfaceId,
        },
        expression::AmplitudeOrientationID,
    },
    gammaloop_integrand::{
        evaluators::SingleOrAllOrientations, param_builder, GenericEvaluator,
        GenericEvaluatorFloat, ParamBuilder, ThresholdParams,
    },
    graph::{FeynmanGraph, Graph, GraphGroupPosition, LoopMomentumBasis},
    momentum::Rotation,
    momentum_sample::{LoopMomenta, MomentumSample},
    settings::{
        runtime::{
            IntegratedCounterTermRange, IntegratedCounterTermSettings, UVLocalisationSettings,
        },
        GlobalSettings, RuntimeSettings,
    },
    subtraction::overlap::{OverlapGroup, OverlapStructure},
    utils::{
        self,
        newton_solver::{newton_iteration_and_derivative, NewtonIterationResult},
        FloatLike, F,
    },
    GammaLoopContext,
};

const MAX_ITERATIONS: usize = 40;
const TOLERANCE: f64 = 1.0;

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermData {
    pub overlap: OverlapStructure,
    pub evaluators: TiVec<EsurfaceID, AmplitudeCountertermEvaluator>,
    pub esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
    pub own_group_position: GraphGroupPosition,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermAtom {
    pub parametric_local: Atom,
    pub parametric_integrated: Atom,
    pub concrete_local: Option<Atom>,
    pub concrete_integrated: Option<Atom>,
}

impl AmplitudeCountertermAtom {
    pub(crate) fn to_evaluator(
        &self,
        param_builder: &ParamBuilder,
        optimiation_settings: OptimizationSettings,
    ) -> AmplitudeCountertermEvaluator {
        let parametric = GenericEvaluator::new_from_builder(
            [
                self.parametric_local.clone(),
                self.parametric_integrated.clone(),
            ],
            param_builder,
            optimiation_settings,
        )
        .unwrap();

        AmplitudeCountertermEvaluator {
            parametric,
            concrete: None,
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermEvaluator {
    pub parametric: GenericEvaluator,
    pub concrete: Option<GenericEvaluator>,
}

impl AmplitudeCountertermData {
    pub fn new_empty(own_group_position: GraphGroupPosition) -> Self {
        Self {
            overlap: OverlapStructure::new_empty(),
            evaluators: TiVec::new(),
            esurface_map: TiVec::new(),
            own_group_position,
        }
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
    ) {
        for (i, e) in self.evaluators.iter_mut_enumerated() {
            e.parametric.compile(
                &path.as_ref().join(format!("esurface_{}", i.0)),
                format!("esurface_{}", i.0),
                &path.as_ref().join(format!("esurface_{}", i.0)),
                settings.generation.gammaloop_compile_options.inline_asm(),
            );
        }
    }

    pub fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
        param_builder: &mut ParamBuilder<f64>,
        orientation: SingleOrAllOrientations<'_, AmplitudeOrientationID>,
    ) -> Complex<F<T>> {
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

        for group in self.overlap.overlap_groups.iter() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            for existing_esurface_id in group.existing_esurfaces.iter() {
                let single_result = overlap_builder
                    .new_esurface_builder(*existing_esurface_id)
                    .map(|esurface_builder| {
                        esurface_builder
                            .solve_rstar()
                            .rstar_samples()
                            .evaluate(param_builder, orientation)
                    })
                    .unwrap_or_else(|| Complex::new_re(momentum_sample.zero()));

                debug!(
                    "Param Builder for {}:\n{}",
                    existing_esurface_id, param_builder
                );

                result += single_result;
            }
        }

        result
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
    fn new(
        graph: &'a Graph,
        rotation_for_overlap: &'a Rotation,
        settings: &'a RuntimeSettings,
        esurface_collection: &'a EsurfaceCollection,
        sample: &'a MomentumSample<T>,
        overlap_structure: &'a OverlapStructure,
        evaluators: &'a TiVec<EsurfaceID, AmplitudeCountertermEvaluator>,
        own_group_position: GraphGroupPosition,
        esurface_map: &'a TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
    ) -> Self {
        let real_mass_vector = graph.underlying.get_real_mass_vector();
        let e_cm = F::from_ff64(settings.kinematics.e_cm);

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
            unrotated_center,
            unit_shifted_momenta,
            radius,
        }
    }
}

struct OverlapBuilder<'a, T: FloatLike> {
    counterterm_builder: &'a CounterTermBuilder<'a, T>,
    overlap_group: &'a OverlapGroup,
    rotated_center: LoopMomenta<F<T>>,
    unrotated_center: LoopMomenta<F<T>>,
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
            existing_esurface_id,
            esurface: &self.counterterm_builder.esurface_collection[esurface_id],
            esurface_id,
        })
    }
}

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    existing_esurface_id: ExistingEsurfaceId,
    esurface: &'a Esurface,
    esurface_id: EsurfaceID,
}

impl<'a, T: FloatLike> EsurfaceCTBuilder<'a, T> {
    fn solve_rstar(self) -> RstarSolution<'a, T> {
        let settings = self.overlap_builder.counterterm_builder.settings;
        let (radius_guess, _) = self.esurface.get_radius_guess(
            &self.overlap_builder.unit_shifted_momenta,
            &self
                .overlap_builder
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
                &self
                    .overlap_builder
                    .counterterm_builder
                    .graph
                    .loop_momentum_basis,
                &self.overlap_builder.counterterm_builder.real_mass_vector,
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

        // some gymnastics to get the sample right
        // do not touch!
        let mut sample_to_modify = self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .sample
            .clone();

        let new_sample = match self
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .sample
            .rotated_sample
        {
            None => {
                sample_to_modify.sample.loop_moms = rstar_loop_momenta;
                sample_to_modify
            }
            Some(_) => {
                sample_to_modify.rotated_sample.as_mut().unwrap().loop_moms = rstar_loop_momenta;

                let new_unrotated_momenta = &(&self
                    .esurface_ct_builder
                    .overlap_builder
                    .counterterm_builder
                    .sample
                    .sample
                    .loop_moms
                    - &self.esurface_ct_builder.overlap_builder.unrotated_center)
                    .rescale(&self.esurface_ct_builder.overlap_builder.radius.inv(), None)
                    .rescale(&self.solution.solution, None)
                    + &self.esurface_ct_builder.overlap_builder.unrotated_center;

                sample_to_modify.sample.loop_moms = new_unrotated_momenta;
                sample_to_modify
            }
        };

        RstarSample {
            rstar_solution: self,
            rstar_sample: new_sample,
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
        orientations: SingleOrAllOrientations<'a, AmplitudeOrientationID>,
    ) -> Complex<F<T>> {
        let esurface_ct_builder = &self.rstar_solution.esurface_ct_builder;
        let ct_builder = esurface_ct_builder.overlap_builder.counterterm_builder;

        let existing_esurface_id = self.rstar_solution.esurface_ct_builder.existing_esurface_id;

        let esurface_id = esurface_ct_builder.esurface_id;

        let esurfaces = ct_builder.esurface_collection;
        let masses = &ct_builder.real_mass_vector;
        let lmb = &ct_builder.graph.loop_momentum_basis;

        let model_params = todo!();

        let prefactor = self.evaluate_multichanneling_prefactor(&self.rstar_sample, model_params);

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

        let uv_damp_plus = evaluate_uv_damper(&radius, &radius_star, &e_cm, &settings);
        let uv_damp_minus = evaluate_uv_damper(&-&radius, &radius_star, &e_cm, &settings);

        debug!("uv_damp_plus: {:?}", uv_damp_plus);
        debug!("uv_damp_minus: {:?}", uv_damp_minus);

        let h_function = evaluate_integrated_ct_normalisation(
            &radius,
            &radius_star,
            &e_cm,
            &integrated_settings,
        );

        let threshold_params = ThresholdParams {
            radius,
            radius_star,
            esurface_derivative: self.rstar_solution.solution.derivative_at_solution,

            uv_damp_plus,
            uv_damp_minus,
            h_function,
        };

        let evaluator = &ct_builder.evaluators[esurface_id].parametric;

        let mut local_ct = Complex::new_re(F::from_f64(0.0));
        let mut integrated_ct = Complex::new_re(F::from_f64(0.0));

        for (i, orientation) in orientations.iter() {
            param_builder.orientation_value(orientation);

            let params = T::get_parameters(
                param_builder,
                ct_builder.graph,
                &self.rstar_sample,
                ct_builder.settings.kinematics.externals.get_helicities(),
                Some(&threshold_params),
            );
            let result = <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params);
            local_ct += &result[0];
            integrated_ct += &result[1];
        }

        debug!(
            "Evaluating ct for esurface {}\nwith params:\n {}",
            esurface_id.0,
            param_builder.clone()
        );

        debug!(
            "results\nlocal ct:      {:+16e}\nintegrated ct: {:+16e}\nprefactor:     {:+16e}",
            local_ct, integrated_ct, prefactor
        );

        let final_result = (local_ct + integrated_ct) * prefactor;

        debug!(
            "sum of local and integrated ct (with prefactor): {}",
            final_result
        );

        final_result
    }

    fn evaluate_multichanneling_prefactor(
        &self,
        momentum_sample: &MomentumSample<T>,
        model_params: &[Complex<F<T>>],
    ) -> F<T> {
        let overlap_builder = self.rstar_solution.esurface_ct_builder.overlap_builder;
        let overlap = overlap_builder.counterterm_builder.overlap_structure;

        if overlap.overlap_groups.len() < 2 {
            return momentum_sample.one();
        }

        let evaluator = overlap_builder
            .overlap_group
            .prefactor_evaluator
            .as_ref()
            .unwrap();

        todo!()
    }
}

fn evaluate_uv_damper<T: FloatLike>(
    radius: &F<T>,
    radius_star: &F<T>,
    e_cm: &F<T>,
    settings: &UVLocalisationSettings,
) -> F<T> {
    let normalizing_scale = match settings.dynamic_width {
        true => radius_star,
        false => e_cm,
    };

    let delta_r = radius - radius_star;

    if delta_r.abs() > F::from_f64(settings.sliver_width) * normalizing_scale {
        return radius.zero();
    }

    let delta_r_sq = &delta_r * &delta_r;
    let width = F::from_f64(settings.gaussian_width) * normalizing_scale;
    let width_sq = &width * &width;

    (-delta_r_sq / width_sq).exp()
}

fn evaluate_integrated_ct_normalisation<T: FloatLike>(
    radius: &F<T>,
    radius_star: &F<T>,
    _e_cm: &F<T>,
    settings: &IntegratedCounterTermSettings,
) -> F<T> {
    match &settings.range {
        IntegratedCounterTermRange::Infinite {
            h_function_settings,
        } => {
            let h = utils::h(&(radius_star / radius), None, None, h_function_settings);
            h * (radius_star).inv()
        }
        IntegratedCounterTermRange::Compact => {
            todo!();
        }
    }
}
