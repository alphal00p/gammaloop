use bincode_trait_derive::{Decode, Encode};
use dot_parser::canonical::Edge;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use spenso::algebra::complex::Complex;
use symbolica::domains::float::{NumericalFloatLike, Real};
use typed_index_collections::TiVec;
use vakint::Momentum;

use crate::{
    cff::esurface::{self, Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId},
    momentum::{ExternalMomenta, Rotation},
    momentum_sample::{self, ExternalFourMomenta, LoopMomenta, MomentumSample},
    new_gammaloop_integrand::{
        GenericEvaluator, GenericEvaluatorFloat, ParamBuilder, ThresholdParams,
    },
    new_graph::{FeynmanGraph, Graph, LoopMomentumBasis},
    numerator,
    subtraction::overlap::{OverlapGroup, OverlapStructure},
    utils::{
        newton_solver::{newton_iteration_and_derivative, NewtonIterationResult},
        FloatLike, F,
    },
    GammaLoopContext, RuntimeSettings, UVLocalisationSettings,
};

const MAX_ITERATIONS: usize = 40;
const TOLERANCE: f64 = 1.0;

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeCountertermData {
    pub overlap: OverlapStructure,
    pub evaluators: TiVec<EsurfaceID, GenericEvaluator>,
    pub param_builder: ParamBuilder<f64>,
}

impl AmplitudeCountertermData {
    pub fn new_empty() -> Self {
        Self {
            overlap: OverlapStructure::new_empty(),
            evaluators: TiVec::new(),
            param_builder: ParamBuilder::new(),
        }
    }

    pub fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &RuntimeSettings,
    ) -> Complex<F<T>> {
        let counter_term_builder = CounterTermBuilder::new(
            graph,
            rotation,
            settings,
            esurfaces,
            momentum_sample,
            &self.overlap,
            &self.evaluators,
        );

        let mut result = Complex::new_re(momentum_sample.zero());

        for group in self.overlap.overlap_groups.iter() {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);

            for existing_esurface_id in group.existing_esurfaces.iter() {
                let single_result = overlap_builder
                    .new_esurface_builder(*existing_esurface_id)
                    .solve_rstar()
                    .rstar_samples()
                    .evaluate(&mut self.param_builder);

                result += single_result;
            }
        }

        result
    }
}

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    evaluators: &'a TiVec<EsurfaceID, GenericEvaluator>,
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
        evaluators: &'a TiVec<EsurfaceID, GenericEvaluator>,
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
    ) -> EsurfaceCTBuilder<'a, T> {
        let esurface_id = self
            .counterterm_builder
            .overlap_structure
            .existing_esurfaces[existing_esurface_id];

        EsurfaceCTBuilder {
            overlap_builder: self,
            existing_esurface_id,
            esurface: &self.counterterm_builder.esurface_collection[esurface_id],
        }
    }
}

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    existing_esurface_id: ExistingEsurfaceId,
    esurface: &'a Esurface,
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
    fn evaluate(self, param_builder: &mut ParamBuilder<f64>) -> Complex<F<T>> {
        let ct_builder = self
            .rstar_solution
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder;

        let existing_esurface_id = self.rstar_solution.esurface_ct_builder.existing_esurface_id;

        let esurface_id = ct_builder.overlap_structure.existing_esurfaces[existing_esurface_id];

        let esurfaces = ct_builder.esurface_collection;
        let masses = &ct_builder.real_mass_vector;
        let lmb = &ct_builder.graph.loop_momentum_basis;

        let prefactor =
            self.evaluate_multichanneling_prefactor(&self.rstar_sample, lmb, masses, esurfaces);

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

        let uv_damp_plus = evaluate_uv_damper(&radius, &radius_star, &e_cm, &settings);
        let uv_damp_minus = evaluate_uv_damper(&-&radius, &radius_star, &e_cm, &settings);

        let threshold_params = ThresholdParams {
            radius,
            radius_star,
            esurface_derivative: self.rstar_solution.solution.derivative_at_solution,

            uv_damp_plus,
            uv_damp_minus,
        };

        let params = T::get_parameters(
            param_builder,
            ct_builder.graph,
            &self.rstar_sample,
            ct_builder.settings.kinematics.externals.get_helicities(),
            Some(&threshold_params),
        );

        let evaluator = &ct_builder.evaluators[esurface_id];
        <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params) * &prefactor
    }

    fn evaluate_multichanneling_prefactor(
        &self,
        momentum_sample: &MomentumSample<T>,
        lmb: &LoopMomentumBasis,
        masses: &EdgeVec<F<T>>,
        esurfaces: &EsurfaceCollection,
    ) -> F<T> {
        let overlap_builder = self.rstar_solution.esurface_ct_builder.overlap_builder;
        let overlap = overlap_builder.counterterm_builder.overlap_structure;

        if overlap.overlap_groups.len() < 2 {
            return momentum_sample.one();
        }

        let denominator = overlap
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .complement
                    .iter()
                    .map(|existing_esurface_id| {
                        let esurface_id = overlap.existing_esurfaces[*existing_esurface_id];
                        let esurface = &esurfaces[esurface_id];
                        let esurface_val = esurface.compute_from_momenta(
                            lmb,
                            masses,
                            momentum_sample.loop_moms(),
                            momentum_sample.external_moms(),
                        );

                        &esurface_val * &esurface_val
                    })
                    .reduce(|product, x| &product * &x)
                    .unwrap_or_else(|| momentum_sample.one())
            })
            .reduce(|sum, x| &sum + &x)
            .unwrap_or_else(|| momentum_sample.zero());

        let numerator = overlap_builder
            .overlap_group
            .complement
            .iter()
            .map(|existing_esurface_id| {
                let esurface_id = overlap.existing_esurfaces[*existing_esurface_id];
                let esurface = &esurfaces[esurface_id];
                let esurface_val = esurface.compute_from_momenta(
                    lmb,
                    masses,
                    momentum_sample.loop_moms(),
                    momentum_sample.external_moms(),
                );

                &esurface_val * &esurface_val
            })
            .reduce(|product, x| &product * &x)
            .unwrap_or_else(|| momentum_sample.one());

        numerator / denominator
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
