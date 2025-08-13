use dot_parser::canonical::Edge;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use spenso::algebra::complex::Complex;
use symbolica::domains::float::{NumericalFloatLike, Real};
use typed_index_collections::TiVec;

use crate::{
    cff::esurface::{self, Esurface, EsurfaceCollection, EsurfaceID, ExistingEsurfaceId},
    momentum::{ExternalMomenta, Rotation},
    momentum_sample::{self, ExternalFourMomenta, LoopMomenta, MomentumSample},
    new_gammaloop_integrand::GenericEvaluator,
    new_graph::{FeynmanGraph, Graph, LoopMomentumBasis},
    subtraction::overlap::{OverlapGroup, OverlapStructure},
    utils::{newton_solver::NewtonIterationResult, FloatLike, F},
    Settings,
};

pub struct AmplitudeCountertermData {
    pub overlap: OverlapStructure,
    pub evaluators: TiVec<EsurfaceID, GenericEvaluator>,
}

impl AmplitudeCountertermData {
    pub fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        esurfaces: &EsurfaceCollection,
        rotation: &Rotation,
        settings: &Settings,
    ) -> Complex<T> {
        let counter_term_builder = CounterTermBuilder::new(
            graph,
            rotation,
            settings,
            esurfaces,
            momentum_sample,
            &self.overlap,
        );

        self.overlap.overlap_groups.iter().map(|group| {
            let overlap_builder = counter_term_builder.new_overlap_builder(group);
        });
        todo!()
    }

    pub fn evaluate_multichanneling_denominator<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        lmb: &LoopMomentumBasis,
        masses: &EdgeVec<F<T>>,
        esurfaces: &EsurfaceCollection,
    ) -> F<T> {
        self.overlap
            .overlap_groups
            .iter()
            .map(|group| {
                group
                    .complement
                    .iter()
                    .map(|existing_esurface_id| {
                        let esurface_id = self.overlap.existing_esurfaces[*existing_esurface_id];
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
            .unwrap_or_else(|| momentum_sample.zero())
    }
}

struct CounterTermBuilder<'a, T: FloatLike> {
    overlap_structure: &'a OverlapStructure,
    real_mass_vector: EdgeVec<F<T>>,
    e_cm: F<T>,
    graph: &'a Graph,
    rotation_for_overlap: &'a Rotation,
    settings: &'a Settings,
    esurface_collection: &'a EsurfaceCollection,
    sample: &'a MomentumSample<T>,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    fn new(
        graph: &'a Graph,
        rotation_for_overlap: &'a Rotation,
        settings: &'a Settings,
        esurface_collection: &'a EsurfaceCollection,
        sample: &'a MomentumSample<T>,
        overlap_structure: &'a OverlapStructure,
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

struct RstarSolution<'a, T: FloatLike> {
    esurface_ct_builder: EsurfaceCTBuilder<'a, T>,
    solution: NewtonIterationResult<T>,
}
