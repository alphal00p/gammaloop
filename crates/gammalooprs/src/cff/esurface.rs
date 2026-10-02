use std::collections::VecDeque;
use std::fmt::Display;
use std::mem;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Flow, HedgePair};
use linnet::half_edge::subgraph::{Inclusion, ModifySubSet, OrientedCut, SuBitGraph};
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};

use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::dual::HyperDual;
use symbolica::domains::float::{FloatLike as SymFloatLike, Real};
use symbolica::id::Replacement;
use symbolica::{function, parse};
use thiserror::Error;
use tracing::debug;
use typed_index_collections::TiVec;

use crate::cff::VertexSet;

use crate::cff::expression::{
    CFFExpression, OrientationID, RaisedEsurfaceDataView, RaisedEsurfaceGroupView,
};
pub use crate::cff::surface::EsurfaceID;
use crate::graph::{Graph, GraphGroupPosition, LmbIndex, LoopMomentumBasis};
use crate::{GammaLoopContext, define_index};

use crate::graph::FinalizedCut;
use crate::integrands::process::GenericEvaluator;
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta, SubspaceData,
};
use crate::utils::hyperdual_utils::new_constant;
use crate::utils::{
    DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, ESURFACE_SHIFT_THRESHOLD, F, FloatLike, GS,
    compute_loop_part, compute_loop_part_subspace, compute_shift_part, compute_shift_part_subspace,
    compute_t_part_of_shift_part, cut_energy, external_energy_atom_from_index, ose_atom_from_index,
};
use crate::uv::uv_graph::UVE;
use color_eyre::Result;

use super::generation::ShiftRewrite;

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone, bincode::Encode, bincode::Decode)]
pub struct Esurface {
    pub energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
    //#[bincode(with_serde)]
    //pub subspace_graph: InternalSubGraph,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NonExistingEsurfaceReason {
    NoExternalShift,
    ShiftNotNegative,
    NoRadialDependence,
    NoRealZero,
}

#[derive(Debug, Clone)]
pub(crate) enum EsurfaceExistence<T: FloatLike> {
    NonExisting {
        normalized_margin: Option<F<T>>,
        reason: NonExistingEsurfaceReason,
    },
    Pinched {
        normalized_margin: F<T>,
    },
    Existing {
        normalized_margin: F<T>,
    },
}

/// Pointwise existence classification exposed without the internal diagnostic margin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EsurfaceExistenceStatus {
    NonExisting,
    Pinched,
    Existing,
}

impl<T: FloatLike> EsurfaceExistence<T> {
    fn status(&self) -> EsurfaceExistenceStatus {
        match self {
            Self::NonExisting { .. } => EsurfaceExistenceStatus::NonExisting,
            Self::Pinched { .. } => EsurfaceExistenceStatus::Pinched,
            Self::Existing { .. } => EsurfaceExistenceStatus::Existing,
        }
    }

    pub(crate) fn is_existing(&self) -> bool {
        matches!(self, Self::Existing { .. })
    }

    pub(crate) fn normalized_margin(&self) -> Option<&F<T>> {
        match self {
            Self::NonExisting {
                normalized_margin, ..
            } => normalized_margin.as_ref(),
            Self::Pinched { normalized_margin } | Self::Existing { normalized_margin } => {
                Some(normalized_margin)
            }
        }
    }

    pub(crate) fn label(&self) -> &'static str {
        match self {
            Self::NonExisting { .. } => "non_existing",
            Self::Pinched { .. } => "pinched",
            Self::Existing { .. } => "existing",
        }
    }

    pub(crate) fn non_existing_reason(&self) -> Option<NonExistingEsurfaceReason> {
        match self {
            Self::NonExisting { reason, .. } => Some(*reason),
            Self::Pinched { .. } | Self::Existing { .. } => None,
        }
    }
}

pub(crate) fn esurface_value_is_strictly_inside<T: FloatLike>(value: &F<T>, e_cm: &F<T>) -> bool {
    let interior_tolerance = value.epsilon() * value.from_i64(8) * e_cm;
    !value.is_nan() && !value.is_infinite() && value < &(-interior_tolerance)
}

impl PartialEq for Esurface {
    fn eq(&self, other: &Self) -> bool {
        self.energies == other.energies && self.external_shift == other.external_shift
    }
}

impl Eq for Esurface {}

impl Esurface {
    pub(crate) fn has_radial_dependence_in_subspace(
        &self,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
    ) -> bool {
        let lmb = subspace.get_lmb(all_lmbs);
        subspace.contains(&self.energies, graph).any(|index| {
            subspace
                .project_loop_signature(&lmb.edge_signatures[index].internal)
                .any(|sign| sign.is_sign())
        })
    }

    pub(crate) fn external_shift_is_strictly_negative_for_positive_energies(
        &self,
        incoming_edges: &[EdgeIndex],
        outgoing_edges: &[EdgeIndex],
    ) -> bool {
        if incoming_edges.is_empty()
            || outgoing_edges.is_empty()
            || incoming_edges
                .iter()
                .any(|edge| outgoing_edges.contains(edge))
            || self
                .external_shift
                .iter()
                .any(|(edge, _)| !incoming_edges.contains(edge) && !outgoing_edges.contains(edge))
        {
            return false;
        }

        let coefficient = |edge: &EdgeIndex| {
            self.external_shift
                .iter()
                .filter(|(shift_edge, _)| shift_edge == edge)
                .map(|(_, coefficient)| i128::from(*coefficient))
                .sum::<i128>()
        };

        // Energy conservation makes external-energy coefficient vectors `c` and
        // `c + lambda * sigma` equivalent, where sigma is +1 for incoming and -1
        // for outgoing momenta. The shift is strictly negative for all positive
        // external energies if one representative is component-wise non-positive
        // and not identically zero.
        let lambda_lower_bound = outgoing_edges
            .iter()
            .map(&coefficient)
            .max()
            .expect("outgoing external edges were checked to be non-empty");
        let lambda_upper_bound = incoming_edges
            .iter()
            .map(|edge| -coefficient(edge))
            .min()
            .expect("incoming external edges were checked to be non-empty");

        if lambda_lower_bound > lambda_upper_bound {
            return false;
        }

        let lambda = lambda_lower_bound;
        let adjusted_coefficients = incoming_edges
            .iter()
            .map(|edge| coefficient(edge) + lambda)
            .chain(outgoing_edges.iter().map(|edge| coefficient(edge) - lambda));
        let mut has_strictly_negative_coefficient = false;
        for adjusted_coefficient in adjusted_coefficients {
            if adjusted_coefficient > 0 {
                return false;
            }
            has_strictly_negative_coefficient |= adjusted_coefficient < 0;
        }

        has_strictly_negative_coefficient
    }

    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let symbolic_energies = self
            .energies
            .iter()
            .map(|i| {
                if cut_edges.contains(i) {
                    cut_energy(*i)
                } else {
                    ose_atom_from_index(*i)
                }
            })
            .collect_vec();

        let symbolic_shift = self
            .external_shift
            .iter()
            .fold(Atom::new(), |sum, (i, sign)| {
                external_energy_atom_from_index(*i) * &Atom::num(*sign) + &sum
            });

        let builder_atom = Atom::new();
        let energy_sum = symbolic_energies
            .iter()
            .fold(builder_atom, |acc, energy| acc + energy);

        energy_sum + &symbolic_shift
    }

    #[inline]
    pub(crate) fn compute_from_dual_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        dual_loop_moms: &LoopMomenta<HyperDual<F<T>>>,
        dual_external_moms: &ExternalFourMomenta<HyperDual<F<T>>>,
    ) -> HyperDual<F<T>> {
        let spatial_part_of_externals = dual_external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                let momentum = signature
                    .try_compute_momentum(&dual_loop_moms.0, &spatial_part_of_externals.raw)
                    .unwrap_or_else(|| unreachable!());
                let mass = &real_mass_vector[*index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| dual_loop_moms[LoopIndex(0)].px.zero());

        let shift_part = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                let sign = energy_sum.values[0].from_i64(*sign);
                new_constant(&energy_sum, &sign)
                    * external_signature
                        .try_apply(&dual_external_moms.raw)
                        .map(|mom| mom.temporal.value)
                        .unwrap_or_else(|| energy_sum.zero())
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| energy_sum.zero());

        energy_sum + shift_part
    }

    /// Compute the value of the esurface from the momenta, needed to check if an arbitrary point
    /// is inside the esurface
    #[inline]
    pub(crate) fn compute_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) -> F<T> {
        let spatial_part_of_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                let momentum = signature.compute_momentum(loop_moms, &spatial_part_of_externals);
                let mass = &real_mass_vector[*index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| loop_moms[LoopIndex(0)].px.zero());

        let shift_part = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                energy_sum.from_i64(*sign)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| energy_sum.zero());

        energy_sum + shift_part
    }

    fn classify_invariant_margin<T: FloatLike>(
        shift_part: &F<T>,
        invariant_margin: F<T>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistence<T> {
        // /!\ In alphaLoop this boundary needed special care so that `mul_unit` could not
        // create a discontinuous non-pinched-to-pinched transition for a massless 2->2
        // E-surface sandwich. Such a surface is `Existing` at non-collinear kinematics, where
        // its invariant margin is positive, and `Pinched` only at the collinear null boundary.
        // Here the normalized tolerance band around that boundary is classified as `Pinched`,
        // so perturbations within the band cannot toggle threshold subtraction. GammaLoop's
        // statuses are exclusive and only `Existing` surfaces are subtraction targets.
        // Keep the decision in the original, dimensionful scale. Besides avoiding an
        // unnecessary division, this preserves the previous existence boundary exactly.
        // Deserialization rejects invalid configured tolerances. Keep this defensive fallback for
        // settings mutated programmatically (for example through a language binding), so an
        // invalid sign can never invert the existence test.
        let normalized_margin_tolerance =
            if normalized_margin_tolerance.is_nan() || normalized_margin_tolerance.is_infinite() {
                F::from_f64(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD)
            } else {
                normalized_margin_tolerance.abs()
            };
        let invariant_tolerance = &normalized_margin_tolerance * e_cm * e_cm;
        let normalized_margin = &invariant_margin / (e_cm * e_cm);
        let shift_tolerance = F::from_f64(ESURFACE_SHIFT_THRESHOLD) * e_cm;

        if invariant_margin.abs() <= invariant_tolerance && shift_part <= &shift_tolerance {
            return EsurfaceExistence::Pinched { normalized_margin };
        }

        if shift_part >= &(-shift_tolerance) {
            return EsurfaceExistence::NonExisting {
                normalized_margin: Some(normalized_margin),
                reason: NonExistingEsurfaceReason::ShiftNotNegative,
            };
        }

        if invariant_margin > invariant_tolerance {
            EsurfaceExistence::Existing { normalized_margin }
        } else if invariant_margin >= -invariant_tolerance {
            EsurfaceExistence::Pinched { normalized_margin }
        } else {
            EsurfaceExistence::NonExisting {
                normalized_margin: Some(normalized_margin),
                reason: NonExistingEsurfaceReason::NoRealZero,
            }
        }
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    /// Classify a surface in an active loop-momentum subspace. Pinched and
    /// non-existing surfaces remain distinct and are never subtraction targets.
    ///
    /// # Errors
    /// Returns a structural error when fixed energies require a boundary that
    /// cannot be certified from this surface's topology and exact routing.
    pub(crate) fn classify_existence_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> Result<EsurfaceExistence<T>, EsurfaceBoundaryError> {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            debug!("esurface has no external shift, cannot exist");
            return Ok(EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoExternalShift,
            });
        }

        let shift_part = self.compute_shift_part_from_momenta_in_subspace(
            loop_moms,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            real_mass_vector,
        );

        let subspace_energy_indices = subspace.contains(&self.energies, graph).collect_vec();

        if !self.has_radial_dependence_in_subspace(subspace, all_lmbs, graph) {
            debug!(
                "esurface has no radial energy in this subspace, cannot bound a threshold region"
            );
            return Ok(EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoRadialDependence,
            });
        }

        let lmb = subspace.get_lmb(all_lmbs);
        let mass_sum: F<T> = subspace_energy_indices
            .iter()
            .map(|&index| &real_mass_vector[index])
            .fold(e_cm.zero(), |acc, x| acc + x);

        let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

        let graph_vector = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                compute_shift_part(external_signature, external_moms).spatial * e_cm.from_i64(*sign)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

        let fixed_energies = subspace
            .does_not_contain(&self.energies, graph)
            .collect_vec();
        let other_part = if fixed_energies.is_empty() {
            zero_vector.clone()
        } else {
            let boundary = self.oriented_boundary(graph)?;
            let spatial_externals: ExternalThreeMomenta<F<T>> = external_moms
                .iter()
                .map(|mom| mom.spatial.clone())
                .collect();
            fixed_energies
                .into_iter()
                .map(|index| {
                    let signature = &lmb.edge_signatures[index];
                    let momentum = signature.compute_momentum(loop_moms, &spatial_externals);
                    let (flow, _, _) = boundary
                        .iter_edges_flow(&graph.underlying)
                        .find(|(_, edge, _)| *edge == index)
                        .expect("every surface energy belongs to its certified boundary");
                    match flow {
                        Flow::Source => momentum,
                        Flow::Sink => -momentum,
                    }
                })
                .reduce(|acc, x| acc + x)
                .expect("the fixed-energy list is nonempty")
        };

        let shift_vector_sq = (&graph_vector + &other_part).norm_squared();
        let invariant_margin = &shift_part * &shift_part - &shift_vector_sq - &mass_sum * &mass_sum;
        let classification = Self::classify_invariant_margin(
            &shift_part,
            invariant_margin,
            e_cm,
            normalized_margin_tolerance,
        );

        if !classification.is_existing() {
            debug!(
                "subspace esurface classified as {}: shift_part^2: {}, shift_vector_sq: {}, mass_sum^2: {}, normalized_margin: {:?}",
                classification.label(),
                &shift_part * &shift_part,
                shift_vector_sq,
                &mass_sum * &mass_sum,
                classification.normalized_margin(),
            );
        }

        Ok(classification)
    }

    #[inline]
    /// Classify a full-space surface. Pinched and non-existing surfaces remain
    /// distinct and are never subtraction targets.
    pub(crate) fn classify_existence<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistence<T> {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            return EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoExternalShift,
            };
        }

        let shift_part = self.compute_shift_part_from_momenta(external_moms, lmb);
        let mass_sum: F<T> = self
            .energies
            .iter()
            .map(|index| &real_mass_vector[*index])
            .fold(F::from_f64(0.0), |acc, x| acc + x);

        let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

        let shift_vector = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                compute_shift_part(external_signature, external_moms).spatial
                    * F::from_f64(*sign as f64)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

        let shift_vector_sq = shift_vector.norm_squared();
        let invariant_margin = &shift_part * &shift_part - shift_vector_sq - &mass_sum * &mass_sum;

        Self::classify_invariant_margin(
            &shift_part,
            invariant_margin,
            e_cm,
            normalized_margin_tolerance,
        )
    }

    /// Classify a full-space surface while keeping normalized margins and rejection reasons
    /// internal to the subtraction implementation.
    pub fn existence_status<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistenceStatus {
        self.classify_existence(
            external_moms,
            lmb,
            real_mass_vector,
            e_cm,
            normalized_margin_tolerance,
        )
        .status()
    }

    /// Only compute the shift part, useful for center finding.
    pub(crate) fn compute_shift_part_from_momenta_in_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
    ) -> F<T> {
        let lmb = subspace.get_lmb(all_lmbs);

        let full_external_shift = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                external_moms[ExternalIndex(0)]
                    .temporal
                    .value
                    .from_i64(*sign)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[ExternalIndex(0)].temporal.value.zero());

        let spatial_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let remaining_shift = subspace
            .does_not_contain(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];
                //panic!("signature: {:?}", signature);
                let momentum = signature.compute_momentum(loop_moms, &spatial_externals);
                let mass = &masses[index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| full_external_shift.zero());

        full_external_shift + remaining_shift
    }

    /// Only compute the shift part, useful for center finding.
    pub(crate) fn compute_shift_part_from_momenta<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> F<T> {
        self.external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                external_moms[ExternalIndex(0)]
                    .temporal
                    .value
                    .from_i64(*sign)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[ExternalIndex(0)].temporal.value.zero())
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute_self_and_r_derivative_subspace<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops_in_subspace: &LoopMomenta<F<T>>,
        center_in_subspace: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        real_mass_vector: &EdgeVec<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = shifted_unit_loops_in_subspace
            .iter_enumerated()
            .map(|(loop_index, shifted_unit_momenta)| {
                if subspace.contains_loop_index(loop_index) {
                    shifted_unit_momenta * radius + &center_in_subspace[loop_index]
                } else {
                    shifted_unit_momenta.clone()
                }
            })
            .collect();

        let shift = self.compute_shift_part_from_momenta_in_subspace(
            shifted_unit_loops_in_subspace,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            real_mass_vector,
        );

        let lmb = subspace.get_lmb(all_lmbs);
        let (derivative, energy_sum) = subspace
            .contains(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];

                let momentum = signature.compute_momentum(&loops, &spatial_part_of_externals);
                let unit_loop_part = compute_loop_part_subspace(
                    &signature.internal,
                    shifted_unit_loops_in_subspace,
                    subspace,
                );

                let energy = (momentum.norm_squared()
                    + &real_mass_vector[index] * &real_mass_vector[index])
                    .sqrt();

                let numerator = momentum * &unit_loop_part;

                (numerator / &energy, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    #[inline]
    pub(crate) fn compute_self_and_r_derivative<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        real_mass_vector: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = shifted_unit_loops
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum * radius + center)
            .collect();

        let shift = self.compute_shift_part_from_momenta(external_moms, lmb);

        let (derivative, energy_sum) = self
            .energies
            .iter()
            .map(|&index| {
                let signature = &lmb.edge_signatures[index];

                let momentum = signature.compute_momentum(&loops, &spatial_part_of_externals);
                let unit_loop_part = compute_loop_part(&signature.internal, shifted_unit_loops);

                let energy = (momentum.norm_squared()
                    + &real_mass_vector[index] * &real_mass_vector[index])
                    .sqrt();

                let numerator = momentum * &unit_loop_part;

                (numerator / &energy, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    // #[inline]
    /// the "loops_unit_in_subspace" means that the loop momenta that are part of the subspace are jointly normalized to unit length
    pub(crate) fn get_radius_guess_subspace<T: FloatLike>(
        &self,
        loops_unit_in_subspace: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
    ) -> (F<T>, F<T>) {
        let const_builder = &loops_unit_in_subspace[LoopIndex(0)].px;

        let esurface_shift = self.compute_shift_part_from_momenta_in_subspace(
            loops_unit_in_subspace,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            masses,
        );

        debug!("shift part for radius guess: {}", esurface_shift);

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        let lmb = subspace.get_lmb(all_lmbs);

        debug!("unit loops in subspace: {:?}", loops_unit_in_subspace);
        let external_3_momenta = external_moms.iter().map(|x| x.spatial.clone()).collect();

        //println!("got to energy loop");
        for energy in subspace.contains(&self.energies, graph) {
            debug!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);

            let unit_loop_part =
                compute_loop_part_subspace(&signature.internal, loops_unit_in_subspace, subspace);
            //println!("computed_loop_part {:?}", unit_loop_part);

            let three_shift = compute_shift_part_subspace(
                &signature.internal,
                &signature.external,
                loops_unit_in_subspace,
                &external_3_momenta,
                subspace,
            );
            //./bprintln!("computed_shift {:?}", shift);

            let norm_unit_loop_part_squared = unit_loop_part.norm_squared();
            let loop_dot_shift = &unit_loop_part * three_shift;

            debug!(
                "norm unit loop part squared: {}",
                norm_unit_loop_part_squared
            );

            radius_guess += loop_dot_shift.abs() / &norm_unit_loop_part_squared;
            debug!("current radius guess: {}", radius_guess);
            denominator += norm_unit_loop_part_squared.sqrt();
            debug!("current denominator: {}", denominator);
        }

        if denominator != const_builder.zero() {
            radius_guess += esurface_shift.abs() / denominator;
        }
        debug!("final radius guess: {}", radius_guess);
        let negative_radius = radius_guess.ref_neg();
        (radius_guess, negative_radius)
    }

    pub(crate) fn get_radius_guess<T: FloatLike>(
        &self,
        unit_loops: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let const_builder = &unit_loops[LoopIndex(0)].px;

        let esurface_shift = self.compute_shift_part_from_momenta(external_moms, lmb);

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        //println!("got to energy loop");
        for &energy in &self.energies {
            //println!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);

            let unit_loop_part = compute_loop_part(&signature.internal, unit_loops);
            //println!("computed_loop_part {:?}", unit_loop_part);

            let three_shift = compute_shift_part(&signature.external, external_moms).spatial;
            //./bprintln!("computed_shift {:?}", shift);

            let norm_unit_loop_part_squared = unit_loop_part.norm_squared();
            let loop_dot_shift = &unit_loop_part * three_shift;

            radius_guess += loop_dot_shift.abs() / &norm_unit_loop_part_squared;
            denominator += norm_unit_loop_part_squared.sqrt();
        }

        radius_guess += esurface_shift.abs() / denominator;
        let negative_radius = radius_guess.ref_neg();
        (radius_guess, negative_radius)
    }

    pub(crate) fn canonicalize_shift(&mut self, shift_rewrite: &ShiftRewrite) {
        if let Some(dep_mom_pos) = self
            .external_shift
            .iter()
            .position(|(index, _)| *index == shift_rewrite.dependent_momentum)
        {
            let (_, dep_mom_sign) = self.external_shift.remove(dep_mom_pos);

            let external_shift = shift_rewrite
                .dependent_momentum_expr
                .iter()
                .map(|(index, sign)| (*index, dep_mom_sign * sign))
                .collect();

            self.external_shift = add_external_shifts(&self.external_shift, &external_shift);
        }
    }

    pub(crate) fn new_from_cut_left<E, V, H>(
        graph: &HedgeGraph<E, V, H>,
        cut: &FinalizedCut,
        initial_state_cut: Option<&OrientedCut>,
    ) -> Self {
        let edges = graph
            .iter_edges_of(&cut.cut)
            .map(|(_, id, _)| id)
            .sorted()
            .collect();

        let external_shift = if let Some(is_cut) = initial_state_cut {
            graph
                .iter_edges_of(is_cut)
                .map(|(_, edge_index, __)| (edge_index, -1))
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect()
        } else {
            graph
                .iter_edges_of(&cut.left)
                .filter_map(|(hedge_pair, edge_index, _)| match hedge_pair {
                    HedgePair::Unpaired { flow, .. } => match flow {
                        Flow::Sink => Some((edge_index, -1)),
                        Flow::Source => Some((edge_index, 1)),
                    },
                    _ => None,
                })
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect()
        };

        let vertex_set = graph
            .iter_nodes_of(&cut.left)
            .map(|(node_id, _, _)| VertexSet::from_usize(node_id.into()))
            .reduce(|acc, v| acc.join(&v))
            .unwrap();

        Self {
            energies: edges,
            external_shift,
            vertex_set,
            //subspace_graph: graph.full_graph(),
        }
    }

    pub(crate) fn lmb_atom(
        &self,
        graph: &Graph,
        model: &crate::model::Model,
        lmb_reps: &[Replacement],
    ) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom(model);
                let emr_symbols = (0..3)
                    .map(|i| function!(GS.emr_mom, usize::from(*index), i + 1))
                    .collect_vec();

                (&emr_symbols[0] * &emr_symbols[0]
                    + &emr_symbols[1] * &emr_symbols[1]
                    + &emr_symbols[2] * &emr_symbols[2]
                    + &mass_symbol * &mass_symbol)
                    .sqrt()
            })
            .chain(self.external_shift.iter().map(|(index, sign)| {
                function!(GS.emr_mom, usize::from(*index), 0) * Atom::num(*sign)
            }))
            .reduce(|sum, atom| sum + atom)
            .unwrap_or_else(Atom::new)
            .replace_multiple(lmb_reps)
            .replace(parse!("ZERO"))
            .with(Atom::new())
            .expand() // ensure canonical form
    }

    // more readable version for debugging, because it doesn't write out components
    pub(crate) fn lmb_atom_simplified(
        &self,
        graph: &Graph,
        model: &crate::model::Model,
        lmb_reps: &[Replacement],
    ) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom(model);
                let emr_symbol = function!(GS.emr_mom, usize::from(*index));

                (&emr_symbol * &emr_symbol + &mass_symbol * &mass_symbol).sqrt()
            })
            .chain(self.external_shift.iter().map(|(index, sign)| {
                function!(GS.emr_mom, usize::from(*index), 0) * Atom::num(*sign)
            }))
            .reduce(|sum, atom| sum + atom)
            .unwrap_or_else(Atom::new)
            .replace_multiple(lmb_reps)
            .replace(parse!("ZERO"))
            .with(Atom::new())
            .expand() // ensure canonical form
    }
}

define_index! {pub struct GroupEsurfaceId;}

pub type EsurfaceCollection = TiVec<EsurfaceID, Esurface>;

pub type EsurfaceCache<T> = TiVec<EsurfaceID, T>;

/// Container for esurfaces that exist at a given point in the phase space
pub type ExistingEsurfaces = TiVec<ExistingEsurfaceId, GroupEsurfaceId>;
pub type ExistingThresholds = TiVec<ExistingEsurfaceId, EsurfaceID>;

pub(crate) fn get_representative<T: Copy>(
    esurface_map: &TiVec<GraphGroupPosition, Option<T>>,
) -> Result<(GraphGroupPosition, T)> {
    for (group_pos, esurface_option) in esurface_map.iter_enumerated() {
        if let Some(esurface_id) = esurface_option {
            return Ok((group_pos, *esurface_id));
        }
    }

    Err(eyre!(
        "No representative esurface found, esurface map corrupted"
    ))
}

/// Index in the list of all existing esurfaces, essentially a pointer to a pointer to an esurface
#[derive(
    Debug,
    From,
    Into,
    Copy,
    Clone,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Encode,
    Decode,
)]
pub struct ExistingEsurfaceId(usize);

impl Display for ExistingEsurfaceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ExistingEsurfaceID({})", self.0)
    }
}

pub type ExternalShift = Vec<(EdgeIndex, i64)>;

/// add two external shifts, eliminates zero signs and sorts
pub(crate) fn add_external_shifts(lhs: &ExternalShift, rhs: &ExternalShift) -> ExternalShift {
    let mut res = lhs.clone();

    for rhs_element in rhs.iter() {
        if let Some(lhs_element) = res
            .iter_mut()
            .find(|lhs_element| rhs_element.0 == lhs_element.0)
        {
            lhs_element.1 += rhs_element.1;
        } else {
            res.push(*rhs_element)
        }
    }

    res.retain(|(_index, sign)| *sign != 0);
    res.sort_by_key(|(index, _)| *index);
    res
}

define_index!(
    pub struct RaisedEsurfaceId;
);

#[derive(Debug, Clone, Default, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct RaisedEsurfaceData {
    pub raised_groups: TiVec<RaisedEsurfaceId, RaisedEsurfaceGroup>,
    pub pass_two_evaluator: Option<Vec<GenericEvaluator>>,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct RaisedEsurfaceGroup {
    pub esurface_ids: Vec<EsurfaceID>,
    pub max_occurence: usize,
}

impl RaisedEsurfaceGroupView for RaisedEsurfaceGroup {
    fn esurface_ids(&self) -> &[EsurfaceID] {
        &self.esurface_ids
    }

    fn max_occurrence(&self) -> usize {
        self.max_occurence
    }
}

impl RaisedEsurfaceDataView for RaisedEsurfaceData {
    fn for_each_raised_group(&self, mut f: impl FnMut(&dyn RaisedEsurfaceGroupView)) {
        for group in &self.raised_groups {
            f(group);
        }
    }
}

impl Graph {
    pub(crate) fn normalize_esurface_with_raised_edge_groups(
        esurface: &Esurface,
        raised_edges: &[Vec<EdgeIndex>],
    ) -> Esurface {
        let mut normalized = esurface.clone();
        for energy in normalized.energies.iter_mut() {
            if let Some(representative) = raised_edges
                .iter()
                .find(|group| group.contains(energy))
                .and_then(|group| group.first())
            {
                *energy = *representative;
            }
        }
        normalized.energies.sort();
        normalized
    }

    pub(crate) fn determine_raised_esurfaces_from_expression(
        &self,
        expr: &CFFExpression<OrientationID>,
    ) -> RaisedEsurfaceData {
        let raised_edges = self.get_raised_edge_groups();
        let normalized_cut_esurfaces = expr
            .surfaces
            .esurface_cache
            .iter()
            .map(|esurface| {
                Self::normalize_esurface_with_raised_edge_groups(esurface, &raised_edges)
            })
            .collect::<TiVec<EsurfaceID, _>>();

        let mut raised_groups = TiVec::<RaisedEsurfaceId, RaisedEsurfaceGroup>::new();

        for (esurface_id, normalized_cut_esurface) in normalized_cut_esurfaces.iter_enumerated() {
            let raised_esurface_group_id = raised_groups.iter_enumerated().find_map(
                |(raised_esurface_group_id, esurface_group)| {
                    if esurface_group
                        .esurface_ids
                        .iter()
                        .all(|esurface_id_in_group| {
                            normalized_cut_esurfaces[*esurface_id_in_group].energies
                                == normalized_cut_esurface.energies
                                && normalized_cut_esurfaces[*esurface_id_in_group].external_shift
                                    == normalized_cut_esurface.external_shift
                        })
                    {
                        Some(raised_esurface_group_id)
                    } else {
                        None
                    }
                },
            );

            if let Some(found_group_id) = raised_esurface_group_id {
                raised_groups[found_group_id].esurface_ids.push(esurface_id);
            } else {
                raised_groups.push(RaisedEsurfaceGroup {
                    esurface_ids: vec![esurface_id],
                    max_occurence: 0,
                });
            }
        }

        let mut result = RaisedEsurfaceData {
            raised_groups,
            pass_two_evaluator: None,
        };

        let mut expression_copy = expr.clone();
        expression_copy.normalize_wrt_all_raisings(&result);

        for cut_group in result.raised_groups.iter_mut() {
            let representative_esurface_id = cut_group.esurface_ids[0];

            let max_occurence_for_this_id = expression_copy
                .orientations
                .iter()
                .map(|orientation_expression| {
                    orientation_expression.max_effective_denominator_value_count_on_branch(
                        &crate::cff::surface::HybridSurfaceID::Esurface(representative_esurface_id),
                    )
                })
                .max()
                .unwrap_or(0);

            cut_group.max_occurence = max_occurence_for_this_id;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;

    use eyre::Report;
    use itertools::Itertools;
    use linnet::half_edge::HedgeGraph;
    use linnet::half_edge::builder::HedgeGraphBuilder;
    use linnet::half_edge::involution::{EdgeIndex, Flow, HedgePair, Orientation};
    use linnet::half_edge::nodestore::NodeStorageVec;
    use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubGraphLike, SubSetLike};
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::parse;

    use crate::cff::VertexSet;
    use crate::graph::FinalizedCut;
    use crate::graph::LoopMomentumBasis;
    use crate::momentum::{FourMomentum, sample::ExternalFourMomenta, signature::LoopExtSignature};
    use crate::{
        cff::{esurface::Esurface, generation::ShiftRewrite},
        utils::{
            DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, ESURFACE_SHIFT_THRESHOLD, F,
            test_utils::dummy_hedge_graph,
        },
    };

    use super::{
        EsurfaceBoundaryError, EsurfaceExistence, SamplingEvaluationError, add_external_shifts,
    };

    #[test]
    fn classification_preserves_the_previous_existing_predicate() {
        let e_cm = F(10.0);
        let shift_tolerance = F(ESURFACE_SHIFT_THRESHOLD) * e_cm;
        let normalized_margin_tolerance = F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD);
        let invariant_tolerance = normalized_margin_tolerance * e_cm * e_cm;

        for shift_factor in [-2, -1, 0, 1, 2] {
            let shift_part = shift_tolerance * shift_tolerance.from_i64(shift_factor);
            for margin_factor in [-2, -1, 0, 1, 2] {
                let invariant_margin =
                    invariant_tolerance * invariant_tolerance.from_i64(margin_factor);
                let was_existing =
                    shift_part < -&shift_tolerance && invariant_margin > invariant_tolerance;
                let classification = Esurface::classify_invariant_margin(
                    &shift_part,
                    invariant_margin,
                    &e_cm,
                    &normalized_margin_tolerance,
                );

                assert_eq!(
                    classification.is_existing(),
                    was_existing,
                    "existence changed for shift factor {shift_factor} and margin factor {margin_factor}",
                );
            }
        }
    }

    /// Run the compiled test binary from the repository root:
    /// ```text
    /// <gammalooprs-test-binary> --ignored --exact \
    ///   cff::esurface::tests::nnlo_example_topological_boundaries_survive_initial_state_sewing \
    ///   --nocapture
    /// ```
    #[test]
    #[ignore = "audit the complete NNLO example; run this binary from the repository root"]
    fn nnlo_example_topological_boundaries_survive_initial_state_sewing() {
        test_initialise().unwrap();
        let model = utils::load_generic_model("sm");
        let directory = Path::new("examples/cli/epem_a_ttxh/NNLO/graphs");
        let mut paths = fs::read_dir(directory)
            .expect("run the NNLO boundary audit from the repository root")
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().is_some_and(|extension| extension == "dot"))
            .collect_vec();
        paths.sort();
        assert!(
            !paths.is_empty(),
            "the NNLO graph inventory must be present"
        );
        let mut boundary_count = 0;
        for path in &paths {
            let graph: Graph = fs::read_to_string(path)
                .unwrap()
                .into_graph(&model)
                .unwrap();
            let (source, target) = graph.get_source_and_target();
            for (left, cut, right) in
                graph.all_st_cuts_for_cs(source, target, &graph.get_initial_state_tree())
            {
                if cut.nedges(&graph) <= 1 {
                    continue;
                }
                let surface = Esurface::new_from_cut_left(
                    &graph.underlying,
                    &CrossSectionCut { left, cut, right },
                    Some(&graph.initial_state_cut),
                );
                surface.oriented_boundary(&graph).unwrap_or_else(|error| {
                    panic!("{}: {error}", path.display());
                });
                boundary_count += 1;
            }
        }
        println!(
            "Certified {boundary_count} topological boundaries on {} sewn NNLO graphs",
            paths.len(),
        );
    }

    #[test]
    fn surface_boundary_excludes_both_halves_of_sewn_initial_edges() {
        test_initialise().unwrap();
        let model = utils::load_generic_model("sm");
        let graph: Graph =
            include_str!("../../../../tests/resources/graphs/ir_safe_thresholds/GL562.dot")
                .into_graph(&model)
                .unwrap();
        let initial_edges = graph.initial_state_cut.as_subgraph();
        let mut sewn_initial_count = 0;
        for (pair, _, _) in graph.underlying.iter_edges() {
            if !initial_edges.includes(&pair) {
                continue;
            }
            let HedgePair::Paired { source, sink } = pair else {
                panic!("initial-state legs must have been sewn into paired edges");
            };
            assert_ne!(
                graph.initial_state_cut.includes(&source),
                graph.initial_state_cut.includes(&sink),
                "an oriented cut includes only one half of a sewn edge",
            );
            assert!(initial_edges.includes(&source) && initial_edges.includes(&sink));
            sewn_initial_count += 1;
        }
        assert_eq!(sewn_initial_count, 2);

        let surface = Esurface {
            energies: [4, 10, 12, 14].map(EdgeIndex).into(),
            external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let boundary = surface.oriented_boundary(&graph).unwrap();
        for (flow, edge, _) in boundary.iter_edges_flow(&graph.underlying) {
            assert_eq!(
                flow == Flow::Source,
                [EdgeIndex(4), EdgeIndex(12)].contains(&edge),
                "incorrect sign on GL562 threshold edge {edge}",
            );
        }
    }

    #[test]
    fn surface_boundary_rejects_independent_component_orientations() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph disconnected_boundary {
            incoming_a [style=invis]
            outgoing_a [style=invis]
            incoming_b [style=invis]
            outgoing_b [style=invis]
            node [num=1]
            edge [num=1 mass=1]
            incoming_a -> A [id=0 mass=0]
            B -> outgoing_a [id=1 mass=0]
            A -> B [id=2 lmb_id=0]
            A -> B [id=3]
            incoming_b -> C [id=4 mass=0]
            D -> outgoing_b [id=5 mass=0]
            C -> D [id=6 lmb_id=1]
            C -> D [id=7]
        })
        .unwrap();
        let surface = Esurface {
            energies: [2, 3, 6, 7].map(EdgeIndex).into(),
            external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(4), -1)],
            vertex_set: VertexSet::dummy(),
        };
        assert!(matches!(
            surface.oriented_boundary(&graph),
            Err(EsurfaceBoundaryError::DisconnectedSupport { .. }),
        ));
    }

    #[test]
    fn projected_gl280_threshold_uses_its_own_signed_boundary() {
        test_initialise().unwrap();
        let model = utils::load_generic_model("sm");
        let graph: Graph =
            include_str!("../../../../tests/resources/graphs/ir_safe_thresholds/GL280.dot")
                .into_graph(&model)
                .unwrap();
        let mut cross_section = CrossSectionGraph::new(graph.clone());
        cross_section.build_lmbs().unwrap();
        let production_lmbs = cross_section.derived_data.lmbs.as_ref().unwrap();
        let masses = graph.get_real_mass_vector::<f64>(&model);
        let invalid = Esurface {
            energies: [5, 9, 13].map(EdgeIndex).into(),
            external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let error = invalid.oriented_boundary(&graph).unwrap_err();
        assert!(matches!(
            &error,
            EsurfaceBoundaryError::NonSeparating { .. }
        ));
        let error: Report = error.into();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_none());
        assert!(error.downcast_ref::<EsurfaceBoundaryError>().is_some());
        let point = [
            [-4.731473194719612, 6.308630926292807, 257.1380514480527],
            [57.0, -31.0, 89.0],
            [-140.22095511751874, 186.961273490025, 311.60212248337507],
            [72.97175339348193, -80.64567119130925, -82.36278531884878],
        ]
        .map(|[x, y, z]| ThreeMomentum::new(F(x), F(y), F(z)));
        // Directed energies are fixed by the simultaneous physical host and
        // threshold cuts. Their linear routing therefore survives a Lorentz boost.
        let energies = [
            (point[0].norm_squared() + F(173.0).square()).sqrt(),
            (point[1].norm_squared() + F(173.0).square()).sqrt(),
            -(point[2].norm_squared() + F(173.0).square()).sqrt(),
            -point[3].norm(),
        ];
        let base_loops = LoopMomenta::from_iter(point);
        let base_externals = ExternalFourMomenta::from_iter([
            FourMomentum::from_args(F(500.0), F(0.0), F(0.0), F(500.0)),
            FourMomentum::from_args(F(500.0), F(0.0), F(0.0), F(-500.0)),
        ]);
        for edges in [vec![2, 4, 9], vec![5, 9, 13, 14]] {
            let surface = Esurface {
                energies: edges.iter().copied().map(EdgeIndex).collect(),
                external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(1), -1)],
                vertex_set: VertexSet::dummy(),
            };
            let eta = surface.compute_from_momenta(
                &graph.loop_momentum_basis,
                &masses,
                &base_loops,
                &base_externals,
            );
            assert!(eta.abs() < F(1.0e-9), "GL280 shell {edges:?}: eta={eta}");
        }
        let signature = &graph.loop_momentum_basis.edge_signatures[EdgeIndex(14)];
        let directed_time = signature
            .internal
            .iter()
            .zip(energies)
            .map(|(sign, energy)| energy * energy.from_i64(i64::from(i8::from(*sign))))
            .chain(
                signature
                    .external
                    .iter()
                    .map(|sign| F(500.0) * F(500.0).from_i64(i64::from(i8::from(*sign)))),
            )
            .fold(F(0.0), |sum, value| sum + value);
        let spatial_externals: ExternalThreeMomenta<F<f64>> =
            base_externals.iter().map(|beam| beam.spatial).collect();
        let q14 = signature.compute_momentum(&base_loops, &spatial_externals);
        assert!(
            (directed_time + q14.norm()).abs() < F(1.0e-9),
            "GL280 q14 must have its negative on-shell temporal component: time={directed_time}, norm={}",
            q14.norm(),
        );
        let boost = |energy: F<f64>, momentum: ThreeMomentum<F<f64>>| {
            FourMomentum::from_args(
                (F(5.0) * energy + F(3.0) * momentum.px) / F(4.0),
                (F(3.0) * energy + F(5.0) * momentum.px) / F(4.0),
                momentum.py,
                momentum.pz,
            )
        };
        let e_cm = F(1000.0);
        let mut reference_margin = None;
        for boosted in [false, true] {
            let native = point.iter().zip(energies).map(|(&momentum, energy)| {
                if boosted {
                    boost(energy, momentum).spatial
                } else {
                    momentum
                }
            });
            let native: LoopMomenta<F<f64>> = LoopMomenta::from_iter(native);
            let beams = [F(500.0), F(-500.0)].map(|z| {
                let momentum = ThreeMomentum::new(F(0.0), F(0.0), z);
                if boosted {
                    boost(F(500.0), momentum)
                } else {
                    FourMomentum::from_args(F(500.0), momentum.px, momentum.py, momentum.pz)
                }
            });
            for reversed in [false, true] {
                let external_moms = ExternalFourMomenta::from_iter(
                    beams
                        .iter()
                        .map(|beam| if reversed { -beam } else { *beam }),
                );
                let spatial_externals: ExternalThreeMomenta<F<f64>> =
                    external_moms.iter().map(|beam| beam.spatial).collect();
                let native = LoopMomenta::from_iter(
                    native
                        .iter()
                        .map(|momentum| if reversed { -momentum } else { *momentum }),
                );
                let surface = Esurface {
                    energies: [5, 9, 13, 14].map(EdgeIndex).into(),
                    external_shift: vec![
                        (EdgeIndex(0), if reversed { 1 } else { -1 }),
                        (EdgeIndex(1), if reversed { 1 } else { -1 }),
                    ],
                    // General CFF surfaces need not carry a topology vertex set.
                    vertex_set: VertexSet::dummy(),
                };
                let boundary = surface.oriented_boundary(&graph).unwrap();
                for (flow, edge, _) in boundary.iter_edges_flow(&graph.underlying) {
                    let expected = (edge == EdgeIndex(5)) != reversed;
                    assert_eq!(flow == Flow::Source, expected);
                }
                for edges in [vec![2, 4, 9], vec![5, 9, 13, 14]] {
                    let shell = Esurface {
                        energies: edges.iter().copied().map(EdgeIndex).collect(),
                        external_shift: surface.external_shift.clone(),
                        vertex_set: VertexSet::dummy(),
                    };
                    let eta = shell.compute_from_momenta(
                        &graph.loop_momentum_basis,
                        &masses,
                        &native,
                        &external_moms,
                    );
                    assert!(
                        eta.abs() < F(1.0e-9),
                        "GL280 transformed shell {edges:?}: boost={boosted}, reversed={reversed}, eta={eta}",
                    );
                }
                for parent in [[2, 3, 5, 12], [2, 5, 9, 12]] {
                    let parent = parent.map(EdgeIndex);
                    let lmbs = ti_vec![
                        production_lmbs
                            .iter()
                            .find(|lmb| {
                                lmb.loop_edges.len() == parent.len()
                                    && parent.iter().all(|edge| lmb.loop_edges.contains(edge))
                            })
                            .expect("production must generate the requested physical parent")
                            .clone()
                    ];
                    let lmb_id = LmbIndex::from(0);
                    let lmb = &lmbs[lmb_id];
                    assert_eq!(lmb.loop_edges.len(), 4, "parent={parent:?}");
                    assert_eq!(lmb.ext_edges.len(), 2, "parent={parent:?}");
                    assert_eq!(lmb.ext_edges, graph.loop_momentum_basis.ext_edges);
                    let subspace = SubspaceData::new_from_parent_basis_edges(
                        &[EdgeIndex(5)],
                        &graph.full_filter(),
                        lmb_id,
                        &graph,
                        &lmbs,
                    )
                    .unwrap();
                    let loops = LoopMomenta::from_iter(lmb.loop_edges.iter().map(|edge| {
                        graph.loop_momentum_basis.edge_signatures[*edge]
                            .compute_momentum(&native, &spatial_externals)
                    }));
                    let classification = surface
                        .classify_existence_subspace(
                            &loops,
                            &external_moms,
                            &subspace,
                            &lmbs,
                            &graph,
                            &masses,
                            &e_cm,
                            &F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
                        )
                        .unwrap();
                    let active = subspace.iter_lmb_indices().next().unwrap();
                    let mut center = loops.clone();
                    center[active] = ThreeMomentum::new(F(0.0), F(0.0), F(0.0));
                    let offset = -lmb.edge_signatures[EdgeIndex(13)]
                        .compute_momentum(&center, &spatial_externals);
                    center[active] = offset;
                    let available = -surface.compute_shift_part_from_momenta_in_subspace(
                        &loops,
                        &external_moms,
                        &subspace,
                        &lmbs,
                        &graph,
                        &masses,
                    );
                    let invariant =
                        available.square() - offset.norm_squared() - masses[EdgeIndex(5)].square();
                    println!(
                        "GL280 parent={parent:?} boost={boosted} reversed={reversed} loops={:?} external={:?} mass={} shift={} fixed_vector={offset:?} invariant={} status={}",
                        lmb.loop_edges,
                        lmb.ext_edges,
                        masses[EdgeIndex(5)],
                        -available,
                        invariant,
                        classification.label(),
                    );
                    assert!(
                        classification.is_existing(),
                        "parent={parent:?}, boost={boosted}, reversed={reversed}, shift={}, offset={offset:?}, invariant={invariant}, classification={classification:?}",
                        -available,
                    );
                    let margin = *classification.normalized_margin().unwrap();
                    assert!(margin > F(0.05) && margin < F(0.06));
                    let reference = *reference_margin.get_or_insert(margin);
                    assert!((margin - reference).abs() < F(1.0e-10));

                    assert!(
                        surface.compute_from_momenta(lmb, &masses, &center, &external_moms)
                            < F(-1.0),
                    );
                    if !boosted && !reversed {
                        let fixed_energy = e_cm - available;
                        let minimum =
                            (offset.norm_squared() + masses[EdgeIndex(5)].square()).sqrt();
                        for (delta, expected) in [
                            (F(-1.0), "non_existing"),
                            (F(0.0), "pinched"),
                            (F(1.0), "existing"),
                        ] {
                            let energy = (fixed_energy + minimum + delta) / F(2.0);
                            let beams = ExternalFourMomenta::from_iter([
                                FourMomentum::from_args(energy, F(0.0), F(0.0), energy),
                                FourMomentum::from_args(energy, F(0.0), F(0.0), -energy),
                            ]);
                            let status = surface
                                .classify_existence_subspace(
                                    &loops,
                                    &beams,
                                    &subspace,
                                    &lmbs,
                                    &graph,
                                    &masses,
                                    &e_cm,
                                    &F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
                                )
                                .unwrap();
                            assert_eq!(status.label(), expected);
                        }
                    }
                    for direction in [F(-1.0), F(1.0)] {
                        let radius = invariant / (F(2.0) * (available + direction * offset.px));
                        assert!(radius > F(0.0));
                        let mut root = center.clone();
                        root[active].px += direction * radius;
                        assert!(
                            surface
                                .compute_from_momenta(lmb, &masses, &root, &external_moms)
                                .abs()
                                < F(1.0e-9),
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn projected_gl358_threshold_distinguishes_open_pinched_and_closed_fibres() {
        test_initialise().unwrap();
        let model = utils::load_generic_model("sm");
        let graph: Graph =
            include_str!("../../../../tests/resources/graphs/ir_safe_thresholds/GL358.dot")
                .into_graph(&model)
                .unwrap();
        let mut cross_section = CrossSectionGraph::new(graph.clone());
        cross_section.build_lmbs().unwrap();
        let masses = graph.get_real_mass_vector::<f64>(&model);
        let parent = [3, 6, 7, 11].map(EdgeIndex);
        let lmbs = ti_vec![
            cross_section
                .derived_data
                .lmbs
                .as_ref()
                .unwrap()
                .iter()
                .find(|lmb| {
                    lmb.loop_edges.len() == parent.len()
                        && parent.iter().all(|edge| lmb.loop_edges.contains(edge))
                })
                .expect("production must generate the requested physical parent")
                .clone()
        ];
        let lmb_id = LmbIndex::from(0);
        let lmb = &lmbs[lmb_id];
        assert_eq!(lmb.loop_edges.len(), 4);
        assert_eq!(lmb.ext_edges.len(), 2);
        assert_eq!(lmb.ext_edges, graph.loop_momentum_basis.ext_edges);
        let subspace = SubspaceData::new_from_parent_basis_edges(
            &[EdgeIndex(3)],
            &graph.full_filter(),
            lmb_id,
            &graph,
            &lmbs,
        )
        .unwrap();
        let external_moms = ExternalFourMomenta::from_iter([
            FourMomentum::from_args(F(500.0), F(0.0), F(0.0), F(500.0)),
            FourMomentum::from_args(F(500.0), F(0.0), F(0.0), F(-500.0)),
        ]);
        let spatial_externals: ExternalThreeMomenta<F<f64>> =
            external_moms.iter().map(|p| p.spatial).collect();
        let surface = Esurface {
            energies: [3, 6, 12, 13].map(EdgeIndex).into(),
            external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        // Three points already rescaled on physical host [2,6,7]. The old
        // host-oriented invariant admitted the closed side and failed center finding.
        for (point, expected) in [
            (
                [
                    239.41681001042383,
                    -114.46145845228274,
                    129.29264510269138,
                    239.4167273843055,
                    -114.46178575901338,
                    129.2926253365453,
                    95.59905078888387,
                    -225.95843361500184,
                    -121.37915228406406,
                    342.0091274491519,
                    101.42159085418247,
                    134.8239582852383,
                ],
                "existing",
            ),
            (
                [
                    239.41676577049503,
                    -114.46119915654612,
                    129.2927988212799,
                    239.41676577049503,
                    -114.46119915654612,
                    129.2927988212799,
                    95.59919031992825,
                    -225.9582666603062,
                    -121.37910999326795,
                    342.0093104427367,
                    101.42231992498913,
                    134.82413188309832,
                ],
                "pinched",
            ),
            (
                [
                    239.41672153007048,
                    -114.46093986071739,
                    129.2929525394926,
                    239.41680415613854,
                    -114.4606125541859,
                    129.29297230562662,
                    95.59932985067904,
                    -225.95809970521887,
                    -121.37906770223259,
                    342.00949343546347,
                    101.42304899513076,
                    134.82430548055834,
                ],
                "non_existing",
            ),
        ] {
            let loops = LoopMomenta::from_iter(lmb.loop_edges.iter().map(|edge| {
                let index = parent.iter().position(|parent| parent == edge).unwrap();
                let [x, y, z] = point.as_chunks::<3>().0[index];
                ThreeMomentum::new(F(x), F(y), F(z))
            }));
            let status = surface
                .classify_existence_subspace(
                    &loops,
                    &external_moms,
                    &subspace,
                    &lmbs,
                    &graph,
                    &masses,
                    &F(1000.0),
                    &F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
                )
                .unwrap();
            assert_eq!(
                status.label(),
                expected,
                "parent={parent:?}, loops={:?}, externals={:?}, classification={status:?}",
                lmb.loop_edges,
                lmb.ext_edges,
            );
            let x = lmb.edge_signatures[EdgeIndex(6)].compute_momentum(&loops, &spatial_externals);
            let d = lmb.edge_signatures[EdgeIndex(12)].compute_momentum(&loops, &spatial_externals);
            let minimum = ((x + d).norm_squared() + F(173.0).square()).sqrt()
                + (x.norm_squared() + F(173.0).square()).sqrt()
                + d.norm()
                - F(1000.0);
            let mut center = loops.clone();
            center[subspace.iter_lmb_indices().next().unwrap()] = x + d;
            assert!(
                (surface.compute_from_momenta(lmb, &masses, &center, &external_moms) - minimum)
                    .abs()
                    < F(1.0e-10),
            );
        }
    }

    #[test]
    fn invalid_programmatic_tolerances_cannot_invert_classification() {
        let e_cm = F(10.0);
        let shift_part = F(-1.0);
        let invariant_margin = F(0.0);

        for tolerance in [
            F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            F(-DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            F(f64::NAN),
            F(f64::INFINITY),
        ] {
            assert!(matches!(
                Esurface::classify_invariant_margin(
                    &shift_part,
                    invariant_margin,
                    &e_cm,
                    &tolerance,
                ),
                EsurfaceExistence::Pinched { .. }
            ));
        }
    }

    #[test]
    fn positive_external_energy_filter_uses_energy_conservation() {
        let incoming_edges = (0..4).map(EdgeIndex::from).collect_vec();
        let outgoing_edges = (4..9).map(EdgeIndex::from).collect_vec();
        let shift_is_negative = |external_shift: &[(usize, i64)]| {
            Esurface {
                energies: vec![],
                external_shift: external_shift
                    .iter()
                    .map(|(edge, coefficient)| (EdgeIndex::from(*edge), *coefficient))
                    .collect(),
                vertex_set: VertexSet::dummy(),
            }
            .external_shift_is_strictly_negative_for_positive_energies(
                &incoming_edges,
                &outgoing_edges,
            )
        };

        assert!(
            shift_is_negative(&[(0, -1), (1, -1)]),
            "a negative proper subset of incoming energies must be retained"
        );
        assert!(
            shift_is_negative(&[(4, -1)]),
            "a negative proper subset of outgoing energies must be retained"
        );
        assert!(
            shift_is_negative(&[(0, -1), (1, -1), (2, -1), (3, -1), (4, 1)]),
            "energy conservation turns minus all incoming plus one outgoing into minus the remaining outgoing energies"
        );
        assert!(
            !shift_is_negative(&[
                (0, -1),
                (1, -1),
                (2, -1),
                (3, -1),
                (4, 1),
                (5, 1),
                (6, 1),
                (7, 1),
                (8, 1),
            ]),
            "the energy-conservation identity is zero rather than strictly negative"
        );
        assert!(
            !shift_is_negative(&[(0, -1), (4, 1)]),
            "a sign-indefinite difference must not be accepted from positivity alone"
        );
        assert!(
            !shift_is_negative(&[(9, -1)]),
            "a shift outside the external-energy partition must be rejected conservatively"
        );
    }

    #[test]
    fn massless_two_to_two_surface_is_existing_away_from_collinear_pinch() {
        let dummy_graph = dummy_hedge_graph(4);
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![EdgeIndex::from(2)].into(),
            ext_edges: vec![].into(),
            edge_signatures: dummy_graph
                .new_edgevec_from_iter(vec![
                    LoopExtSignature::from((vec![0], vec![1, 0])),
                    LoopExtSignature::from((vec![0], vec![0, 1])),
                    LoopExtSignature::from((vec![1], vec![0, 0])),
                    LoopExtSignature::from((vec![-1], vec![-1, -1])),
                ])
                .unwrap(),
        };
        let esurface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift: vec![(EdgeIndex::from(0), -1), (EdgeIndex::from(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = dummy_graph.new_edgevec_from_iter(vec![F(0.0); 4]).unwrap();
        // Either side of the 2->2 sandwich determines the same total four-momentum by
        // momentum conservation, so the incoming pair is sufficient for this classification.
        let classify_pair = |second_spatial_momentum: (f64, f64), threshold: f64| {
            let external_momenta = ExternalFourMomenta::from_iter([
                FourMomentum::from_args(F(5.0), F(5.0), F(0.0), F(0.0)),
                FourMomentum::from_args(
                    F(5.0),
                    F(second_spatial_momentum.0),
                    F(second_spatial_momentum.1),
                    F(0.0),
                ),
            ]);
            esurface.classify_existence(&external_momenta, &lmb, &masses, &F(10.0), &F(threshold))
        };

        assert!(matches!(
            classify_pair((0.0, 5.0), DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            EsurfaceExistence::Existing { .. }
        ));
        assert!(matches!(
            classify_pair((5.0, 0.0), DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classify_pair((0.0, 5.0), 1.0),
            EsurfaceExistence::Pinched { .. }
        ));
    }

    #[test]
    fn classifies_existing_pinched_and_non_existing_surfaces() {
        let dummy_graph = dummy_hedge_graph(5);
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![EdgeIndex::from(2), EdgeIndex::from(3)].into(),
            ext_edges: vec![].into(),
            edge_signatures: dummy_graph
                .new_edgevec_from_iter(vec![
                    LoopExtSignature::from((vec![0, 0], vec![1])),
                    LoopExtSignature::from((vec![0, 0], vec![-1])),
                    LoopExtSignature::from((vec![1, 0], vec![0])),
                    LoopExtSignature::from((vec![0, 1], vec![0])),
                    LoopExtSignature::from((vec![1, 1], vec![-1])),
                ])
                .unwrap(),
        };
        let esurface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3), EdgeIndex::from(4)],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = dummy_graph.new_edgevec_from_iter(vec![F(0.0); 5]).unwrap();

        let classification = |energy| {
            let external_momenta = ExternalFourMomenta::from_iter([FourMomentum::from_args(
                F(energy),
                F(-10.0),
                F(0.0),
                F(0.0),
            )]);
            esurface.classify_existence(
                &external_momenta,
                &lmb,
                &masses,
                &F(10.0),
                &F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            )
        };

        assert!(matches!(
            classification(11.0),
            EsurfaceExistence::Existing { .. }
        ));
        assert!(matches!(
            classification(10.0),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classification(10.0 + 1.0e-8),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classification(9.0),
            EsurfaceExistence::NonExisting { .. }
        ));
    }

    #[test]
    fn test_esurface() {
        let dummy_graph = dummy_hedge_graph(5);

        let _energies_cache = dummy_graph
            .new_edgevec_from_iter([F(1.), F(2.), F(3.), F(4.), F(5.)])
            .unwrap();

        let energies = vec![EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(2)];

        let external_shift = vec![(EdgeIndex::from(3), 1), (EdgeIndex::from(4), 1)];

        let mut esurface = Esurface {
            energies,
            external_shift,
            vertex_set: VertexSet::dummy(),
            //subspace_graph: dummy_graph.full_graph(),
        };

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: EdgeIndex::from(4),
            dependent_momentum_expr: vec![
                (EdgeIndex::from(1), -1),
                (EdgeIndex::from(2), -1),
                (EdgeIndex::from(3), -1),
            ],
        };

        esurface.canonicalize_shift(&shift_rewrite);

        assert_eq!(
            esurface.external_shift,
            vec![(EdgeIndex::from(1), -1), (EdgeIndex::from(2), -1)]
        );

        let energies = vec![EdgeIndex::from(0), EdgeIndex::from(2)];

        let external_shift = vec![(EdgeIndex::from(1), -1)];

        let _esurface = Esurface {
            energies,
            external_shift,
            vertex_set: VertexSet::dummy(),
            //subspace_graph: dummy_graph.full_graph(),
        };
    }

    #[test]
    fn test_add_external_shifts() {
        let shift_1 = vec![
            (EdgeIndex::from(0), 1),
            (EdgeIndex::from(1), 1),
            (EdgeIndex::from(2), -1),
        ];
        let shift_2 = vec![(EdgeIndex::from(1), -1), (EdgeIndex::from(2), 1)];

        let add = add_external_shifts(&shift_1, &shift_2);

        assert_eq!(add, vec![(EdgeIndex::from(0), 1)]);

        let shift_3 = vec![(EdgeIndex::from(3), 1), (EdgeIndex::from(4), -1)];
        let shift_4 = vec![
            (EdgeIndex::from(0), 1),
            (EdgeIndex::from(1), 1),
            (EdgeIndex::from(2), 1),
            (EdgeIndex::from(4), 1),
        ];

        let add = add_external_shifts(&shift_3, &shift_4);

        assert_eq!(
            add,
            vec![
                (EdgeIndex::from(0), 1),
                (EdgeIndex::from(1), 1),
                (EdgeIndex::from(2), 1),
                (EdgeIndex::from(3), 1)
            ]
        );
    }

    #[test]
    fn test_esurface_equality() {
        let esurface_1 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
            //subspace_graph: unsafe { InternalSubGraph::new_unchecked(SuBitGraph::new()) },
        };

        let esurface_2 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
            //subspace_graph: unsafe { InternalSubGraph::new_unchecked(SuBitGraph::new()) },
        };

        assert_eq!(esurface_1, esurface_2);
    }

    mod failing {
        use super::*;

        #[test]
        fn test_to_atom() {
            let external_shift = vec![(EdgeIndex::from(1), -1)];

            let esurface = Esurface {
                energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
                external_shift,
                vertex_set: VertexSet::dummy(),
                // subspace_graph: unsafe { InternalSubGraph::new_unchecked(SuBitGraph::new()) },
            };

            let esurface_atom = esurface.to_atom(&[]);
            let expected_atom = parse!("Q(2, cind(0)) + Q(3, cind(0)) - P(1, cind(0))");

            let diff = esurface_atom - &expected_atom;
            let diff = diff.expand();
            assert_eq!(diff, Atom::new());
        }

        #[test]
        fn test_from_cut_left_dt() {
            let mut hedge_graph_builder = HedgeGraphBuilder::new();
            let nodes = (0..4)
                .map(|_| hedge_graph_builder.add_node(()))
                .collect_vec();

            hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

            hedge_graph_builder.add_external_edge(
                nodes[0],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
            hedge_graph_builder.add_external_edge(
                nodes[3],
                (),
                Orientation::Undirected,
                Flow::Source,
            );

            let double_triangle: HedgeGraph<(), (), ()> =
                hedge_graph_builder.build::<NodeStorageVec<()>>();
            let node_0 = double_triangle.iter_crown(nodes[0]).into();
            let node_3 = double_triangle.iter_crown(nodes[3]).into();

            let cuts = double_triangle.all_cuts(node_0, node_3);

            let cross_section_cuts = cuts
                .into_iter()
                .map(|(node_l, cut, node_r)| FinalizedCut {
                    cut,
                    left: node_l,
                    right: node_r,
                })
                .map(|cut| Esurface::new_from_cut_left(&double_triangle, &cut, None))
                .collect_vec();

            let expected_esurfaces = vec![
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(2), EdgeIndex::from(4)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(3), EdgeIndex::from(4)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(1), EdgeIndex::from(2), EdgeIndex::from(3)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
            ];

            for expected_esurface in expected_esurfaces {
                assert!(cross_section_cuts.contains(&expected_esurface));
            }
        }

        #[test]
        fn test_from_cut_left_box() {
            let mut hedge_graph_builder = HedgeGraphBuilder::new();
            let nodes = (0..4)
                .map(|_| hedge_graph_builder.add_node(()))
                .collect_vec();

            hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[3], nodes[0], (), Orientation::Undirected);

            hedge_graph_builder.add_external_edge(
                nodes[0],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
            hedge_graph_builder.add_external_edge(
                nodes[1],
                (),
                Orientation::Undirected,
                Flow::Source,
            );
            hedge_graph_builder.add_external_edge(
                nodes[2],
                (),
                Orientation::Undirected,
                Flow::Source,
            );
            hedge_graph_builder.add_external_edge(
                nodes[3],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );

            let box_graph: HedgeGraph<(), (), ()> =
                hedge_graph_builder.build::<NodeStorageVec<()>>();

            let node_0 = box_graph.iter_crown(nodes[0]).into();
            let node_2 = box_graph.iter_crown(nodes[2]).into();

            let cuts = box_graph.all_cuts(node_0, node_2);
            assert_eq!(cuts.len(), 4);

            let cross_section_cuts = cuts
                .into_iter()
                .map(|(node_l, cut, node_r)| FinalizedCut {
                    cut,
                    left: node_l,
                    right: node_r,
                })
                .map(|cut| Esurface::new_from_cut_left(&box_graph, &cut, None))
                .collect_vec();

            let expected_esurfaces = vec![
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(3)],
                    external_shift: vec![(EdgeIndex::from(4), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(2)],
                    external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(7), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(1), EdgeIndex::from(3)],
                    external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(1), EdgeIndex::from(2)],
                    external_shift: vec![
                        (EdgeIndex::from(4), -1),
                        (EdgeIndex::from(5), 1),
                        (EdgeIndex::from(7), -1),
                    ],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
            ];

            for expected_esurface in expected_esurfaces {
                assert!(cross_section_cuts.contains(&expected_esurface));
            }
        }
    }
}
