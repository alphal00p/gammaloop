use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Flow, HedgePair};
use linnet::half_edge::subgraph::OrientedCut;
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};

use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::dual::HyperDual;
use symbolica::domains::float::{FloatLike as SymFloatLike, Real};
use symbolica::id::Replacement;
use symbolica::{function, parse};
use tracing::debug;
use typed_index_collections::TiVec;

use crate::cff::VertexSet;

use crate::cff::expression::{
    CFFExpression, OrientationID, RaisedEsurfaceDataView, RaisedEsurfaceGroupView,
};
pub use crate::cff::surface::EsurfaceID;
use crate::graph::{Graph, GraphGroupPosition, LmbIndex, LoopMomentumBasis};
use crate::{GammaLoopContext, define_index};

use crate::integrands::process::{GenericEvaluator, ImplicitSurfaceRadialMap};
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta, SubspaceData,
};
use crate::processes::CrossSectionCut;
use crate::utils::hyperdual_utils::new_constant;
use crate::utils::newton_solver::{
    NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity, SafeguardedNewtonError,
};
use crate::utils::{
    DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, ESURFACE_SHIFT_THRESHOLD, F, FloatLike, GS, Length,
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

// Edge identity, radial velocity, constant spatial offset, and mass.
type EsurfaceRayEnergy<T> = (EdgeIndex, ThreeMomentum<F<T>>, ThreeMomentum<F<T>>, F<T>);

/// One represented affine energy equation, retaining its ordered edge occurrences.
/// This is native geometry, not a cross-precision or selected-host transport record.
#[derive(Debug, Clone)]
pub(crate) struct EsurfaceRay<T: FloatLike> {
    // A Vec deliberately preserves repeated energy occurrences; an edge-keyed
    // map would not. The tuple's fields are documented on EsurfaceRayEnergy.
    energies: Vec<EsurfaceRayEnergy<T>>,
    shift: F<T>,
}

impl<T: FloatLike> EsurfaceRay<T> {
    #[inline]
    pub(crate) fn evaluate(&self, radius: &F<T>) -> (F<T>, F<T>) {
        let zero = radius.zero();
        let (derivative, energy_sum) = self
            .energies
            .iter()
            .map(|(_, velocity, offset, mass)| {
                let momentum = velocity * radius + offset;
                let energy = (momentum.norm_squared() + mass * mass).sqrt();
                // At a massless endpoint E(r)=r|v| the radial right derivative is
                // |v|, whereas the two-sided formula q.v/E would evaluate 0/0.
                // Share this convention with sampling charts, using exact source
                // zeros: a computed E=0 can instead be numerical underflow.
                let derivative = if radius == &zero
                    && mass == &zero
                    && momentum.px == zero
                    && momentum.py == zero
                    && momentum.pz == zero
                {
                    velocity.norm_squared().sqrt()
                } else {
                    momentum * velocity / &energy
                };
                (derivative, energy)
            })
            .fold(
                (zero.clone(), zero.clone()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );
        (energy_sum + &self.shift, derivative)
    }

    /// Original eta jet at the same represented coefficients used by the root.
    /// The caller retains the existing dual shape and factorial convention.
    #[inline]
    pub(crate) fn evaluate_dual(&self, radius: &HyperDual<F<T>>) -> HyperDual<F<T>> {
        let energy_sum = self
            .energies
            .iter()
            .map(|(_, velocity, offset, mass)| {
                let momentum = velocity.map_ref(&|v| new_constant(radius, v) * radius)
                    + offset.map_ref(&|b| new_constant(radius, b));
                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|sum, energy| sum + energy)
            .unwrap_or_else(|| radius.zero());
        energy_sum + new_constant(radius, &self.shift)
    }
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
    /// Match two energy sums on one active spatial momentum. Their shared
    /// energy is identified by a single sign of the complete affine routing
    /// and the mass expression, never by edge identity or rounded momenta.
    /// The returned common route has active coefficient +1: the kernel uses
    /// x=L+c0, so the existing affine embedding must return L=x-c0.
    // Physical binders must transport one prerequisite-only disk policy
    // through native precision retries alongside this reconstructed geometry.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sampling_joint_geometry_in_subspace<T: FloatLike>(
        &self,
        other: &Self,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        complement: &[LoopIndex],
    ) -> Result<(
        crate::integrands::process::sampling_joint::SharedEnergyJointGeometryEvaluator<T>,
        crate::momentum::signature::LoopExtSignature,
    )> {
        use crate::integrands::process::{
            SharedEnergyJointGeometry, sampling_maps::SamplingEvaluationError,
        };
        use crate::momentum::{SignOrZero, signature::LoopExtSignature};
        use std::sync::Arc;

        let lmb = subspace.get_lmb(all_lmbs);
        let active = subspace.iter_lmb_indices().collect_vec();
        let [active] = active.as_slice() else {
            return Err(eyre!(
                "joint energy sampling requires exactly one active loop cycle"
            ));
        };
        let active = *active;
        let mut covered = complement.iter().copied().chain([active]).collect_vec();
        covered.sort();
        if covered.iter().any(|index| index.0 >= lmb.loop_edges.len())
            || covered.windows(2).any(|pair| pair[0] == pair[1])
        {
            return Err(eyre!(
                "joint energy sampling has repeated or invalid active/prior loop indices"
            ));
        }
        if external_momenta.is_empty() || external_momenta.len() != lmb.ext_edges.len() {
            return Err(eyre!(
                "joint energy sampling requires all {} external ports, received {}",
                lmb.ext_edges.len(),
                external_momenta.len()
            ));
        }
        let surfaces = [self, other];
        for &surface in &surfaces {
            for &edge in &surface.energies {
                let signature = &lmb.edge_signatures[edge];
                for index in (0..lmb.loop_edges.len())
                    .map(LoopIndex)
                    .filter(|index| !covered.contains(index))
                {
                    if signature.internal[index] != SignOrZero::Zero {
                        return Err(eyre!(
                            "joint energy surface {:?} depends on unsampled parent edge {}",
                            surface.energies,
                            lmb.loop_edges[index]
                        ));
                    }
                }
                if !masses[edge].0.is_finite() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "joint energy mass",
                        detail: format!("non-finite mass on edge {edge}"),
                    }
                    .into());
                }
                if masses[edge] < masses[edge].zero() {
                    return Err(eyre!(
                        "joint energy edge {edge} requires a nonnegative mass"
                    ));
                }
            }
        }
        let varying = surfaces.map(|surface| {
            surface
                .energies
                .iter()
                .copied()
                .filter(|edge| lmb.edge_signatures[*edge].internal[active] != SignOrZero::Zero)
                .sorted_by_key(|edge| edge.0)
                .collect_vec()
        });
        if varying.iter().any(|edges| edges.len() != 2) {
            return Err(eyre!(
                "joint energy sampling requires exactly two varying energy occurrences per equation, got {:?}",
                varying
            ));
        }
        for (surface, varying) in surfaces.iter().zip(&varying) {
            let contained = subspace.contains(&surface.energies, graph).collect_vec();
            if varying.iter().any(|edge| !contained.contains(edge)) {
                return Err(eyre!(
                    "joint varying energy routing is inconsistent with its selected cycle subgraph"
                ));
            }
        }
        let shared = varying[0]
            .iter()
            .enumerate()
            .flat_map(|(left, &a)| {
                varying[1]
                    .iter()
                    .enumerate()
                    .filter_map(move |(right, &b)| {
                        (lmb.edge_signatures[a].equality_up_to_sign(&lmb.edge_signatures[b])
                            && graph[a].mass_atom() == graph[b].mass_atom())
                        .then_some((left, right))
                    })
            })
            .collect_vec();
        let [(left, right)] = shared.as_slice() else {
            return Err(eyre!(
                "joint energy sampling requires exactly one common full signed routing and mass expression; found {} candidates for {:?}. Distinct external shifts, including fixed-kinematics spatial specializations, are not supported",
                shared.len(),
                varying
            ));
        };
        let edges = [
            varying[0][*left],
            varying[0][1 - *left],
            varying[1][1 - *right],
        ];
        if masses[edges[0]] != masses[varying[1][*right]] {
            return Err(eyre!(
                "joint shared mass expression has inconsistent native values on edges {} and {}",
                edges[0],
                varying[1][*right]
            ));
        }
        let canonical = |edge| {
            let signature = &lmb.edge_signatures[edge];
            if signature.internal[active] == SignOrZero::Minus {
                LoopExtSignature {
                    internal: signature.internal.iter().map(|sign| -*sign).collect(),
                    external: signature.external.iter().map(|sign| -*sign).collect(),
                }
            } else {
                signature.clone()
            }
        };
        let routes = edges.map(canonical);
        let common = routes[0].clone();
        let fixed = surfaces.map(|surface| {
            surface
                .energies
                .iter()
                .copied()
                .filter(|edge| lmb.edge_signatures[*edge].internal[active] == SignOrZero::Zero)
                .collect_vec()
        });
        // Keep the original global equations. Replacing this shift by a
        // cut-eliminated identity would move the target by a finite LU residual.
        let shifts =
            surfaces.map(|surface| surface.compute_shift_part_from_momenta(external_momenta, lmb));
        let externals: ExternalThreeMomenta<F<T>> =
            external_momenta.iter().map(|p| p.spatial.clone()).collect();
        if shifts
            .iter()
            .chain(externals.iter().flat_map(|p| [&p.px, &p.py, &p.pz]))
            .any(|value| !value.0.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint external data",
                detail: "non-finite external spatial vector or original energy shift".into(),
            }
            .into());
        }
        let lmb = lmb.clone();
        let masses = masses.clone();
        let complement = complement.to_vec();
        let geometry = Arc::new(move |context: &[T]| {
            if context.len() != 3 * complement.len() {
                return Err(eyre!(
                    "joint energy context has {} components, expected {}",
                    context.len(),
                    3 * complement.len()
                ));
            }
            if context.iter().any(|x| !x.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "joint energy context",
                    detail: "non-finite declared prerequisite".into(),
                }
                .into());
            }
            let zero = masses[edges[0]].zero();
            let mut loops = LoopMomenta::from_iter(
                (0..lmb.loop_edges.len())
                    .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            for (&index, point) in complement.iter().zip(context.chunks_exact(3)) {
                loops[index] = ThreeMomentum::new(
                    F(point[0].clone()),
                    F(point[1].clone()),
                    F(point[2].clone()),
                );
            }
            let offsets: [ThreeMomentum<F<T>>; 3] = routes
                .each_ref()
                .map(|route| route.compute_momentum(&loops, &externals));
            let differences = [&offsets[1] - &offsets[0], &offsets[2] - &offsets[0]];
            let sums = fixed.iter().enumerate().map(|(index, fixed)| {
                let fixed_sum = fixed.iter().try_fold(zero.clone(), |sum, &edge| -> Result<F<T>> {
                    let momentum: ThreeMomentum<F<T>> = lmb.edge_signatures[edge].compute_momentum(&loops, &externals);
                    let energy = (momentum.norm_squared()+masses[edge].square()).sqrt();
                    if energy == zero && (masses[edge] != zero || momentum.px != zero || momentum.py != zero || momentum.pz != zero) {
                        return Err(SamplingEvaluationError::Unrepresentable { operation: "joint fixed energy", detail: format!("nonzero mass or momentum on edge {edge} produced zero energy") }.into());
                    }
                    Ok(sum+energy)
                })?;
                Ok((-(&shifts[index]+fixed_sum)).0)
            }).collect::<Result<Vec<T>>>()?;
            let energy_sums = [sums[0].clone(), sums[1].clone()];
            let geometry = SharedEnergyJointGeometry {
                shifts: differences.map(|v| [v.px.0, v.py.0, v.pz.0]),
                masses: edges.map(|edge| masses[edge].0.clone()),
                energy_sums,
            };
            if geometry
                .shifts
                .iter()
                .flatten()
                .chain(&geometry.energy_sums)
                .any(|x| !x.is_finite())
            {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "joint prepared energies",
                    detail: "non-finite routed offset or fixed energy sum".into(),
                }
                .into());
            }
            Ok(geometry)
        });
        Ok((geometry, common))
    }

    /// Compile an exact radial chart around this graph-routed energy surface.
    /// The callback retains the complete parent frame, masses and external
    /// data in the evaluation precision; its ray root is the same cut equation
    /// used by LU, without promoting already rounded lower-precision vectors.
    pub(crate) fn sampling_radial_map<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        masses: &EdgeVec<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        beta: f64,
        power: f64,
    ) -> Result<ImplicitSurfaceRadialMap<T>> {
        let dimension = 3 * lmb.loop_edges.len();
        let surface = self.clone();
        let lmb = lmb.clone();
        let masses = masses.clone();
        let external_momenta = external_momenta.clone();
        let zero = F::<T>::from_f64(0.0);
        let center = LoopMomenta::from_iter(
            (0..dimension / 3)
                .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        let evaluator = std::sync::Arc::new(move |direction: &[T], radius: T| {
            let radius = F(radius);
            let unit_loops = LoopMomenta::from_iter(direction.chunks_exact(3).map(|components| {
                ThreeMomentum::new(
                    F(components[0].clone()),
                    F(components[1].clone()),
                    F(components[2].clone()),
                )
            }));
            let (value, derivative) = surface.compute_self_and_r_derivative(
                &radius,
                &unit_loops,
                &center,
                &external_momenta,
                &masses,
                &lmb,
            );
            Ok((value.0, derivative.0))
        });
        ImplicitSurfaceRadialMap::new(dimension, vec![zero.0; dimension], beta, power, evaluator)
    }

    /// Bind a routed energy surface on active coordinates, conditional on an
    /// already sampled complement in the same complete parent LMB. Later
    /// coordinates may be omitted only when every energy row vanishes on them.
    /// Preparation
    /// never consumes active cube coordinates, so its derivatives occupy only
    /// the off-diagonal block of an ordered composition's Jacobian. Active
    /// outputs follow `subspace.iter_lmb_indices()` (parent-index order), while
    /// complement inputs follow the explicit `complement` order. A compiler
    /// must permute these blocks before attaching user-ordered edge metadata.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sampling_radial_map_in_subspace<T: FloatLike>(
        &self,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        complement: &[LoopIndex],
        settings: &crate::settings::RuntimeSettings,
        beta: f64,
        power: f64,
    ) -> Result<ImplicitSurfaceRadialMap<T>> {
        use crate::integrands::process::PreparedSurfaceStatus;
        use crate::integrands::process::sampling_maps::SamplingEvaluationError;
        use crate::momentum::SignOrZero;
        use crate::subtraction::overlap_subspace::{OverlapInput, find_center};
        use std::sync::Arc;
        use three_dimensional_reps::utils::rank_i64;

        if !settings.kinematics.e_cm.is_finite() || settings.kinematics.e_cm <= 0.0 {
            return Err(eyre!("sampling fiber requires finite positive e_cm"));
        }
        let lmb = subspace.get_lmb(all_lmbs);
        let active = subspace.iter_lmb_indices().collect_vec();
        let mut covered = active.iter().chain(complement).copied().collect_vec();
        covered.sort();
        if active.is_empty()
            || covered.windows(2).any(|pair| pair[0] == pair[1])
            || covered.iter().any(|index| index.0 >= lmb.loop_edges.len())
            || external_momenta.len() != lmb.ext_edges.len()
            || external_momenta.is_empty()
        {
            return Err(eyre!(
                "sampling fiber requires disjoint active/complement slots in its parent frame and complete external momenta: active {active:?}, complement {complement:?}, parent {:?}, external count {} (expected {})",
                lmb.loop_edges,
                external_momenta.len(),
                lmb.ext_edges.len(),
            ));
        }
        // Later blocks may be absent from this prefix only when their exact
        // routed rows vanish. Their zero placeholders then represent proven
        // spectators, never unknown active coordinates of the energy equation.
        for index in (0..lmb.loop_edges.len())
            .map(LoopIndex)
            .filter(|index| !covered.contains(index))
        {
            if self
                .energies
                .iter()
                .any(|edge| lmb.edge_signatures[*edge].internal[index] != SignOrZero::Zero)
            {
                return Err(eyre!(
                    "sampling fiber depends on unsampled parent edge {}",
                    lmb.loop_edges[index]
                ));
            }
        }
        for momentum in external_momenta {
            if !momentum.temporal.value.0.is_finite()
                || [
                    &momentum.spatial.px,
                    &momentum.spatial.py,
                    &momentum.spatial.pz,
                ]
                .iter()
                .any(|value| !value.0.is_finite())
            {
                return Err(eyre!("sampling fiber external momenta must be finite"));
            }
        }
        let rows = self
            .energies
            .iter()
            .map(|edge| {
                active
                    .iter()
                    .map(|&index| match lmb.edge_signatures[*edge].internal[index] {
                        SignOrZero::Minus => -1_i64,
                        SignOrZero::Zero => 0,
                        SignOrZero::Plus => 1,
                    })
                    .collect_vec()
            })
            .collect_vec();
        let rank = rank_i64(&rows);
        if rank != active.len() {
            return Err(eyre!(
                "sampling fiber {:?} has active routing rank {rank}, expected {}; spectator directions require a complement block",
                self.energies,
                active.len()
            ));
        }
        let varying = self
            .energies
            .iter()
            .zip(&rows)
            .filter_map(|(&edge, row)| row.iter().any(|&sign| sign != 0).then_some(edge))
            .collect_vec();
        let subspace_energies = subspace.contains(&self.energies, graph).collect_vec();
        if varying.iter().any(|edge| !subspace_energies.contains(edge)) {
            return Err(eyre!(
                "sampling fiber routing is inconsistent with its selected cycle subgraph"
            ));
        }
        for &edge in &self.energies {
            let mass = &masses[edge];
            if !mass.0.is_finite() || mass < &mass.zero() {
                return Err(eyre!(
                    "sampling fiber edge {edge} requires a finite nonnegative mass"
                ));
            }
        }
        let surface = Arc::new(self.clone());
        let lmbs = Arc::new(all_lmbs.clone());
        let graph = Arc::new(graph.clone());
        let subspace = Arc::new(subspace.clone());
        let masses = Arc::new(masses.clone());
        let externals = Arc::new(external_momenta.clone());
        let settings = Arc::new(settings.clone());
        let complement = Arc::new(complement.to_vec());
        let active = Arc::new(active);
        let zero = external_momenta[ExternalIndex(0)].temporal.value.zero();
        let dimension = active.len() * 3;
        let freeze_context = complement.is_empty();
        let n_loops = lmb.loop_edges.len();
        // One shared routed-ray evaluator serves both full and conditional
        // charts. Complement entries have zero radial velocity.
        let evaluator = {
            let (surface, lmbs, subspace, masses, externals, complement, active) = (
                surface.clone(),
                lmbs.clone(),
                subspace.clone(),
                masses.clone(),
                externals.clone(),
                complement.clone(),
                active.clone(),
            );
            Arc::new(
                move |direction: &[T], radius: T, center: &[T], context: &[T]| {
                    let zero = F(radius.zero());
                    let mut loops = LoopMomenta::from_iter(
                        (0..n_loops)
                            .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                    );
                    let mut velocity = loops.clone();
                    for (&index, components) in complement.iter().zip(context.chunks_exact(3)) {
                        loops[index] = ThreeMomentum::new(
                            F(components[0].clone()),
                            F(components[1].clone()),
                            F(components[2].clone()),
                        );
                    }
                    for ((&index, components), unit) in active
                        .iter()
                        .zip(center.chunks_exact(3))
                        .zip(direction.chunks_exact(3))
                    {
                        loops[index] = ThreeMomentum::new(
                            F(components[0].clone()),
                            F(components[1].clone()),
                            F(components[2].clone()),
                        );
                        velocity[index] = ThreeMomentum::new(
                            F(unit[0].clone()),
                            F(unit[1].clone()),
                            F(unit[2].clone()),
                        );
                    }
                    let (value, derivative) = surface.compute_self_and_r_derivative(
                        &F(radius),
                        &velocity,
                        &loops,
                        &externals,
                        &masses,
                        subspace.get_lmb(&lmbs),
                    );
                    Ok((value.0, derivative.0))
                },
            )
        };
        let preparer = Arc::new(move |context: &[T]| {
            let zero = externals[ExternalIndex(0)].temporal.value.zero();
            if context.len() != complement.len() * 3
                || context.iter().any(|value| !value.is_finite())
            {
                return Err(eyre!(
                    "sampling fiber expected {} finite complement components, got {}",
                    complement.len() * 3,
                    context.len()
                ));
            }
            let lmb = subspace.get_lmb(&lmbs);
            let spatial: ExternalThreeMomenta<F<T>> =
                externals.iter().map(|p| p.spatial.clone()).collect();
            let mut center = LoopMomenta::from_iter(
                (0..n_loops).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            for (&index, components) in complement.iter().zip(context.chunks_exact(3)) {
                center[index] = ThreeMomentum::new(
                    F(components[0].clone()),
                    F(components[1].clone()),
                    F(components[2].clone()),
                );
            }
            let zero_velocity = LoopMomenta::from_iter(
                (0..n_loops).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            let mut constant = surface.compute_shift_part_from_momenta(&externals, lmb);
            // Include the inputs to routed differences, not just their small
            // residuals. Large common spatial shifts may otherwise round an
            // actually nonzero separation to zero before minimum testing.
            let mut scale = center
                .iter()
                .flat_map(|p| [&p.px, &p.py, &p.pz])
                .chain(externals.iter().flat_map(|p| {
                    [
                        &p.temporal.value,
                        &p.spatial.px,
                        &p.spatial.py,
                        &p.spatial.pz,
                    ]
                }))
                .fold(constant.abs(), |sum, component| sum + component.abs());
            for &edge in &surface.energies {
                if !varying.contains(&edge) {
                    let momentum = lmb.edge_signatures[edge].compute_momentum(&center, &spatial);
                    let energy = (momentum.norm_squared() + masses[edge].square()).sqrt();
                    constant += &energy;
                    scale += energy;
                }
            }
            let minimum = if active.len() == 1 && varying.len() == 2 {
                // For q_i = s_i k + b_i, s_i = +/-1, the exact convex minimum
                // is sqrt((b1-s1*s2*b2)^2 + (m1+m2)^2). With zero total mass
                // the minimizing set is a segment; its midpoint is a valid
                // proposal center, not a claim that the pinch is isolated.
                let [first, second] = [varying[0], varying[1]];
                let index = active[0];
                let sign = |edge| match lmb.edge_signatures[edge].internal[index] {
                    SignOrZero::Plus => zero.one(),
                    SignOrZero::Minus => -zero.one(),
                    SignOrZero::Zero => unreachable!("two varying rank-one energies"),
                };
                let first_shift = lmb.edge_signatures[first].compute_momentum(&center, &spatial);
                let second_shift = lmb.edge_signatures[second].compute_momentum(&center, &spatial);
                let separation = &first_shift - &(second_shift * (sign(first) * sign(second)));
                let total_mass = &masses[first] + &masses[second];
                let fraction = if total_mass > zero {
                    &masses[first] / &total_mass
                } else {
                    zero.one() / zero.from_i64(2)
                };
                center[index] = (&separation * fraction - first_shift) * sign(first);
                let energy = (separation.norm_squared() + total_mass.square()).sqrt();
                scale += &energy;
                Some(energy + &constant)
            } else {
                None
            };
            let tolerance = zero.epsilon()
                * zero.from_usize(64 * (surface.energies.len() + 1))
                * scale.max(zero.one());
            let normalized = |value: &F<T>| -> Result<T> {
                let margin = -value / F::<T>::from_f64(settings.kinematics.e_cm);
                if !margin.0.is_finite() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "normalized sampling surface margin",
                        detail: "finite surface value and energy scale produce a nonfinite margin"
                            .to_owned(),
                    }
                    .into());
                }
                Ok(margin.0)
            };
            let status = if let Some(minimum) = minimum {
                if !minimum.0.is_finite() || minimum.abs() <= tolerance {
                    return Err(SamplingEvaluationError::UncertainGeometry { detail: format!("two-energy fiber minimum is not sign-certified: minimum={minimum}, tolerance={tolerance}, edges={:?}", surface.energies) }.into());
                }
                if minimum > zero {
                    PreparedSurfaceStatus::absent(format!(
                        "positive exact two-energy minimum {minimum}"
                    ))?
                } else {
                    PreparedSurfaceStatus::existing(Some(normalized(&minimum)?))?
                }
            } else {
                // The sum of masses is a necessary lower bound for every
                // routed topology. A positive bound certifies absence; a
                // failed approximate optimizer by itself never does.
                let lower = varying
                    .iter()
                    .fold(constant, |sum, &edge| sum + &masses[edge]);
                if lower > tolerance {
                    PreparedSurfaceStatus::absent(format!("positive energy lower bound {lower}"))?
                } else {
                    let value = surface
                        .compute_self_and_r_derivative(
                            &zero,
                            &zero_velocity,
                            &center,
                            &externals,
                            &masses,
                            lmb,
                        )
                        .0;
                    if !esurface_value_is_strictly_inside(
                        &value,
                        &F::<T>::from_f64(settings.kinematics.e_cm),
                    ) {
                        let thresholds = EsurfaceCollection::from(vec![surface.as_ref().clone()]);
                        let existing = ExistingThresholds::from(vec![EsurfaceID(0)]);
                        let input = OverlapInput {
                            graph: &graph,
                            settings: &settings,
                            subspace: &subspace,
                            threshold_subspaces: None,
                            lmbs: &lmbs,
                            thresholds: &thresholds,
                            edge_masses: masses.iter().map(|(_, mass)| mass.into_ff64()).collect(),
                            surface_kinematics: None,
                        };
                        let loops: LoopMomenta<F<f64>> =
                            center.iter().map(ThreeMomentum::to_f64).collect();
                        let external: ExternalFourMomenta<F<f64>> =
                            externals.iter().map(|p| p.to_f64()).collect();
                        if loops
                            .iter()
                            .flat_map(|p| [&p.px, &p.py, &p.pz])
                            .any(|v| !v.0.is_finite())
                            || external
                                .iter()
                                .flat_map(|p| {
                                    [
                                        &p.temporal.value,
                                        &p.spatial.px,
                                        &p.spatial.py,
                                        &p.spatial.pz,
                                    ]
                                })
                                .any(|v| !v.0.is_finite())
                            || input
                                .edge_masses
                                .iter()
                                .any(|(_, mass)| !mass.0.is_finite())
                        {
                            return Err(SamplingEvaluationError::UncertainGeometry {
                                detail: "native fiber data exceed the finite range of the approximate f64 SOCP center seed".to_owned(),
                            }.into());
                        }
                        let candidate = find_center(&input, &[ExistingEsurfaceId(0)], &existing, &loops, &external, false)
                            .map_err(|error| SamplingEvaluationError::UncertainGeometry { detail: format!("fiber center seed failed without an absence certificate: {error}") })?
                            .ok_or_else(|| SamplingEvaluationError::UncertainGeometry { detail: "fiber SOCP reported infeasibility without a native absence certificate".to_owned() })?;
                        for &index in active.iter() {
                            center[index] = ThreeMomentum::from_ff64(candidate[index]);
                        }
                    }
                    let value = surface
                        .compute_self_and_r_derivative(
                            &zero,
                            &zero_velocity,
                            &center,
                            &externals,
                            &masses,
                            lmb,
                        )
                        .0;
                    if !esurface_value_is_strictly_inside(
                        &value,
                        &F::<T>::from_f64(settings.kinematics.e_cm),
                    ) {
                        return Err(SamplingEvaluationError::UncertainGeometry {
                            detail: format!(
                                "fiber center is not natively certified inside: value={value}"
                            ),
                        }
                        .into());
                    }
                    PreparedSurfaceStatus::existing(Some(normalized(&value)?))?
                }
            };
            if status.is_existing() {
                // The analytic minimum does not excuse a rounded minimizing
                // point: the actual bound center must also be strictly inside.
                let value = surface
                    .compute_self_and_r_derivative(
                        &zero,
                        &zero_velocity,
                        &center,
                        &externals,
                        &masses,
                        lmb,
                    )
                    .0;
                if !esurface_value_is_strictly_inside(
                    &value,
                    &F::<T>::from_f64(settings.kinematics.e_cm),
                ) || value >= -&tolerance
                {
                    return Err(SamplingEvaluationError::UncertainGeometry {
                        detail: format!(
                            "prepared fiber center is not natively certified inside: value={value}, tolerance={tolerance}"
                        ),
                    }
                    .into());
                }
            }
            let center = active
                .iter()
                .flat_map(|&index| {
                    [
                        center[index].px.0.clone(),
                        center[index].py.0.clone(),
                        center[index].pz.0.clone(),
                    ]
                })
                .collect();
            Ok((center, status))
        });
        ImplicitSurfaceRadialMap::new(
            dimension,
            vec![zero.0; dimension],
            beta,
            power,
            Arc::new(|_, _| Err(eyre!("sampling fiber requires prepared complement data"))),
        )
        .and_then(|map| {
            let map = map
                .with_context_evaluator(evaluator)
                .with_context_preparer(preparer);
            if freeze_context {
                map.freeze_context(&[])
            } else {
                Ok(map)
            }
        })
    }

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
    pub(crate) fn classify_existence_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        real_mass_vector: &EdgeVec<F<T>>,
        reversed_edges: &[EdgeIndex],
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistence<T> {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            debug!("esurface has no external shift, cannot exist");
            return EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoExternalShift,
            };
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
            return EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoRadialDependence,
            };
        }

        let lmb = subspace.get_lmb(all_lmbs);
        let mass_sum: F<T> = subspace_energy_indices
            .iter()
            .map(|&index| &real_mass_vector[index])
            .fold(F::from_f64(0.0), |acc, x| acc + x);

        let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

        let graph_vector = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                compute_shift_part(external_signature, external_moms).spatial
                    * F::from_f64(*sign as f64)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

        let other_part = subspace
            .does_not_contain(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];
                let sign = if reversed_edges.contains(&index) {
                    -F::from_f64(1.0)
                } else {
                    F::from_f64(1.0)
                };

                signature.compute_momentum(
                    loop_moms,
                    &external_moms
                        .iter()
                        .map(|mom| mom.spatial.clone())
                        .collect::<TiVec<ExternalIndex, _>>(),
                ) * sign
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

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

        classification
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

    /// Shared routed ray evaluation for physical roots and full-space/proper-fiber charts.
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

                // At a massless endpoint E(r)=r|v| the radial right derivative is
                // |v|, whereas the two-sided formula q.v/E would evaluate 0/0.
                // Share this convention with sampling charts, using exact source
                // zeros: a computed E=0 can instead be numerical underflow.
                let zero = radius.zero();
                let derivative = if radius == &zero
                    && real_mass_vector[index] == zero
                    && momentum.px == zero
                    && momentum.py == zero
                    && momentum.pz == zero
                {
                    unit_loop_part.norm_squared().sqrt()
                } else {
                    momentum * &unit_loop_part / &energy
                };

                (derivative, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    /// Explicit prepared-ray boundary for LU root and jet evaluation. Routing
    /// before radial scaling changes finite-precision association. Existing
    /// amplitude/static/fiber evaluation keeps scale-then-route semantics.
    /// The center is supplied explicitly; physical LU uses generation zero.
    pub(crate) fn routed_ray<T: FloatLike>(
        &self,
        velocity: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> EsurfaceRay<T> {
        let spatial: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|momentum| momentum.spatial.clone())
            .collect();
        EsurfaceRay {
            energies: self
                .energies
                .iter()
                .map(|&edge| {
                    let signature = &lmb.edge_signatures[edge];
                    (
                        edge,
                        compute_loop_part(&signature.internal, velocity),
                        signature.compute_momentum(center, &spatial),
                        masses[edge].clone(),
                    )
                })
                .collect(),
            shift: self.compute_shift_part_from_momenta(external_moms, lmb),
        }
    }

    /// Solve the physical LU scaling about zero generation momentum. Sampling
    /// hosts use the same equation and policy; the caller owns the identity and
    /// precision history, so sharing this entry does not imply shared root data.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_lu_cut<T: FloatLike>(
        &self,
        loop_momenta: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> std::result::Result<(EsurfaceRay<T>, NewtonIterationResult<T>), SafeguardedNewtonError<T>>
    {
        let zero = loop_momenta[LoopIndex(0)].px.zero();
        let center = LoopMomenta::from_iter(
            (0..loop_momenta.len())
                .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        let ray = self.routed_ray(loop_momenta, &center, external_momenta, masses, lmb);
        let guess = Self::radius_guess_from_terms(
            &ray.shift,
            ray.energies
                .iter()
                .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
        );
        crate::debug_tags!(#integration, #cut, #solver;
            radial_root = %identity,
            initial_guess = %guess,
            residual_tolerance = %(e_cm * guess.epsilon()),
            "LU radial root setup"
        );
        let solution = diagnostics.solve(
            identity,
            &zero,
            &guess,
            |t| ray.evaluate(t),
            &guess.one(),
            2000,
            64,
            e_cm,
        )?;
        Ok((ray, solution))
    }

    // #[inline]
    /// the "loops_unit_in_subspace" means that the loop momenta that are part of the subspace are jointly normalized to unit length
    /// A ray with no radial dependence has no isolated root and returns zero guesses.
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
            // An energy with zero radial momentum is constant along this ray. It contributes
            // no large-radius growth or directional shift, even when its mass is nonzero.
            if norm_unit_loop_part_squared == const_builder.zero() {
                continue;
            }
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

    /// A ray with no radial dependence has no isolated root and returns zero guesses.
    pub(crate) fn get_radius_guess<T: FloatLike>(
        &self,
        unit_loops: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let shift = self.compute_shift_part_from_momenta(external_moms, lmb);
        //println!("got to energy loop");
        let terms = self.energies.iter().map(|&energy| {
            //println!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);
            let unit_loop_part = compute_loop_part(&signature.internal, unit_loops);
            //println!("computed_loop_part {:?}", unit_loop_part);
            let shift = compute_shift_part(&signature.external, external_moms).spatial;
            //./bprintln!("computed_shift {:?}", shift);
            (unit_loop_part.norm_squared(), &unit_loop_part * shift)
        });
        // Existing CT callers intentionally estimate with true external offsets,
        // even when their later root evaluation uses a nonzero overlap center.
        let guess = Self::radius_guess_from_terms(&shift, terms);
        let negative = guess.ref_neg();
        (guess, negative)
    }

    fn radius_guess_from_terms<T: FloatLike>(
        shift: &F<T>,
        terms: impl Iterator<Item = (F<T>, F<T>)>,
    ) -> F<T> {
        let zero = shift.zero();
        let mut guess = zero.clone();
        let mut denominator = zero.clone();
        // Each pair is (|v|^2, v.b); LU supplies its represented ray, while
        // non-LU callers preserve the existing external-only seed convention.
        for (norm_squared, dot_offset) in terms {
            // Constant energies have no directional contribution to the radius estimate.
            if norm_squared == zero {
                continue;
            }
            guess += dot_offset.abs() / &norm_squared;
            denominator += norm_squared.sqrt();
        }
        if denominator != zero {
            guess += shift.abs() / denominator;
        }
        guess
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
        cut: &CrossSectionCut,
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

    pub(crate) fn lmb_atom(&self, graph: &Graph, lmb_reps: &[Replacement]) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom();
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
    pub(crate) fn lmb_atom_simplified(&self, graph: &Graph, lmb_reps: &[Replacement]) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom();
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

#[derive(Debug, Clone, Encode, Decode)]
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
    use itertools::Itertools;
    use linnet::half_edge::HedgeGraph;
    use linnet::half_edge::builder::HedgeGraphBuilder;
    use linnet::half_edge::involution::{EdgeIndex, Flow, Orientation};
    use linnet::half_edge::nodestore::NodeStorageVec;
    use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::parse;

    use crate::cff::VertexSet;
    use crate::graph::{Graph, LmbIndex, LoopMomentumBasis, parse::from_dot::IntoGraph};
    use crate::initialisation::test_initialise;
    use crate::integrands::process::SamplingMapAffine;
    use crate::momentum::{
        FourMomentum, ThreeMomentum,
        sample::{ExternalFourMomenta, LoopMomenta, SubspaceData},
        signature::LoopExtSignature,
    };
    use crate::processes::CrossSectionCut;
    use crate::{
        cff::{esurface::Esurface, generation::ShiftRewrite},
        dot,
        utils::{
            ArbPrec, DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, ESURFACE_SHIFT_THRESHOLD, F, FloatLike,
            QuadFloat,
            newton_solver::{
                RadialRootDiagnostics, RadialRootIdentity, SafeguardedNewtonError,
                safeguarded_newton_iteration_and_derivative,
            },
            test_utils::dummy_hedge_graph,
        },
    };
    use typed_index_collections::ti_vec;

    use super::{EsurfaceExistence, add_external_shifts};

    #[test]
    fn joint_sampling_matches_distinct_reversed_edges_and_unequal_masses() {
        test_initialise().unwrap();
        // Splitting the kite's common line through E gives distinct edge IDs
        // with opposite complete affine routes, without altering its two cycles.
        let graph: Graph = dot!(digraph joint_serial_energy {
            ext_in [style=invis]
            ext_out [style=invis]
            node [num=1]
            edge [num=1 mass=1]
            ext_in -> A:0 [id=0 mass=0]
            C:1 -> ext_out [id=1 mass=0]
            A -> B [id=2 mass=2]
            B -> C [id=3]
            C -> D [id=4 lmb_id=0]
            D -> A [id=5 mass=3]
            B -> E [id=6 lmb_id=1]
            D -> E [id=7]
        })
        .unwrap();
        let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
        let id = LmbIndex::from(0);
        let lmb = &lmbs[id];
        let subspace = SubspaceData::new_from_parent_basis_edges(
            &[EdgeIndex(6)],
            &graph.full_filter(),
            id,
            &graph,
            &lmbs,
        )
        .unwrap();
        let active = subspace.iter_lmb_indices().next().unwrap();
        let prior = lmb
            .loop_edges
            .iter_enumerated()
            .find_map(|(index, edge)| (*edge == EdgeIndex(4)).then_some(index))
            .unwrap();
        let left = Esurface {
            energies: vec![EdgeIndex(2), EdgeIndex(4), EdgeIndex(6)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let right = Esurface {
            energies: vec![EdgeIndex(3), EdgeIndex(5), EdgeIndex(7)],
            ..left.clone()
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter([0., 0., 2., 1., 1., 3., 1., 1.].map(F))
            .unwrap();
        let externals = ExternalFourMomenta::from_iter(
            [FourMomentum::from_args(F(26_f64.sqrt()), F(0.), F(0.), F(1.)); 2],
        );
        let (prepare, common) = left
            .sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &graph,
                &masses,
                &externals,
                &[prior],
            )
            .unwrap();
        let geometry = prepare(&[1., 0., -0.5]).unwrap();
        assert_eq!(geometry.masses, [1., 2., 3.]);
        assert_eq!(geometry.shifts, [[1., 0., 0.5], [1., 0., -0.5]]);
        assert_ne!(
            lmb.edge_signatures[EdgeIndex(6)],
            lmb.edge_signatures[EdgeIndex(7)]
        );
        assert!(lmb.edges_are_raised(EdgeIndex(6), EdgeIndex(7)));
        assert_eq!(common, lmb.edge_signatures[EdgeIndex(6)]);
        let mut loops = LoopMomenta::from_iter([ThreeMomentum::new(F(0.), F(0.), F(0.)); 2]);
        loops[prior] = ThreeMomentum::new(F(1.), F(0.), F(-0.5));
        for components in [[0.5, -1., 1.], [-0.7, 0.3, -0.4]] {
            let x = ThreeMomentum::new(F(components[0]), F(components[1]), F(components[2]));
            loops[active] = x;
            let common_energy = (x.norm_squared() + F(1.)).sqrt();
            for (i, surface) in [&left, &right].into_iter().enumerate() {
                let [a, b, c] = geometry.shifts[i];
                let partner = x + ThreeMomentum::new(F(a), F(b), F(c));
                let expected = common_energy
                    + (partner.norm_squared() + F(geometry.masses[i + 1]).square()).sqrt()
                    - F(geometry.energy_sums[i]);
                assert!(
                    (surface.compute_from_momenta(lmb, &masses, &loops, &externals) - expected)
                        .abs()
                        < F(1e-12)
                );
            }
        }
        // Fixed energy multiplicities remain part of the original equation.
        let mut repeated = left.clone();
        repeated.energies.push(EdgeIndex(4));
        let (prepare_repeated, _) = repeated
            .sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &graph,
                &masses,
                &externals,
                &[prior],
            )
            .unwrap();
        assert!(
            (prepare_repeated(&[1., 0., -0.5]).unwrap().energy_sums[0]
                - (geometry.energy_sums[0] - 1.5))
                .abs()
                < 1e-12
        );

        let mut inconsistent_masses = masses.clone();
        inconsistent_masses[EdgeIndex(7)] = F(2.);
        assert!(
            left.sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &graph,
                &inconsistent_masses,
                &externals,
                &[prior],
            )
            .err()
            .unwrap()
            .to_string()
            .contains("inconsistent native values")
        );
        // Changing only the symbolic mass distinguishes formal identity from
        // accidentally equal numeric values in a supplied mass cache.
        let mut different_mass_graph = graph.clone();
        different_mass_graph.underlying[EdgeIndex(7)].particle =
            crate::graph::edge::PossibleParticle::JustMass { expr: parse!("2") };
        assert!(
            left.sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &different_mass_graph,
                &masses,
                &externals,
                &[prior],
            )
            .err()
            .unwrap()
            .to_string()
            .contains("found 0 candidates")
        );
    }

    #[test]
    fn radial_guesses_and_lu_roots_handle_constant_massive_energies() {
        test_initialise().unwrap();
        // A ttH cut with back-to-back tops and a Higgs at rest has the analytic root
        // sqrt(((Q-mH)/2)^2-mt^2)/|p|. The Higgs momentum vanishes by cancellation of
        // two nonzero basis momenta, not because its graph signature is constant.
        let graph: Graph = dot!(digraph massive_energy_at_rest {
            ext [style=invis]
            node [num=1]
            edge [num=1 mass=0]
            ext -> a:0 [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            a -> b [id=3]
            b:1 -> ext [id=4]
        })
        .unwrap();
        let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
        let lmb = &lmbs[LmbIndex::from(0)];
        let subspace = SubspaceData::new_from_parent_basis_edges(
            &[EdgeIndex(1), EdgeIndex(2)],
            &graph.underlying.full_filter(),
            LmbIndex::from(0),
            &graph,
            &lmbs,
        )
        .unwrap();
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter([F(0.0), F(173.0), F(173.0), F(125.0), F(0.0)])
            .unwrap();
        // Supply both external ports as future-directed momenta; the graph's incoming and
        // outgoing edge flows provide their relative sign in momentum conservation.
        let externals = ExternalFourMomenta::from_iter(
            [FourMomentum::from_args(F(1000.0), F(0.0), F(0.0), F(0.0)); 2],
        );
        let center = LoopMomenta::from_iter([ThreeMomentum::new(F(0.0), F(0.0), F(0.0)); 2]);

        for momentum in [0.0, 100.0] {
            let loops = LoopMomenta::from_iter([
                ThreeMomentum::new(F(momentum), F(0.0), F(0.0)),
                ThreeMomentum::new(F(-momentum), F(0.0), F(0.0)),
            ]);
            let guesses = [
                surface.get_radius_guess(&loops, &externals, lmb),
                surface.get_radius_guess_subspace(
                    &loops, &externals, &subspace, &lmbs, &graph, &masses,
                ),
            ];
            for (positive, negative) in guesses {
                assert!(positive.0.is_finite());
                assert_eq!(negative, -positive);
                let result = safeguarded_newton_iteration_and_derivative(
                    &F(0.0),
                    &positive,
                    |radius| {
                        surface.compute_self_and_r_derivative(
                            radius, &loops, &center, &externals, &masses, lmb,
                        )
                    },
                    &F(1.0),
                    2000,
                    64,
                    &F(1000.0),
                );
                if momentum == 0.0 {
                    assert_eq!(positive, F(0.0));
                    // The constant surface is negative everywhere: no finite radial root
                    // exists. Reject the ray instead of propagating a NaN or inventing a root.
                    assert!(matches!(
                        result,
                        Err(SafeguardedNewtonError::InvalidOutside { .. })
                    ));
                } else {
                    let result = result.unwrap();
                    let top_energy = (1000.0_f64 - 125.0) / 2.0;
                    let expected = (top_energy.powi(2) - 173.0_f64.powi(2)).sqrt() / momentum;
                    assert!((result.solution.0 - expected).abs() < 1.0e-14);
                    assert!(result.derivative_at_solution.0 > 0.0);
                    assert!(result.error_of_function.0.abs() < 1.0e-12);
                }
            }
        }

        fn check_native<T: FloatLike>(graph: &Graph, surface: &Esurface) {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let f = |value| one.from_i64(value);
            let lmb = &graph.loop_momentum_basis;
            let mut masses = graph
                .underlying
                .new_edgevec_from_iter([0, 173, 173, 125, 0].map(f))
                .unwrap();
            let e_cm = f(1000);
            let externals = ExternalFourMomenta::from_iter((0..2).map(|_| {
                FourMomentum::from_args(e_cm.clone(), zero.clone(), zero.clone(), zero.clone())
            }));
            let zero_vector = ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone());
            let center = LoopMomenta::from_iter([zero_vector.clone(), zero_vector.clone()]);
            let perturbation = &one / f(10).powi(25);
            let p = f(100) + &perturbation;
            let loops = LoopMomenta::from_iter([
                ThreeMomentum::new(p.clone(), zero.clone(), zero.clone()),
                ThreeMomentum::new(-&p, zero.clone(), zero.clone()),
            ]);
            let identity = RadialRootIdentity::new("native shared LU host".into());
            let mut diagnostics = RadialRootDiagnostics::default();
            let mut prior_policy = diagnostics.clone();
            let (ray, result) = surface
                .solve_lu_cut(
                    &loops,
                    &externals,
                    &masses,
                    lmb,
                    &e_cm,
                    &mut diagnostics,
                    &identity,
                )
                .unwrap();
            // The extracted entry must preserve the existing physical policy on
            // identical routed inputs and histories, including diagnostic data.
            let guess = surface.get_radius_guess(&loops, &externals, lmb).0;
            assert_eq!(
                Esurface::radius_guess_from_terms(
                    &ray.shift,
                    ray.energies
                        .iter()
                        .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
                ),
                guess
            );
            let prior = prior_policy
                .solve(
                    &identity,
                    &zero,
                    &guess,
                    |t| {
                        surface.compute_self_and_r_derivative(
                            t, &loops, &center, &externals, &masses, lmb,
                        )
                    },
                    &one,
                    2000,
                    64,
                    &e_cm,
                )
                .unwrap();
            assert_eq!(result.solution, prior.solution);
            assert_eq!(result.derivative_at_solution, prior.derivative_at_solution);
            assert_eq!(result.error_of_function, prior.error_of_function);
            assert_eq!(result.num_iterations_used, prior.num_iterations_used);
            let top_energy = (&e_cm - f(125)) / f(2);
            let radial_momentum = (top_energy.square() - f(173).square()).sqrt();
            assert!((&result.solution - &radial_momentum / &p).abs() < one.epsilon() * f(128));
            if perturbation > one.epsilon() * f(1000) {
                assert!(
                    (&result.solution - radial_momentum / f(100)).abs() > &one / f(10).powi(28)
                );
            }

            // A retained ray preserves energy multiplicity and a nonzero affine
            // center. Independently eta(t)=2*sqrt((3*t+1)^2+9)+5-e_cm.
            use crate::utils::hyperdual_utils::{
                extract_t_derivatives, new_constant, simple_n_deriv_shape,
            };
            use symbolica::domains::dual::HyperDual;
            let mut jet_surface = surface.clone();
            jet_surface.energies = vec![EdgeIndex(1), EdgeIndex(1), EdgeIndex(3)];
            let mut jet_masses = masses.clone();
            jet_masses[EdgeIndex(1)] = f(3);
            jet_masses[EdgeIndex(3)] = f(4);
            let jet_velocity = LoopMomenta::from_iter([
                ThreeMomentum::new(f(3), zero.clone(), zero.clone()),
                ThreeMomentum::new(-f(3), zero.clone(), zero.clone()),
            ]);
            let jet_center = LoopMomenta::from_iter([
                ThreeMomentum::new(one.clone(), zero.clone(), zero.clone()),
                ThreeMomentum::new(-&one, zero.clone(), zero.clone()),
            ]);
            let jet_externals =
                ExternalFourMomenta::from_iter((0..2).map(|_| {
                    FourMomentum::from_args(e_cm.clone(), zero.clone(), f(3), zero.clone())
                }));
            let jet_ray = jet_surface.routed_ray(
                &jet_velocity,
                &jet_center,
                &jet_externals,
                &jet_masses,
                lmb,
            );
            assert_eq!(
                jet_ray.energies.iter().map(|term| term.0).collect_vec(),
                jet_surface.energies
            );
            assert_eq!(jet_ray.energies[2].1.norm_squared(), zero);
            assert_eq!(jet_ray.energies[2].2.norm_squared(), f(9));
            // A CT seed still ignores the center, whereas the explicit LU ray
            // seed consumes its actual affine offsets. Here the difference is2/3.
            let old_seed = jet_surface
                .get_radius_guess(&jet_velocity, &jet_externals, lmb)
                .0;
            assert_eq!(old_seed, &e_cm / f(6));
            let ray_seed = Esurface::radius_guess_from_terms(
                &jet_ray.shift,
                jet_ray
                    .energies
                    .iter()
                    .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
            );
            assert!((ray_seed - old_seed - f(2) / f(3)).abs() < one.epsilon() * f(1024));
            let expected = [
                f(15) - &e_cm,
                f(24) / f(5),
                f(162) / f(125),
                -f(5832) / f(3125),
            ];
            let scalar = jet_ray.evaluate(&one);
            assert_eq!(scalar.0, expected[0]);
            assert!((&scalar.1 - &expected[1]).abs() < one.epsilon() * f(128));
            for order in 1..=3 {
                let t =
                    HyperDual::<F<T>>::new(simple_n_deriv_shape(order)).variable(0, one.clone());
                let actual = extract_t_derivatives(jet_ray.evaluate_dual(&t));
                // Independent legacy path: first scale complete loop jets, then
                // route original Esurface energies and true external ports.
                let dual_loops = LoopMomenta::from_iter(
                    jet_velocity.iter().zip(jet_center.iter()).map(|(v, b)| {
                        v.map_ref(&|x| new_constant(&t, x) * &t)
                            + b.map_ref(&|x| new_constant(&t, x))
                    }),
                );
                let dual_externals = jet_externals
                    .iter()
                    .map(|p| p.map_ref(&|x| new_constant(&t, x)))
                    .collect();
                let original = extract_t_derivatives(jet_surface.compute_from_dual_momenta(
                    lmb,
                    &jet_masses,
                    &dual_loops,
                    &dual_externals,
                ));
                for ((actual, original), expected) in actual.iter().zip(&original).zip(&expected) {
                    let tolerance = one.epsilon() * f(2048) * (one.clone() + expected.abs());
                    assert!((actual - expected).abs() < tolerance);
                    assert!((actual - original).abs() < tolerance);
                }
            }

            // Reassociation is deliberately confined to explicit LU preparation.
            // With exact binary64 t=0.1 and k=[10^16,-10^16+2], the old
            // scale-then-route momentum is 1/8; route-then-scale gives 2*t.
            // This is a negative equivalence control, not a tolerance failure.
            let mut cancellation = surface.clone();
            cancellation.energies = vec![EdgeIndex(3)];
            cancellation.external_shift.clear();
            let large = f(10).powi(16);
            let velocity = LoopMomenta::from_iter([
                ThreeMomentum::new(large.clone(), zero.clone(), zero.clone()),
                ThreeMomentum::new(-large + f(2), zero.clone(), zero.clone()),
            ]);
            let t = F(T::from_f64_exact_binary(0.1));
            for (offset, mass) in [(0, 0), (3, 0), (3, 4)] {
                let center = LoopMomenta::from_iter([
                    ThreeMomentum::new(f(offset), zero.clone(), zero.clone()),
                    zero_vector.clone(),
                ]);
                let mut case_masses = masses.clone();
                case_masses[EdgeIndex(3)] = f(mass);
                let represented =
                    cancellation.routed_ray(&velocity, &center, &externals, &case_masses, lmb);
                let prepared = represented.evaluate(&t);
                let original = cancellation.compute_self_and_r_derivative(
                    &t,
                    &velocity,
                    &center,
                    &externals,
                    &case_masses,
                    lmb,
                );
                let prepared_p = f(2) * &t + f(offset);
                let prepared_e = (prepared_p.square() + f(mass).square()).sqrt();
                let tolerance = one.epsilon() * f(128);
                assert!((&prepared.0 - &prepared_e).abs() < tolerance);
                assert!((&prepared.1 - f(2) * &prepared_p / &prepared_e).abs() < tolerance);
                if one.epsilon() > &one / f(10).powi(20) {
                    let original_p = &one / f(8) + f(offset);
                    let original_e = (original_p.square() + f(mass).square()).sqrt();
                    assert!((&original.0 - &original_e).abs() < tolerance);
                    assert!((&original.1 - f(2) * original_p / original_e).abs() < tolerance);
                    assert!((&prepared.0 - &original.0).abs() > &one / f(100));
                } else {
                    // Quad/Arb have enough mantissa for these exact represented
                    // products; they recover the same point by either ordering.
                    assert!((&original.0 - &prepared.0).abs() < tolerance);
                    assert!((&original.1 - &prepared.1).abs() < tolerance);
                }
            }

            // At the origin combine a massless moving energy with a shifted
            // massive energy E=sqrt(4^2+3^2). Their right derivative is exactly5.
            let mut endpoint = surface.clone();
            endpoint.energies = vec![EdgeIndex(1), EdgeIndex(2)];
            masses[EdgeIndex(1)] = zero.clone();
            masses[EdgeIndex(2)] = f(3);
            let velocity = LoopMomenta::from_iter([
                ThreeMomentum::new(f(3), f(4), zero.clone()),
                zero_vector.clone(),
            ]);
            let mut shifted_center = LoopMomenta::from_iter([
                zero_vector.clone(),
                ThreeMomentum::new(f(4), zero.clone(), zero.clone()),
            ]);
            let (value, derivative) = endpoint.compute_self_and_r_derivative(
                &zero,
                &velocity,
                &shifted_center,
                &externals,
                &masses,
                lmb,
            );
            assert_eq!(value, f(5) - &e_cm);
            assert_eq!(derivative, f(5));
            let prepared_endpoint = endpoint
                .routed_ray(&velocity, &shifted_center, &externals, &masses, lmb)
                .evaluate(&zero);
            assert_eq!(prepared_endpoint, (f(5) - &e_cm, f(5)));
            let stationary = LoopMomenta::from_iter([zero_vector.clone(), zero_vector]);
            assert_eq!(
                endpoint
                    .compute_self_and_r_derivative(
                        &zero,
                        &stationary,
                        &shifted_center,
                        &externals,
                        &masses,
                        lmb,
                    )
                    .1,
                zero
            );
            assert!(matches!(
                endpoint.solve_lu_cut(
                    &stationary,
                    &externals,
                    &masses,
                    lmb,
                    &e_cm,
                    &mut diagnostics,
                    &identity,
                ),
                Err(SafeguardedNewtonError::InvalidOutside { .. })
            ));

            // Squaring a nonzero source can underflow in Double/Quad. It must
            // never activate the exact massless/zero-momentum endpoint rule.
            let tiny = &one / f(10).powi(200);
            masses[EdgeIndex(1)] = tiny.clone();
            let derivative = endpoint
                .compute_self_and_r_derivative(
                    &zero,
                    &velocity,
                    &shifted_center,
                    &externals,
                    &masses,
                    lmb,
                )
                .1;
            let prepared_derivative = endpoint
                .routed_ray(&velocity, &shifted_center, &externals, &masses, lmb)
                .evaluate(&zero)
                .1;
            for derivative in [derivative, prepared_derivative] {
                if tiny.square() == zero {
                    assert!(derivative.is_nan() || derivative.is_infinite());
                } else {
                    assert_eq!(derivative, zero);
                }
            }
            masses[EdgeIndex(1)] = zero.clone();
            for (axis, expected) in [3, 4, 0].into_iter().enumerate() {
                shifted_center[crate::momentum::sample::LoopIndex(0)] = ThreeMomentum::new(
                    if axis == 0 {
                        tiny.clone()
                    } else {
                        zero.clone()
                    },
                    if axis == 1 {
                        tiny.clone()
                    } else {
                        zero.clone()
                    },
                    if axis == 2 {
                        tiny.clone()
                    } else {
                        zero.clone()
                    },
                );
                let derivative = endpoint
                    .compute_self_and_r_derivative(
                        &zero,
                        &velocity,
                        &shifted_center,
                        &externals,
                        &masses,
                        lmb,
                    )
                    .1;
                let prepared_derivative = endpoint
                    .routed_ray(&velocity, &shifted_center, &externals, &masses, lmb)
                    .evaluate(&zero)
                    .1;
                for derivative in [derivative, prepared_derivative] {
                    if tiny.square() == zero {
                        assert!(derivative.is_nan() || derivative.is_infinite());
                    } else {
                        assert!((derivative - f(expected)).abs() < one.epsilon() * f(128));
                    }
                }
            }
        }
        check_native::<f64>(&graph, &surface);
        check_native::<QuadFloat>(&graph, &surface);
        check_native::<ArbPrec>(&graph, &surface);
    }

    #[test]
    fn implicit_sampling_map_uses_graph_esurface_root_and_roundtrips() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph implicit_sampling_esurface {
            ext [style=invis]
            node [num=1]
            edge [num=1 mass=0]
            ext -> a:0 [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            a -> b [id=3]
            b:1 -> ext [id=4]
        })
        .unwrap();
        let lmb = graph.loop_momentum_basis.clone();
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter([F(0.0), F(173.0), F(173.0), F(125.0), F(0.0)])
            .unwrap();
        let externals = ExternalFourMomenta::from_iter(
            [FourMomentum::from_args(F(1000.0), F(0.0), F(0.0), F(0.0)); 2],
        );
        let map = surface
            .sampling_radial_map(&lmb, &masses, &externals, 400.0, 2.0)
            .unwrap();
        // Exercise both sides of the cut shell with a sizeable branch
        // split; the outer compactification carries its (1-split) factor.
        for radial_coordinate in [0.1, 0.9] {
            let coordinates = [radial_coordinate, 0.27, 0.61, 0.39, 0.72, 0.58];
            let forward = map.forward(&coordinates).unwrap();
            let step = 1.0e-5;
            let mut derivative_matrix = vec![vec![0.0; 6]; 6];
            for axis in 0..6 {
                let mut plus = coordinates;
                let mut minus = coordinates;
                plus[axis] += step;
                minus[axis] -= step;
                let plus = map.forward(&plus).unwrap();
                let minus = map.forward(&minus).unwrap();
                for component in 0..6 {
                    derivative_matrix[component][axis] =
                        (plus.point[component] - minus.point[component]) / (2.0 * step);
                }
            }
            let numerical_jacobian = SamplingMapAffine::new(derivative_matrix, vec![0.0; 6])
                .unwrap()
                .determinant();
            assert!((numerical_jacobian / forward.jacobian - 1.0).abs() < 1.0e-4);

            assert!(
                forward
                    .diagnostics
                    .iter()
                    .any(|diagnostic| diagnostic == "implicit_surface:regular_root")
            );
            assert!(forward.jacobian.is_finite() && forward.jacobian > 0.0);
            let inverse = map.inverse(&forward.point).unwrap();
            assert!(inverse.residual < 1.0e-9, "{}", inverse.residual);
            for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
                assert!((actual - expected).abs() < 1.0e-9);
            }
        }

        // The graph factory must bind source data in the evaluation precision.
        // These two energies are indistinguishable in f64 but define different
        // physical shells and therefore different forward maps at Quad precision.
        let one = F::<crate::utils::f128>::from_f64(1.0);
        let energy = one.from_i64(1000);
        let displacement = one.from_i64(10).powi(-20);
        assert_eq!(energy.into_f64(), (energy + displacement).into_f64());
        let native_masses = graph
            .underlying
            .new_edgevec_from_iter(
                masses
                    .iter()
                    .map(|(_, mass)| F::<crate::utils::f128>::from_ff64(*mass)),
            )
            .unwrap();
        let coordinates = [0.19, 0.27, 0.61, 0.39, 0.72, 0.58]
            .map(|value| F::<crate::utils::f128>::from_f64(value).0);
        let mut points = Vec::new();
        for energy in [energy, energy + displacement] {
            let externals = ExternalFourMomenta::from_iter(
                [FourMomentum::from_args(energy, one.zero(), one.zero(), one.zero()); 2],
            );
            let map = surface
                .sampling_radial_map(&lmb, &native_masses, &externals, 400.0, 2.0)
                .unwrap();
            let forward = map.forward(&coordinates).unwrap();
            let inverse = map.inverse(&forward.point).unwrap();
            for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
                assert!((F(*actual) - F(expected)).abs() < one.from_i64(10).powi(-25));
            }
            points.push(forward.point);
        }
        assert_ne!(points[0], points[1]);
        assert_eq!(
            points[0]
                .iter()
                .map(|value| F(*value).into_f64())
                .collect::<Vec<_>>(),
            points[1]
                .iter()
                .map(|value| F(*value).into_f64())
                .collect::<Vec<_>>(),
        );
    }

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
                .map(|(node_l, cut, node_r)| CrossSectionCut {
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
                .map(|(node_l, cut, node_r)| CrossSectionCut {
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
