//! Prepared graph and cut data consumed by graph-aware sampling maps.
//!
//! A sampling map must not solve a cut lazily from a partially parameterized
//! point.  In particular, left and right threshold maps may be combined in one
//! channel only after the cut host has supplied its solved `t*` rescaling and
//! external data.  The small, serializable records in this module are the
//! boundary between that kinematic preparation and the graph-independent map
//! kernels in [`super::sampling_maps`].

use crate::{
    momentum::{
        FourMomentum, Rotatable, Rotation, ThreeMomentum,
        sample::{ExternalFourMomenta, LoopMomenta, MomentumSample},
    },
    utils::{F, FloatLike},
};
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};

use super::SamplingChannelId;

/// Which amplitude side owns a prepared threshold map.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingCutSide {
    Left,
    Right,
}

/// Exact result of a cut-dependent sampling-map evaluation.
///
/// This is the typed hand-off a cross-section host can use once it has solved
/// the cut LU root and constructed a conditional map.  The mapped momentum
/// sample, determinant and the channel's runtime context stay together so a
/// stability rotation or precision rescue cannot accidentally pair a map with
/// context from another canonical channel.  Construction is deliberately
/// independent of the physical map implementation; the host supplies the
/// already evaluated values.
#[derive(Clone, Debug)]
pub struct PreparedCrossSectionMapEvaluation<T: FloatLike> {
    channel_id: SamplingChannelId,
    mapped_sample: MomentumSample<T>,
    jacobian: F<T>,
    inverse_jacobian: F<T>,
    runtime_context: Vec<F<T>>,
    prepared_cut_context: PreparedCutSamplingContext<T>,
}

impl<T: FloatLike> PreparedCrossSectionMapEvaluation<T> {
    pub fn new(
        channel_id: SamplingChannelId,
        mapped_sample: MomentumSample<T>,
        jacobian: F<T>,
        inverse_jacobian: F<T>,
        runtime_context: Vec<F<T>>,
        prepared_cut_context: PreparedCutSamplingContext<T>,
    ) -> Result<Self> {
        prepared_cut_context.validate_loop_dimension(mapped_sample.loop_moms().0.len())?;
        for (name, value) in [
            ("map Jacobian", &jacobian),
            ("inverse map Jacobian", &inverse_jacobian),
        ] {
            if !value.0.is_finite() || value <= &value.zero() {
                return Err(eyre!("{name} must be finite and positive, got {value}"));
            }
        }
        if runtime_context.iter().any(|value| !value.0.is_finite()) {
            return Err(eyre!(
                "deferred map runtime context contains a non-finite value"
            ));
        }
        Ok(Self {
            channel_id,
            mapped_sample,
            jacobian,
            inverse_jacobian,
            runtime_context,
            prepared_cut_context,
        })
    }

    pub fn channel_id(&self) -> SamplingChannelId {
        self.channel_id
    }
    pub fn mapped_sample(&self) -> &MomentumSample<T> {
        &self.mapped_sample
    }
    pub fn jacobian(&self) -> &F<T> {
        &self.jacobian
    }
    pub fn inverse_jacobian(&self) -> &F<T> {
        &self.inverse_jacobian
    }
    pub fn runtime_context(&self) -> &[F<T>] {
        &self.runtime_context
    }
    pub fn prepared_cut_context(&self) -> &PreparedCutSamplingContext<T> {
        &self.prepared_cut_context
    }

    pub fn rotate(
        &self,
        rotation: &Rotation,
        loop_mom_cache_id: usize,
        external_mom_cache_id: usize,
    ) -> Self {
        Self {
            channel_id: self.channel_id,
            mapped_sample: self.mapped_sample.rotate(
                rotation,
                loop_mom_cache_id,
                external_mom_cache_id,
            ),
            jacobian: self.jacobian.clone(),
            inverse_jacobian: self.inverse_jacobian.clone(),
            runtime_context: self.runtime_context.clone(),
            prepared_cut_context: self.prepared_cut_context.rotated(rotation),
        }
    }

    pub fn cast_sample<T2: FloatLike>(&self) -> PreparedCrossSectionMapEvaluation<T2>
    where
        F<T2>: From<F<T>>,
    {
        PreparedCrossSectionMapEvaluation {
            channel_id: self.channel_id,
            mapped_sample: self.mapped_sample.cast_sample(),
            jacobian: self.jacobian.clone().into(),
            inverse_jacobian: self.inverse_jacobian.clone().into(),
            runtime_context: self
                .runtime_context
                .iter()
                .cloned()
                .map(F::<T2>::from)
                .collect(),
            prepared_cut_context: self
                .prepared_cut_context
                .cast_sample(|value| F::<T2>::from(F(value.clone())).0),
        }
    }
}

/// Pointwise classification of an energy surface for the prepared cut data.
///
/// `Absent` is a normal branch for cross-section kinematics.  Callers should
/// use the full-support fallback map in that branch rather than dropping the
/// channel.  `Pinched` remains explicit because it has different diagnostics
/// and may require a dedicated local map in a later implementation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum PreparedSurfaceStatus<T: FloatLike = f64> {
    Existing {
        threshold_radius: T,
        normalized_margin: Option<T>,
    },
    Pinched {
        normalized_margin: T,
    },
    Absent {
        reason: String,
    },
}

impl<T: FloatLike> PreparedSurfaceStatus<T> {
    pub fn existing(threshold_radius: T, normalized_margin: Option<T>) -> Result<Self> {
        validate_finite_non_negative(&threshold_radius, "threshold radius")?;
        if let Some(margin) = &normalized_margin {
            validate_finite(margin, "normalized surface margin")?;
        }
        Ok(Self::Existing {
            threshold_radius,
            normalized_margin,
        })
    }

    pub fn pinched(normalized_margin: T) -> Result<Self> {
        validate_finite(&normalized_margin, "normalized surface margin")?;
        Ok(Self::Pinched { normalized_margin })
    }

    pub fn absent(reason: impl Into<String>) -> Result<Self> {
        let reason = reason.into();
        if reason.trim().is_empty() {
            return Err(eyre!(
                "an absent sampling surface needs a diagnostic reason"
            ));
        }
        Ok(Self::Absent { reason })
    }

    /// Radius to use for a regular threshold shell, if one exists.
    pub fn threshold_radius(&self) -> Option<T> {
        match self {
            Self::Existing {
                threshold_radius, ..
            } => Some(threshold_radius.clone()),
            Self::Pinched { .. } | Self::Absent { .. } => None,
        }
    }

    pub fn is_existing(&self) -> bool {
        matches!(self, Self::Existing { .. })
    }

    pub fn uses_full_support_fallback(&self) -> bool {
        matches!(self, Self::Absent { .. })
    }
}

/// One surface prepared in a cut host's native loop-momentum frame.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PreparedSamplingSurface<T: FloatLike = f64> {
    /// Edges defining the physical E-surface.
    pub edge_ids: Vec<usize>,
    /// Edges spanning the coordinates on which this surface is solved.
    pub subspace_edges: Vec<usize>,
    pub status: PreparedSurfaceStatus<T>,
}

impl<T: FloatLike> PreparedSamplingSurface<T> {
    pub fn new(
        edge_ids: Vec<usize>,
        subspace_edges: Vec<usize>,
        status: PreparedSurfaceStatus<T>,
    ) -> Result<Self> {
        validate_edges(&edge_ids, "surface")?;
        validate_edges(&subspace_edges, "surface subspace")?;
        Ok(Self {
            edge_ids,
            subspace_edges,
            status,
        })
    }
}

/// Kinematic host prepared for one side of one Cutkosky cut.
///
/// `rescaling_t_star`, `loop_momenta`, and `external_momenta` must all come
/// from the same solved cut point.  Keeping them together prevents a channel
/// which combines thresholds from different cuts from accidentally reusing a
/// stale `t*` or a complement frame from another cut.  The parent LMB is
/// intentionally represented by its complete ordered edge list; no implicit
/// default is permitted at this boundary.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PreparedCutSamplingContext<T: FloatLike = f64> {
    pub graph_name: String,
    pub graph_id: usize,
    pub cut_id: usize,
    pub orientation: Option<usize>,
    pub side: SamplingCutSide,
    /// Complete parent LMB edge list used to interpret all surface data.
    pub parent_lmb: Vec<usize>,
    pub rescaling_t_star: T,
    pub loop_momenta: Vec<[T; 3]>,
    pub external_momenta: Vec<[T; 4]>,
    pub surfaces: Vec<PreparedSamplingSurface<T>>,
}

impl<T: FloatLike> PreparedCutSamplingContext<T> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph_name: impl Into<String>,
        graph_id: usize,
        cut_id: usize,
        orientation: Option<usize>,
        side: SamplingCutSide,
        parent_lmb: Vec<usize>,
        rescaling_t_star: T,
        loop_momenta: Vec<[T; 3]>,
        external_momenta: Vec<[T; 4]>,
        surfaces: Vec<PreparedSamplingSurface<T>>,
    ) -> Result<Self> {
        let graph_name = graph_name.into();
        if graph_name.trim().is_empty() {
            return Err(eyre!("prepared sampling context needs a graph name"));
        }
        validate_edges(&parent_lmb, "parent LMB")?;
        validate_finite(&rescaling_t_star, "cut rescaling t*")?;
        if rescaling_t_star <= rescaling_t_star.zero() {
            return Err(eyre!(
                "cut rescaling t* must be strictly positive, got {rescaling_t_star}"
            ));
        }
        for momentum in loop_momenta.iter().flatten() {
            validate_finite(momentum, "loop momentum")?;
        }
        for momentum in external_momenta.iter().flatten() {
            validate_finite(momentum, "external momentum")?;
        }
        for (index, surface) in surfaces.iter().enumerate() {
            if surface
                .subspace_edges
                .iter()
                .any(|edge| !parent_lmb.contains(edge))
            {
                return Err(eyre!(
                    "surface {index} subspace {:?} is not contained in parent LMB {:?}",
                    surface.subspace_edges,
                    parent_lmb
                ));
            }
        }
        Ok(Self {
            graph_name,
            graph_id,
            cut_id,
            orientation,
            side,
            parent_lmb,
            rescaling_t_star,
            loop_momenta,
            external_momenta,
            surfaces,
        })
    }

    /// Build the kinematic part of a prepared context directly from one
    /// solved LU sample. The input loop momenta are the unrescaled vectors
    /// used by that root solve; applying `t*` here keeps the resulting cut
    /// point and its rescaling together, without any process-global cut cache.
    #[allow(clippy::too_many_arguments)]
    pub fn from_lu_sample(
        graph_name: impl Into<String>,
        graph_id: usize,
        cut_id: usize,
        orientation: Option<usize>,
        side: SamplingCutSide,
        parent_lmb: Vec<usize>,
        rescaling_t_star: &F<T>,
        unrescaled_loop_momenta: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        surfaces: Vec<PreparedSamplingSurface<T>>,
    ) -> Result<Self> {
        let zero = rescaling_t_star.zero();
        if !rescaling_t_star.0.is_finite() || rescaling_t_star.0 <= zero.0 {
            return Err(eyre!(
                "a prepared LU sample requires a finite positive t*, got {rescaling_t_star}"
            ));
        }
        let n_loop_momenta = parent_lmb.len();
        let loop_momenta = unrescaled_loop_momenta
            .rescale(rescaling_t_star, None)
            .iter()
            .map(|momentum| {
                [
                    momentum.px.0.clone(),
                    momentum.py.0.clone(),
                    momentum.pz.0.clone(),
                ]
            })
            .collect();
        let external_momenta = external_momenta
            .iter()
            .map(|momentum: &FourMomentum<F<T>>| {
                [
                    momentum.temporal.value.0.clone(),
                    momentum.spatial.px.0.clone(),
                    momentum.spatial.py.0.clone(),
                    momentum.spatial.pz.0.clone(),
                ]
            })
            .collect();
        let prepared = Self::new(
            graph_name,
            graph_id,
            cut_id,
            orientation,
            side,
            parent_lmb,
            rescaling_t_star.0.clone(),
            loop_momenta,
            external_momenta,
            surfaces,
        )?;
        prepared.validate_loop_dimension(n_loop_momenta)?;
        Ok(prepared)
    }

    /// Rotate the kinematic payload together with a mapped momentum sample.
    /// Surface identities and the parent frame stay unchanged; the numerical
    /// loop and external vectors must remain in the same frame as the sample
    /// after stability rotations.
    pub(crate) fn rotated(&self, rotation: &Rotation) -> Self {
        let rotate_three = |momentum: &[T; 3]| {
            let rotated = ThreeMomentum::new(
                F(momentum[0].clone()),
                F(momentum[1].clone()),
                F(momentum[2].clone()),
            )
            .rotate(rotation);
            [rotated.px.0, rotated.py.0, rotated.pz.0]
        };
        let rotate_four = |momentum: &[T; 4]| {
            let [px, py, pz] = rotate_three(&[
                momentum[1].clone(),
                momentum[2].clone(),
                momentum[3].clone(),
            ]);
            [momentum[0].clone(), px, py, pz]
        };
        Self {
            loop_momenta: self.loop_momenta.iter().map(rotate_three).collect(),
            external_momenta: self.external_momenta.iter().map(rotate_four).collect(),
            ..self.clone()
        }
    }

    /// Convert the whole prepared payload without a lower-precision intermediate.
    /// Precision rescue must still rebuild it from the original source point;
    /// this conversion is only for an already prepared sample's representation.
    pub fn cast_sample<T2: FloatLike>(
        &self,
        cast: impl Fn(&T) -> T2 + Copy,
    ) -> PreparedCutSamplingContext<T2> {
        PreparedCutSamplingContext {
            graph_name: self.graph_name.clone(),
            graph_id: self.graph_id,
            cut_id: self.cut_id,
            orientation: self.orientation,
            side: self.side,
            parent_lmb: self.parent_lmb.clone(),
            rescaling_t_star: cast(&self.rescaling_t_star),
            loop_momenta: self
                .loop_momenta
                .iter()
                .map(|p| p.each_ref().map(cast))
                .collect(),
            external_momenta: self
                .external_momenta
                .iter()
                .map(|p| p.each_ref().map(cast))
                .collect(),
            surfaces: self
                .surfaces
                .iter()
                .map(|surface| PreparedSamplingSurface {
                    edge_ids: surface.edge_ids.clone(),
                    subspace_edges: surface.subspace_edges.clone(),
                    status: match &surface.status {
                        PreparedSurfaceStatus::Existing {
                            threshold_radius,
                            normalized_margin,
                        } => PreparedSurfaceStatus::Existing {
                            threshold_radius: cast(threshold_radius),
                            normalized_margin: normalized_margin.as_ref().map(cast),
                        },
                        PreparedSurfaceStatus::Pinched { normalized_margin } => {
                            PreparedSurfaceStatus::Pinched {
                                normalized_margin: cast(normalized_margin),
                            }
                        }
                        PreparedSurfaceStatus::Absent { reason } => PreparedSurfaceStatus::Absent {
                            reason: reason.clone(),
                        },
                    },
                })
                .collect(),
        }
    }

    /// Check that this context is suitable for a parent frame with the given
    /// loop dimension.  This is intentionally a separate check from
    /// construction: a serialized context can be valid in isolation while
    /// still belonging to a different graph or sample frame.
    pub fn validate_loop_dimension(&self, n_loop_momenta: usize) -> Result<()> {
        if self.parent_lmb.len() != n_loop_momenta {
            return Err(eyre!(
                "prepared parent LMB has {} edges, expected {} loop momenta",
                self.parent_lmb.len(),
                n_loop_momenta
            ));
        }
        if self.loop_momenta.len() != n_loop_momenta {
            return Err(eyre!(
                "prepared LU sample has {} loop momenta, expected {}",
                self.loop_momenta.len(),
                n_loop_momenta
            ));
        }
        Ok(())
    }

    pub fn surface(&self, index: usize) -> Option<&PreparedSamplingSurface<T>> {
        self.surfaces.get(index)
    }

    pub fn regular_surface_radii(&self) -> impl Iterator<Item = T> + '_ {
        self.surfaces
            .iter()
            .filter_map(|surface| surface.status.threshold_radius())
    }

    pub fn has_absent_surface(&self) -> bool {
        self.surfaces
            .iter()
            .any(|surface| surface.status.uses_full_support_fallback())
    }
}

fn validate_edges(edges: &[usize], label: &str) -> Result<()> {
    if edges.is_empty() {
        return Err(eyre!("{label} edge list must not be empty"));
    }
    let mut unique = std::collections::HashSet::with_capacity(edges.len());
    if edges.iter().any(|edge| !unique.insert(*edge)) {
        return Err(eyre!("{label} edge list contains duplicates: {edges:?}"));
    }
    Ok(())
}

fn validate_finite<T: FloatLike>(value: &T, label: &str) -> Result<()> {
    if !value.is_finite() {
        return Err(eyre!("{label} must be finite, got {value}"));
    }
    Ok(())
}

fn validate_finite_non_negative<T: FloatLike>(value: &T, label: &str) -> Result<()> {
    validate_finite(value, label)?;
    if value < &value.zero() {
        return Err(eyre!("{label} must be non-negative, got {value}"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DependentMomentaConstructor, settings::runtime::kinematic::Externals};

    fn surface(status: PreparedSurfaceStatus) -> PreparedSamplingSurface {
        PreparedSamplingSurface::new(vec![3, 7], vec![3], status).unwrap()
    }

    fn context(side: SamplingCutSide) -> PreparedCutSamplingContext {
        PreparedCutSamplingContext::new(
            "GL638",
            12,
            4,
            Some(2),
            side,
            vec![3, 7, 11],
            0.83,
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            vec![[10.0, 0.0, 0.0, 10.0]],
            vec![surface(
                PreparedSurfaceStatus::existing(2.5, Some(0.1)).unwrap(),
            )],
        )
        .unwrap()
    }

    #[test]
    fn prepared_context_keeps_cut_host_identity_and_t_star() {
        let left = context(SamplingCutSide::Left);
        let right = context(SamplingCutSide::Right);
        assert_eq!(left.parent_lmb, vec![3, 7, 11]);
        assert_eq!(left.rescaling_t_star, right.rescaling_t_star);
        assert_eq!(left.external_momenta, right.external_momenta);
        assert_ne!(left.side, right.side);
        assert_eq!(
            left.surface(0).unwrap().status.threshold_radius(),
            Some(2.5)
        );
    }

    #[test]
    fn absent_surface_requests_full_support_fallback() {
        let absent = PreparedSurfaceStatus::absent("external shift is not negative").unwrap();
        assert!(!absent.is_existing());
        assert!(absent.uses_full_support_fallback());
        let context = PreparedCutSamplingContext::new(
            "g",
            0,
            0,
            None,
            SamplingCutSide::Left,
            vec![0],
            1.0,
            vec![[0.0, 0.0, 0.0]],
            vec![],
            vec![PreparedSamplingSurface::new(vec![3, 7], vec![0], absent).unwrap()],
        )
        .unwrap();
        assert!(context.has_absent_surface());
        assert_eq!(context.regular_surface_radii().count(), 0);
    }

    #[test]
    fn parent_lmb_and_surface_data_are_validated() {
        assert!(PreparedSurfaceStatus::existing(-1.0, None).is_err());
        assert!(
            PreparedSamplingSurface::new(
                vec![1, 2, 1],
                vec![1],
                PreparedSurfaceStatus::<f64>::absent("missing").unwrap()
            )
            .is_err()
        );
        assert!(
            PreparedCutSamplingContext::new(
                "g",
                0,
                0,
                None,
                SamplingCutSide::Left,
                vec![],
                1.0,
                vec![],
                vec![],
                vec![],
            )
            .is_err()
        );
    }

    #[test]
    fn lu_sample_constructor_rejects_nonpositive_root_and_checks_frame_dimension() {
        let loop_momenta =
            LoopMomenta::from_iter([crate::momentum::ThreeMomentum::new(F(1.0), F(2.0), F(3.0))]);
        let external_momenta = ExternalFourMomenta::from_iter([FourMomentum::new(
            crate::momentum::Energy::new(F(10.0)),
            crate::momentum::ThreeMomentum::new(F(0.0), F(0.0), F(10.0)),
        )]);
        for t_star in [0.0, -0.5, f64::INFINITY, f64::NAN] {
            assert!(
                PreparedCutSamplingContext::from_lu_sample(
                    "g",
                    0,
                    0,
                    Some(0),
                    SamplingCutSide::Left,
                    vec![0],
                    &F(t_star),
                    &loop_momenta,
                    &external_momenta,
                    vec![],
                )
                .is_err()
            );
        }

        let context = PreparedCutSamplingContext::from_lu_sample(
            "g",
            0,
            0,
            Some(0),
            SamplingCutSide::Left,
            vec![0],
            &F(0.5),
            &loop_momenta,
            &external_momenta,
            vec![],
        )
        .unwrap();
        context.validate_loop_dimension(1).unwrap();
        assert!(context.validate_loop_dimension(2).is_err());
        assert_eq!(context.rescaling_t_star, 0.5);
        assert_eq!(context.loop_momenta, vec![[0.5, 1.0, 1.5]]);
        assert_eq!(context.external_momenta, vec![[10.0, 0.0, 0.0, 10.0]]);
        let wrong_dimension = PreparedCutSamplingContext::from_lu_sample(
            "g",
            0,
            0,
            Some(0),
            SamplingCutSide::Left,
            vec![0, 1],
            &F(0.5),
            &loop_momenta,
            &external_momenta,
            vec![],
        )
        .unwrap_err();
        assert!(
            wrong_dimension
                .to_string()
                .contains("1 loop momenta, expected 2")
        );
    }

    #[test]
    fn prepared_lu_context_preserves_native_precision_through_rotation_and_cast() {
        use crate::utils::{PrecisionUpgradable, f128};

        let one = F::<f128>::from_f64(1.0);
        let t_star = one + one.from_i64(2).powi(-80);
        assert_eq!(t_star.into_f64(), 1.0);
        assert_ne!(t_star, one);
        let loops = LoopMomenta::from_iter([ThreeMomentum::new(t_star, one, -t_star)]);
        let externals = ExternalFourMomenta::from_iter([FourMomentum::from_args(
            t_star,
            one.zero(),
            t_star,
            one.zero(),
        )]);
        let prepared = PreparedCutSamplingContext::from_lu_sample(
            "native",
            2,
            3,
            Some(4),
            SamplingCutSide::Right,
            vec![7],
            &t_star,
            &loops,
            &externals,
            vec![
                PreparedSamplingSurface::new(
                    vec![2, 7],
                    vec![7],
                    PreparedSurfaceStatus::existing(t_star.0, Some((-t_star).0)).unwrap(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(prepared.rescaling_t_star, t_star.0);
        assert_eq!(prepared.loop_momenta[0][0], t_star.square().0);
        assert_eq!(prepared.external_momenta[0][0], t_star.0);
        assert_eq!(
            prepared.regular_surface_radii().collect::<Vec<_>>(),
            vec![t_star.0]
        );
        let restored: PreparedCutSamplingContext<f128> =
            serde_json::from_str(&serde_json::to_string(&prepared).unwrap()).unwrap();
        assert_eq!(restored, prepared);
        let promoted = prepared.cast_sample(|value| value.higher());
        assert_eq!(promoted.cast_sample(|value| value.lower()), prepared);

        let rotation = Rotation::new(crate::momentum::RotationMethod::Pi2X);
        let rotated = promoted.rotated(&rotation);
        assert_eq!(rotated.rescaling_t_star, promoted.rescaling_t_star);
        assert_eq!(
            rotated.loop_momenta[0][1],
            (-F(promoted.loop_momenta[0][2].clone())).0
        );
        assert_eq!(rotated.loop_momenta[0][2], promoted.loop_momenta[0][1]);
        assert_eq!(
            rotated.external_momenta[0][3],
            promoted.external_momenta[0][2]
        );
        assert_eq!(rotated.surfaces, promoted.surfaces);
    }

    #[test]
    fn prepared_lu_context_accepts_finite_native_values_outside_f64_range() {
        let one = F::<crate::utils::ArbPrec>::from_f64(1.0);
        let large = one.from_i64(10).powi(400);
        let small = &one / &large;
        assert!(large.into_f64().is_infinite());
        assert_eq!(small.into_f64(), 0.0);
        let loops =
            LoopMomenta::from_iter([ThreeMomentum::new(large.clone(), one.zero(), one.zero())]);
        let externals = ExternalFourMomenta::from_iter([FourMomentum::from_args(
            large.clone(),
            one.zero(),
            one.zero(),
            large.clone(),
        )]);
        let prepared = PreparedCutSamplingContext::from_lu_sample(
            "native",
            0,
            0,
            None,
            SamplingCutSide::Left,
            vec![0],
            &small,
            &loops,
            &externals,
            vec![
                PreparedSamplingSurface::new(
                    vec![0, 1],
                    vec![0],
                    PreparedSurfaceStatus::existing(large.0.clone(), Some(small.0.clone()))
                        .unwrap(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(prepared.rescaling_t_star, small.0);
        assert_eq!(prepared.loop_momenta[0][0], (&large * &small).0);
        assert_eq!(prepared.external_momenta[0][0], large.0);
        assert_eq!(
            prepared.surface(0).unwrap().status.threshold_radius(),
            Some(large.0)
        );
    }

    #[test]
    fn prepared_map_evaluation_keeps_channel_context_and_rotates_together() {
        let sample = MomentumSample::new(
            LoopMomenta::from_iter([
                crate::momentum::ThreeMomentum::new(F(1.0), F(2.0), F(3.0)),
                crate::momentum::ThreeMomentum::new(F(4.0), F(5.0), F(6.0)),
                crate::momentum::ThreeMomentum::new(F(7.0), F(8.0), F(9.0)),
            ]),
            0,
            &Externals::default(),
            0,
            F(1.0),
            DependentMomentaConstructor::CrossSection,
            Some(0),
        )
        .unwrap();
        let evaluation = PreparedCrossSectionMapEvaluation::new(
            SamplingChannelId::from(7),
            sample,
            F(2.0),
            F(0.5),
            vec![F(0.25), F(0.75)],
            context(SamplingCutSide::Left),
        )
        .unwrap();
        assert_eq!(evaluation.channel_id(), SamplingChannelId::from(7));
        assert_eq!(evaluation.runtime_context(), &[F(0.25), F(0.75)]);

        let rotated = evaluation.rotate(
            &Rotation::new(crate::momentum::RotationMethod::Pi2X),
            11,
            13,
        );
        assert_eq!(rotated.channel_id(), evaluation.channel_id());
        assert_eq!(
            rotated.prepared_cut_context().loop_momenta[0],
            [1.0, -3.0, 2.0]
        );
        assert_eq!(
            rotated.prepared_cut_context().external_momenta[0],
            [10.0, 0.0, -10.0, 0.0]
        );
        assert_eq!(rotated.mapped_sample().sample.loop_mom_cache_id, 11);
        assert_eq!(rotated.mapped_sample().sample.external_mom_cache_id, 13);
    }

    #[test]
    fn prepared_map_evaluation_rejects_invalid_jacobians_and_context() {
        let sample = MomentumSample::new(
            LoopMomenta::from_iter(vec![ThreeMomentum::new(F(0.0), F(0.0), F(0.0)); 3]),
            0,
            &Externals::default(),
            0,
            F(1.0),
            DependentMomentaConstructor::CrossSection,
            None,
        )
        .unwrap();
        for (jacobian, inverse_jacobian) in [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0)] {
            let error = PreparedCrossSectionMapEvaluation::new(
                SamplingChannelId::from(0),
                sample.clone(),
                F(jacobian),
                F(inverse_jacobian),
                vec![],
                context(SamplingCutSide::Right),
            )
            .unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains("Jacobian must be finite and positive")
            );
        }
        let error = PreparedCrossSectionMapEvaluation::new(
            SamplingChannelId::from(0),
            sample.clone(),
            F(1.0),
            F(1.0),
            vec![F(f64::NAN)],
            context(SamplingCutSide::Right),
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("runtime context contains a non-finite value")
        );

        // A finite value at rescue precision need not fit in f64. Validation
        // must not classify an underflowed or overflowed conversion as invalid.
        let one = F::<crate::utils::ArbPrec>::default().one();
        let large = one.from_usize(10).powi(400);
        let small = &one / &large;
        assert!(large.into_f64().is_infinite());
        assert_eq!(small.into_f64(), 0.0);
        PreparedCrossSectionMapEvaluation::new(
            SamplingChannelId::from(0),
            MomentumSample::new(
                LoopMomenta::from_iter(vec![
                    ThreeMomentum::new(one.zero(), one.zero(), one.zero());
                    3
                ]),
                0,
                &Externals::default(),
                0,
                one,
                DependentMomentaConstructor::CrossSection,
                None,
            )
            .unwrap(),
            large.clone(),
            small,
            vec![large],
            context(SamplingCutSide::Right)
                .cast_sample(|value| crate::utils::ArbPrec::from_f64(*value)),
        )
        .unwrap();
    }
}
