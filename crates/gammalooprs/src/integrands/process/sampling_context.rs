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
        FourMomentum,
        sample::{ExternalFourMomenta, LoopMomenta},
    },
    utils::{F, FloatLike},
};
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};

/// Which amplitude side owns a prepared threshold map.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingCutSide {
    Left,
    Right,
}

/// Pointwise classification of an energy surface for the prepared cut data.
///
/// `Absent` is a normal branch for cross-section kinematics.  Callers should
/// use the full-support fallback map in that branch rather than dropping the
/// channel.  `Pinched` remains explicit because it has different diagnostics
/// and may require a dedicated local map in a later implementation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum PreparedSurfaceStatus {
    Existing {
        threshold_radius: f64,
        normalized_margin: Option<f64>,
    },
    Pinched {
        normalized_margin: f64,
    },
    Absent {
        reason: String,
    },
}

impl PreparedSurfaceStatus {
    pub fn existing(threshold_radius: f64, normalized_margin: Option<f64>) -> Result<Self> {
        validate_finite_non_negative(threshold_radius, "threshold radius")?;
        if let Some(margin) = normalized_margin {
            validate_finite(margin, "normalized surface margin")?;
        }
        Ok(Self::Existing {
            threshold_radius,
            normalized_margin,
        })
    }

    pub fn pinched(normalized_margin: f64) -> Result<Self> {
        validate_finite(normalized_margin, "normalized surface margin")?;
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
    pub fn threshold_radius(&self) -> Option<f64> {
        match self {
            Self::Existing {
                threshold_radius, ..
            } => Some(*threshold_radius),
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
pub struct PreparedSamplingSurface {
    /// Edges defining the physical E-surface.
    pub edge_ids: Vec<usize>,
    /// Edges spanning the coordinates on which this surface is solved.
    pub subspace_edges: Vec<usize>,
    pub status: PreparedSurfaceStatus,
}

impl PreparedSamplingSurface {
    pub fn new(
        edge_ids: Vec<usize>,
        subspace_edges: Vec<usize>,
        status: PreparedSurfaceStatus,
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
pub struct PreparedCutSamplingContext {
    pub graph_name: String,
    pub graph_id: usize,
    pub cut_id: usize,
    pub orientation: Option<usize>,
    pub side: SamplingCutSide,
    /// Complete parent LMB edge list used to interpret all surface data.
    pub parent_lmb: Vec<usize>,
    pub rescaling_t_star: f64,
    pub loop_momenta: Vec<[f64; 3]>,
    pub external_momenta: Vec<[f64; 4]>,
    pub surfaces: Vec<PreparedSamplingSurface>,
}

impl PreparedCutSamplingContext {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph_name: impl Into<String>,
        graph_id: usize,
        cut_id: usize,
        orientation: Option<usize>,
        side: SamplingCutSide,
        parent_lmb: Vec<usize>,
        rescaling_t_star: f64,
        loop_momenta: Vec<[f64; 3]>,
        external_momenta: Vec<[f64; 4]>,
        surfaces: Vec<PreparedSamplingSurface>,
    ) -> Result<Self> {
        let graph_name = graph_name.into();
        if graph_name.trim().is_empty() {
            return Err(eyre!("prepared sampling context needs a graph name"));
        }
        validate_edges(&parent_lmb, "parent LMB")?;
        validate_finite(rescaling_t_star, "cut rescaling t*")?;
        for momentum in loop_momenta.iter().flatten() {
            validate_finite(*momentum, "loop momentum")?;
        }
        for momentum in external_momenta.iter().flatten() {
            validate_finite(*momentum, "external momentum")?;
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
    pub fn from_lu_sample<T: FloatLike>(
        graph_name: impl Into<String>,
        graph_id: usize,
        cut_id: usize,
        orientation: Option<usize>,
        side: SamplingCutSide,
        parent_lmb: Vec<usize>,
        rescaling_t_star: &F<T>,
        unrescaled_loop_momenta: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        surfaces: Vec<PreparedSamplingSurface>,
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
                    momentum.px.into_f64(),
                    momentum.py.into_f64(),
                    momentum.pz.into_f64(),
                ]
            })
            .collect();
        let external_momenta = external_momenta
            .iter()
            .map(|momentum: &FourMomentum<F<T>>| {
                [
                    momentum.temporal.value.into_f64(),
                    momentum.spatial.px.into_f64(),
                    momentum.spatial.py.into_f64(),
                    momentum.spatial.pz.into_f64(),
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
            rescaling_t_star.into_f64(),
            loop_momenta,
            external_momenta,
            surfaces,
        )?;
        prepared.validate_loop_dimension(n_loop_momenta)?;
        Ok(prepared)
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

    pub fn surface(&self, index: usize) -> Option<&PreparedSamplingSurface> {
        self.surfaces.get(index)
    }

    pub fn regular_surface_radii(&self) -> impl Iterator<Item = f64> + '_ {
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

fn validate_finite(value: f64, label: &str) -> Result<()> {
    if !value.is_finite() {
        return Err(eyre!("{label} must be finite, got {value}"));
    }
    Ok(())
}

fn validate_finite_non_negative(value: f64, label: &str) -> Result<()> {
    validate_finite(value, label)?;
    if value < 0.0 {
        return Err(eyre!("{label} must be non-negative, got {value}"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
            vec![[1.0, 2.0, 3.0]],
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
                PreparedSurfaceStatus::absent("missing").unwrap()
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
}
