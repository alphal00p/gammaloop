//! Native classification for graph-aware sampling maps.
//!
//! Graph, host cut, side and complete parent identity belong to the immutable
//! resolved block binding. Each conditional map prepares its actual native t*
//! and fixed physical coordinates together from the declared raw prerequisites.
//! The same preparation is repeated for each foreign inverse at its supplied
//! raw point. No partially sampled data are labelled a complete MomentumSample,
//! and precision rescue reconstructs the original cube rather than promoting a
//! previously mapped point or pairing it with stale cut data.

use crate::utils::FloatLike;
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};

/// Which amplitude side owns a prepared threshold map; a cut host has no side.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingCutSide {
    Left,
    Right,
}

/// Pointwise classification of an energy surface for prepared complement or cut data.
///
/// `Absent` is a normal branch for cross-section kinematics.  Callers should
/// use the full-support fallback map in that branch rather than dropping the
/// channel.  `Pinched` remains explicit because it has different diagnostics
/// and may require a dedicated local map in a later implementation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum PreparedSurfaceStatus<T: FloatLike = f64> {
    Existing { normalized_margin: Option<T> },
    Pinched { normalized_margin: T },
    Absent { reason: String },
}

impl<T: FloatLike> PreparedSurfaceStatus<T> {
    pub fn existing(normalized_margin: Option<T>) -> Result<Self> {
        if let Some(margin) = &normalized_margin {
            validate_finite(margin, "normalized surface margin")?;
        }
        Ok(Self::Existing { normalized_margin })
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

    /// A regular threshold radius is direction dependent and belongs to the
    /// radial map evaluation, not to this complement-only classification.
    pub fn is_existing(&self) -> bool {
        matches!(self, Self::Existing { .. })
    }

    pub fn uses_full_support_fallback(&self) -> bool {
        matches!(self, Self::Absent { .. } | Self::Pinched { .. })
    }
}

fn validate_finite<T: FloatLike>(value: &T, label: &str) -> Result<()> {
    if !value.is_finite() {
        return Err(eyre!("{label} must be finite, got {value}"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classified_sampling_surface_preserves_normalized_fallback_and_validates_margins() {
        let absent = PreparedSurfaceStatus::<f64>::absent("external shift is nonnegative").unwrap();
        assert!(!absent.is_existing());
        assert!(absent.uses_full_support_fallback());
        assert!(
            PreparedSurfaceStatus::pinched(0.0)
                .unwrap()
                .uses_full_support_fallback()
        );
        assert!(
            PreparedSurfaceStatus::existing(Some(0.1))
                .unwrap()
                .is_existing()
        );
        assert!(PreparedSurfaceStatus::existing(Some(f64::NAN)).is_err());
        assert!(PreparedSurfaceStatus::pinched(f64::INFINITY).is_err());
        assert!(PreparedSurfaceStatus::<f64>::absent(" ").is_err());
    }
}
