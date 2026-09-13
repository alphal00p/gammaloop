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
use std::collections::BTreeMap;

use super::{sampling_maps::SamplingEvaluationError, sampling_selection::SamplingChannelId};

/// Only the discrete proposal recipe survives native retries. Geometry, roots
/// and the resulting radius are rebuilt from the original draw in every lane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SamplingProposalDecision {
    Compact { dyadic_exponent: u8 },
    Ordinary,
}

/// Immutable compiled identities distinguish generating rows and foreign maps.
/// Neither the generator nor this storage path is an input to the proposal law.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct SamplingProposalKey {
    pub(crate) graph_id: usize,
    pub(crate) generating_channel: SamplingChannelId,
    pub(crate) target_channel: SamplingChannelId,
    pub(crate) block_path: Vec<usize>,
}

/// The original-draw metadata owns these records. Default is sealed and empty:
/// a production decision cannot be invented by an incidental native attempt.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SamplingProposalPolicies {
    records: BTreeMap<SamplingProposalKey, SamplingProposalDecision>,
    collecting: bool,
}

impl SamplingProposalPolicies {
    pub(crate) fn begin_collection(&mut self) {
        self.records.clear();
        self.collecting = true;
    }

    pub(crate) fn seal(&mut self) {
        self.collecting = false;
    }

    pub(crate) fn is_collecting(&self) -> bool {
        self.collecting
    }

    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.records.len()
    }

    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Borrowed view of a channel's current numerical prerequisites and the sole
/// original-draw policy owner. Detached access is explicit for represented-
/// geometry fixtures; production uses collecting or sealed metadata records.
#[derive(Debug)]
pub struct SamplingMapContext<'a, T: FloatLike = f64> {
    pub previous: &'a [T],
    pub(crate) policies: Option<&'a mut SamplingProposalPolicies>,
    pub(crate) key: SamplingProposalKey,
}

impl<'a, T: FloatLike> SamplingMapContext<'a, T> {
    pub fn detached(previous: &'a [T]) -> Self {
        Self {
            previous,
            policies: None,
            key: SamplingProposalKey {
                graph_id: 0,
                generating_channel: SamplingChannelId(0),
                target_channel: SamplingChannelId(0),
                block_path: Vec::new(),
            },
        }
    }

    /// Wrappers keep lineage while replacing prerequisites. A real compiled
    /// child appends its immutable index, identically in forward and inverse.
    pub(crate) fn reborrow<'b>(
        &'b mut self,
        previous: &'b [T],
        child: Option<usize>,
    ) -> SamplingMapContext<'b, T> {
        let mut key = self.key.clone();
        key.block_path.extend(child);
        SamplingMapContext {
            previous,
            policies: self.policies.as_deref_mut(),
            key,
        }
    }

    /// Compare the complete native policy, not merely validity of an older
    /// disk. Two valid radii can induce different densities for the same draw.
    pub(crate) fn validate_decision(&mut self, decision: SamplingProposalDecision) -> Result<()> {
        let Some(policies) = self.policies.as_deref_mut() else {
            return Ok(());
        };
        let canonical = policies.records.get(&self.key).copied();
        if canonical == Some(decision) {
            return Ok(());
        }
        if canonical.is_none() && policies.collecting {
            policies.records.insert(self.key.clone(), decision);
            return Ok(());
        }
        Err(SamplingEvaluationError::UncertainGeometry {
            detail: format!(
                "proposal policy mismatch for {:?}: canonical Arb(1000-bit) decision {canonical:?}, native {:?} decision {decision:?}; phase {}",
                self.key,
                T::sampling_precision(),
                if policies.collecting { "collecting" } else { "sealed" }
            ),
        }.into())
    }
}

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
