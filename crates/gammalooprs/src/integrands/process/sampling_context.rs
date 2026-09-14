//! Native classification for graph-aware sampling maps.
//!
//! Graph, host cut, side and complete parent identity belong to the immutable
//! resolved block binding. Each conditional map prepares its actual native t*
//! and fixed physical coordinates together from the declared raw prerequisites.
//! Each foreign inverse uses its supplied raw point, reusing a native preparation
//! only for an identical plan and source. No partially sampled data are labelled a complete MomentumSample,
//! and precision rescue reconstructs the original cube rather than promoting a
//! previously mapped point or pairing it with stale cut data.

use crate::{
    cff::esurface::EsurfaceRay,
    processes::{CutGroupId, CutId},
    utils::{
        FloatLike,
        newton_solver::{NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity},
    },
};
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, sync::Arc};

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

/// One bound coordinate-subset plan. The graph label is diagnostic; actual graph
/// identity is the existing row/source graph ID, checked at physical adoption.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SamplingLUHostPlan {
    pub(crate) graph_name: String,
    pub(crate) cut_group_id: CutGroupId,
    pub(crate) representative_cut_id: CutId,
    pub(crate) parent_lmb: Vec<usize>,
    pub(crate) required_prior_lmb: Vec<usize>,
}

/// Native source authority, never promoted between precision attempts or
/// replaced by the rounded prerequisites reconstructed during inverse checks.
#[derive(Clone, Debug)]
pub(crate) struct PreparedLUHost<T: FloatLike> {
    pub(crate) plan: Arc<SamplingLUHostPlan>,
    pub(crate) source: SamplingProposalKey,
    pub(crate) prior: Vec<T>,
    pub(crate) ray: EsurfaceRay<T>,
    pub(crate) solution: NewtonIterationResult<T>,
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
    pub(crate) prepared_lu_hosts: Option<&'a mut Vec<PreparedLUHost<T>>>,
    pub(crate) radial_root_diagnostics: Option<&'a mut RadialRootDiagnostics>,
    pub(crate) selecting_lu_hosts: bool,
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
            prepared_lu_hosts: None,
            radial_root_diagnostics: None,
            selecting_lu_hosts: false,
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
            prepared_lu_hosts: self.prepared_lu_hosts.as_deref_mut(),
            radial_root_diagnostics: self.radial_root_diagnostics.as_deref_mut(),
            selecting_lu_hosts: self.selecting_lu_hosts,
        }
    }

    /// Reuse only exact represented prerequisites. During the initial selected
    /// operation one plan has one source; inverse work may prepare a different
    /// source but cannot extend the authoritative prefix frozen by the bridge.
    pub(crate) fn prepare_lu_host(
        &mut self,
        plan: &Arc<SamplingLUHostPlan>,
        prior: Vec<T>,
        prepare: impl FnOnce(
            &mut RadialRootDiagnostics,
            &RadialRootIdentity,
        ) -> Result<(EsurfaceRay<T>, NewtonIterationResult<T>)>,
    ) -> Result<&PreparedLUHost<T>> {
        if prior.len() != 3 * plan.required_prior_lmb.len() {
            return Err(eyre!(
                "LU host plan {plan:?} needs {} canonical prior components, received {}",
                3 * plan.required_prior_lmb.len(),
                prior.len()
            ));
        }
        if prior.iter().any(|value| !value.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "LU host prerequisites",
                detail: format!("nonfinite native prior for {plan:?}"),
            }
            .into());
        }
        let hosts = self
            .prepared_lu_hosts
            .as_deref_mut()
            .ok_or_else(|| eyre!("LU host preparation requires a runtime row context"))?;
        if self.selecting_lu_hosts
            && hosts
                .iter()
                .any(|host| host.plan == *plan && host.prior != prior)
        {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: format!(
                    "selected LU host plan {plan:?} received inconsistent native prerequisites"
                ),
            }
            .into());
        }
        if let Some(index) = hosts
            .iter()
            .position(|host| host.plan == *plan && host.prior == prior)
        {
            return Ok(&hosts[index]);
        }
        let identity = RadialRootIdentity::new(format!(
            "sampling LU host {plan:?}, source {:?}, {} {}",
            self.key,
            if self
                .policies
                .as_deref()
                .is_some_and(SamplingProposalPolicies::is_collecting)
            {
                "canonical"
            } else {
                "native"
            },
            if self.selecting_lu_hosts {
                "selected"
            } else {
                "inverse"
            },
        ));
        let diagnostics = self
            .radial_root_diagnostics
            .as_deref_mut()
            .ok_or_else(|| eyre!("LU host preparation requires the runtime root diagnostics"))?;
        let (ray, solution) = prepare(diagnostics, &identity)?;
        hosts.push(PreparedLUHost {
            plan: Arc::clone(plan),
            source: self.key.clone(),
            prior,
            ray,
            solution,
        });
        Ok(hosts.last().expect("prepared native host was inserted"))
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
