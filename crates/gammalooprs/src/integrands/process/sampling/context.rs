//! Native classification for graph-aware sampling maps.
//!
//! Graph, host cut, side and complete parent identity belong to the immutable
//! resolved block binding. Each conditional map prepares its actual native t*
//! and fixed physical coordinates together from the declared raw prerequisites.
//! Each foreign inverse uses its supplied raw point, reusing a native preparation
//! only for an identical plan and source. No partially sampled data are labelled
//! a complete MomentumSample. Production exactly embeds one completed source
//! draw in canonical Arb storage; physical precision attempts materialize
//! directly from that immutable authority, never promoting a previous physical
//! retry or pairing its point with stale cut data.

use crate::{
    cff::esurface::EsurfaceRay,
    processes::{CutGroupId, CutId},
    utils::{
        ArbPrec, F, FloatLike,
        newton_solver::{NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity},
    },
};
use color_eyre::eyre::{Result, WrapErr, eyre};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, sync::Arc};

use super::{maps::SamplingEvaluationError, selection::SamplingChannelId};

/// A component's discrete recipe is authenticated during canonical preparation.
/// Physical retries consume the completed draw without selecting another recipe.
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

/// Native source authority. Physical attempts materialize directly from its
/// canonical Arb instance; rounded inverse prerequisites cannot replace it.
#[derive(Clone, Debug)]
pub(crate) struct PreparedLUHost<T: FloatLike> {
    pub(crate) plan: Arc<SamplingLUHostPlan>,
    pub(crate) source: SamplingProposalKey,
    pub(crate) prior: Vec<T>,
    pub(crate) ray: EsurfaceRay<T>,
    pub(crate) solution: NewtonIterationResult<T>,
}

impl<T: FloatLike> PreparedLUHost<T> {
    /// Freeze the completed source authority without reconstructing its ray or
    /// revisiting any root, prerequisite or discrete source decision.
    pub(crate) fn to_arb_exact(&self) -> Result<PreparedLUHost<ArbPrec>> {
        Ok(PreparedLUHost {
            plan: Arc::clone(&self.plan),
            source: self.source.clone(),
            prior: self
                .prior
                .iter()
                .map(|value| F(value.clone()).to_arb_exact().map(|value| value.0))
                .collect::<Result<_>>()?,
            ray: self.ray.to_arb_exact()?,
            solution: NewtonIterationResult {
                solution: self.solution.solution.to_arb_exact()?,
                derivative_at_solution: self.solution.derivative_at_solution.to_arb_exact()?,
                error_of_function: self.solution.error_of_function.to_arb_exact()?,
                num_iterations_used: self.solution.num_iterations_used,
            },
        })
    }
}

impl PreparedLUHost<ArbPrec> {
    /// Convert the canonical authority, retaining its identities and root
    /// history. This is not a new solve or a native map-preparation observation.
    pub(crate) fn materialize<T: FloatLike>(&self, tolerance: f64) -> Result<PreparedLUHost<T>> {
        let positive = |value: &F<ArbPrec>, label: &str| -> Result<F<T>> {
            let native = F::<T>::from_arb(&value.0)?;
            if value <= &value.zero() || native <= native.zero() {
                return Err(SamplingEvaluationError::UncertainGeometry {
                    detail: format!("LU host materialization requires positive {label}"),
                }
                .into());
            }
            // Compare to the original represented source, including separated
            // native limbs. The existing physical adoption subsequently checks
            // the completed ray, residual, slope and inverse-volume budget.
            native
                .verify_arb_materialization(&value.0, tolerance)
                .wrap_err_with(|| format!("LU host {label} for {:?}", self.source))?;
            Ok(native)
        };
        Ok(PreparedLUHost {
            plan: Arc::clone(&self.plan),
            source: self.source.clone(),
            prior: self
                .prior
                .iter()
                .map(|value| F::<T>::from_arb(value).map(|value| value.0))
                .collect::<Result<_>>()?,
            ray: self.ray.materialize()?,
            solution: NewtonIterationResult {
                solution: positive(&self.solution.solution, "root")?,
                derivative_at_solution: positive(
                    &self.solution.derivative_at_solution,
                    "root derivative",
                )?,
                error_of_function: F::<T>::from_arb(&self.solution.error_of_function.0)?,
                num_iterations_used: self.solution.num_iterations_used,
            },
        })
    }
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
    fn canonical_lu_host_materializes_directly_without_new_root_history() -> Result<()> {
        use crate::{
            DependentMomentaConstructor,
            cff::{VertexSet, esurface::Esurface},
            dot,
            graph::{Graph, GroupId, parse::from_dot::IntoGraph},
            initialisation::test_initialise,
            integrands::process::gammaloop_sample::{DiscreteGraphSample, GammaLoopSample},
            momentum::{
                ExternalMomenta, FourMomentum, ThreeMomentum,
                sample::{ExternalFourMomenta, LoopMomenta, MomentumSample},
            },
            settings::runtime::kinematic::Externals,
            utils::{QuadFloat, SamplingFloat},
        };
        use linnet::half_edge::involution::EdgeIndex;

        test_initialise()?;
        let graph: Graph = dot!(digraph canonical_host {
            ext [style=invis]
            node [num=1]
            edge [num=1 mass=1]
            ext -> a:0 [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b:1 -> ext [id=3]
        })?;
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let one = F::<ArbPrec>::default().one();
        let zero = one.zero();
        let prior = vec![
            (&one + &one / one.from_usize(10).powi(25)).0,
            (&one / one.from_usize(7)).0,
            (&one / one.from_usize(11)).0,
        ];
        let loops = LoopMomenta::from_iter([ThreeMomentum::new(
            F(prior[0].clone()),
            F(prior[1].clone()),
            F(prior[2].clone()),
        )]);
        let externals = ExternalFourMomenta::from_iter((0..2).map(|_| {
            FourMomentum::from_args(one.from_usize(4), zero.clone(), zero.clone(), zero.clone())
        }));
        let masses = graph.underlying.new_edgevec_from_iter([
            zero.clone(),
            one.clone(),
            one.clone(),
            zero.clone(),
        ])?;
        let mut history = RadialRootDiagnostics::default();
        let (ray, solution) = surface
            .solve_lu_cut(
                &loops,
                &externals,
                &masses,
                &graph.loop_momentum_basis,
                &one.from_usize(4),
                &mut history,
                &RadialRootIdentity::new("canonical fixture".into()),
            )
            .map_err(|error| eyre!("{error:?}"))?;
        let original_history = format!("{history:?}");
        let host = PreparedLUHost {
            plan: Arc::new(SamplingLUHostPlan {
                graph_name: graph.name.clone(),
                cut_group_id: CutGroupId::from(0),
                representative_cut_id: CutId(0),
                parent_lmb: vec![1],
                required_prior_lmb: vec![1],
            }),
            source: SamplingProposalKey {
                graph_id: 2,
                generating_channel: SamplingChannelId(3),
                target_channel: SamplingChannelId(3),
                block_path: vec![1, 0],
            },
            prior,
            ray,
            solution,
        };
        let double = host.materialize::<f64>(1e-10)?;
        let quad = host.materialize::<QuadFloat>(1e-28)?;
        let arb = host.materialize::<ArbPrec>(1e-90)?;
        assert_eq!(double.prior[0], 1.0);
        assert_ne!(
            quad.prior[0],
            QuadFloat::from_f64_exact_binary(double.prior[0])
        );
        assert_eq!(arb.prior, host.prior);
        assert_eq!(arb.solution.solution, host.solution.solution);
        assert_eq!(
            arb.solution.error_of_function,
            host.solution.error_of_function
        );
        assert!(Arc::ptr_eq(&quad.plan, &host.plan));
        assert_eq!(quad.source, host.source);
        assert_eq!(
            quad.solution.num_iterations_used,
            host.solution.num_iterations_used
        );
        assert_eq!(format!("{history:?}"), original_history);
        let mut source = host.materialize::<SamplingFloat>(1e-60)?;
        source.solution.error_of_function = source.ray.evaluate(&source.solution.solution).0;
        let promoted = source.to_arb_exact()?;
        assert!(Arc::ptr_eq(&promoted.plan, &source.plan));
        assert_eq!(promoted.source, source.source);
        assert_eq!(
            promoted.prior,
            source
                .prior
                .iter()
                .map(|value| F(value.clone()).to_arb_exact().unwrap().0)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            promoted.solution.solution,
            source.solution.solution.to_arb_exact()?
        );
        assert_eq!(
            promoted.solution.derivative_at_solution,
            source.solution.derivative_at_solution.to_arb_exact()?
        );
        assert_eq!(
            promoted.solution.error_of_function,
            source.solution.error_of_function.to_arb_exact()?
        );
        assert_eq!(
            promoted.solution.num_iterations_used,
            source.solution.num_iterations_used
        );
        let roundtrip = promoted.materialize::<SamplingFloat>(1e-60)?;
        assert_eq!(format!("{roundtrip:?}"), format!("{source:?}"));

        // Exercise the completed row boundary too: a selected host survives
        // exact promotion and every physical retry without a second root solve.
        let mut external_settings = Externals::default();
        let Externals::Constant { momenta, .. } = &mut external_settings;
        *momenta = vec![ExternalMomenta::Independent([4.0, 0.0, 0.0, 0.0].map(F)); 2];
        let source_one = F::<SamplingFloat>::default().one();
        let sample = MomentumSample::new(
            LoopMomenta::from_iter([ThreeMomentum::new(
                F(source.prior[0].clone()),
                F(source.prior[1].clone()),
                F(source.prior[2].clone()),
            )]),
            33,
            &external_settings,
            35,
            &source_one / source_one.from_usize(7),
            DependentMomentaConstructor::CrossSection,
            Some(4),
        )?;
        let draw = GammaLoopSample {
            groups: vec![(
                Some(GroupId(1)),
                vec![DiscreteGraphSample {
                    graph_id: source.source.graph_id,
                    channel_id: Some(source.source.generating_channel),
                    prepared_lu_hosts: vec![source.clone()],
                    physical_overlaps: None,
                    sample,
                    integrand_prefactor: source_one,
                }],
            )],
        };
        let canonical = draw.clone().into_canonical(
            &external_settings,
            DependentMomentaConstructor::CrossSection,
        )?;
        let retained = &canonical.groups[0].1[0].prepared_lu_hosts;
        assert_eq!(retained.len(), 1);
        assert!(Arc::ptr_eq(&retained[0].plan, &source.plan));
        assert_eq!(format!("{:?}", retained[0]), format!("{promoted:?}"));
        let frozen = format!("{canonical:?}");
        canonical.materialize::<f64>(1e-10)?;
        canonical.materialize::<QuadFloat>(1e-28)?;
        canonical.materialize::<ArbPrec>(1e-60)?;
        assert_eq!(format!("{canonical:?}"), frozen);
        assert_eq!(format!("{history:?}"), original_history);
        let mut invalid_row = draw;
        invalid_row.groups[0].1[0].prepared_lu_hosts[0].prior[0] =
            SamplingFloat::from_f64_exact_binary(f64::INFINITY);
        assert!(
            invalid_row
                .into_canonical(
                    &external_settings,
                    DependentMomentaConstructor::CrossSection
                )
                .is_err()
        );
        for field in 0..3 {
            let mut invalid = source.clone();
            let nonfinite = F(SamplingFloat::from_f64_exact_binary(f64::INFINITY));
            match field {
                0 => invalid.solution.solution = nonfinite,
                1 => invalid.solution.derivative_at_solution = nonfinite,
                _ => invalid.solution.error_of_function = nonfinite,
            }
            assert!(invalid.to_arb_exact().is_err());
        }
        assert!(matches!(
            host.materialize::<f64>(1e-25)
                .unwrap_err()
                .downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::UncertainGeometry { .. })
        ));
        let mut invalid = host.clone();
        invalid.solution.solution = -one.clone();
        assert!(invalid.materialize::<f64>(1e-10).is_err());
        invalid.solution.solution = &one / one.from_usize(10).powi(400);
        assert!(matches!(
            invalid
                .materialize::<f64>(1e-10)
                .unwrap_err()
                .downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::Unrepresentable { .. })
        ));
        Ok(())
    }

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
