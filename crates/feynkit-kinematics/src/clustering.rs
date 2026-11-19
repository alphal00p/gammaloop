use std::{cmp::Ordering, collections::BTreeSet};

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{FourMomentum, KinematicScalar, momentum::wrapped_delta_phi};

/// Generalized-kt exponent choices supported by FastJet's standard algorithms.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum JetAlgorithm {
    /// `p = 1`: soft particles cluster first.
    Kt,
    /// `p = 0`: clustering is ordered only by angular separation.
    CambridgeAachen,
    /// `p = -1`: hard particles act as stable clustering seeds.
    #[default]
    AntiKt,
}

/// A generalized-kt clustering definition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JetDefinition<T> {
    algorithm: JetAlgorithm,
    radius: T,
    minimum_pt: T,
}

impl<T: KinematicScalar> JetDefinition<T> {
    /// Construct a definition with zero transverse-momentum threshold.
    pub fn new(algorithm: JetAlgorithm, radius: T) -> Self {
        let minimum_pt = radius.zero();
        Self {
            algorithm,
            radius,
            minimum_pt,
        }
    }

    /// Construct a kt definition.
    pub fn kt(radius: T) -> Self {
        Self::new(JetAlgorithm::Kt, radius)
    }

    /// Construct a Cambridge/Aachen definition.
    pub fn cambridge_aachen(radius: T) -> Self {
        Self::new(JetAlgorithm::CambridgeAachen, radius)
    }

    /// Construct an anti-kt definition.
    pub fn anti_kt(radius: T) -> Self {
        Self::new(JetAlgorithm::AntiKt, radius)
    }

    /// Set the inclusive-jet transverse-momentum threshold.
    pub fn with_minimum_pt(mut self, minimum_pt: T) -> Self {
        self.minimum_pt = minimum_pt;
        self
    }

    /// Selected clustering algorithm.
    pub fn algorithm(&self) -> JetAlgorithm {
        self.algorithm
    }

    /// Jet-radius parameter.
    pub fn radius(&self) -> &T {
        &self.radius
    }

    /// Inclusive-jet transverse-momentum threshold.
    pub fn minimum_pt(&self) -> &T {
        &self.minimum_pt
    }

    /// Cluster input momenta using E-scheme recombination.
    ///
    /// Constituents are identified by their positions in `momenta`. Particle
    /// filtering intentionally remains the caller's responsibility.
    pub fn cluster(
        &self,
        momenta: &[FourMomentum<T>],
    ) -> Result<ClusteringResult<T>, ClusteringError> {
        self.cluster_indexed(momenta.iter().cloned().enumerate())
    }

    /// Cluster momenta carrying caller-assigned constituent indices.
    ///
    /// This lets an event layer filter particles before clustering while keeping
    /// indices relative to the unfiltered event record.
    pub fn cluster_indexed(
        &self,
        momenta: impl IntoIterator<Item = (usize, FourMomentum<T>)>,
    ) -> Result<ClusteringResult<T>, ClusteringError> {
        self.validate()?;

        let mut indices = BTreeSet::new();
        let mut active = momenta
            .into_iter()
            .map(|(index, momentum)| {
                if !indices.insert(index) {
                    return Err(ClusteringError::DuplicateInputIndex { index });
                }
                PseudoJet::from_momentum(index, momentum)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut jets = Vec::new();
        let inverse_radius_squared = (self.radius.clone() * self.radius.clone()).inv();

        while !active.is_empty() {
            let mut best_action = BestAction::Beam(0);
            let mut best_distance: Option<Distance<T>> = None;

            for left in 0..active.len() {
                update_best(
                    &mut best_action,
                    &mut best_distance,
                    BestAction::Beam(left),
                    beam_distance(&active[left], self.algorithm),
                );

                for right in (left + 1)..active.len() {
                    update_best(
                        &mut best_action,
                        &mut best_distance,
                        BestAction::Pair(left, right),
                        pair_distance(
                            &active[left],
                            &active[right],
                            &inverse_radius_squared,
                            self.algorithm,
                        ),
                    );
                }
            }

            match best_action {
                BestAction::Beam(index) => {
                    let jet = active.remove(index).into_jet();
                    if jet.pt() >= self.minimum_pt {
                        jets.push(jet);
                    }
                }
                BestAction::Pair(left, right) => {
                    let merged = PseudoJet::merged(&active[left], &active[right])?;
                    active[left] = merged;
                    active.remove(right);
                }
            }
        }

        jets.sort_unstable_by(compare_jets_descending_pt);
        Ok(ClusteringResult { jets })
    }

    fn validate(&self) -> Result<(), ClusteringError> {
        let zero = self.radius.zero();
        if !self.radius.is_finite() || self.radius <= zero {
            return Err(ClusteringError::InvalidRadius);
        }
        if !self.minimum_pt.is_finite() || self.minimum_pt < self.minimum_pt.zero() {
            return Err(ClusteringError::InvalidMinimumPt);
        }
        Ok(())
    }
}

/// A clustered jet and the input positions contributing to it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Jet<T> {
    /// E-scheme sum of constituent four-momenta.
    pub momentum: FourMomentum<T>,
    constituents: Vec<usize>,
}

impl<T: KinematicScalar> Jet<T> {
    /// Transverse momentum.
    pub fn pt(&self) -> T {
        self.momentum.pt()
    }

    /// Transverse momentum squared.
    pub fn pt_squared(&self) -> T {
        self.momentum.pt_squared()
    }

    /// Rapidity.
    pub fn rapidity(&self) -> T {
        self.momentum.rapidity()
    }

    /// Azimuthal angle in `[0, 2*pi)`.
    pub fn phi(&self) -> T {
        self.momentum.phi()
    }

    /// Sorted positions of contributing input momenta.
    pub fn constituent_indices(&self) -> &[usize] {
        &self.constituents
    }
}

/// Inclusive jets sorted by decreasing transverse momentum.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ClusteringResult<T> {
    /// Clustered inclusive jets.
    pub jets: Vec<Jet<T>>,
}

impl<T> ClusteringResult<T> {
    /// Number of inclusive jets.
    pub fn len(&self) -> usize {
        self.jets.len()
    }

    /// Whether no jet passed the transverse-momentum threshold.
    pub fn is_empty(&self) -> bool {
        self.jets.is_empty()
    }

    /// Highest-transverse-momentum jet.
    pub fn leading_jet(&self) -> Option<&Jet<T>> {
        self.jets.first()
    }
}

/// Invalid clustering configuration or input kinematics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum ClusteringError {
    /// Radius must be finite and positive.
    #[error("jet radius must be finite and positive")]
    InvalidRadius,
    /// Minimum transverse momentum must be finite and nonnegative.
    #[error("minimum jet transverse momentum must be finite and nonnegative")]
    InvalidMinimumPt,
    /// An input four-momentum has a non-finite component.
    #[error("input momentum {index} contains a non-finite component")]
    NonFiniteMomentum {
        /// Position in the input slice.
        index: usize,
    },
    /// Constituent indices must identify inputs unambiguously.
    #[error("input constituent index {index} occurs more than once")]
    DuplicateInputIndex {
        /// Repeated caller-assigned index.
        index: usize,
    },
    /// Rapidity or azimuth is undefined for an input or merged momentum.
    #[error("jet kinematics are non-finite for constituents containing input {index}")]
    UndefinedKinematics {
        /// First input constituent identifying the failed pseudojet.
        index: usize,
    },
}

#[derive(Clone, Copy)]
enum BestAction {
    Beam(usize),
    Pair(usize, usize),
}

#[derive(Clone)]
enum Distance<T> {
    Finite(T),
    Infinite,
}

impl<T: PartialOrd> Distance<T> {
    fn is_less_than(&self, rhs: &Self) -> bool {
        match (self, rhs) {
            (Self::Finite(left), Self::Finite(right)) => left < right,
            (Self::Finite(_), Self::Infinite) => true,
            (Self::Infinite, _) => false,
        }
    }

    fn minimum(self, rhs: Self) -> Self {
        if rhs.is_less_than(&self) { rhs } else { self }
    }
}

#[derive(Debug, Clone)]
struct PseudoJet<T> {
    momentum: FourMomentum<T>,
    constituents: Vec<usize>,
    pt_squared: T,
    rapidity: T,
    phi: T,
}

impl<T: KinematicScalar> PseudoJet<T> {
    fn from_momentum(index: usize, momentum: FourMomentum<T>) -> Result<Self, ClusteringError> {
        if momentum
            .components()
            .into_iter()
            .any(|component| !component.is_finite())
        {
            return Err(ClusteringError::NonFiniteMomentum { index });
        }
        Self::new(momentum, vec![index])
    }

    fn new(
        momentum: FourMomentum<T>,
        mut constituents: Vec<usize>,
    ) -> Result<Self, ClusteringError> {
        constituents.sort_unstable();
        let pt_squared = momentum.pt_squared();
        let rapidity = momentum.rapidity_with_pt_squared(&pt_squared);
        let phi = momentum.phi();
        if !pt_squared.is_finite() || !rapidity.is_finite() || !phi.is_finite() {
            return Err(ClusteringError::UndefinedKinematics {
                index: constituents[0],
            });
        }
        Ok(Self {
            momentum,
            constituents,
            pt_squared,
            rapidity,
            phi,
        })
    }

    fn merged(lhs: &Self, rhs: &Self) -> Result<Self, ClusteringError> {
        let mut constituents = lhs.constituents.clone();
        constituents.extend_from_slice(&rhs.constituents);
        Self::new(lhs.momentum.clone() + rhs.momentum.clone(), constituents)
    }

    fn into_jet(self) -> Jet<T> {
        Jet {
            momentum: self.momentum,
            constituents: self.constituents,
        }
    }
}

fn update_best<T: KinematicScalar>(
    best_action: &mut BestAction,
    best_distance: &mut Option<Distance<T>>,
    action: BestAction,
    candidate: Distance<T>,
) {
    if best_distance
        .as_ref()
        .is_none_or(|current| candidate.is_less_than(current))
    {
        *best_action = action;
        *best_distance = Some(candidate);
    }
}

fn beam_distance<T: KinematicScalar>(jet: &PseudoJet<T>, algorithm: JetAlgorithm) -> Distance<T> {
    scale_for_algorithm(&jet.pt_squared, algorithm)
}

fn pair_distance<T: KinematicScalar>(
    left: &PseudoJet<T>,
    right: &PseudoJet<T>,
    inverse_radius_squared: &T,
    algorithm: JetAlgorithm,
) -> Distance<T> {
    let scale = scale_for_algorithm(&left.pt_squared, algorithm)
        .minimum(scale_for_algorithm(&right.pt_squared, algorithm));
    let delta_rapidity = left.rapidity.clone() - right.rapidity.clone();
    let delta_phi = wrapped_delta_phi(left.phi.clone(), right.phi.clone());
    let angular_distance = delta_rapidity.clone() * delta_rapidity + delta_phi.clone() * delta_phi;
    match scale {
        Distance::Finite(scale) => {
            Distance::Finite(scale * angular_distance * inverse_radius_squared.clone())
        }
        Distance::Infinite => Distance::Infinite,
    }
}

fn scale_for_algorithm<T: KinematicScalar>(pt_squared: &T, algorithm: JetAlgorithm) -> Distance<T> {
    match algorithm {
        JetAlgorithm::Kt => Distance::Finite(pt_squared.clone()),
        JetAlgorithm::CambridgeAachen => Distance::Finite(pt_squared.one()),
        JetAlgorithm::AntiKt if *pt_squared == pt_squared.zero() => Distance::Infinite,
        JetAlgorithm::AntiKt => Distance::Finite(pt_squared.inv()),
    }
}

fn compare_jets_descending_pt<T: KinematicScalar>(lhs: &Jet<T>, rhs: &Jet<T>) -> Ordering {
    rhs.pt()
        .partial_cmp(&lhs.pt())
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            lhs.rapidity()
                .partial_cmp(&rhs.rapidity())
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| lhs.phi().partial_cmp(&rhs.phi()).unwrap_or(Ordering::Equal))
        .then_with(|| lhs.constituents.cmp(&rhs.constituents))
}

#[cfg(test)]
mod tests {
    use numerica::domains::float::Real;

    use super::*;

    fn massless_momentum(pt: f64, rapidity: f64, phi: f64) -> FourMomentum<f64> {
        FourMomentum::from_args(
            pt * rapidity.cosh(),
            pt * phi.cos(),
            pt * phi.sin(),
            pt * rapidity.sinh(),
        )
    }

    fn assert_close(lhs: f64, rhs: f64) {
        let scale = lhs.abs().max(rhs.abs()).max(1.0);
        assert!((lhs - rhs).abs() <= 1.0e-11 * scale, "{lhs} != {rhs}");
    }

    #[test]
    fn all_generalized_kt_algorithms_merge_nearby_particles() {
        let momenta = [
            massless_momentum(80.0, 0.2, 0.1),
            massless_momentum(25.0, 0.25, 0.18),
            massless_momentum(35.0, -1.1, 2.4),
        ];

        for algorithm in [
            JetAlgorithm::Kt,
            JetAlgorithm::CambridgeAachen,
            JetAlgorithm::AntiKt,
        ] {
            let result = JetDefinition::new(algorithm, 0.6)
                .with_minimum_pt(5.0)
                .cluster(&momenta)
                .unwrap();
            assert_eq!(result.len(), 2, "algorithm: {algorithm:?}");
            assert_eq!(
                result.leading_jet().unwrap().constituent_indices(),
                &[0, 1],
                "algorithm: {algorithm:?}"
            );
            assert_eq!(result.jets[1].constituent_indices(), &[2]);
            assert_close(
                result.leading_jet().unwrap().momentum.temporal.value,
                momenta[0].temporal.value + momenta[1].temporal.value,
            );
        }
    }

    #[test]
    fn generalized_kt_exponent_changes_the_clustering_history() {
        let momenta = [
            massless_momentum(100.0, 0.0, 0.0),
            massless_momentum(10.0, 0.5, 0.0),
            massless_momentum(1.0, 0.9, 0.0),
        ];

        let kt = JetDefinition::kt(0.6).cluster(&momenta).unwrap();
        let cambridge = JetDefinition::cambridge_aachen(0.6)
            .cluster(&momenta)
            .unwrap();
        let anti_kt = JetDefinition::anti_kt(0.6).cluster(&momenta).unwrap();

        assert_eq!(kt.jets[0].constituent_indices(), &[0, 1, 2]);
        assert_eq!(cambridge.jets[0].constituent_indices(), &[0, 1, 2]);
        assert_eq!(anti_kt.len(), 2);
        assert_eq!(anti_kt.jets[0].constituent_indices(), &[0, 1]);
        assert_eq!(anti_kt.jets[1].constituent_indices(), &[2]);
    }

    #[test]
    fn wraps_azimuth_when_deciding_a_merge() {
        let epsilon = 0.02;
        let two_pi = 0.0.pi() * 2.0;
        let momenta = [
            massless_momentum(20.0, 0.0, epsilon),
            massless_momentum(15.0, 0.0, two_pi - epsilon),
        ];

        let result = JetDefinition::anti_kt(0.2).cluster(&momenta).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result.jets[0].constituent_indices(), &[0, 1]);
    }

    #[test]
    fn minimum_pt_is_applied_after_recombination() {
        let momenta = [
            massless_momentum(6.0, 0.0, 0.0),
            massless_momentum(6.0, 0.0, 0.05),
            massless_momentum(4.0, 2.0, 2.0),
        ];

        let result = JetDefinition::kt(0.4)
            .with_minimum_pt(10.0)
            .cluster(&momenta)
            .unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result.jets[0].constituent_indices(), &[0, 1]);
    }

    #[test]
    fn indexed_clustering_preserves_pre_filter_event_positions() {
        let inputs = [
            (4, massless_momentum(20.0, 0.0, 0.0)),
            (9, massless_momentum(15.0, 0.05, 0.05)),
        ];

        let result = JetDefinition::anti_kt(0.4).cluster_indexed(inputs).unwrap();
        assert_eq!(result.jets[0].constituent_indices(), &[4, 9]);
    }

    #[test]
    fn rejects_invalid_configuration_and_input() {
        let momentum = massless_momentum(10.0, 0.0, 0.0);
        assert_eq!(
            JetDefinition::anti_kt(0.0).cluster(&[momentum]),
            Err(ClusteringError::InvalidRadius)
        );
        assert_eq!(
            JetDefinition::anti_kt(0.4)
                .with_minimum_pt(-1.0)
                .cluster(&[momentum]),
            Err(ClusteringError::InvalidMinimumPt)
        );

        let invalid = FourMomentum::from_args(f64::NAN, 1.0, 0.0, 0.0);
        assert_eq!(
            JetDefinition::anti_kt(0.4).cluster(&[invalid]),
            Err(ClusteringError::NonFiniteMomentum { index: 0 })
        );
        assert_eq!(
            JetDefinition::anti_kt(0.4).cluster_indexed([
                (3, massless_momentum(10.0, 0.0, 0.0)),
                (3, massless_momentum(5.0, 1.0, 1.0)),
            ]),
            Err(ClusteringError::DuplicateInputIndex { index: 3 })
        );
    }
}
