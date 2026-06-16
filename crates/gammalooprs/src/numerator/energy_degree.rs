use std::{
    borrow::Borrow,
    collections::{BTreeMap, BTreeSet},
};

use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski},
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::domains::rational::Rational;
use thiserror::Error;
use three_dimensional_reps::utils::{rank_i64, solve_rational_system};

use crate::{
    graph::Graph,
    momentum::SignOrZero,
    utils::{GS, symbols::UvMomentumProvenanceRole},
    uv::UltravioletGraph,
};

// Bound native CFF trials while sharing the same budget for hard and soft dispatch.
const ENERGY_ASSIGNMENT_PROPOSAL_BUDGET: usize = 3;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct EnergyPowerCapMap {
    degrees: BTreeMap<EdgeIndex, usize>,
}

impl EnergyPowerCapMap {
    pub fn is_empty(&self) -> bool {
        self.degrees.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (EdgeIndex, usize)> + '_ {
        self.degrees.iter().map(|(edge, degree)| (*edge, *degree))
    }

    pub fn into_generation_bounds(self) -> Vec<(usize, usize)> {
        self.degrees
            .into_iter()
            .map(|(edge, degree)| (usize::from(edge), degree))
            .collect()
    }

    fn unit(edge: EdgeIndex) -> Self {
        let mut degrees = BTreeMap::new();
        degrees.insert(edge, 1);
        Self { degrees }
    }

    fn add_assign(&mut self, other: Self) -> Result<(), EnergyPowerAnalysisError> {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = slot
                .checked_add(degree)
                .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
        }
        Ok(())
    }

    fn max_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = (*slot).max(degree);
        }
    }

    fn scale(&mut self, exponent: usize) -> Result<(), EnergyPowerAnalysisError> {
        for degree in self.degrees.values_mut() {
            *degree = degree
                .checked_mul(exponent)
                .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
        }
        Ok(())
    }
}

/// Serial exact-energy occurrences certified for one immutable source owner.
///
/// Untagged and fixed numerator factors stay on the explicitly designated base
/// occurrence. Denominator-derived factors are load-balanced only over the
/// derivative-created copies of that owner. Momentum-routing signs are restored
/// while mapping the numerator, and candidate sets belonging to distinct source
/// owners remain disjoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum EnergyCandidateFamily {
    Unprovenanced(EdgeIndex),
    Fixed(EdgeIndex),
    DenominatorDerived(EdgeIndex),
}

impl EnergyCandidateFamily {
    fn edge(self) -> EdgeIndex {
        match self {
            Self::Unprovenanced(edge) | Self::Fixed(edge) | Self::DenominatorDerived(edge) => edge,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EquivalentEnergyCandidates {
    by_family: BTreeMap<EnergyCandidateFamily, Vec<usize>>,
}

impl EquivalentEnergyCandidates {
    /// Build candidate sets from the serial denominator occurrences certified
    /// by source-edge topology reconstruction. Algebraic energy equality is
    /// not rediscovered here: members already share one immutable source owner.
    pub(crate) fn try_from_source_occurrences(
        groups: impl IntoIterator<Item = (EdgeIndex, Vec<usize>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        let mut by_family = BTreeMap::new();
        let mut candidate_owners = BTreeMap::<usize, EdgeIndex>::new();

        for (edge, mut candidates) in groups {
            if by_family.contains_key(&EnergyCandidateFamily::Fixed(edge)) {
                return Err(
                    EnergyPowerAnalysisError::DuplicatePhysicalEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            candidates.sort_unstable();
            if candidates.is_empty() {
                return Err(
                    EnergyPowerAnalysisError::EmptyEquivalentEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            for pair in candidates.windows(2) {
                if pair[0] == pair[1] {
                    return Err(
                        EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate {
                            edge: edge.into(),
                            candidate: pair[0],
                        },
                    );
                }
            }
            for candidate in &candidates {
                if let Some(first_edge) = candidate_owners.insert(*candidate, edge) {
                    return Err(
                        EnergyPowerAnalysisError::OverlappingEquivalentEnergyCandidates {
                            candidate: *candidate,
                            first_edge: first_edge.into(),
                            second_edge: edge.into(),
                        },
                    );
                }
            }

            // Serial copies can reduce the maximal denominator-derived degree.
            // The lowest canonical occurrence supplies the base when
            // the caller has not explicitly partitioned the serial copies.
            by_family.insert(
                EnergyCandidateFamily::Unprovenanced(edge),
                vec![candidates[0]],
            );
            by_family.insert(EnergyCandidateFamily::Fixed(edge), vec![candidates[0]]);
            by_family.insert(EnergyCandidateFamily::DenominatorDerived(edge), candidates);
        }

        Ok(Self { by_family })
    }

    /// Build production candidate sets from one retained source occurrence and
    /// its derivative-created serial copies. The base is explicit because
    /// exact-topology canonicalization may place it after a generated copy in
    /// the occurrence-local energy namespace.
    pub(crate) fn try_from_partitioned_source_occurrences(
        groups: impl IntoIterator<Item = (EdgeIndex, usize, Vec<usize>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        let mut groups = groups.into_iter().collect::<Vec<_>>();
        let mut candidates = Self::try_from_source_occurrences(groups.iter().map(
            |(edge, base, derivative_copies)| {
                (
                    *edge,
                    std::iter::once(*base)
                        .chain(derivative_copies.iter().copied())
                        .collect(),
                )
            },
        ))?;
        for (edge, base, derivative_copies) in &mut groups {
            derivative_copies.sort_unstable();
            candidates
                .by_family
                .insert(EnergyCandidateFamily::Unprovenanced(*edge), vec![*base]);
            candidates
                .by_family
                .insert(EnergyCandidateFamily::Fixed(*edge), vec![*base]);
            if derivative_copies.is_empty() {
                candidates
                    .by_family
                    .remove(&EnergyCandidateFamily::DenominatorDerived(*edge));
            } else {
                candidates.by_family.insert(
                    EnergyCandidateFamily::DenominatorDerived(*edge),
                    derivative_copies.clone(),
                );
            }
        }
        Ok(candidates)
    }

    fn get(&self, family: EnergyCandidateFamily) -> Option<&[usize]> {
        self.by_family.get(&family).map(Vec::as_slice)
    }
}

// Factor ID -> original edge owner -> certified exact energy occurrence.
type FactorEnergyAssignments = BTreeMap<usize, BTreeMap<EdgeIndex, usize>>;

/// One immutable assignment of factor-local physical energy dependencies to
/// certified equivalent exact energy variables.
///
/// The same plan owns both the bounds passed to CFF generation and the
/// factor-local assignments used when mapping the numerator. This prevents a
/// balanced bound from understating a numerator which is still evaluated
/// wholly through one canonical exact energy.
#[derive(Debug, Clone)]
pub(crate) struct EnergyPowerAssignmentPlan {
    expression: PlannedEnergyExpression,
    factor_assignments: FactorEnergyAssignments,
    energy_degree_bounds: Vec<(usize, usize)>,
}

impl EnergyPowerAssignmentPlan {
    pub(crate) fn energy_degree_bounds(&self) -> &[(usize, usize)] {
        &self.energy_degree_bounds
    }

    pub(crate) fn map_factors<E>(
        &self,
        mut map: impl FnMut(&Atom, &BTreeMap<EdgeIndex, usize>) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        self.expression.map(
            &self.factor_assignments,
            &mut |_, expression, assignments| map(expression, assignments),
        )
    }
}

#[derive(Debug, Clone)]
enum PlannedEnergyExpression {
    Factor {
        id: usize,
        expression: Atom,
        degrees: EnergyPowerCapMap,
        families: BTreeMap<EdgeIndex, EnergyCandidateFamily>,
    },
    Add(Vec<Self>),
    Mul(Vec<Self>),
    MultilinearFunction {
        symbol: Symbol,
        arguments: Vec<Self>,
    },
}

impl PlannedEnergyExpression {
    /// Select independently at factor leaves, adding ranks in products and
    /// taking their componentwise maximum across sums. The rank Pareto frontier
    /// is exact; no numerator polynomial or Cartesian product of its monomials is
    /// ever constructed.
    fn soft_momentum_assignments(
        &self,
        alternatives: &BTreeMap<Atom, Vec<(Atom, EnergyPowerCapMap)>>,
    ) -> Result<Vec<(EnergyPowerCapMap, BTreeMap<usize, usize>)>, EnergyPowerAnalysisError> {
        match self {
            Self::Factor {
                id,
                expression,
                degrees,
                ..
            } => Ok(alternatives
                .get(expression)
                .map(|choices| {
                    choices
                        .iter()
                        .enumerate()
                        .map(|(choice, (_, degrees))| {
                            (degrees.clone(), BTreeMap::from([(*id, choice)]))
                        })
                        .collect()
                })
                .unwrap_or_else(|| vec![(degrees.clone(), BTreeMap::new())])),
            Self::Add(children)
            | Self::Mul(children)
            | Self::MultilinearFunction {
                arguments: children,
                ..
            } => {
                let additive = matches!(self, Self::Add(_));
                children.iter().try_fold(
                    vec![(EnergyPowerCapMap::default(), BTreeMap::new())],
                    |states, child| {
                        let right = child.soft_momentum_assignments(alternatives)?;
                        let mut unique = BTreeMap::new();
                        for (left_degrees, left_assignments) in &states {
                            for (right_degrees, right_assignments) in &right {
                                let mut degrees = left_degrees.clone();
                                if additive {
                                    degrees.max_assign(right_degrees.clone());
                                } else {
                                    degrees.add_assign(right_degrees.clone())?;
                                }
                                unique.entry(degrees.degrees).or_insert_with(|| {
                                    let mut assignments = left_assignments.clone();
                                    assignments.extend(right_assignments);
                                    assignments
                                });
                            }
                        }
                        Ok(unique
                            .iter()
                            .filter(|(degrees, _)| {
                                !unique.keys().any(|other| {
                                    other != *degrees
                                        && other.iter().all(|(edge, degree)| {
                                            *degree <= degrees.get(edge).copied().unwrap_or(0)
                                        })
                                })
                            })
                            .map(|(degrees, assignments)| {
                                (
                                    EnergyPowerCapMap {
                                        degrees: degrees.clone(),
                                    },
                                    assignments.clone(),
                                )
                            })
                            .collect())
                    },
                )
            }
        }
    }

    fn degrees(&self) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        match self {
            Self::Factor { degrees, .. } => Ok(degrees.clone()),
            Self::Add(terms) => {
                terms
                    .iter()
                    .try_fold(EnergyPowerCapMap::default(), |mut sum, term| {
                        sum.max_assign(term.degrees()?);
                        Ok(sum)
                    })
            }
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => factors
                .iter()
                .try_fold(EnergyPowerCapMap::default(), |mut product, factor| {
                    product.add_assign(factor.degrees()?)?;
                    Ok(product)
                }),
        }
    }

    fn assignment_proposals(
        &self,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<Vec<FactorEnergyAssignments>, EnergyPowerAnalysisError> {
        #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
        struct OwnerState {
            energy_loads: Vec<usize>,
            derived_unit_loads: Vec<usize>,
            assignments: BTreeMap<usize, usize>,
        }

        fn prune(mut states: Vec<OwnerState>) -> Vec<OwnerState> {
            states.sort();
            states.dedup_by(|left, right| {
                left.energy_loads == right.energy_loads
                    && left.derived_unit_loads == right.derived_unit_loads
            });
            let mut pareto = Vec::new();
            'candidate: for (index, state) in states.iter().enumerate() {
                for (other_index, other) in states.iter().enumerate() {
                    if index != other_index
                        && other
                            .energy_loads
                            .iter()
                            .zip(&state.energy_loads)
                            .all(|(left, right)| left <= right)
                        && other
                            .derived_unit_loads
                            .iter()
                            .zip(&state.derived_unit_loads)
                            .all(|(left, right)| left <= right)
                    {
                        continue 'candidate;
                    }
                }
                pareto.push(state.clone());
            }
            pareto
        }

        fn combine(
            left: Vec<OwnerState>,
            right: Vec<OwnerState>,
            additive: bool,
        ) -> Result<Vec<OwnerState>, EnergyPowerAnalysisError> {
            let mut combined = Vec::new();
            for left in &left {
                for right in &right {
                    let combine_loads = |left: &[usize], right: &[usize]| {
                        left.iter()
                            .zip(right)
                            .map(|(left, right)| {
                                if additive {
                                    Ok((*left).max(*right))
                                } else {
                                    left.checked_add(*right)
                                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)
                                }
                            })
                            .collect::<Result<Vec<_>, _>>()
                    };
                    let energy_loads = combine_loads(&left.energy_loads, &right.energy_loads)?;
                    let derived_unit_loads =
                        combine_loads(&left.derived_unit_loads, &right.derived_unit_loads)?;
                    let mut assignments = left.assignments.clone();
                    for (factor, candidate) in &right.assignments {
                        // This insertion is part of the algorithm: keeping it
                        // inside `debug_assert!` would erase it in release builds.
                        let previous = assignments.insert(*factor, *candidate);
                        debug_assert!(previous.is_none());
                    }
                    combined.push(OwnerState {
                        energy_loads,
                        derived_unit_loads,
                        assignments,
                    });
                }
            }
            Ok(prune(combined))
        }

        fn frontier(
            expression: &PlannedEnergyExpression,
            owner: EdgeIndex,
            owner_candidates: &[usize],
            candidate_positions: &BTreeMap<usize, usize>,
            candidates: &EquivalentEnergyCandidates,
        ) -> Result<Vec<OwnerState>, EnergyPowerAnalysisError> {
            let zero = || {
                vec![OwnerState {
                    energy_loads: vec![0; owner_candidates.len()],
                    derived_unit_loads: vec![0; owner_candidates.len()],
                    assignments: BTreeMap::new(),
                }]
            };
            match expression {
                PlannedEnergyExpression::Factor {
                    id,
                    degrees,
                    families,
                    ..
                } => {
                    let Some(degree) = degrees.degrees.get(&owner).copied() else {
                        return Ok(zero());
                    };
                    let family = families[&owner];
                    let eligible = candidates.get(family).ok_or(
                        EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates {
                            edge: owner.into(),
                            degree,
                        },
                    )?;
                    Ok(eligible
                        .iter()
                        .map(|candidate| {
                            let mut energy_loads = vec![0; owner_candidates.len()];
                            let position = candidate_positions[candidate];
                            energy_loads[position] = degree;
                            let mut derived_unit_loads = vec![0; owner_candidates.len()];
                            if degree == 1
                                && family == EnergyCandidateFamily::DenominatorDerived(owner)
                            {
                                derived_unit_loads[position] = 1;
                            }
                            OwnerState {
                                energy_loads,
                                derived_unit_loads,
                                assignments: BTreeMap::from([(*id, *candidate)]),
                            }
                        })
                        .collect())
                }
                PlannedEnergyExpression::Add(terms) => {
                    terms.iter().try_fold(zero(), |left, right| {
                        combine(
                            left,
                            frontier(
                                right,
                                owner,
                                owner_candidates,
                                candidate_positions,
                                candidates,
                            )?,
                            true,
                        )
                    })
                }
                PlannedEnergyExpression::Mul(factors)
                | PlannedEnergyExpression::MultilinearFunction {
                    arguments: factors, ..
                } => factors.iter().try_fold(zero(), |left, right| {
                    combine(
                        left,
                        frontier(
                            right,
                            owner,
                            owner_candidates,
                            candidate_positions,
                            candidates,
                        )?,
                        false,
                    )
                }),
            }
        }

        let mut factor_assignments = FactorEnergyAssignments::new();
        let mut owner_alternatives = Vec::new();
        for (owner, _) in self.degrees()?.iter() {
            let mut owner_candidates = candidates
                .get(EnergyCandidateFamily::DenominatorDerived(owner))
                .unwrap_or_default()
                .to_vec();
            if let Some(fixed_candidates) = candidates.get(EnergyCandidateFamily::Fixed(owner)) {
                owner_candidates.extend_from_slice(fixed_candidates);
            }
            owner_candidates.sort_unstable();
            owner_candidates.dedup();
            if owner_candidates.is_empty() {
                return Err(
                    EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates {
                        edge: owner.into(),
                        degree: self.degrees()?.degrees[&owner],
                    },
                );
            }
            let candidate_positions = owner_candidates
                .iter()
                .enumerate()
                .map(|(position, candidate)| (*candidate, position))
                .collect::<BTreeMap<_, _>>();
            let mut states = frontier(
                self,
                owner,
                &owner_candidates,
                &candidate_positions,
                candidates,
            )?;
            states.sort_by_cached_key(|state| {
                // Compare the complete descending envelope: an unrelated
                // saturated occurrence must not hide an avoidable (3, 1)
                // distribution when the same factors fit in (2, 2).
                let mut envelope = state.energy_loads.clone();
                envelope.sort_unstable_by(|left, right| right.cmp(left));
                (
                    envelope,
                    state.derived_unit_loads.iter().copied().max().unwrap_or(0),
                    state
                        .energy_loads
                        .iter()
                        .map(|degree| *degree as u128)
                        .sum::<u128>(),
                    state
                        .derived_unit_loads
                        .iter()
                        .map(|degree| *degree as u128)
                        .sum::<u128>(),
                    state.assignments.clone(),
                    state.energy_loads.clone(),
                    state.derived_unit_loads.clone(),
                )
            });
            let baseline = states.remove(0);
            for (factor, candidate) in baseline.assignments {
                factor_assignments
                    .entry(factor)
                    .or_default()
                    .insert(owner, candidate);
            }
            for alternative in states
                .into_iter()
                .filter(|state| state.energy_loads != baseline.energy_loads)
                .unique_by(|state| state.energy_loads.clone())
                .take(ENERGY_ASSIGNMENT_PROPOSAL_BUDGET - 1)
            {
                owner_alternatives.push((owner, alternative.assignments));
            }
        }

        // Rank supplies the baseline and up to two deterministic one-owner
        // deviations. The caller compares their actual CFF maps.
        // The Pareto frontiers are only a proposal heuristic for that metric;
        // no count-optimality or independent-owner factorization is assumed.
        let baseline_bounds = self.exact_degrees(&factor_assignments)?;
        let occurrence_ids = candidates
            .by_family
            .values()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>();
        let mut alternatives = Vec::new();
        for (owner, owner_assignments) in owner_alternatives {
            let mut assignments = factor_assignments.clone();
            for (factor, candidate) in owner_assignments {
                assignments
                    .entry(factor)
                    .or_default()
                    .insert(owner, candidate);
            }
            let bounds = self.exact_degrees(&assignments)?;
            if bounds == baseline_bounds {
                continue;
            }
            let mut envelope = occurrence_ids
                .iter()
                .map(|edge| bounds.get(edge).copied().unwrap_or(0))
                .collect::<Vec<_>>();
            envelope.sort_unstable_by(|left, right| right.cmp(left));
            alternatives.push((envelope, assignments, bounds));
        }
        alternatives.sort();
        Ok(std::iter::once(factor_assignments)
            .chain(
                alternatives
                    .into_iter()
                    .unique_by(|(_, _, bounds)| bounds.clone())
                    .take(ENERGY_ASSIGNMENT_PROPOSAL_BUDGET - 1)
                    .map(|(_, assignments, _)| assignments),
            )
            .collect())
    }

    fn exact_degrees(
        &self,
        assignments: &FactorEnergyAssignments,
    ) -> Result<BTreeMap<usize, usize>, EnergyPowerAnalysisError> {
        match self {
            Self::Factor {
                id,
                expression,
                degrees,
                ..
            } => {
                let mut exact = BTreeMap::new();
                for (edge, degree) in degrees.iter() {
                    let candidate = assignments
                        .get(id)
                        .and_then(|factor| factor.get(&edge))
                        .copied()
                        .ok_or_else(|| EnergyPowerAnalysisError::MissingFactorEnergyAssignment {
                            factor_id: *id,
                            expression: expression.log_print(None),
                            edge: edge.into(),
                        })?;
                    let slot = exact.entry(candidate).or_insert(0usize);
                    *slot = slot
                        .checked_add(degree)
                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
                }
                Ok(exact)
            }
            Self::Add(terms) => {
                let mut envelope = BTreeMap::new();
                for term in terms {
                    for (candidate, degree) in term.exact_degrees(assignments)? {
                        let slot = envelope.entry(candidate).or_insert(0);
                        *slot = (*slot).max(degree);
                    }
                }
                Ok(envelope)
            }
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => {
                let mut product = BTreeMap::new();
                for factor in factors {
                    for (candidate, degree) in factor.exact_degrees(assignments)? {
                        let slot = product.entry(candidate).or_insert(0usize);
                        *slot = slot
                            .checked_add(degree)
                            .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
                    }
                }
                Ok(product)
            }
        }
    }

    fn map<E>(
        &self,
        assignments: &FactorEnergyAssignments,
        map: &mut impl FnMut(usize, &Atom, &BTreeMap<EdgeIndex, usize>) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        match self {
            Self::Factor { id, expression, .. } => {
                let empty = BTreeMap::new();
                map(*id, expression, assignments.get(id).unwrap_or(&empty))
            }
            Self::Add(terms) => {
                terms.iter().try_fold(
                    Atom::Zero,
                    |sum, term| Ok(sum + term.map(assignments, map)?),
                )
            }
            Self::Mul(factors) => factors.iter().try_fold(Atom::one(), |product, factor| {
                Ok(product * factor.map(assignments, map)?)
            }),
            Self::MultilinearFunction { symbol, arguments } => {
                let mut builder = FunctionBuilder::new(*symbol);
                for argument in arguments {
                    let mapped = argument.map(assignments, map)?;
                    builder = builder.add_arg(mapped.as_view());
                }
                Ok(builder.finish())
            }
        }
    }
}

#[derive(Debug, Error)]
pub enum EnergyPowerAnalysisError {
    #[error("numerator energy degree exceeds the supported integer range")]
    EnergyDegreeOverflow,
    #[error("physical EMR edge {edge} has more than one exact-energy candidate set")]
    DuplicatePhysicalEnergyCandidateSet { edge: usize },
    #[error("physical EMR edge {edge} has an empty exact-energy candidate set")]
    EmptyEquivalentEnergyCandidateSet { edge: usize },
    #[error(
        "physical EMR edge {edge} repeats exact-energy candidate {candidate} in one candidate set"
    )]
    DuplicateEquivalentEnergyCandidate { edge: usize, candidate: usize },
    #[error(
        "exact-energy candidate {candidate} is shared by physical EMR edges {first_edge} and {second_edge}; minimax assignment requires disjoint certified candidate sets"
    )]
    OverlappingEquivalentEnergyCandidates {
        candidate: usize,
        first_edge: usize,
        second_edge: usize,
    },
    #[error(
        "physical EMR energy-power degree {degree} for edge {edge} has no certified equivalent exact-energy candidates"
    )]
    MissingEquivalentEnergyCandidates { edge: usize, degree: usize },
    #[error(
        "factor {factor_id} (`{expression}`) depends on physical EMR edge {edge} but has no exact-energy assignment"
    )]
    MissingFactorEnergyAssignment {
        factor_id: usize,
        expression: String,
        edge: usize,
    },
    #[error(
        "positive denominator wrapper `{expression}` is owned by edge {owner} but its active energy depends on edge {dependent}"
    )]
    NonlocalPositiveDenominator {
        expression: String,
        owner: usize,
        dependent: usize,
    },
    #[error(
        "positive denominator wrapper `{expression}` mixes incompatible energy-provenance families for owner {owner}"
    )]
    MixedPositiveDenominatorProvenance { expression: String, owner: usize },
    #[error(
        "physical-EMR energy-power analysis cannot assign loop-basis expression `{expression}` to a source edge; normalize it to Q(edge, index) before CFF analysis"
    )]
    LoopMomentumInPhysicalEmrAnalysis { expression: String },
    #[error("loop momentum K({loop_index}, ...) has no carrier edge in this loop momentum basis")]
    MissingLoopCarrierEdge { loop_index: usize },
    #[error(
        "energy-power analysis expected an integer edge id in Q(edge, index), found `{argument}`"
    )]
    InvalidEmrEdgeArgument { argument: String },
    #[error("energy-power analysis found malformed UV momentum provenance `{argument}`")]
    InvalidUvMomentumProvenance { argument: String },
    #[error(
        "energy-power analysis expected an integer loop id in K(loop, index), found `{argument}`"
    )]
    InvalidLoopMomentumArgument { argument: String },
    #[error(
        "non-polynomial energy dependence encountered in numerator subexpression `{expression}` with exponent `{exponent}`"
    )]
    NonPolynomialEnergyPower {
        expression: String,
        exponent: String,
    },
    #[error(
        "energy dependence encountered in function without declared multilinearity: `{expression}`"
    )]
    OpaqueEnergyFunction { expression: String },
}

pub struct EnergyPowerAnalyzer {
    loop_edges: Option<Vec<EdgeIndex>>,
    internal_edges: Option<BTreeSet<EdgeIndex>>,
    minkowski_symbol: Symbol,
}

#[derive(Clone, Copy)]
enum EnergyPowerAnalysisMode {
    StrictPolynomial,
    ConservativeUpperBound,
}

#[derive(Clone, Copy)]
enum LorentzIndexKind {
    Temporal,
    Spatial,
    Abstract,
}

impl EnergyPowerAnalyzer {
    pub fn new(loop_edges: impl IntoIterator<Item = EdgeIndex>) -> Self {
        Self {
            loop_edges: Some(loop_edges.into_iter().collect()),
            internal_edges: None,
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn with_internal_edges(
        loop_edges: impl IntoIterator<Item = EdgeIndex>,
        internal_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Self {
        Self {
            loop_edges: Some(loop_edges.into_iter().collect()),
            internal_edges: Some(internal_edges.into_iter().collect()),
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn for_physical_emr_edges(internal_edges: impl IntoIterator<Item = EdgeIndex>) -> Self {
        Self {
            loop_edges: None,
            internal_edges: Some(internal_edges.into_iter().collect()),
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn analyze_atom(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::StrictPolynomial,
        )
    }

    pub fn analyze_atom_upper_bound(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::ConservativeUpperBound,
        )
    }

    #[cfg(test)]
    pub(crate) fn plan_atom_assignment(
        &self,
        expression: &Atom,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        Ok(self
            .plan_atom_assignment_proposals(expression, candidates)?
            .remove(0))
    }

    pub(crate) fn plan_atom_assignment_proposals(
        &self,
        expression: &Atom,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<Vec<EnergyPowerAssignmentPlan>, EnergyPowerAnalysisError> {
        let expected_degrees = self.analyze_atom(expression)?;
        let mut next_factor_id = 0;
        let planned = self.plan_view(expression.as_view(), &mut next_factor_id)?;
        let planned_degrees = planned.degrees()?;
        debug_assert_eq!(planned_degrees, expected_degrees);

        planned
            .assignment_proposals(candidates)?
            .into_iter()
            .map(|factor_assignments| {
                let exact_bounds = planned.exact_degrees(&factor_assignments)?;
                Ok(EnergyPowerAssignmentPlan {
                    expression: planned.clone(),
                    factor_assignments,
                    energy_degree_bounds: exact_bounds.into_iter().collect(),
                })
            })
            .collect()
    }

    fn plan_view(
        &self,
        expression: AtomView<'_>,
        next_factor_id: &mut usize,
    ) -> Result<PlannedEnergyExpression, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Add(add) => Ok(PlannedEnergyExpression::Add(
                add.iter()
                    .map(|term| self.plan_view(term, next_factor_id))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
            AtomView::Mul(mul) => Ok(PlannedEnergyExpression::Mul(
                mul.iter()
                    .map(|factor| self.plan_view(factor, next_factor_id))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let base_degrees =
                    self.analyze_view(base, EnergyPowerAnalysisMode::StrictPolynomial)?;
                if base_degrees.is_empty() {
                    return self.plan_factor(expression, next_factor_id);
                }
                if !self
                    .analyze_view(exponent, EnergyPowerAnalysisMode::StrictPolynomial)?
                    .is_empty()
                {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                }
                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                };
                if exponent < 0 {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_string(),
                    });
                }
                if exponent == 0 {
                    return self.plan_factor(expression, next_factor_id);
                }
                Ok(PlannedEnergyExpression::Mul(
                    (0..exponent)
                        .map(|_| self.plan_view(base, next_factor_id))
                        .collect::<Result<Vec<_>, _>>()?,
                ))
            }
            AtomView::Fun(function)
                if function.get_nargs() == 4 && function.get_symbol() == GS.den =>
            {
                let owner = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
                    EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                        argument: function.get(0).to_owned().log_print(None),
                    }
                })?);
                let expression = expression.to_owned();
                let degrees =
                    self.analyze_view(function.get(3), EnergyPowerAnalysisMode::StrictPolynomial)?;
                let mut value_families = BTreeSet::new();
                let mut family_error = None;
                let _ = function.get(3).to_owned().replace_map(|view, _, _| {
                    if family_error.is_some() {
                        return;
                    }
                    let AtomView::Fun(momentum) = view else {
                        return;
                    };
                    if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
                        return;
                    }
                    let family = match self.emr_candidate_family(momentum.get(0)) {
                        Ok(family) => family,
                        Err(error) => {
                            family_error = Some(error);
                            return;
                        }
                    };
                    if self.is_internal_edge(family.edge())
                        && !self
                            .component_degree(family.edge(), momentum.get(1))
                            .is_empty()
                    {
                        value_families.insert(family);
                    }
                });
                if let Some(error) = family_error {
                    return Err(error);
                }
                if value_families.len() > 1 {
                    return Err(
                        EnergyPowerAnalysisError::MixedPositiveDenominatorProvenance {
                            expression: expression.log_print(None),
                            owner: owner.into(),
                        },
                    );
                }
                let value_family = value_families.into_iter().next();
                let mut families = BTreeMap::new();
                for (dependent, _) in degrees.iter() {
                    if dependent != owner {
                        return Err(EnergyPowerAnalysisError::NonlocalPositiveDenominator {
                            expression: expression.log_print(None),
                            owner: owner.into(),
                            dependent: dependent.into(),
                        });
                    }
                    let family = value_family.ok_or_else(|| {
                        EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                            argument: expression.log_print(None),
                        }
                    })?;
                    families.insert(dependent, family);
                }
                let id = *next_factor_id;
                *next_factor_id += 1;
                // Conservatively keep a positive common-denominator
                // completion as one assignment unit. This gives its temporal
                // and spatial/mass pieces one occurrence frame without
                // expanding the wrapper. Separate wrappers may still be
                // load-balanced over serial copies of the same line, and the
                // mapper restores each selected routing sign.
                Ok(PlannedEnergyExpression::Factor {
                    id,
                    expression,
                    degrees,
                    families,
                })
            }
            AtomView::Fun(function)
                if function.get_symbol() == GS.dot
                    || (function.get_symbol().is_linear()
                        && function.get_symbol() != GS.emr_mom
                        && function.get_symbol() != GS.loop_mom
                        && function.get_symbol() != GS.emr_vec) =>
            {
                Ok(PlannedEnergyExpression::MultilinearFunction {
                    symbol: function.get_symbol(),
                    arguments: function
                        .iter()
                        .map(|argument| self.plan_view(argument, next_factor_id))
                        .collect::<Result<Vec<_>, _>>()?,
                })
            }
            _ => self.plan_factor(expression, next_factor_id),
        }
    }

    fn plan_factor(
        &self,
        expression: AtomView<'_>,
        next_factor_id: &mut usize,
    ) -> Result<PlannedEnergyExpression, EnergyPowerAnalysisError> {
        let id = *next_factor_id;
        *next_factor_id += 1;
        let degrees = self.analyze_view(expression, EnergyPowerAnalysisMode::StrictPolynomial)?;
        let mut families = BTreeMap::new();
        if !degrees.is_empty() {
            let AtomView::Fun(momentum) = expression else {
                unreachable!("factorized energy dependence must be one momentum component");
            };
            debug_assert_eq!(momentum.get_symbol(), GS.emr_mom);
            let family = self.emr_candidate_family(momentum.get(0))?;
            families.insert(family.edge(), family);
        }
        Ok(PlannedEnergyExpression::Factor {
            id,
            expression: expression.to_owned(),
            degrees,
            families,
        })
    }

    fn analyze_view(
        &self,
        expression: AtomView<'_>,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Num(_) | AtomView::Var(_) => Ok(EnergyPowerCapMap::default()),
            AtomView::Add(add) => {
                let mut max_degree = EnergyPowerCapMap::default();
                for term in add.iter() {
                    max_degree.max_assign(self.analyze_view(term, mode)?);
                }
                Ok(max_degree)
            }
            AtomView::Mul(mul) => {
                let mut product_degree = EnergyPowerCapMap::default();
                for factor in mul.iter() {
                    product_degree.add_assign(self.analyze_view(factor, mode)?)?;
                }
                Ok(product_degree)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut base_degree = self.analyze_view(base, mode)?;
                if !self.analyze_view(exponent, mode)?.is_empty() {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                }
                if base_degree.is_empty() {
                    return Ok(base_degree);
                }

                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                };
                if exponent < 0 {
                    if matches!(mode, EnergyPowerAnalysisMode::ConservativeUpperBound) {
                        return Ok(EnergyPowerCapMap::default());
                    }
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_string(),
                    });
                }
                base_degree.scale(
                    usize::try_from(exponent)
                        .map_err(|_| EnergyPowerAnalysisError::EnergyDegreeOverflow)?,
                )?;
                Ok(base_degree)
            }
            AtomView::Fun(function) => {
                if function.get_symbol() == GS.emr_vec {
                    return Ok(EnergyPowerCapMap::default());
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.emr_mom {
                    let edge = self.emr_candidate_family(function.get(0))?.edge();
                    return Ok(if self.is_internal_edge(edge) {
                        self.component_degree(edge, function.get(1))
                    } else {
                        EnergyPowerCapMap::default()
                    });
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.loop_mom {
                    let Some(loop_edges) = &self.loop_edges else {
                        return Err(
                            EnergyPowerAnalysisError::LoopMomentumInPhysicalEmrAnalysis {
                                expression: Atom::from(function.to_owned()).log_print(None),
                            },
                        );
                    };
                    let loop_index = usize::try_from(function.get(0)).map_err(|_| {
                        EnergyPowerAnalysisError::InvalidLoopMomentumArgument {
                            argument: function.get(0).to_owned().log_print(None),
                        }
                    })?;
                    let edge = loop_edges
                        .get(loop_index)
                        .copied()
                        .ok_or(EnergyPowerAnalysisError::MissingLoopCarrierEdge { loop_index })?;
                    return Ok(if self.is_internal_edge(edge) {
                        self.component_degree(edge, function.get(1))
                    } else {
                        EnergyPowerCapMap::default()
                    });
                }

                // `den(edge, momentum, mass, full_expr)` is a provenance
                // wrapper whose algebraic value is its final argument. A
                // positive wrapper can remain in a factorized numerator after
                // a local 4D Taylor expansion, so analyze that polynomial
                // without also counting its edge/momentum metadata.
                if function.get_nargs() == 4 && function.get_symbol() == GS.den {
                    return self.analyze_view(function.get(3), mode);
                }

                // Dot products are bilinear. Keep that composition explicit:
                // the Symbolica symbol is intentionally not globally declared
                // linear, while numerator energy degrees still add between its
                // two vector slots.
                if function.get_nargs() == 2 && function.get_symbol() == GS.dot {
                    let mut degree = self.analyze_view(function.get(0), mode)?;
                    degree.add_assign(self.analyze_view(function.get(1), mode)?)?;
                    return Ok(degree);
                }

                let argument_degrees = function
                    .iter()
                    .map(|argument| self.analyze_view(argument, mode))
                    .collect::<Result<Vec<_>, _>>()?;
                if argument_degrees.iter().any(|degree| !degree.is_empty())
                    && !function.get_symbol().is_linear()
                {
                    return Err(EnergyPowerAnalysisError::OpaqueEnergyFunction {
                        expression: Atom::from(function.to_owned()).log_print(None),
                    });
                }
                if function.get_symbol().is_linear() {
                    let mut degree = EnergyPowerCapMap::default();
                    for argument_degree in argument_degrees {
                        degree.add_assign(argument_degree)?;
                    }
                    Ok(degree)
                } else {
                    Ok(EnergyPowerCapMap::default())
                }
            }
        }
    }

    fn component_degree(&self, edge: EdgeIndex, index: AtomView<'_>) -> EnergyPowerCapMap {
        match self.lorentz_index_kind(index) {
            LorentzIndexKind::Temporal | LorentzIndexKind::Abstract => {
                EnergyPowerCapMap::unit(edge)
            }
            LorentzIndexKind::Spatial => EnergyPowerCapMap::default(),
        }
    }

    fn emr_candidate_family(
        &self,
        argument: AtomView<'_>,
    ) -> Result<EnergyCandidateFamily, EnergyPowerAnalysisError> {
        if let Ok(edge) = usize::try_from(argument) {
            return Ok(EnergyCandidateFamily::Unprovenanced(EdgeIndex(edge)));
        }

        if let Some((edge, role, _)) = GS.uv_momentum_provenance_data(argument) {
            return Ok(match role {
                UvMomentumProvenanceRole::DenominatorDerived => {
                    EnergyCandidateFamily::DenominatorDerived(edge)
                }
                UvMomentumProvenanceRole::TaylorFixed
                | UvMomentumProvenanceRole::PhysicalSourceFixed => {
                    EnergyCandidateFamily::Fixed(edge)
                }
                UvMomentumProvenanceRole::DenominatorDerivedSoft => {
                    EnergyCandidateFamily::Unprovenanced(edge)
                }
            });
        }
        match argument {
            AtomView::Fun(provenance) if provenance.get_symbol() == GS.uv_momentum_provenance => {
                Err(EnergyPowerAnalysisError::InvalidUvMomentumProvenance {
                    argument: argument.to_owned().log_print(None),
                })
            }
            _ => Err(EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                argument: argument.to_owned().log_print(None),
            }),
        }
    }

    fn lorentz_index_kind(&self, index: AtomView<'_>) -> LorentzIndexKind {
        if let Ok(index) = usize::try_from(index) {
            return if index == 0 {
                LorentzIndexKind::Temporal
            } else {
                LorentzIndexKind::Spatial
            };
        }

        let AtomView::Fun(function) = index else {
            return LorentzIndexKind::Spatial;
        };

        if function.get_symbol() == AIND_SYMBOLS.cind {
            return if function.get_nargs() == 1
                && usize::try_from(function.get(0)).is_ok_and(|index| index == 0)
            {
                LorentzIndexKind::Temporal
            } else {
                LorentzIndexKind::Spatial
            };
        }

        if function.get_symbol() == self.minkowski_symbol && matches!(function.get_nargs(), 1 | 2) {
            // `mink(4)` is the compact full-vector slot and `mink(4, n)` is an
            // abstract Lorentz-index label in the numerator algebra. Neither
            // selects a concrete spatial component; those use `cind(n)` and are
            // handled above. Both forms can contract onto the temporal
            // component and therefore contribute to the energy-power bound.
            return LorentzIndexKind::Abstract;
        }

        LorentzIndexKind::Spatial
    }

    fn is_internal_edge(&self, edge: EdgeIndex) -> bool {
        self.internal_edges
            .as_ref()
            .is_none_or(|internal_edges| internal_edges.contains(&edge))
    }
}

impl Graph {
    #[cfg(test)]
    pub(crate) fn optimize_soft_momentum_routing(
        &self,
        numerator: &Atom,
        active_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> color_eyre::Result<Atom> {
        Ok(self
            .soft_momentum_routing_proposals(numerator, active_edges)?
            .remove(0))
    }

    /// Route only denominator-derived soft momenta once the untouched outer
    /// numerator is present. Original factors retain their literal owners.
    /// Every candidate is an exact off-shell linear identity in the production
    /// LMB, including the fixed external shift; on-shell energies play no role.
    pub(crate) fn soft_momentum_routing_proposals(
        &self,
        numerator: &Atom,
        active_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> color_eyre::Result<Vec<Atom>> {
        if !numerator.contains_symbol(GS.uv_momentum_provenance) {
            return Ok(vec![numerator.clone()]);
        }
        let mut soft_components = BTreeSet::new();
        let _ = numerator.replace_map(|view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() == 2
                && let Some((_, UvMomentumProvenanceRole::DenominatorDerivedSoft, _)) =
                    GS.uv_momentum_provenance_data(momentum.get(0))
            {
                soft_components.insert(view.to_owned());
                **output = view.to_owned();
            }
        });
        if soft_components.is_empty() {
            return Ok(vec![numerator.clone()]);
        }

        let active_edges = active_edges.into_iter().collect::<BTreeSet<_>>();
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(active_edges.iter().copied());
        let lmb = &self.loop_momentum_basis;
        let row = |edge: EdgeIndex| {
            lmb.edge_signatures[edge]
                .internal
                .iter()
                .map(|sign| match sign {
                    SignOrZero::Minus => -1,
                    SignOrZero::Zero => 0,
                    SignOrZero::Plus => 1,
                })
                .collect::<Vec<i64>>()
        };
        let rows = active_edges
            .iter()
            .map(|edge| (*edge, row(*edge)))
            .filter(|(_, row)| row.iter().any(|coefficient| *coefficient != 0))
            .collect::<Vec<_>>();
        let rank = rank_i64(&rows.iter().map(|(_, row)| row.clone()).collect::<Vec<_>>());
        let mut alternatives = BTreeMap::new();
        let mut planning_edges = active_edges.clone();
        for component in soft_components {
            let AtomView::Fun(momentum) = component.as_view() else {
                unreachable!()
            };
            let (owner, _, payload) = GS.uv_momentum_provenance_data(momentum.get(0)).unwrap();
            planning_edges.insert(owner);
            let original = GS.erase_uv_momentum_provenance(&component);
            // Soft tags name actual crown carriers. A nested Taylor operation
            // consumes its new hard part separately before reaching this stage.
            if payload
                != FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(owner))
                    .finish()
            {
                return Err(eyre::eyre!(
                    "soft momentum provenance does not retain its literal crown carrier"
                ));
            }
            let target = row(owner);
            let indices = [momentum.get(1).to_owned()];
            let spatial = matches!(
                analyzer.lorentz_index_kind(momentum.get(1)),
                LorentzIndexKind::Spatial
            );
            let fixed_external = target.iter().all(|coefficient| *coefficient == 0);
            let mut candidates = BTreeSet::new();
            if spatial || fixed_external || active_edges.contains(&owner) {
                candidates.insert(original);
            }
            if spatial {
                // Concrete spatial components carry no energy rank and keep
                // their original routing.
            } else if fixed_external {
                candidates.insert(lmb.ext_atom(owner, GS.emr_mom, &indices, true));
            } else {
                // Every minimal support is linearly independent and extends to
                // a basis of the available row space. Enumerating those bases
                // therefore includes every nondominated energy support.
                for selected in rows.iter().combinations(rank) {
                    for columns in (0..target.len()).combinations(rank) {
                        let matrix = columns
                            .iter()
                            .map(|axis| {
                                selected
                                    .iter()
                                    .map(|(_, row)| Rational::from(row[*axis]))
                                    .collect()
                            })
                            .collect();
                        let rhs = columns
                            .iter()
                            .map(|axis| Rational::from(target[*axis]))
                            .collect();
                        let Some(coefficients) = solve_rational_system(matrix, rhs) else {
                            continue;
                        };
                        if !(0..target.len()).all(|axis| {
                            selected.iter().zip(&coefficients).fold(
                                Rational::from(0),
                                |sum, ((_, row), coefficient)| {
                                    sum + coefficient * &Rational::from(row[axis])
                                },
                            ) == target[axis]
                        }) {
                            continue;
                        }
                        let mut candidate = lmb.ext_atom(owner, GS.emr_mom, &indices, true);
                        for ((edge, _), coefficient) in selected.iter().zip(coefficients) {
                            if coefficient != 0 {
                                candidate += Atom::num(coefficient)
                                    * (GS.emr_mom(*edge, &indices[0])
                                        - lmb.ext_atom(*edge, GS.emr_mom, &indices, true));
                            }
                        }
                        candidates.insert(candidate);
                        // Different invertible minors of one basis give the
                        // same unique exact solution.
                        break;
                    }
                }
            }
            if candidates.is_empty() {
                return Err(eyre::eyre!(
                    "soft energy on crown edge {owner} has no exact routing through the active outer CFF edges {active_edges:?}"
                ));
            }
            alternatives.insert(
                component,
                candidates
                    .into_iter()
                    .map(|candidate| Ok((candidate.clone(), analyzer.analyze_atom(&candidate)?)))
                    .collect::<Result<Vec<_>, EnergyPowerAnalysisError>>()?,
            );
        }

        // An unavailable crown carrier must still expose its power factors to
        // the planner; treating it as constant would silently lose the ranks
        // of the active edges onto which its energy has to be routed. Fixed
        // loads on these extra axes are identical in every assignment and do
        // not enter the final active-edge objective.
        let planned = EnergyPowerAnalyzer::for_physical_emr_edges(planning_edges)
            .plan_view(numerator.as_view(), &mut 0)?;
        let mut proposals = planned.soft_momentum_assignments(&alternatives)?;
        proposals.sort_by_cached_key(|(degrees, assignments)| {
            let mut envelope = active_edges
                .iter()
                .map(|edge| degrees.degrees.get(edge).copied().unwrap_or(0))
                .collect::<Vec<_>>();
            envelope.sort_unstable_by(|left, right| right.cmp(left));
            (envelope, assignments.clone())
        });
        // The bounded search uses the certified basis/Pareto class only to
        // propose candidates. Actual native map count selects between them before
        // graph surface conversion; no global count optimum is claimed.
        proposals
            .into_iter()
            .unique_by(|(degrees, _)| {
                active_edges
                    .iter()
                    .map(|edge| degrees.degrees.get(edge).copied().unwrap_or(0))
                    .collect::<Vec<_>>()
            })
            .take(ENERGY_ASSIGNMENT_PROPOSAL_BUDGET)
            .map(|(_, assignments)| {
                let optimized = planned.map(&BTreeMap::new(), &mut |id, factor, _| {
                    Ok::<_, EnergyPowerAnalysisError>(assignments.get(&id).map_or_else(
                        || factor.clone(),
                        |choice| alternatives[factor][*choice].0.clone(),
                    ))
                })?;
                // Energy-independent powers/functions need no assignment traversal.
                // Erase their remaining soft metadata without changing their routing.
                Ok(optimized.replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() == 2
                        && let Some((_, UvMomentumProvenanceRole::DenominatorDerivedSoft, _)) =
                            GS.uv_momentum_provenance_data(momentum.get(0))
                    {
                        **output = GS.erase_uv_momentum_provenance(&view.to_owned());
                    }
                }))
            })
            .collect()
    }

    pub fn full_numerator_atom(&self) -> Atom {
        self.numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .expect("Graph numerator should be available")
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor
    }

    pub fn numerator_energy_power_caps(
        &self,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_numerator_energy_powers(EnergyPowerAnalysisMode::StrictPolynomial)
    }

    pub(crate) fn numerator_energy_power_upper_bounds(
        &self,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_numerator_energy_powers(EnergyPowerAnalysisMode::ConservativeUpperBound)
    }

    fn analyze_numerator_energy_powers(
        &self,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        // Keep edge-energy variables attached to their source denominator.
        // Rewriting them to a loop-momentum basis loses cancellation
        // information in non-coordinate UV directions, e.g. when q_i-q_j is
        // kept fixed while q_i and q_j are both large.
        let numerator = self.full_numerator_atom();
        self.analyze_numerator_atom_energy_powers(&numerator, mode)
    }

    fn analyze_numerator_atom_energy_powers(
        &self,
        numerator: &Atom,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        let internal_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy).then_some(edge)
            });

        EnergyPowerAnalyzer::for_physical_emr_edges(internal_edges)
            .analyze_view(numerator.as_view(), mode)
    }

    pub fn automatic_numerator_energy_degree_bounds(
        &self,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        Ok(self
            .numerator_energy_power_upper_bounds()?
            .into_generation_bounds())
    }

    pub(crate) fn automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
        &self,
        numerators: impl IntoIterator<Item = impl Borrow<Atom>>,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        min_degree: usize,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(active_edges);
        // Branches are evaluated independently. Take the maximum of their
        // ranks without constructing a sum in which leading powers could cancel.
        let mut bounds = EnergyPowerCapMap::default();
        for numerator in numerators {
            bounds.max_assign(analyzer.analyze_atom(numerator.borrow())?);
        }
        Ok(bounds
            .iter()
            .filter_map(|(edge, degree)| (degree >= min_degree).then_some((edge.into(), degree)))
            .collect())
    }

    #[cfg(test)]
    pub(crate) fn plan_numerator_energy_assignment_in_atom_excluding(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        Ok(self
            .plan_numerator_energy_assignment_proposals_in_atom_excluding(
                numerator,
                excluded_edges,
                candidates,
            )?
            .remove(0))
    }

    pub(crate) fn plan_numerator_energy_assignment_proposals_in_atom_excluding(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<Vec<EnergyPowerAssignmentPlan>, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        EnergyPowerAnalyzer::for_physical_emr_edges(active_edges)
            .plan_atom_assignment_proposals(numerator, candidates)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::{
        dot,
        graph::{Graph, LMBext, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::energy_degree::{
            EnergyCandidateFamily, EnergyPowerAnalysisError, EnergyPowerAnalyzer,
            EquivalentEnergyCandidates,
        },
        utils::{GS, W_, symbols::UvMomentumProvenanceRole},
    };
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::shadowing::symbolica_utils::ReplaceBuilderExt;
    use spenso::structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, Minkowski},
    };
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        function, symbol,
    };

    fn mink_index(label: &str) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::var(symbol!(label)).as_view())
            .finish()
    }

    fn mink_component(index: i64) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::num(index).as_view())
            .finish()
    }

    fn compact_minkowski_vector() -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .finish()
    }

    fn bounds(expression: Atom) -> Vec<(usize, usize)> {
        EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap()
            .into_generation_bounds()
    }

    fn bounds_with_internal_edges(expression: Atom) -> Vec<(usize, usize)> {
        EnergyPowerAnalyzer::with_internal_edges([EdgeIndex(7), EdgeIndex(9)], [EdgeIndex(3)])
            .analyze_atom(&expression)
            .unwrap()
            .into_generation_bounds()
    }

    #[test]
    fn emr_minkowski_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_index("mu"));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_numeric_minkowski_label_counts_as_abstract_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_component(1));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_numeric_zero_minkowski_label_counts_as_abstract_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_component(0));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn compact_emr_vector_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, compact_minkowski_vector());
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn dot_product_composes_compact_vector_energy_powers() {
        let q3 = function!(GS.emr_mom, 3, compact_minkowski_vector());
        let q7 = function!(GS.emr_mom, 7, compact_minkowski_vector());
        assert_eq!(bounds(function!(GS.dot, q3.clone(), q3)), vec![(3, 2)]);
        assert_eq!(bounds(function!(GS.dot, q7, q3)), vec![(3, 1), (7, 1)]);
    }

    #[test]
    fn positive_denominator_wrapper_analyzes_only_its_polynomial_value() {
        let q3 = function!(GS.emr_mom, 3, compact_minkowski_vector());
        let full_expr = function!(GS.dot, q3.clone(), q3.clone()) - Atom::one();
        let denominator = GS.den(3, q3, Atom::one(), full_expr);
        assert_eq!(bounds(denominator), vec![(3, 2)]);
    }

    #[test]
    fn emr_external_edge_minkowski_index_does_not_count_when_internal_edges_are_known() {
        let expression = function!(GS.emr_mom, 0, mink_index("mu"));
        assert!(bounds_with_internal_edges(expression).is_empty());
    }

    #[test]
    fn emr_temporal_concrete_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.call(0));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_spatial_concrete_index_does_not_count() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.call(2));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn emr_direct_temporal_component_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, 0);
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_direct_spatial_component_does_not_count() {
        let expression = function!(GS.emr_mom, 3, 1);
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn spatial_emr_vector_does_not_count() {
        let expression = function!(GS.emr_vec, 3, mink_index("mu"));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn lmb_energy_maps_to_carrier_edge() {
        let expression = function!(GS.loop_mom, 1, mink_index("mu"));
        assert_eq!(bounds(expression), vec![(9, 1)]);
    }

    #[test]
    fn physical_emr_analysis_rejects_loop_basis_energy() {
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(7)]);
        let loop_momentum = function!(GS.loop_mom, 0, mink_index("mu"));
        assert!(matches!(
            analyzer.analyze_atom(&loop_momentum),
            Err(EnergyPowerAnalysisError::LoopMomentumInPhysicalEmrAnalysis { .. })
        ));

        let emr_momentum = function!(GS.emr_mom, 7, mink_index("mu"));
        assert_eq!(
            analyzer
                .analyze_atom(&emr_momentum)
                .unwrap()
                .into_generation_bounds(),
            vec![(7, 1)]
        );
    }

    #[test]
    fn sums_take_edgewise_max_and_products_add() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let k1 = function!(GS.loop_mom, 1, mink_index("nu"));
        let expression = q3.clone() * k1.clone() + q3.pow(Atom::num(3));
        assert_eq!(bounds(expression), vec![(3, 3), (9, 1)]);
    }

    #[test]
    fn powers_of_affine_atoms_produce_edgewise_bounds() {
        let k0 = function!(GS.loop_mom, 0, mink_index("mu"));
        let k1 = function!(GS.loop_mom, 1, mink_index("nu"));
        let expression = (k0 + k1).pow(Atom::num(2));
        assert_eq!(bounds(expression), vec![(7, 2), (9, 2)]);
    }

    #[test]
    fn factorized_products_compose_energy_bounds_without_expansion() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression = (q3.clone() + Atom::num(1)) * (q3 + Atom::num(2));
        assert_eq!(bounds(expression), vec![(3, 2)]);
    }

    fn two_equivalent_candidates() -> EquivalentEnergyCandidates {
        EquivalentEnergyCandidates::try_from_source_occurrences([(EdgeIndex(3), vec![10, 11])])
            .unwrap()
    }

    fn uv_owned_q3(denominator_derived: bool, index: Atom) -> Atom {
        let hard_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_momentum_provenance_tag(
                Atom::num(3).as_view(),
                denominator_derived,
                hard_momentum.as_view(),
            ))
            .add_arg(index.as_view())
            .finish()
    }

    fn production_candidates_with_two_derived_copies() -> EquivalentEnergyCandidates {
        EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            EdgeIndex(3),
            10,
            vec![11, 12],
        )])
        .unwrap()
    }

    fn planned_candidates(expression: &Atom, candidates: &EquivalentEnergyCandidates) {
        let edge = EdgeIndex(3);
        let plans = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment_proposals(expression, candidates)
            .unwrap();
        assert!(!plans.is_empty());
        for plan in plans {
            let mapped = plan.map_factors(|factor, assignments| {
                let Some(candidate) = assignments.get(&edge) else {
                    return Ok::<_, ()>(factor.clone());
                };
                let _ = factor.replace_map(|view, _, _| {
                    let AtomView::Fun(momentum) = view else { return; };
                    if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 2 {
                        return;
                    }
                    let family = match GS.uv_momentum_provenance_data(momentum.get(0)) {
                        Some((owner, UvMomentumProvenanceRole::DenominatorDerived, _)) =>
                            EnergyCandidateFamily::DenominatorDerived(owner),
                        Some((owner, _, _)) => EnergyCandidateFamily::Fixed(owner),
                        None => EnergyCandidateFamily::Unprovenanced(edge),
                    };
                    assert!(candidates.get(family).unwrap().contains(candidate),
                        "every original or derivative-owned factor must stay in its certified occurrence domain");
                });
                Ok(GS.erase_uv_momentum_provenance(factor)
                    .replace(function!(GS.emr_mom, 3, W_.x___))
                    .with(function!(GS.emr_mom, *candidate, W_.x___)))
            }).unwrap();
            let bounds = plan
                .energy_degree_bounds()
                .iter()
                .copied()
                .collect::<BTreeMap<_, _>>();
            let actual = EnergyPowerAnalyzer::new([])
                .analyze_atom(&mapped)
                .unwrap()
                .into_generation_bounds();
            assert!(
                actual
                    .iter()
                    .all(|(edge, degree)| bounds.get(edge).is_some_and(|bound| bound >= degree)),
                "the returned capacity must certify the fully mapped numerator"
            );
            let mut reconstructed = mapped;
            for candidate in bounds.keys() {
                reconstructed = reconstructed
                    .replace(function!(GS.emr_mom, *candidate, W_.x___))
                    .with(function!(GS.emr_mom, 3, W_.x___));
            }
            assert!(
                (reconstructed - GS.erase_uv_momentum_provenance(expression))
                    .expand()
                    .is_zero(),
                "every proposed assignment must reconstruct the complete numerator value: {expression}"
            );
        }
    }

    #[test]
    fn assignment_plans_preserve_owner_domains_and_complete_factorized_values() {
        let fixed = uv_owned_q3(false, GS.cind(0));
        let derived = uv_owned_q3(true, GS.cind(0));
        let plain = function!(GS.emr_mom, 3, GS.cind(0));
        let vector = uv_owned_q3(true, compact_minkowski_vector());
        let a = Atom::var(symbol!("energy_assignment_test::a"));
        let b = Atom::var(symbol!("energy_assignment_test::b"));
        let linear = symbol!("energy_assignment_test::linear"; Linear);
        let denominator = function!(GS.den, 3, function!(GS.emr_mom, 3), 1, derived.pow(2));
        // The quadratic wrapper stays indivisible even below a fixed quartic.
        // The complete certificate permits any valid distribution of the
        // derivative-created factors, including additive capacity reuse.
        let expressions = [
            (&plain + Atom::one()) * (&plain + Atom::num(2)),
            function!(GS.dot, &vector, &vector),
            function!(linear, &derived) + function!(linear, &derived * &a),
            (&derived + &a).pow(2),
            denominator.pow(2) + derived.pow(2),
            fixed.pow(4) * &denominator * derived.pow(2),
            fixed.pow(4) * derived.pow(2),
            (&derived + Atom::one()) * (&derived + Atom::num(2)) * (&derived + Atom::num(3)),
            &fixed * &derived,
            &fixed * (&a * &derived + &b * &derived),
            &a * &fixed * &derived + &b * derived.pow(2),
            derived.pow(3) + derived.pow(2),
            derived.pow(2),
        ];
        for candidates in [
            two_equivalent_candidates(),
            production_candidates_with_two_derived_copies(),
            EquivalentEnergyCandidates::try_from_source_occurrences([(
                EdgeIndex(3),
                vec![10, 11, 12],
            )])
            .unwrap(),
        ] {
            for expression in &expressions {
                planned_candidates(expression, &candidates);
            }
        }
    }

    #[test]
    fn soft_taylor_routing_balances_the_outer_envelope_without_moving_originals()
    -> color_eyre::Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph soft_taylor_cograph {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b -> c [id=3 lmb_id=1]
            c -> d [id=4]
            d -> a [id=5]
            d -> outgoing [id=6]
        })?;
        let soft_component = |index: Atom| {
            let owner = Atom::num(4);
            let payload = function!(GS.emr_mom, 4);
            function!(
                GS.emr_mom,
                GS.uv_momentum_provenance_tag(
                    owner.as_view(),
                    UvMomentumProvenanceRole::DenominatorDerivedSoft,
                    payload.as_view(),
                ),
                index
            )
        };
        let soft = |index| soft_component(GS.cind(index));
        let fixed_quartic = GS.emr_mom(EdgeIndex(3), GS.cind(0)).pow(4);
        let fixed_linear = GS.emr_mom(EdgeIndex(4), GS.cind(0));
        let numerator = &fixed_quartic * &fixed_linear * soft(0).pow(3);
        let active = [EdgeIndex(3), EdgeIndex(4), EdgeIndex(5)];
        let optimized = graph.optimize_soft_momentum_routing(&numerator, active)?;
        assert_ne!(
            optimized.replace(fixed_quartic.clone()).with(Atom::one()),
            optimized,
            "the original quartic remains on its original edge",
        );
        let neutral = graph.normal_emr_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[W_.x___],
            |_| true,
        );
        assert!(
            (optimized.replace_multiple(&neutral).collect_factors()
                - GS.erase_uv_momentum_provenance(&numerator)
                    .replace_multiple(&neutral)
                    .collect_factors())
            .collect_factors()
            .is_zero(),
            "soft routing must be an exact off-shell identity"
        );
        assert!(!optimized.contains_symbol(GS.uv_momentum_provenance));
        assert_eq!(
            graph.optimize_soft_momentum_routing(&soft(1).pow(2), active)?,
            GS.emr_mom(EdgeIndex(4), GS.cind(1)).pow(2),
            "spatial powers lose metadata without being expanded or rerouted",
        );
        let surviving = [EdgeIndex(3), EdgeIndex(5)];
        let rerouted = graph.optimize_soft_momentum_routing(&soft(0).pow(2), surviving)?;
        assert!(
            !rerouted
                .replace(GS.emr_mom(EdgeIndex(4), GS.cind(0)))
                .matches(),
            "an unavailable crown carrier must be routed through active edges"
        );
        assert!(
            (rerouted.replace_multiple(&neutral).collect_factors()
                - GS.erase_uv_momentum_provenance(&soft(0).pow(2))
                    .replace_multiple(&neutral)
                    .collect_factors())
            .collect_factors()
            .is_zero()
        );
        assert!(
            graph
                .optimize_soft_momentum_routing(&soft(0), [EdgeIndex(1)])
                .unwrap_err()
                .to_string()
                .contains("no exact routing")
        );
        let vector_index = compact_minkowski_vector();
        let dot = function!(
            GS.dot,
            soft_component(vector_index.clone()),
            GS.emr_mom(EdgeIndex(0), vector_index)
        );
        // This diagnostic routing polynomial exercises abstract vector slots
        // and a factorized additive power without using numerator expansion.
        let factorized = &fixed_quartic * (dot + Atom::one()).pow(2);
        let routed = graph.optimize_soft_momentum_routing(&factorized, active)?;
        assert!(
            (routed.replace_multiple(&neutral).collect_factors()
                - GS.erase_uv_momentum_provenance(&factorized)
                    .replace_multiple(&neutral)
                    .collect_factors())
            .collect_factors()
            .is_zero(),
            "abstract-slot routing must obey the same exact linear identity"
        );
        // A dependent external carrier has two literal exact representations,
        // but both request the same all-zero active CFF capacity. A lone factor
        // bypasses composite-node collection; every retained proposal must still
        // preserve its complete external value.
        let external_owner = [EdgeIndex(0), EdgeIndex(6)]
            .into_iter()
            .find(|edge| {
                graph
                    .loop_momentum_basis
                    .ext_atom(*edge, GS.emr_mom, &[GS.cind(0)], true)
                    != GS.emr_mom(*edge, GS.cind(0))
            })
            .expect("one of the two external momenta is dependent");
        let owner = Atom::num(usize::from(external_owner));
        let payload = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(external_owner))
            .finish();
        let external_soft = function!(
            GS.emr_mom,
            GS.uv_momentum_provenance_tag(
                owner.as_view(),
                UvMomentumProvenanceRole::DenominatorDerivedSoft,
                payload.as_view(),
            ),
            GS.cind(0)
        );
        let literal = GS.erase_uv_momentum_provenance(&external_soft);
        let external =
            graph
                .loop_momentum_basis
                .ext_atom(external_owner, GS.emr_mom, &[GS.cind(0)], true);
        assert_ne!(
            literal, external,
            "the fixture must expose duplicate capacity through different atoms"
        );
        for source in [external_soft.clone(), (external_soft + Atom::one()).pow(8)] {
            let proposals = graph.soft_momentum_routing_proposals(&source, active)?;
            assert!(!proposals.is_empty());
            for proposal in proposals {
                assert!(
                    EnergyPowerAnalyzer::for_physical_emr_edges(active)
                        .analyze_atom(&proposal)?
                        .is_empty()
                );
                assert!(
                    (proposal.replace_multiple(&neutral).collect_factors()
                        - GS.erase_uv_momentum_provenance(&source)
                            .replace_multiple(&neutral)
                            .collect_factors())
                    .is_zero(),
                    "all zero-capacity routes must retain the complete factorized external coefficient"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn positive_denominator_completion_stays_in_one_occurrence_frame() {
        let q3 = uv_owned_q3(true, GS.cind(0));
        let spatial_mass = Atom::var(symbol!("energy_assignment_test::spatial_mass"));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg((q3.pow(2) - &spatial_mass).as_view())
            .finish();
        planned_candidates(&denominator, &two_equivalent_candidates());
    }

    #[test]
    fn positive_denominator_completion_obeys_embedded_emr_provenance() {
        let edge = EdgeIndex(3);
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = |energy: &Atom| {
            FunctionBuilder::new(GS.den)
                .add_arg(3)
                .add_arg(source_momentum.as_view())
                .add_arg(Atom::one().as_view())
                .add_arg(energy.clone().pow(2).as_view())
                .finish()
        };
        let assert_assignment = |energy: Atom, candidates: EquivalentEnergyCandidates| {
            planned_candidates(&denominator(&energy), &candidates);
        };

        let base_only = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            edge,
            10,
            Vec::new(),
        )])
        .unwrap();
        assert_assignment(function!(GS.emr_mom, 3, GS.cind(0)), base_only);

        assert_assignment(
            uv_owned_q3(false, GS.cind(0)),
            production_candidates_with_two_derived_copies(),
        );
        assert_assignment(
            uv_owned_q3(true, GS.cind(0)),
            production_candidates_with_two_derived_copies(),
        );
    }

    #[test]
    fn positive_denominator_completion_rejects_mixed_emr_provenance() {
        let edge = EdgeIndex(3);
        let fixed = uv_owned_q3(false, GS.cind(0));
        let derived = uv_owned_q3(true, GS.cind(0));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg((fixed * derived).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(
                &denominator,
                &production_candidates_with_two_derived_copies(),
            )
            .unwrap_err();

        assert!(matches!(
            error,
            EnergyPowerAnalysisError::MixedPositiveDenominatorProvenance { owner: 3, .. }
        ));
    }

    #[test]
    fn separate_positive_denominator_completions_may_use_separate_occurrences() {
        let q3 = uv_owned_q3(true, GS.cind(0));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg(q3.pow(2).as_view())
            .finish();
        let expression = denominator.pow(2);
        planned_candidates(&expression, &two_equivalent_candidates());
    }

    #[test]
    fn positive_denominator_completion_rejects_a_nonlocal_energy_owner() {
        let q4 = FunctionBuilder::new(GS.emr_mom)
            .add_arg(
                GS.uv_momentum_provenance_tag(
                    Atom::num(4).as_view(),
                    true,
                    FunctionBuilder::new(GS.emr_mom)
                        .add_arg(4)
                        .finish()
                        .as_view(),
                ),
            )
            .add_arg(GS.cind(0))
            .finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(
                FunctionBuilder::new(GS.emr_mom)
                    .add_arg(3)
                    .finish()
                    .as_view(),
            )
            .add_arg(Atom::one().as_view())
            .add_arg(q4.pow(2).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(3), EdgeIndex(4)])
            .plan_atom_assignment(&denominator, &two_equivalent_candidates())
            .unwrap_err();

        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonlocalPositiveDenominator {
                owner: 3,
                dependent: 4,
                ..
            }
        ));
    }

    #[test]
    fn source_energy_candidates_fix_originals_and_offer_all_derived_occurrences() {
        let candidates = EquivalentEnergyCandidates::try_from_source_occurrences([
            (EdgeIndex(3), vec![12, 10, 11]),
            (EdgeIndex(7), vec![22, 21, 20]),
            (EdgeIndex(9), vec![33, 31, 32, 30]),
        ])
        .unwrap();

        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(3))),
            Some([10, 11, 12].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(7))),
            Some([20, 21, 22].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(9))),
            Some([30, 31, 32, 33].as_slice())
        );
        for (edge, canonical) in [(EdgeIndex(3), 10), (EdgeIndex(7), 20), (EdgeIndex(9), 30)] {
            assert_eq!(
                candidates.get(EnergyCandidateFamily::Unprovenanced(edge)),
                Some([canonical].as_slice())
            );
            assert_eq!(
                candidates.get(EnergyCandidateFamily::Fixed(edge)),
                Some([canonical].as_slice())
            );
        }
    }

    #[test]
    fn source_energy_candidates_keep_occurrences_disjoint_between_owners() {
        let overlapping = EquivalentEnergyCandidates::try_from_source_occurrences([
            (EdgeIndex(3), vec![10, 11, 12]),
            (EdgeIndex(7), vec![12]),
        ])
        .unwrap_err();
        assert!(matches!(
            overlapping,
            EnergyPowerAnalysisError::OverlappingEquivalentEnergyCandidates { .. }
        ));

        let duplicate = EquivalentEnergyCandidates::try_from_source_occurrences([(
            EdgeIndex(3),
            vec![10, 11, 11],
        )])
        .unwrap_err();
        assert!(matches!(
            duplicate,
            EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate { .. }
        ));
    }

    #[test]
    fn partitioned_source_candidates_do_not_infer_the_base_from_energy_id_order() {
        let candidates = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([
            (EdgeIndex(3), 12, vec![10, 11]),
            (EdgeIndex(7), 20, Vec::new()),
        ])
        .unwrap();

        assert_eq!(
            candidates.get(EnergyCandidateFamily::Unprovenanced(EdgeIndex(3))),
            Some([12].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::Fixed(EdgeIndex(3))),
            Some([12].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(3))),
            Some([10, 11].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::Fixed(EdgeIndex(7))),
            Some([20].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(7))),
            None,
            "an owner without a derivative-created serial copy must not accept a derivative-owned numerator factor"
        );
    }

    #[test]
    fn partitioned_source_candidates_reject_derived_factors_without_a_copy() {
        let edge = EdgeIndex(3);
        let candidates = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            edge,
            12,
            Vec::new(),
        )])
        .unwrap();
        let derived = uv_owned_q3(true, mink_index("mu"));
        assert!(matches!(
            EnergyPowerAnalyzer::for_physical_emr_edges([edge])
                .plan_atom_assignment(&derived, &candidates),
            Err(EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates { .. })
        ));
    }

    #[test]
    fn assignment_rejects_a_missing_certified_candidate_set() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let candidates =
            EquivalentEnergyCandidates::try_from_source_occurrences([(EdgeIndex(7), vec![12])])
                .unwrap();
        let error = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(3)])
            .plan_atom_assignment(&q3, &candidates)
            .unwrap_err();

        assert!(matches!(
            error,
            EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates { edge: 3, degree: 1 }
        ));
    }

    #[test]
    fn overflowing_polynomial_degrees_return_an_error() {
        let energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let large_degree = usize::MAX / 3 + 1;
        let exponent = Atom::num(large_degree as u64);
        let powered = energy.clone().pow(exponent.clone());
        let overflow_product = (1..=3).fold(Atom::one(), |product, shift| {
            product * (energy.clone() + Atom::num(shift)).pow(exponent.clone())
        });
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(0)]);
        for expression in [(powered.clone() + Atom::one()).pow(3), overflow_product] {
            assert!(matches!(
                analyzer.analyze_atom(&expression),
                Err(EnergyPowerAnalysisError::EnergyDegreeOverflow)
            ));
            assert!(matches!(
                analyzer.analyze_atom_upper_bound(&expression),
                Err(EnergyPowerAnalysisError::EnergyDegreeOverflow)
            ));
        }
        assert_eq!(
            analyzer
                .analyze_atom(&(powered + Atom::one()).pow(2))
                .unwrap()
                .into_generation_bounds(),
            vec![(0, 2 * large_degree)]
        );
    }

    #[test]
    fn strict_analysis_rejects_negative_energy_powers() {
        let expression = function!(GS.emr_mom, 3, mink_index("mu")).pow(Atom::num(-1));
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonPolynomialEnergyPower { .. }
        ));
    }

    #[test]
    fn strict_analysis_ignores_energy_powers_outside_the_active_edge_set() {
        let analyzer =
            EnergyPowerAnalyzer::with_internal_edges([EdgeIndex(7), EdgeIndex(9)], [EdgeIndex(5)]);
        for expression in [
            function!(GS.emr_mom, 3, mink_index("mu")).pow(Atom::num(-1)),
            function!(GS.loop_mom, 0, mink_index("mu")).pow(Atom::num(-1)),
        ] {
            assert!(analyzer.analyze_atom(&expression).unwrap().is_empty());
        }
    }

    #[test]
    fn energy_dependent_exponents_are_rejected_even_for_constant_bases() {
        let exponent = function!(GS.emr_mom, 3, 0);
        let expression = Atom::num(2).pow(exponent);
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonPolynomialEnergyPower { .. }
        ));
    }

    #[test]
    fn strict_analysis_rejects_opaque_energy_functions() {
        let expression = FunctionBuilder::new(symbol!("opaque_energy_function"))
            .add_arg(function!(GS.emr_mom, 3, mink_index("mu")).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::OpaqueEnergyFunction { .. }
        ));
    }

    #[test]
    fn zero_arguments_do_not_hide_energy_in_opaque_functions() {
        let expression = FunctionBuilder::new(symbol!("opaque_energy_function_with_zero"))
            .add_arg(Atom::zero().as_view())
            .add_arg(function!(GS.emr_mom, 3, 0).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom_upper_bound(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::OpaqueEnergyFunction { .. }
        ));
    }

    #[test]
    fn declared_multilinear_functions_preserve_energy_bounds() {
        let expression = FunctionBuilder::new(symbol!("multilinear_energy_function"; Linear))
            .add_arg(function!(GS.emr_mom, 3, mink_index("mu")).as_view())
            .finish();
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn malformed_emr_edge_argument_is_reported() {
        let expression = FunctionBuilder::new(GS.emr_mom)
            .add_arg(Atom::var(symbol!("edge")).as_view())
            .add_arg(mink_index("mu").as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::InvalidEmrEdgeArgument { .. }
        ));
    }

    #[test]
    fn malformed_loop_argument_is_reported() {
        let expression = FunctionBuilder::new(GS.loop_mom)
            .add_arg(Atom::var(symbol!("loop_id")).as_view())
            .add_arg(mink_index("mu").as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::InvalidLoopMomentumArgument { .. }
        ));
    }
}
