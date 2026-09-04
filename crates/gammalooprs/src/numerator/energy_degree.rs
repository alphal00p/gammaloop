use std::{
    collections::{BTreeMap, BTreeSet},
    ops::MulAssign,
};

use linnet::half_edge::involution::EdgeIndex;
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski},
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use thiserror::Error;

use crate::{
    graph::Graph,
    utils::{GS, symbols::UvMomentumProvenanceRole},
    uv::UltravioletGraph,
};

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

    fn add_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            *self.degrees.entry(edge).or_insert(0) += degree;
        }
    }

    fn max_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = (*slot).max(degree);
        }
    }

    fn scale(&mut self, exponent: usize) {
        for degree in self.degrees.values_mut() {
            degree.mul_assign(exponent);
        }
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
    fn try_new<K: Eq>(
        groups: impl IntoIterator<Item = (EdgeIndex, Vec<(usize, K)>)>,
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
            candidates.sort_by_key(|(candidate, _)| *candidate);
            if candidates.is_empty() {
                return Err(
                    EnergyPowerAnalysisError::EmptyEquivalentEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            for pair in candidates.windows(2) {
                if pair[0].0 == pair[1].0 {
                    return Err(
                        EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate {
                            edge: edge.into(),
                            candidate: pair[0].0,
                        },
                    );
                }
            }
            for (candidate, _) in &candidates {
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

            let mut equivalent_classes = Vec::<Vec<(usize, K)>>::new();
            for candidate in candidates {
                if let Some(class) = equivalent_classes
                    .iter_mut()
                    .find(|class| class[0].1 == candidate.1)
                {
                    class.push(candidate);
                } else {
                    equivalent_classes.push(vec![candidate]);
                }
            }
            // More equivalent occurrences minimize the maximal assigned
            // degree. The lowest canonical occurrence fixes equal-size ties.
            equivalent_classes.sort_by(|left, right| {
                right
                    .len()
                    .cmp(&left.len())
                    .then_with(|| left[0].0.cmp(&right[0].0))
            });
            let candidate_ids: Vec<usize> = equivalent_classes
                .remove(0)
                .into_iter()
                .map(|(candidate, _)| candidate)
                .collect();
            by_family.insert(
                EnergyCandidateFamily::Unprovenanced(edge),
                vec![candidate_ids[0]],
            );
            by_family.insert(EnergyCandidateFamily::Fixed(edge), vec![candidate_ids[0]]);
            by_family.insert(
                EnergyCandidateFamily::DenominatorDerived(edge),
                candidate_ids,
            );
        }

        Ok(Self { by_family })
    }

    /// Build candidate sets from the serial denominator occurrences certified
    /// by source-edge topology reconstruction. Algebraic energy equality is
    /// not rediscovered here: members already share one immutable source owner.
    pub(crate) fn try_from_source_occurrences(
        groups: impl IntoIterator<Item = (EdgeIndex, Vec<usize>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        Self::try_new(groups.into_iter().map(|(edge, occurrences)| {
            (
                edge,
                occurrences
                    .into_iter()
                    .map(|occurrence| (occurrence, ()))
                    .collect(),
            )
        }))
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
    factor_assignments: BTreeMap<usize, BTreeMap<EdgeIndex, usize>>,
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
        self.expression.map(&self.factor_assignments, &mut map)
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
    fn degrees(&self) -> EnergyPowerCapMap {
        match self {
            Self::Factor { degrees, .. } => degrees.clone(),
            Self::Add(terms) => terms
                .iter()
                .fold(EnergyPowerCapMap::default(), |mut sum, term| {
                    sum.max_assign(term.degrees());
                    sum
                }),
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => factors
                .iter()
                .fold(EnergyPowerCapMap::default(), |mut product, factor| {
                    product.add_assign(factor.degrees());
                    product
                }),
        }
    }

    fn minimax_assignments(
        &self,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<BTreeMap<usize, BTreeMap<EdgeIndex, usize>>, EnergyPowerAnalysisError> {
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
        ) -> Vec<OwnerState> {
            let mut combined = Vec::new();
            for left in &left {
                for right in &right {
                    let combine_loads = |left: &[usize], right: &[usize]| {
                        left.iter()
                            .zip(right)
                            .map(|(left, right)| {
                                if additive {
                                    (*left).max(*right)
                                } else {
                                    left + right
                                }
                            })
                            .collect()
                    };
                    let energy_loads = combine_loads(&left.energy_loads, &right.energy_loads);
                    let derived_unit_loads =
                        combine_loads(&left.derived_unit_loads, &right.derived_unit_loads);
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
            prune(combined)
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
                        Ok(combine(
                            left,
                            frontier(
                                right,
                                owner,
                                owner_candidates,
                                candidate_positions,
                                candidates,
                            )?,
                            true,
                        ))
                    })
                }
                PlannedEnergyExpression::Mul(factors)
                | PlannedEnergyExpression::MultilinearFunction {
                    arguments: factors, ..
                } => factors.iter().try_fold(zero(), |left, right| {
                    Ok(combine(
                        left,
                        frontier(
                            right,
                            owner,
                            owner_candidates,
                            candidate_positions,
                            candidates,
                        )?,
                        false,
                    ))
                }),
            }
        }

        let mut factor_assignments = BTreeMap::<usize, BTreeMap<EdgeIndex, usize>>::new();
        for (owner, _) in self.degrees().iter() {
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
                        degree: self.degrees().degrees[&owner],
                    },
                );
            }
            let candidate_positions = owner_candidates
                .iter()
                .enumerate()
                .map(|(position, candidate)| (*candidate, position))
                .collect::<BTreeMap<_, _>>();
            let states = frontier(
                self,
                owner,
                &owner_candidates,
                &candidate_positions,
                candidates,
            )?;
            let assignments = states
                .into_iter()
                .min_by_key(|state| {
                    (
                        state.energy_loads.iter().copied().max().unwrap_or(0),
                        state.derived_unit_loads.iter().copied().max().unwrap_or(0),
                        state.energy_loads.iter().sum::<usize>(),
                        state.derived_unit_loads.iter().sum::<usize>(),
                        state.assignments.clone(),
                        state.energy_loads.clone(),
                        state.derived_unit_loads.clone(),
                    )
                })
                .expect("every planned expression has one owner-assignment state")
                .assignments;
            for (factor, candidate) in assignments {
                factor_assignments
                    .entry(factor)
                    .or_default()
                    .insert(owner, candidate);
            }
        }
        Ok(factor_assignments)
    }

    fn exact_degrees(
        &self,
        assignments: &BTreeMap<usize, BTreeMap<EdgeIndex, usize>>,
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
                    *exact.entry(candidate).or_insert(0) += degree;
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
                        *product.entry(candidate).or_insert(0) += degree;
                    }
                }
                Ok(product)
            }
        }
    }

    fn map<E>(
        &self,
        assignments: &BTreeMap<usize, BTreeMap<EdgeIndex, usize>>,
        map: &mut impl FnMut(&Atom, &BTreeMap<EdgeIndex, usize>) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        match self {
            Self::Factor { id, expression, .. } => {
                let empty = BTreeMap::new();
                map(expression, assignments.get(id).unwrap_or(&empty))
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

    pub(crate) fn plan_atom_assignment(
        &self,
        expression: &Atom,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        let expected_degrees = self.analyze_atom(expression)?;
        let mut next_factor_id = 0;
        let planned = self.plan_view(expression.as_view(), &mut next_factor_id)?;
        debug_assert_eq!(planned.degrees(), expected_degrees);

        let factor_assignments = planned.minimax_assignments(candidates)?;
        let exact_bounds = planned.exact_degrees(&factor_assignments)?;

        Ok(EnergyPowerAssignmentPlan {
            expression: planned,
            factor_assignments,
            energy_degree_bounds: exact_bounds.into_iter().collect(),
        })
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
                    product_degree.add_assign(self.analyze_view(factor, mode)?);
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
                base_degree.scale(exponent as usize);
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
                    degree.add_assign(self.analyze_view(function.get(1), mode)?);
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
                        degree.add_assign(argument_degree);
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

    pub(crate) fn automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
        &self,
        numerator: &Atom,
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
        Ok(EnergyPowerAnalyzer::for_physical_emr_edges(active_edges)
            .analyze_atom(numerator)?
            .iter()
            .filter_map(|(edge, degree)| (degree >= min_degree).then_some((edge.into(), degree)))
            .collect())
    }

    pub(crate) fn plan_numerator_energy_assignment_in_atom_excluding(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        EnergyPowerAnalyzer::for_physical_emr_edges(active_edges)
            .plan_atom_assignment(numerator, candidates)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::{
        numerator::energy_degree::{
            EnergyCandidateFamily, EnergyPowerAnalysisError, EnergyPowerAnalyzer,
            EquivalentEnergyCandidates,
        },
        utils::GS,
    };
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, Minkowski},
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
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

    fn planned_candidates(expression: &Atom) -> (Vec<(usize, usize)>, Vec<usize>, Atom) {
        let edge = EdgeIndex(3);
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(expression, &two_equivalent_candidates())
            .unwrap();
        let bounds = plan.energy_degree_bounds().to_vec();
        let mut assignments = Vec::new();
        let mapped = plan
            .map_factors(|factor, factor_assignments| {
                if let Some(candidate) = factor_assignments.get(&edge) {
                    assignments.push(*candidate);
                }
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
        (bounds, assignments, mapped)
    }

    #[test]
    fn original_factors_stay_on_the_canonical_owner_occurrence() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression = (q3.clone() + Atom::num(1)) * (q3 + Atom::num(2));
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 2)]);
        assert_eq!(assignments, vec![10, 10]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn denominator_derived_minimax_balances_multilinear_function_slots() {
        let q3 = uv_owned_q3(true, compact_minkowski_vector());
        let expression = function!(GS.dot, q3.clone(), q3);
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1), (11, 1)]);
        assert_eq!(assignments, vec![10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn denominator_derived_minimax_reuses_capacity_across_additive_branches() {
        let q3 = uv_owned_q3(true, mink_index("mu"));
        let left = FunctionBuilder::new(symbol!("energy_assignment_test::left"; Linear))
            .add_arg(q3.as_view())
            .finish();
        let right = FunctionBuilder::new(symbol!("energy_assignment_test::right"; Linear))
            .add_arg(q3.as_view())
            .finish();
        let expression = left + right;
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1)]);
        assert_eq!(assignments, vec![10, 10]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn denominator_derived_minimax_repeats_power_bases_without_expansion() {
        let q3 = uv_owned_q3(true, mink_index("mu"));
        let expression =
            (q3 + Atom::var(symbol!("energy_assignment_test::shift"))).pow(Atom::num(2));
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1), (11, 1)]);
        assert_eq!(assignments, vec![10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn production_minimax_reserves_base_and_splits_units_in_a_saturated_branch() {
        let edge = EdgeIndex(3);
        let q3 = uv_owned_q3(true, GS.cind(0));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg(q3.clone().pow(2).as_view())
            .finish();
        let expression = denominator.clone().pow(2) + q3.clone().pow(2);
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(
                &expression,
                &production_candidates_with_two_derived_copies(),
            )
            .unwrap();

        assert_eq!(plan.energy_degree_bounds(), &[(11, 2), (12, 2)]);
        let mut denominator_assignments = Vec::new();
        let mut unit_assignments = Vec::new();
        let _ = plan
            .map_factors(|factor, assignments| {
                if factor == &denominator {
                    denominator_assignments.push(assignments[&edge]);
                } else if factor == &q3 {
                    unit_assignments.push(assignments[&edge]);
                }
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
        denominator_assignments.sort_unstable();
        unit_assignments.sort_unstable();
        assert_eq!(denominator_assignments, vec![11, 12]);
        assert_eq!(unit_assignments, vec![11, 12]);
    }

    #[test]
    fn positive_denominator_completion_stays_in_one_occurrence_frame() {
        let edge = EdgeIndex(3);
        let q3 = uv_owned_q3(true, GS.cind(0));
        let spatial_mass = Atom::var(symbol!("energy_assignment_test::spatial_mass"));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg((q3.pow(2) - &spatial_mass).as_view())
            .finish();
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(&denominator, &two_equivalent_candidates())
            .unwrap();

        assert_eq!(plan.energy_degree_bounds(), &[(10, 2)]);
        let mut assignments_seen = Vec::new();
        let mapped = plan
            .map_factors(|factor, assignments| {
                assert_eq!(factor, &denominator);
                assignments_seen.push(assignments[&edge]);
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
        assert_eq!(assignments_seen, vec![10]);
        assert_eq!(mapped, denominator);
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
        let assert_assignment = |energy: Atom,
                                 candidates: EquivalentEnergyCandidates,
                                 expected_bound: (usize, usize)| {
            let expression = denominator(&energy);
            let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
                .plan_atom_assignment(&expression, &candidates)
                .unwrap();
            assert_eq!(plan.energy_degree_bounds(), &[expected_bound]);
            let mut assignments = Vec::new();
            let mapped = plan
                .map_factors(|factor, factor_assignments| {
                    assert_eq!(factor, &expression);
                    assignments.push(factor_assignments[&edge]);
                    Ok::<_, ()>(factor.clone())
                })
                .unwrap();
            assert_eq!(assignments, vec![expected_bound.0]);
            assert_eq!(mapped, expression);
        };

        let base_only = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            edge,
            10,
            Vec::new(),
        )])
        .unwrap();
        assert_assignment(function!(GS.emr_mom, 3, GS.cind(0)), base_only, (10, 2));

        assert_assignment(
            uv_owned_q3(false, GS.cind(0)),
            production_candidates_with_two_derived_copies(),
            (10, 2),
        );
        assert_assignment(
            uv_owned_q3(true, GS.cind(0)),
            production_candidates_with_two_derived_copies(),
            (11, 2),
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
        let edge = EdgeIndex(3);
        let q3 = uv_owned_q3(true, GS.cind(0));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg(q3.pow(2).as_view())
            .finish();
        let expression = denominator.pow(2);
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(&expression, &two_equivalent_candidates())
            .unwrap();

        assert_eq!(plan.energy_degree_bounds(), &[(10, 2), (11, 2)]);
        let mut assignments_seen = Vec::new();
        let mapped = plan
            .map_factors(|factor, assignments| {
                assert_eq!(factor, &denominator);
                assignments_seen.push(assignments[&edge]);
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
        assignments_seen.sort_unstable();
        assert_eq!(assignments_seen, vec![10, 11]);
        assert_eq!(mapped, expression);
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
    fn denominator_derived_minimax_uses_canonical_lexicographic_tie_break() {
        let q3 = uv_owned_q3(true, mink_index("mu"));
        let expression =
            (q3.clone() + Atom::num(1)) * (q3.clone() + Atom::num(2)) * (q3 + Atom::num(3));
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 2), (11, 1)]);
        assert_eq!(assignments, vec![10, 10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn owner_local_minimax_accounts_for_fixed_load_across_factorized_sums() {
        let edge = EdgeIndex(3);
        let hard = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let tagged = |denominator_derived| {
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(3).as_view(),
                    denominator_derived,
                    hard.as_view(),
                ))
                .add_arg(mink_index("mu"))
                .finish()
        };
        let fixed = tagged(false);
        let derived = tagged(true);
        let a = Atom::var(symbol!("owner_local_minimax::a"));
        let b = Atom::var(symbol!("owner_local_minimax::b"));
        let candidate_10 = Atom::var(symbol!("owner_local_minimax::candidate_10"));
        let candidate_11 = Atom::var(symbol!("owner_local_minimax::candidate_11"));

        for (expression, expected) in [
            (&fixed * &derived, &candidate_10 * &candidate_11),
            (
                &fixed * (&a * &derived + &b * &derived),
                &candidate_10 * (&a * &candidate_11 + &b * &candidate_11),
            ),
            (
                &a * &fixed * &derived + &b * derived.pow(2),
                &a * &candidate_10 * &candidate_11 + &b * &candidate_10 * &candidate_11,
            ),
        ] {
            let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
                .plan_atom_assignment(&expression, &two_equivalent_candidates())
                .unwrap();
            assert_eq!(plan.energy_degree_bounds(), &[(10, 1), (11, 1)]);
            let mapped = plan
                .map_factors(|factor, assignments| {
                    Ok::<_, ()>(match assignments.get(&edge) {
                        Some(10) => candidate_10.clone(),
                        Some(11) => candidate_11.clone(),
                        Some(candidate) => panic!("unexpected candidate {candidate}"),
                        None => factor.clone(),
                    })
                })
                .unwrap();
            assert_eq!(mapped, expected);
        }
    }

    #[test]
    fn denominator_derived_minimax_reuses_the_additive_envelope() {
        let edge = EdgeIndex(3);
        let q3 = uv_owned_q3(true, mink_index("mu"));
        let quadratic = q3.clone().pow(2);
        let expression = q3.pow(3) + &quadratic;
        let candidate_10 = Atom::var(symbol!("energy_assignment_test::candidate_10"));
        let candidate_11 = Atom::var(symbol!("energy_assignment_test::candidate_11"));
        let candidate_12 = Atom::var(symbol!("energy_assignment_test::candidate_12"));

        for (
            candidate_ids,
            expected_bounds,
            expected_cubic,
            expected_combined_quadratic,
            expected_standalone_quadratic,
        ) in [
            (
                vec![10, 11],
                vec![(10, 2), (11, 1)],
                candidate_10.clone().pow(2) * &candidate_11,
                candidate_10.clone().pow(2),
                &candidate_10 * &candidate_11,
            ),
            (
                vec![10, 11, 12],
                vec![(10, 1), (11, 1), (12, 1)],
                &candidate_10 * &candidate_11 * &candidate_12,
                &candidate_10 * &candidate_11,
                &candidate_10 * &candidate_11,
            ),
        ] {
            let candidates =
                EquivalentEnergyCandidates::try_from_source_occurrences([(edge, candidate_ids)])
                    .unwrap();
            let combined = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
                .plan_atom_assignment(&expression, &candidates)
                .unwrap();
            let standalone = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
                .plan_atom_assignment(&quadratic, &candidates)
                .unwrap();
            assert_eq!(combined.energy_degree_bounds(), expected_bounds);

            let map_candidate = |candidate: usize| match candidate {
                10 => candidate_10.clone(),
                11 => candidate_11.clone(),
                12 => candidate_12.clone(),
                _ => unreachable!("the test uses only three canonical candidates"),
            };
            let mut map =
                |factor: &Atom, assignments: &BTreeMap<EdgeIndex, usize>| -> Result<Atom, ()> {
                    Ok(assignments
                        .get(&edge)
                        .map(|candidate| map_candidate(*candidate))
                        .unwrap_or_else(|| factor.clone()))
                };
            let mapped_combined = combined.map_factors(&mut map).unwrap();
            let mapped_quadratic = standalone.map_factors(&mut map).unwrap();
            assert_eq!(mapped_quadratic, expected_standalone_quadratic);
            assert_eq!(
                mapped_combined,
                expected_cubic + expected_combined_quadratic,
                "an additive branch may reuse capacity already present in the global minimax envelope",
            );
        }
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
